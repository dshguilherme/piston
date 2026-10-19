function [S, reason] = r0_robust_core(ctx, S, Fobj, ops, opt, nits, cfg)
% R0_ROBUST_CORE  Copy of tests/r0_long_core.m (resumable GCMMA loop of the pilot tracker: thesis dB
% objective scaling, windowed stopping rule, full optimiser state in S) with the objective replaced
% by the robust sample-average objective of r0_robust_eval:
%   J(t) = mean_i P(t (1 + sigma xi_i))   [or mean + k std],  xi_i = cfg.R.Xi(:,i) FIXED for the
%   whole run (common random numbers), gradient by the per-sample adjoint and the chain rule.
% Fobj: struct('factor', F [, 'iwZ', fluid loading]) already multiplied by sqrt(S.cstart) so that
% J(start) = Jref (pilot warm-start rule, keeps the dB scaling's p0 > 0).
% Recorded per outer iteration (S.H): the robust objective (unscaled by cstart) and its sample std,
% and the NOMINAL (sigma = 0) exact and sinc powers of the current design (normalised factors
% ops.X, ops.S; one extra LU through r0_long_eval).
% The switching rule of r0_long_core is removed; the move limit cfg.move is kept (0 = none).
% reason: 'nits' | 'converged' | 'stopfile'
if ~isfield(cfg,'move'), cfg.move=0; end
if ~isfield(cfg,'stopfile'), cfg.stopfile=''; end
ne=ctx.ne;
tmin_e = opt.tmin*ones(ne,1); tmax_e = opt.tmax*ones(ne,1);
vw = opt.vweights(:);
V0 = opt.volfrac*sum(vw.*tmax_e);
gopt=opt; gopt.kmax=20; gopt.growth=1.30; gopt.raa0eps=1e-3; gopt.raaeps=1e-3;
reason='nits';
volmax=opt.volfrac*(1+1e-5);
R=cfg.R;
for k=1:nits
  t0=tic;
  t=S.t;
  [J,dJ,inf1]=r0_robust_eval(ctx,t,Fobj,R);
  Jx=r0_long_eval(ctx,t,ops.X); Js=r0_long_eval(ctx,t,ops.S);
  S.git=S.git+1; it=S.git;
  gv=sum(vw.*t)/V0-1; dgv=vw/V0;
  Jdb=opt.Abar+opt.dbfac*log10(max(J,eps));
  if isempty(S.p0), S.p0=Jdb; end
  p0=S.p0;
  f0s=100*Jdb/p0; df0s=100*(opt.dbfac/log(10))*dJ/max(J,eps)/p0;
  S.f0seg(end+1)=f0s;
  gopt.eval=@(tq) deal(obj_scaled(tq,ctx,Fobj,R,opt,p0), sum(vw.*tq)/V0-1);
  if cfg.move>0
    lo=max(tmin_e,t-cfg.move); hi=min(tmax_e,t+cfg.move);
  else
    lo=tmin_e; hi=tmax_e;
  end
  lastwarn('','');
  [tnew,S.st,info]=gcmma_solver(t,lo,hi,f0s,df0s,gv,dgv,S.st,it,gopt);
  [~,wid]=lastwarn; if ~isempty(wid), S.nwarn=S.nwarn+1; end
  chg=max(abs(tnew(:)-t));
  S.it=S.it+1;
  H=S.H;
  H.it(end+1)=S.it; H.Jrob(end+1)=J/S.cstart; H.Jstd(end+1)=inf1.std/S.cstart;
  H.Jx(end+1)=Jx; H.Js(end+1)=Js;
  H.vol(end+1)=sum(vw.*t)/sum(vw.*tmax_e); H.chg(end+1)=chg; H.inner(end+1)=info.innerit;
  H.f0(end+1)=f0s; H.time(end+1)=toc(t0); H.gmres(end+1)=mean(inf1.gmres_it); H.nfb(end+1)=inf1.nfallback;
  S.H=H; S.walltime=S.walltime+H.time(end);
  if H.vol(end)<=volmax && J/S.cstart<S.Jrobbest
    S.Jrobbest=J/S.cstart; S.trobbest=t; S.itrobbest=S.it;
  end
  if mod(S.it,10)==0, S.T(:,end+1)=t; S.Tit(end+1)=S.it; end
  if mod(S.it,25)==0 || S.it==1
    fprintf('  [%c %s N%d s%.3g] it %5d: Jrob=%.5e (std %.2e) f0=%+9.4f  nominal exact %+8.3f dB  sinc %+8.3f dB  robust %+8.3f dB  vol=%.4f chg=%.2e inner=%d %.1fs\n',...
      cfg.rung,cfg.src,size(R.Xi,2),R.sigma,S.it,J/S.cstart,inf1.std/S.cstart,f0s,10*log10(Jx/S.Jref),10*log10(Js/S.Jref),10*log10(J/S.cstart/S.Jref),H.vol(end),chg,info.innerit,H.time(end));
  end
  S.t=tnew(:);
  if opt.stop && numel(S.f0seg)>=2*opt.window+1
    W=opt.window; fo=S.f0seg;
    fstat=abs(mean(fo(end-W+1:end))-mean(fo(end-2*W+1:end-W)))/max(abs(mean(fo(end-2*W+1:end-W))),eps);
    xstat=chg/max(opt.tmax-opt.tmin,eps);
    if xstat<opt.xchgtol && fstat<opt.ftol && gv<=opt.feastol, S.nhold=S.nhold+1; else, S.nhold=0; end
    if S.nhold>=opt.nhold
      fprintf('  [%c %s] STOP at it %d: converged (xstat %.2e, fstat %.2e)\n',cfg.rung,cfg.src,S.it,xstat,fstat);
      reason='converged'; return
    end
  end
  if ~isempty(cfg.stopfile) && mod(S.it,25)==0 && exist(cfg.stopfile,'file')
    reason='stopfile'; return
  end
end
end

function f0q = obj_scaled(tq,ctx,Fobj,R,opt,p0)
Jq=r0_robust_eval(ctx,tq(:),Fobj,R);
f0q=100*(opt.Abar+opt.dbfac*log10(max(Jq,eps)))/p0;
end
