function [S, reason] = r0_long_core(ctx, S, Fobj, ops, opt, nits, cfg)
% R0_LONG_CORE  Resumable version of tests/r0_pilot_track.m (itself a copy of
% r0_opt_one_track / articleB_r0_opt_one: GCMMA, thesis dB objective scaling, windowed stopping
% rule). The optimiser logic is unchanged; what is new:
%   * the whole optimiser state lives in S (current design S.t, GCMMA state S.st = asymptotes,
%     xold1/xold2, raa0/raa; dB reference S.p0; GCMMA iteration counter S.git; stopping-rule
%     counter S.nhold and objective window S.f0seg; segment best), so a run split into chunks
%     with checkpoints is the same computation as an uninterrupted run (checked bitwise by the
%     validate job of r0_long_runs);
%   * at every outer iteration the exact, sinc and ERP powers of the TARGET cell (normalised so
%     that the uniform design has J = Jref) are recorded from the same state solve (LU cache of
%     r0_long_eval, ~1% extra cost);
%   * optional move limit cfg.move (box [t-move, t+move] intersected with the bounds, a trust
%     region; Svanberg's code has no other step cap: opt.move is unused by gcmma_solver);
%   * optional switching rule cfg.switch (multi-step strategy iv, see r0_long_runs).
% Fobj: struct('factor', F) = the objective factor of this segment, already multiplied by
% sqrt(cstart) (pilot warm-start rule: J(start) = Jref). S.cstart holds that factor.
% reason: 'nits' | 'converged' | 'switch' | 'stopfile'
if nargin<7, cfg=struct(); end
if ~isfield(cfg,'move'), cfg.move=0; end
if ~isfield(cfg,'switch'), cfg.switch=[]; end
if ~isfield(cfg,'stopfile'), cfg.stopfile=''; end
ne=ctx.ne;
tmin_e = opt.tmin*ones(ne,1); tmax_e = opt.tmax*ones(ne,1);
vw = opt.vweights(:);
V0 = opt.volfrac*sum(vw.*tmax_e);
gopt=opt; gopt.kmax=20; gopt.growth=1.30; gopt.raa0eps=1e-3; gopt.raaeps=1e-3;
reason='nits';
volmax=opt.volfrac*(1+1e-5);
if ~isfield(S,'Jxfbest'), S.Jxfbest=inf; S.txfbest=[]; S.itxfbest=0; end
if ~isfield(S,'Jfbest_seg'), S.Jfbest_seg=inf; S.tfbest_seg=[]; S.itfbest_seg=0; end
for k=1:nits
  t0=tic;
  t=S.t;
  [J,dJ]=r0_long_eval(ctx,t,Fobj);
  Jx=r0_long_eval(ctx,t,ops.X); Js=r0_long_eval(ctx,t,ops.S); Je=r0_long_eval(ctx,t,ops.E);
  cosSX=NaN;
  if ~isempty(cfg.switch)
    [~,gs]=r0_long_eval(ctx,t,ops.S); [~,gx]=r0_long_eval(ctx,t,ops.X);
    fr=(t>opt.tmin+1e-6)&(t<opt.tmax-1e-6);
    gs=gs/Js; gx=gx/Jx; gs(~fr)=0; gx(~fr)=0; gs(fr)=gs(fr)-mean(gs(fr)); gx(fr)=gx(fr)-mean(gx(fr));
    cosSX=gs'*gx/(norm(gs)*norm(gx));
    sw=cfg.switch;
    if cosSX<sw.cosmin, S.swcount=S.swcount+1; else, S.swcount=0; end
    nH=numel(S.H.Jx); ns=S.stage_it;
    rise = ns>=sw.lag && Jx > S.H.Jx(nH-sw.lag+1);
    if S.swcount>=sw.ncons || rise || ns>=sw.maxits
      if S.swcount>=sw.ncons, S.swreason=sprintf('cos<%.2f for %d its',sw.cosmin,sw.ncons);
      elseif rise, S.swreason=sprintf('exact J rose over %d its',sw.lag);
      else, S.swreason='max sinc iterations'; end
      reason='switch'; return
    end
  end
  S.git=S.git+1; it=S.git;
  gv=sum(vw.*t)/V0-1; dgv=vw/V0;
  Jdb=opt.Abar+opt.dbfac*log10(max(J,eps));
  if isempty(S.p0), S.p0=Jdb; end
  p0=S.p0;
  f0s=100*Jdb/p0; df0s=100*(opt.dbfac/log(10))*dJ/max(J,eps)/p0;
  S.f0seg(end+1)=f0s;
  gopt.eval=@(tq) deal(obj_scaled(tq,ctx,Fobj,opt,p0), sum(vw.*tq)/V0-1);
  if cfg.move>0
    lo=max(tmin_e,t-cfg.move); hi=min(tmax_e,t+cfg.move);
  else
    lo=tmin_e; hi=tmax_e;
  end
  lastwarn('','');
  [tnew,S.st,info]=gcmma_solver(t,lo,hi,f0s,df0s,gv,dgv,S.st,it,gopt);
  [~,wid]=lastwarn; if ~isempty(wid), S.nwarn=S.nwarn+1; end
  chg=max(abs(tnew(:)-t));
  if J<S.Jbest_seg, S.Jbest_seg=J; S.tbest_seg=t; end
  S.it=S.it+1; S.stage_it=S.stage_it+1;
  H=S.H;
  H.it(end+1)=S.it; H.stage(end+1)=S.stage; H.rung(end+1)=cfg.rung;
  H.Jown(end+1)=J/S.cstart; H.Jx(end+1)=Jx; H.Js(end+1)=Js; H.Je(end+1)=Je;
  H.vol(end+1)=sum(vw.*t)/sum(vw.*tmax_e); H.chg(end+1)=chg; H.inner(end+1)=info.innerit;
  H.cos(end+1)=cosSX; H.f0(end+1)=f0s; H.time(end+1)=toc(t0);
  S.H=H; S.walltime=S.walltime+H.time(end);
  if Jx<S.Jxbest, S.Jxbest=Jx; S.txbest=t; S.itxbest=S.it; end
  if H.vol(end)<=volmax
    if Jx<S.Jxfbest, S.Jxfbest=Jx; S.txfbest=t; S.itxfbest=S.it; end
    if J<S.Jfbest_seg, S.Jfbest_seg=J; S.tfbest_seg=t; S.itfbest_seg=S.it; end
  end
  if mod(S.it,10)==0, S.T(:,end+1)=t; S.Tit(end+1)=S.it; end
  if mod(S.stage_it,25)==0 || S.stage_it==1
    fprintf('  [%c %s] it %5d (stage %d, %4d): Jown=%.5e f0=%+9.4f  exact %+8.3f dB  sinc %+8.3f dB  erp %+8.3f dB  vol=%.4f chg=%.2e inner=%d\n',...
      cfg.rung,cfg.src,S.it,S.stage,S.stage_it,J/S.cstart,f0s,10*log10(Jx/S.Jref),10*log10(Js/S.Jref),10*log10(Je/S.Jref),H.vol(end),chg,info.innerit);
  end
  S.t=tnew(:);
  if opt.stop && numel(S.f0seg)>=2*opt.window+1
    W=opt.window; fo=S.f0seg;
    fstat=abs(mean(fo(end-W+1:end))-mean(fo(end-2*W+1:end-W)))/max(abs(mean(fo(end-2*W+1:end-W))),eps);
    xstat=chg/max(opt.tmax-opt.tmin,eps);
    if xstat<opt.xchgtol && fstat<opt.ftol && gv<=opt.feastol, S.nhold=S.nhold+1; else, S.nhold=0; end
    if S.nhold>=opt.nhold
      fprintf('  [%c %s] STOP at it %d (stage it %d): converged (xstat %.2e, fstat %.2e)\n',cfg.rung,cfg.src,S.it,S.stage_it,xstat,fstat);
      reason='converged'; return
    end
  end
  if ~isempty(cfg.stopfile) && mod(S.it,25)==0 && exist(cfg.stopfile,'file')
    reason='stopfile'; return
  end
end
end

function f0q = obj_scaled(tq,ctx,Fobj,opt,p0)
Jq=r0_long_eval(ctx,tq(:),Fobj);
f0q=100*(opt.Abar+opt.dbfac*log10(max(Jq,eps)))/p0;
end
