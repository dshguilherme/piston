function [tbest, hist, meta] = r0_pilot_track(ctx, metric, R0op, R0mode, opt)
% R0_PILOT_TRACK  Copy of tests/r0_opt_one_track.m, itself a verbatim copy of
% tests/articleB_r0_opt_one.m (GCMMA, thesis dB objective scaling, windowed stopping rule)
% that stores the design at every outer iteration in hist.T (ne x niter).
% Pilot change (computational only): the objective is evaluated by r0_pilot_eval, which is
% numerically the same as articleB_r0_eval(...,'quadvel',R0op,'consistent') but uses one LU per
% evaluation and skips the adjoint in the GCMMA inner loop. The optimiser logic is unchanged.
if nargin<5, opt=struct(); end
if ~isfield(opt,'tmin'),    opt.tmin=0.25; end
if ~isfield(opt,'tmax'),    opt.tmax=1.0;  end
if ~isfield(opt,'volfrac'), opt.volfrac=0.5; end
if ~isfield(opt,'nmax'),    opt.nmax=60;   end
if ~isfield(opt,'stop'),    opt.stop=true; end
if ~isfield(opt,'ftol'),    opt.ftol=1e-4; end
if ~isfield(opt,'xchgtol'), opt.xchgtol=1e-4; end
if ~isfield(opt,'nhold'),   opt.nhold=3;   end
if ~isfield(opt,'feastol'), opt.feastol=1e-6; end
if ~isfield(opt,'window'),  opt.window=10;  end
if ~isfield(opt,'Abar'),    opt.Abar=100;   end
if ~isfield(opt,'dbfac'),   opt.dbfac=10;   end
if ~isfield(opt,'vweights'), opt.vweights=ones(ctx.ne,1); end
opt.vweights=opt.vweights(:);
ne=ctx.ne;
tmin_e = opt.tmin*ones(ne,1); tmax_e = opt.tmax*ones(ne,1);
if isfield(opt,'tinit') && ~isempty(opt.tinit)
  t=opt.tinit(:);
  assert(numel(t)==ne && all(t>=tmin_e-1e-12) && all(t<=tmax_e+1e-12),'tinit must satisfy design bounds');
  t=min(max(t,tmin_e),tmax_e);
else
  t=opt.volfrac*opt.tmax*ones(ne,1);
end
V0 = opt.volfrac*sum(opt.vweights.*tmax_e);
st=[]; Jbest=inf; tbest=t; p0=[]; nholdcnt=0;
hist=struct('J',[],'Jbest',[],'f0',[],'vol',[],'chg',[],'times',[],'innerit',[],'T',zeros(ne,0));
gopt=opt; gopt.kmax=20; gopt.growth=1.30; gopt.raa0eps=1e-3; gopt.raaeps=1e-3;
stopReason='nmax'; ts=tic;
for it=1:opt.nmax
  t0=tic;
  [J,dJ]=r0_pilot_eval(ctx,t,R0op);   % = articleB_r0_eval(ctx,t,'quadvel',R0op,'consistent')
  hist.T(:,end+1)=t;
  gv=sum(opt.vweights.*t)/V0-1; dgv=opt.vweights/V0;
  Jdb=opt.Abar+opt.dbfac*log10(max(J,eps));
  if isempty(p0), p0=Jdb; end
  f0s=100*Jdb/p0; df0s=100*(opt.dbfac/log(10))*dJ/max(J,eps)/p0;
  hist.f0(end+1)=f0s;
  gopt.eval=@(tq) deal(obj_scaled(tq,ctx,metric,R0op,R0mode,opt,p0), sum(opt.vweights.*tq)/V0-1);
  [tnew,st,info]=gcmma_solver(t,tmin_e,tmax_e,f0s,df0s,gv,dgv,st,it,gopt);
  hist.innerit(it)=info.innerit;
  chg=max(abs(tnew(:)-t));
  if J<Jbest, Jbest=J; tbest=t; end
  hist.J(end+1)=J; hist.Jbest(end+1)=Jbest;
  hist.vol(end+1)=sum(opt.vweights.*t)/sum(opt.vweights.*tmax_e); hist.chg(end+1)=chg; hist.times(end+1)=toc(t0);
  if mod(it,25)==0 || it==1
    fprintf('  [%-11s|%-10s] it %3d: J=%.5e  f0=%+.4f  Jbest=%.5e  vol=%.4f  chg=%.2e  inner=%d\n',...
      metric,R0mode,it,J,f0s,Jbest,hist.vol(end),chg,info.innerit);
  end
  t=tnew(:);
  if opt.stop && numel(hist.J)>=2*opt.window+1
    W=opt.window; fo=hist.f0;
    fstat=abs(mean(fo(end-W+1:end))-mean(fo(end-2*W+1:end-W)))/max(abs(mean(fo(end-2*W+1:end-W))),eps);
    xstat=chg/max(opt.tmax-opt.tmin,eps);
    if xstat<opt.xchgtol && fstat<opt.ftol && gv<=opt.feastol, nholdcnt=nholdcnt+1; else, nholdcnt=0; end
    if nholdcnt>=opt.nhold
      stopReason='converged';
      fprintf('  [%s|%s] STOP at it %d: converged (xstat %.2e, fstat %.2e)\n',metric,R0mode,it,xstat,fstat);
      break;
    end
  end
end
hist.stopReason=stopReason; hist.walltime=toc(ts); hist.numeval=sum(hist.innerit);
hist.metric=metric; hist.R0mode=R0mode; hist.p0=p0;
[~,ib]=min(hist.J); hist.iters_to_best=ib;
meta=struct('metric',metric,'R0mode',R0mode,'niter',numel(hist.J),'numeval',sum(hist.innerit),...
  'walltime',hist.walltime,'final_vol',sum(opt.vweights.*tbest)/sum(opt.vweights.*tmax_e),...
  'stopReason',stopReason,'Jbest',Jbest);
end

function f0q = obj_scaled(tq,ctx,metric,R0op,R0mode,opt,p0)
Jq=r0_pilot_eval(ctx,tq(:),R0op);    % objective only (no adjoint) for the GCMMA inner loop
f0q=100*(opt.Abar+opt.dbfac*log10(max(Jq,eps)))/p0;
end
