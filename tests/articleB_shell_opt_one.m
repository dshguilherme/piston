function [tbest, hist, meta] = articleB_shell_opt_one(ctx, metric, opt)
% ARTICLEB_SHELL_OPT_ONE  Thickness topology optimization of the curved KL
% shell for one harmonic metric, using the VERIFIED src/gcmma_solver.m
% (Svanberg GCMMA) with the THESIS-EXACT dB objective scaling (Olavo thesis
% eqs. 6.113-6.115) and the driver's windowed stopping criterion.
if nargin<3, opt=struct(); end
if ~isfield(opt,'tmin'),    opt.tmin=0.25; end
if ~isfield(opt,'tmax'),    opt.tmax=1.0;  end
if ~isfield(opt,'volfrac'), opt.volfrac=0.5; end
if ~isfield(opt,'nmax'),    opt.nmax=100;  end
if ~isfield(opt,'stop'),    opt.stop=true; end
if ~isfield(opt,'ftol'),    opt.ftol=1e-4; end
if ~isfield(opt,'xchgtol'), opt.xchgtol=1e-4; end
if ~isfield(opt,'nhold'),   opt.nhold=3;   end
if ~isfield(opt,'feastol'), opt.feastol=1e-6; end
if ~isfield(opt,'window'),  opt.window=10;  end
if ~isfield(opt,'Abar'),    opt.Abar=100;   end
if ~isfield(opt,'dbfac'),   opt.dbfac=10;   end

ne=ctx.ne;
t = opt.volfrac*opt.tmax*ones(ne,1);
tmin_e = opt.tmin*ones(ne,1); tmax_e = opt.tmax*ones(ne,1);
V0 = opt.volfrac*sum(tmax_e);
st=[]; Jbest=inf; tbest=t; p0=[]; nholdcnt=0;
hist=struct('J',[],'Jbest',[],'f0',[],'vol',[],'chg',[],'times',[],'innerit',[]);
gopt=opt; gopt.kmax=20; gopt.growth=1.30; gopt.raa0eps=1e-3; gopt.raaeps=1e-3;
stopReason='nmax'; ts=tic;
for it=1:opt.nmax
  t0=tic;
  [J,dJ]=articleB_shell_eval(ctx,t,metric);
  gv=sum(t)/V0-1; dgv=(1/V0)*ones(ne,1);
  Jdb=opt.Abar+opt.dbfac*log10(max(J,eps));
  if isempty(p0), p0=Jdb; end
  f0s=100*Jdb/p0; df0s=100*(opt.dbfac/log(10))*dJ/max(J,eps)/p0;
  hist.f0(end+1)=f0s;
  gopt.eval=@(tq) deal(obj_scaled(tq,ctx,metric,opt,V0,p0), sum(tq)/V0-1);
  [tnew,st,info]=gcmma_solver(t,tmin_e,tmax_e,f0s,df0s,gv,dgv,st,it,gopt);
  hist.innerit(it)=info.innerit;
  chg=max(abs(tnew(:)-t));
  if J<Jbest, Jbest=J; tbest=t; end
  hist.J(end+1)=J; hist.Jbest(end+1)=Jbest;
  hist.vol(end+1)=sum(t)/sum(tmax_e); hist.chg(end+1)=chg; hist.times(end+1)=toc(t0);
  if mod(it,5)==0 || it==1
    fprintf('  [%-11s] it %3d: J=%.5e  f0(dB)=%+.4f  Jbest=%.5e  vol=%.4f  chg=%.2e  inner=%d\n',...
      metric,it,J,f0s,Jbest,hist.vol(end),chg,info.innerit);
  end
  t=tnew(:);
  if opt.stop && numel(hist.J)>=2*opt.window+1
    W=opt.window; fo=hist.f0;
    fstat=abs(mean(fo(end-W+1:end))-mean(fo(end-2*W+1:end-W)))/max(abs(mean(fo(end-2*W+1:end-W))),eps);
    xstat=chg/max(opt.tmax-opt.tmin,eps);
    if xstat<opt.xchgtol && fstat<opt.ftol && gv<=opt.feastol, nholdcnt=nholdcnt+1; else, nholdcnt=0; end
    if nholdcnt>=opt.nhold
      stopReason='converged';
      fprintf('  [%s] STOP at it %d: converged (xstat %.2e, fstat %.2e)\n',metric,it,xstat,fstat);
      break;
    end
  end
end
hist.stopReason=stopReason; hist.walltime=toc(ts); hist.numeval=sum(hist.innerit);
hist.metric=metric; hist.p0=p0; hist.tmin=opt.tmin; hist.tmax=opt.tmax;
[~,ib]=min(hist.J); hist.iters_to_best=ib;
tl10=hist.J(max(1,end-9):end); hist.tailspread10=max(tl10)/max(min(tl10),eps);
hist.tailspread=NaN; if numel(hist.J)>=50, tl=hist.J(end-49:end); hist.tailspread=max(tl)/max(min(tl),eps); end
meta=struct('metric',metric,'niter',numel(hist.J),'numeval',sum(hist.innerit),...
  'walltime',hist.walltime,'final_vol',sum(tbest)/sum(tmax_e),...
  'tailspread10',hist.tailspread10,'tailspread',hist.tailspread,'stopReason',stopReason,'Jbest',Jbest);
end

function f0q = obj_scaled(tq,ctx,metric,opt,V0,p0)
[Jq,~]=articleB_shell_eval(ctx,tq(:),metric);
f0q=100*(opt.Abar+opt.dbfac*log10(max(Jq,eps)))/p0;
end
