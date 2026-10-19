function [tbest, hist, meta] = articleB_band_opt_one(cts, metric, band, opt)
% ARTICLEB_BAND_OPT_ONE  Thickness topology optimization of the curved KL shell
% against the FREQUENCY-BAND objective (src/harmonic_band.m), using the verified
% src/gcmma_solver.m (Svanberg GCMMA) with the thesis-exact dB scaling and the
% driver's windowed stopping criterion.
%
% This is the band analogue of articleB_shell_opt_one: identical optimizer,
% scaling, stopping and bookkeeping, only the objective is aggregated over
% band.f (band.Nf samples) instead of a single frequency. Because the aggregate
% penalises the WHOLE band, an optimizer cannot "win" by relocating a resonance
% peak out of the excitation frequency -- which is exactly the de-aliasing claim.
%
%   cts   cell of per-frequency contexts (articleB_band_ctxs)
%   band  struct: band.Nf, band.f, band.f_lo, band.f_hi
%   opt.bandmode  'mean' | 'ks'  (default 'ks' = worst-case band maximum)
%   opt.rho       KS parameter   (default 60)
if nargin<4, opt=struct(); end
if ~isfield(opt,'tmin'),     opt.tmin=0.25; end
if ~isfield(opt,'tmax'),     opt.tmax=1.0;  end
if ~isfield(opt,'volfrac'),  opt.volfrac=0.5; end
if ~isfield(opt,'nmax'),     opt.nmax=100;  end
if ~isfield(opt,'stop'),     opt.stop=true; end
if ~isfield(opt,'ftol'),     opt.ftol=1e-4; end
if ~isfield(opt,'xchgtol'),  opt.xchgtol=1e-4; end
if ~isfield(opt,'nhold'),    opt.nhold=3;   end
if ~isfield(opt,'feastol'),  opt.feastol=1e-6; end
if ~isfield(opt,'window'),   opt.window=8;  end
if ~isfield(opt,'Abar'),     opt.Abar=100;   end
if ~isfield(opt,'dbfac'),    opt.dbfac=10;   end
if ~isfield(opt,'bandmode'), opt.bandmode='ks'; end
if ~isfield(opt,'rho'),      opt.rho=60;    end

ne=cts{1}.ne;
t = opt.volfrac*opt.tmax*ones(ne,1);
tmin_e = opt.tmin*ones(ne,1); tmax_e = opt.tmax*ones(ne,1);
V0 = opt.volfrac*sum(tmax_e);
st=[]; Jbest=inf; tbest=t; p0=[]; nholdcnt=0; sxa=0; sfa=0;
hist=struct('J',[],'Jbest',[],'f0',[],'vol',[],'chg',[],'times',[],'innerit',[],...
            'bandmax',[],'bandavg',[],'fmax',[]);
gopt=opt; gopt.kmax=20; gopt.growth=1.30; gopt.raa0eps=1e-3; gopt.raaeps=1e-3;
stopReason='nmax'; ts=tic;
efun=@(tt,k) articleB_band_eval(cts,k,tt,metric);
for it=1:opt.nmax
  t0=tic;
  [J,dJ,bi]=harmonic_band(band,t,efun,opt.bandmode,opt.rho);
  gv=sum(t)/V0-1; dgv=(1/V0)*ones(ne,1);
  Jdb=opt.Abar+opt.dbfac*log10(max(J,eps));
  if isempty(p0), p0=Jdb; end
  f0s=100*Jdb/p0; df0s=100*(opt.dbfac/log(10))*dJ/max(J,eps)/p0;
  hist.f0(end+1)=f0s;
  gopt.eval=@(tq) deal(obj_scaled(tq,cts,metric,band,opt,p0), sum(tq)/V0-1);
  [tnew,st,info]=gcmma_solver(t,tmin_e,tmax_e,f0s,df0s,gv,dgv,st,it,gopt);
  hist.innerit(it)=info.innerit;
  chg=max(abs(tnew(:)-t));
  if J<Jbest, Jbest=J; tbest=t; end
  hist.J(end+1)=J; hist.Jbest(end+1)=Jbest;
  hist.bandmax(end+1)=bi.maxJk; hist.bandavg(end+1)=mean(bi.Jk); hist.fmax(end+1)=bi.argmax_f;
  hist.vol(end+1)=sum(t)/sum(tmax_e); hist.chg(end+1)=chg; hist.times(end+1)=toc(t0);
  if mod(it,5)==0 || it==1
    fprintf('  [band-%-4s %-11s] it %3d: Jb=%.5e maxJk=%.5e avg=%.5e f*[%.4f] chg=%.2e inner=%d\n',...
      opt.bandmode,metric,it,J,bi.maxJk,mean(bi.Jk),bi.argmax_f,chg,info.innerit);
  end
  t=tnew(:);
  if opt.stop && numel(hist.J)>=2*opt.window+1
    W=opt.window; fo=hist.f0;
    fstat=abs(mean(fo(end-W+1:end))-mean(fo(end-2*W+1:end-W)))/max(abs(mean(fo(end-2*W+1:end-W))),eps);
    xstat=chg/max(opt.tmax-opt.tmin,eps);
    if xstat<opt.xchgtol && fstat<opt.ftol && gv<=opt.feastol, nholdcnt=nholdcnt+1; else, nholdcnt=0; end
    if nholdcnt>=opt.nhold
      stopReason='converged'; sxa=xstat; sfa=fstat;
      fprintf('  [band-%s] STOP at it %d: converged (xstat %.2e, fstat %.2e)\n',opt.bandmode,it,xstat,fstat);
      break;
    end
  end
end
hist.stopReason=stopReason; hist.walltime=toc(ts); hist.numeval=sum(hist.innerit);
hist.metric=metric; hist.p0=p0; hist.tmin=opt.tmin; hist.tmax=opt.tmax;
hist.bandmode=opt.bandmode; hist.rho=opt.rho; hist.band=band;
[~,ib]=min(hist.J); hist.iters_to_best=ib;
tl10=hist.J(max(1,end-9):end); hist.tailspread10=max(tl10)/max(min(tl10),eps);
hist.tailspread=NaN; if numel(hist.J)>=50, tl=hist.J(end-49:end); hist.tailspread=max(tl)/max(min(tl),eps); end
meta=struct('metric',metric,'bandmode',opt.bandmode,'rho',opt.rho,...
  'niter',numel(hist.J),'numeval',sum(hist.innerit),...
  'walltime',hist.walltime,'final_vol',sum(tbest)/sum(tmax_e),...
  'tailspread10',hist.tailspread10,'tailspread',hist.tailspread,...
  'stopReason',stopReason,'Jbest',Jbest,...
  'bandmax_best',min(hist.bandmax),'bandavg_best',min(hist.bandavg));
end

function f0q = obj_scaled(tq,cts,metric,band,opt,p0)
Jq = harmonic_band(band,tq,@(tt,k) articleB_band_eval(cts,k,tt,metric),opt.bandmode,opt.rho);
f0q=100*(opt.Abar+opt.dbfac*log10(max(Jq,eps)))/p0;
end
