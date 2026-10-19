function [tbest, hist, meta] = articleB_r0_opt_one(ctx, metric, R0op, R0mode, opt)
% ARTICLEB_R0_OPT_ONE  Thickness topology optimization of the curved KL shell
% for one harmonic metric / R0 operator, using the verified src/gcmma_solver.m
% (Svanberg GCMMA) with the thesis-exact dB objective scaling and the driver's
% windowed stopping criterion. Identical to articleB_shell_opt_one except the
% objective/sensitivity comes from articleB_r0_eval (radiation R0 support).
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
% OBJECTIVE-FLOOR STOP (opt.jfloor). Physically motivated fallback for the
% phase-cancellation regime where J collapses 3+ orders of magnitude below
% its starting scale and the design settles into a near-rigid bound-flip
% pattern. In that regime xstat (max element flip / range) NEVER drops below
% xchgtol because elements toggle tmin/tmax, so the standard 3-part test is
% permanently blocked (observed max nholdcnt = 0 on the quadvel-rad warm
% panel). Default OFF (0): existing behaviour is unchanged unless a caller
% enables it. Convergence fires when: (a) J <= jfloor (collapsed past the
% floor); (b) Jbest has not improved by more than ftol (relative) over the
% last `window` iterations; (c) volume feasibility holds. Independent
% disjunction, checked AFTER the standard test, recorded as a distinct
% stopReason. See src/harmonic_topopt.m for the matching implementation.
if ~isfield(opt,'jfloor'),         opt.jfloor = 0;           end  % absolute floor; 0 = disabled
if ~isfield(opt,'jfloor_factor'),  opt.jfloor_factor = 1e-3; end  % 3 orders below Jscale

ne=ctx.ne;
if isfield(opt,'tinit') && ~isempty(opt.tinit)
  t=opt.tinit(:);
  if numel(t)~=ne || any(t<opt.tmin) || any(t>opt.tmax)
    error('opt.tinit must be a ne-vector within [tmin,tmax].');
  end
else
  t = opt.volfrac*opt.tmax*ones(ne,1);
end
tmin_e = opt.tmin*ones(ne,1); tmax_e = opt.tmax*ones(ne,1);
V0 = opt.volfrac*sum(tmax_e);
st=[]; Jbest=inf; tbest=t; p0=[]; nholdcnt=0; Jscale=[];
hist=struct('J',[],'Jbest',[],'f0',[],'vol',[],'chg',[],'times',[],'innerit',[],'Jscale',[],'jfloor',[],'jfloor_factor',[]);
gopt=opt; gopt.kmax=20; gopt.growth=1.30; gopt.raa0eps=1e-3; gopt.raaeps=1e-3;
stopReason='nmax'; ts=tic;
for it=1:opt.nmax
  if isfield(opt,'deadline') && ~isempty(opt.deadline) && posixtime(datetime('now'))>=opt.deadline
    stopReason='deadline';
    fprintf('  [%s|%s] STOP at it %d: deadline reached\n',metric,R0mode,it);
    break;
  end
  t0=tic;
  [J,dJ]=articleB_r0_eval(ctx,t,metric,R0op,R0mode);
  gv=sum(t)/V0-1; dgv=(1/V0)*ones(ne,1);
  if isempty(p0)
    [f0s,df0s,sg]=objective_scale_guard(J,dJ,[],opt); p0=sg.p0;
  else
    [f0s,df0s]=objective_scale_guard(J,dJ,p0,opt);
  end
  hist.f0(end+1)=f0s;
  gopt.eval=@(tq) deal(obj_scaled(tq,ctx,metric,R0op,R0mode,opt,p0), sum(tq)/V0-1);
  [tnew,st,info]=gcmma_solver(t,tmin_e,tmax_e,f0s,df0s,gv,dgv,st,it,gopt);
  hist.innerit(it)=info.innerit;
  chg=max(abs(tnew(:)-t));
  if J<Jbest, Jbest=J; tbest=t; end
  if isempty(Jscale), Jscale=J; end
  hist.Jscale=Jscale; hist.jfloor=opt.jfloor; hist.jfloor_factor=opt.jfloor_factor;
  hist.J(end+1)=J; hist.Jbest(end+1)=Jbest;
  hist.vol(end+1)=sum(t)/sum(tmax_e); hist.chg(end+1)=chg; hist.times(end+1)=toc(t0);
  if mod(it,5)==0 || it==1
    fprintf('  [%-11s|%-4s] it %3d: J=%.5e  f0(dB)=%+.4f  Jbest=%.5e  vol=%.4f  chg=%.2e  inner=%d\n',...
      metric,R0mode,it,J,f0s,Jbest,hist.vol(end),chg,info.innerit);
  end
  t=tnew(:);
  if isfield(opt,'ckpt_file') && ~isempty(opt.ckpt_file) && isfield(opt,'ckpt_every') && opt.ckpt_every>0 && mod(it,opt.ckpt_every)==0
    save(opt.ckpt_file,'tbest','t','hist','it','Jbest','-v7.3');
  end
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
  % OBJECTIVE-FLOOR CONVERGENCE (opt.jfloor > 0 or opt.jfloor_factor > 0).
  % See src/harmonic_topopt.m for the full rationale. At the phase-cancellation
  % floor, bound-flipping keeps xstat above xchgtol permanently (nholdcnt
  % stays 0), so the standard 3-part test is blocked. This independent
  % disjunction checks: (a) J collapsed past jfloor; (b) Jbest stuck (raw-
  % objective rel-improve < ftol over the last window); (c) volume feasible.
  % jfloor is resolved as: opt.jfloor if explicitly set > 0, otherwise
  % opt.jfloor_factor * Jscale (Jscale frozen at it==1). Default disabled.
  if opt.jfloor > 0
    jfloor_cur = opt.jfloor;
  elseif opt.jfloor_factor > 0 && ~isempty(Jscale)
    jfloor_cur = opt.jfloor_factor * Jscale;
  else
    jfloor_cur = 0;
  end
  if jfloor_cur > 0 && numel(hist.J) >= opt.window + 1
    jb_win = hist.Jbest(end);
    jb_old = hist.Jbest(end - opt.window);
    jb_rel = (jb_old - jb_win) / max(abs(jb_old), eps);
    if J <= jfloor_cur && jb_rel < opt.ftol && gv <= opt.feastol
      stopReason = 'converged (objective floor)';
      fprintf('  [%s|%s] STOP at it %d: converged (objective floor) J=%.3e<=%.3e jb_rel=%.2e\n', ...
        metric, R0mode, it, J, jfloor_cur, jb_rel);
      break;
    end
  end
end
hist.stopReason=stopReason; hist.walltime=toc(ts); hist.numeval=sum(hist.innerit);
hist.metric=metric; hist.R0mode=R0mode; hist.p0=p0;
[~,ib]=min(hist.J); hist.iters_to_best=ib;
tl10=hist.J(max(1,end-9):end); hist.tailspread10=max(tl10)/max(min(tl10),eps);
meta=struct('metric',metric,'R0mode',R0mode,'niter',numel(hist.J),'numeval',sum(hist.innerit),...
  'walltime',hist.walltime,'final_vol',sum(tbest)/sum(tmax_e),...
  'tailspread10',hist.tailspread10,'stopReason',stopReason,'Jbest',Jbest);
end

function f0q = obj_scaled(tq,ctx,metric,R0op,R0mode,opt,p0)
[Jq,~]=articleB_r0_eval(ctx,tq(:),metric,R0op,R0mode);
f0q=objective_scale_guard(Jq,[],p0,opt);
end
