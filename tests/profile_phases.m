function R = profile_phases(nsub, niter, w, alpha, beta)
% PROFILE_PHASES Replicate the harmonic_topopt outer loop, independently timing
% each phase with tic/toc, to find where the time actually goes.
%   R = profile_phases(nsub, niter, w, alpha, beta)
% nsub is [n1 n2] (default [64 32]); niter is #outer iterations (default 60).

if nargin < 1 || isempty(nsub),  nsub  = [64 32]; end
if nargin < 2 || isempty(niter), niter = 60;      end
if nargin < 3 || isempty(w),     w     = 2*pi*12.09; end
if nargin < 4 || isempty(alpha), alpha = 25;      end
if nargin < 5 || isempty(beta),  beta  = 0;       end

E=210e9; RHO=7850; NU=0.3;
pd = cantilever_beam(1.0,1.0);
md = struct('degree',[3 3],'regularity',[2 2],'nsub',nsub,'nquad',[4 4]);
[gm,mh,sp] = buildSpaces(pd,md);
H = harmonic_kd_setup(mh,sp,gm,E,NU,RHO,'cpu','double','plate');
F = zeros(H.N,1); [~,ic] = min(abs((1:H.N).'-round(H.N/2))); F(ic)=1;

opt = struct('tmin',0.4,'tmax',1.0,'volfrac',0.5,'nmax',niter,'stop',false,...
             'objscale','log','optimizer','gcmma');

nel1 = H.nel_dir(1); nel2 = H.nel_dir(2);
t = opt.volfrac * opt.tmax * ones(nel1, nel2);
tmin_e = opt.tmin*ones(nel1,nel2); tmax_e = opt.tmax*ones(nel1,nel2);
V0 = opt.volfrac * sum(tmax_e(:));
Jbest = inf;

P.form=0; P.solve=0; P.sens=0; P.gcmma=0; P.kkt=0; P.other=0;
Tin = struct('kmax',20,'growth',1.30,'raa0eps',1e-3,'raaeps',1e-3);
Tin = opt; Tin.kmax=20; Tin.growth=1.30; Tin.raa0eps=1e-3; Tin.raaeps=1e-3;
st = []; innerits = zeros(1,niter); modes = cell(1,niter); Js = zeros(1,niter);

% one warm-up phase so JIT/allocation noise does not land on phase 1
for it = 1:niter
    tt = tic;
    KD = harmonic_kd_form(H, t, w, alpha, beta, 'cpu_sparse'); P.form = P.form + toc(tt);
    tt = tic; U = KD \ F; P.solve = P.solve + toc(tt);
    g = F.'*U; J = abs(g);
    tt = tic; dJ = harmonic_kd_sens(H, t, F, U, w, alpha, beta, 'dyncomp'); P.sens = P.sens + toc(tt);
    gv = sum(t(:))/V0 - 1; dgv = (1/V0)*ones(nel1,nel2);
    Abar=100; fac=10; Jdb = Abar + fac*log10(max(J,eps));
    if ~isfield(opt,'p0') || isempty(opt.p0), opt.p0 = Jdb; end
    f0s = 100*Jdb/opt.p0; df0s = 100*(fac/log(10))*dJ/max(J,eps)/opt.p0;
    tt = tic;
    dL = df0s(:) + max(0, -(df0s(:).'*dgv(:))/(dgv(:).'*dgv(:)))*dgv(:);
    tflat=t(:); atlo=tflat<=tmin_e(:)+1e-12; athi=tflat>=tmax_e(:)-1e-12;
    keep=~((atlo & dL>0)|(athi & dL<0));
    kkt=max(abs(dL(keep)))/max(max(abs(df0s(:))),eps);
    P.kkt = P.kkt + toc(tt);
    tt = tic;
    [tnew, st, info] = gcmma_solver(t, tmin_e, tmax_e, f0s, df0s, gv, dgv, st, it, Tin);
    P.gcmma = P.gcmma + toc(tt);
    P.other = P.other + (toc(tt) - P.gcmma - P.kkt);  % not meaningful; keep 0 below
    innerits(it) = info.innerit; modes{it} = info.mode;
    if J < Jbest, Jbest = J; end
    Js(it) = J;
    t = reshape(tnew, nel1, nel2);
end

tot = P.form+P.solve+P.sens+P.gcmma+P.kkt;
ph = {'formation','solve','sensitivity','gcmma','kkt+metric'};
vals = [P.form P.solve P.sens P.gcmma P.kkt];
R = struct('nsub',nsub,'niter',niter,'tot',tot,'vals',vals,'ph',{ph}, ...
           'innerits',innerits,'mode',modes{end},'Js',Js,'N',H.N);

fprintf('\n=== PHASE PROFILE  mesh %dx%d  N=%d  iters=%d ===\n', nel1, nel2, H.N, niter);
fprintf('%-14s %12s %10s\n','phase','mean [s]','% iter');
for k=1:numel(ph)
    fprintf('%-14s %12.5f %9.1f%%\n', ph{k}, vals(k)/niter, 100*vals(k)/tot);
end
fprintf('%-14s %12.5f %9.1f%%\n','TOTAL',tot/niter,100);
fprintf('innerit: mean %.2f  max %d   mode: %s\n', mean(innerits), max(innerits), modes{end});
end
