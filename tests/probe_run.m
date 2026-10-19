% PROBE_RUN Replicate the driver loop twice: once with the real gcmma_solver
% (baseline J trajectory) and once with the counter-instrumented probe, then
% print the subsolv iteration statistics.
E=210e9; RHO=7850; NU=0.3;
pd = cantilever_beam(1.0,1.0);
md = struct('degree',[3 3],'regularity',[2 2],'nsub',[64 32],'nquad',[4 4]);
[gm,mh,sp] = buildSpaces(pd,md);
H = harmonic_kd_setup(mh,sp,gm,E,NU,RHO,'cpu','double','plate');
F = zeros(H.N,1); [~,ic] = min(abs((1:H.N).'-round(H.N/2))); F(ic)=1;
opt = struct('tmin',0.4,'tmax',1.0,'volfrac',0.5,'nmax',50,'stop',false,...
             'objscale','log','optimizer','gcmma');

Jbase = run_loop(H,F,opt,@gcmma_solver);

global GCMMA_PROBE %#ok<GVMIS>
GCMMA_PROBE = [];
Jprobe = run_loop(H,F,opt,@gcmma_solver_probe);

fprintf('\n=== SUBSOLV WORK (probe, 50 outer iters) ===\n');
disp(GCMMA_PROBE);
fprintf('newton iters per subsolv call : %.1f\n', GCMMA_PROBE.newton_iters/GCMMA_PROBE.subsolv_calls);
fprintf('barrier levels per subsolv call: %.2f\n', GCMMA_PROBE.bar_levels/GCMMA_PROBE.subsolv_calls);
fprintf('linesearch iters per newton   : %.2f\n', GCMMA_PROBE.ls_iters/GCMMA_PROBE.newton_iters);
fprintf('trajectory match: max rel diff J = %.3e\n', max(abs(Jbase-Jprobe)./max(abs(Jbase),eps)));

function J = run_loop(H,F,opt,gfun)
nel1=H.nel_dir(1); nel2=H.nel_dir(2);
t = opt.volfrac*opt.tmax*ones(nel1,nel2);
tmin_e=opt.tmin*ones(nel1,nel2); tmax_e=opt.tmax*ones(nel1,nel2);
V0 = opt.volfrac*sum(tmax_e(:)); st=[]; w=2*pi*12.09; alpha=25; beta=0;
J = zeros(1,opt.nmax); gopt=opt; gopt.kmax=20; gopt.growth=1.30;
gopt.raa0eps=1e-3; gopt.raaeps=1e-3; p0=[];
for it=1:opt.nmax
    KD = harmonic_kd_form(H,t,w,alpha,beta,'cpu_sparse');
    U = KD\F; g=F.'*U; Jv=abs(g);
    dJ = harmonic_kd_sens(H,t,F,U,w,alpha,beta,'dyncomp');
    gv = sum(t(:))/V0 - 1; dgv=(1/V0)*ones(nel1,nel2);
    Jdb = 100 + 10*log10(max(Jv,eps)); if isempty(p0), p0=Jdb; end
    f0s=100*Jdb/p0; df0s=100*(10/log(10))*dJ/max(Jv,eps)/p0;
    [tnew,st] = gfun(t,tmin_e,tmax_e,f0s,df0s,gv,dgv,st,it,gopt);
    J(it)=Jv; t=reshape(tnew,nel1,nel2);
end
end
