function bench_impl(meshes, niter)
% BENCH_IMPL Before/after timing for gcmma_impl='std' vs 'fast', with the
% objective-trajectory correctness check.  bench_impl(meshes, niter)
if nargin < 1 || isempty(meshes), meshes = {[64 32],[128 64]}; end
if nargin < 2 || isempty(niter),  niter  = 50; end

for mi = 1:numel(meshes)
    nsub = meshes{mi};
    [H,F] = setup_case(nsub);
    fprintf('\n==================== mesh %dx%d  (N=%d) ====================\n', nsub(1), nsub(2), H.N);
    [h_s, J_s] = run_one(H,F,nsub,niter,'std');
    [h_f, J_f] = run_one(H,F,nsub,niter,'fast');
    print_split('BEFORE  (gcmma_impl=std )', h_s);
    print_split('AFTER   (gcmma_impl=fast)', h_f);
    fprintf('speedup (sum hist.times): %.2fx  (%.3f s -> %.3f s)\n', ...
        sum(h_s.times)/sum(h_f.times), sum(h_s.times), sum(h_f.times));
    d = max(abs(J_s(:)-J_f(:))./max(abs(J_s(:)),eps));
    fprintf('CORRECTNESS: max rel diff on J per outer iter = %.3e\n', d);
end
end

function [H,F] = setup_case(nsub)
E=210e9; RHO=7850; NU=0.3;
pd = cantilever_beam(1.0,1.0);
md = struct('degree',[3 3],'regularity',[2 2],'nsub',nsub,'nquad',[4 4]);
[gm,mh,sp] = buildSpaces(pd,md);
H = harmonic_kd_setup(mh,sp,gm,E,NU,RHO,'cpu','double','plate');
F = zeros(H.N,1); [~,ic] = min(abs((1:H.N).'-round(H.N/2))); F(ic)=1;
end

function [h,J] = run_one(H,F,nsub,niter,impl)
opt = struct('tmin',0.4,'tmax',1.0,'volfrac',0.5,'nmax',niter,'stop',false,...
             'objscale','log','optimizer','gcmma','gcmma_impl',impl);
[~,h] = harmonic_topopt(H,F,2*pi*12.09,25,0,opt);
J = h.J;
end

function print_split(label, h)
m = h.tphase_mean; p = h.tphase_pct;
fprintf('\n%s :  [phase split over %d outer iters]\n', label, numel(h.times));
fprintf('  %-12s %12s %8s   %12s\n','phase','mean [s]','% ','total [s]');
fprintf('  %-12s %12.5f %7.1f%%   %12.4f\n','formation',m.form, p.form, h.tphase.form);
fprintf('  %-12s %12.5f %7.1f%%   %12.4f\n','solve',    m.solve,p.solve,h.tphase.solve);
fprintf('  %-12s %12.5f %7.1f%%   %12.4f\n','sens',     m.sens, p.sens, h.tphase.sens);
fprintf('  %-12s %12.5f %7.1f%%   %12.4f\n','kkt+metric',m.kkt,p.kkt,  h.tphase.kkt);
fprintf('  %-12s %12.5f %7.1f%%   %12.4f\n','gcmma',    m.gcmma,p.gcmma,h.tphase.gcmma);
fprintf('  %-12s %12.5f %7.1f%%   %12.4f\n','TOTAL',    m.total,100,    h.tphase_tot);
end