function diag_stop(km, gr, re, nmax)
% DIAG_STOP Rerun one GCMMA config with stop=false and recompute the three
% stopping-criterion components (fstat, xstat, feasibility) exactly as the
% driver does, iteration by iteration, so we can see which condition blocks
% convergence. Usage: diag_stop(8,1.3,1e-3,60)
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
if nargin < 1, km=8; end
if nargin < 2, gr=1.3; end
if nargin < 3, re=1e-3; end
if nargin < 4, nmax=60; end
E=210e9; RHO=7850; NU=0.3;
pd = cantilever_beam(1.0,1.0);
md = struct('degree',[3 3],'regularity',[2 2],'nsub',[64 32],'nquad',[4 4]);
[gm,mh,sp] = buildSpaces(pd,md);
H = harmonic_kd_setup(mh,sp,gm,E,NU,RHO,'cpu','double','plate');
F = zeros(H.N,1); [~,ic] = min(abs((1:H.N).'-round(H.N/2))); F(ic)=1;
opt = struct('tmin',0.4,'tmax',1.0,'volfrac',0.5,'nmax',nmax,'stop',false, ...
             'objscale','log','optimizer','gcmma','kmax',km,'growth',gr, ...
             'raa0eps',re,'raaeps',re,'gcmma_impl','fast');
[~,h] = harmonic_topopt(H,F,2*pi*12.09,25,0,opt);
W=10; ftol=1e-4; xchgtol=1e-4; feastol=1e-6;
fo=h.f0;
fprintf('\nkmax=%d growth=%.2f raa0eps=%.0e  iters=%d\n',km,gr,re,numel(h.J));
fprintf('%4s %12s %12s %12s %8s %8s %8s\n','it','fstat','xstat','gv','okF','okX','okFeas');
for k=2*W+1:numel(fo)
    fstat = abs(mean(fo(k-W+1:k)) - mean(fo(k-2*W+1:k-W))) / max(abs(mean(fo(k-2*W+1:k-W))),eps);
    xstat = h.chg(k)/max(1.0-0.4,eps);
    gv    = h.vol(k)/0.5 - 1;
    fprintf('%4d %12.3e %12.3e %12.3e %8d %8d %8d\n', k, fstat, xstat, gv, ...
        fstat<ftol, xstat<xchgtol, gv<=feastol);
end
end
