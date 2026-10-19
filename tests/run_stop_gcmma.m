function run_stop_gcmma(configfile, outfile)
% RUN_STOP_GCMMA Test the stopping criterion (opt.stop=true) across the best
% GCMMA configurations, now that the inner loop is live. Reports whether it
% fires, at which outer iteration, and the resulting tail spread. Longer budget
% (nmax=150) than the sweep so the criterion has room to fire.
%
% Usage: run_stop_gcmma(cfg.csv, out.csv)

addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
if nargin < 1 || isempty(configfile)
    configfile = 'C:/Users/dshgu/AppData/Local/hermes/cache/scratch/stop_cfg.csv';
end
if nargin < 2 || isempty(outfile)
    outfile = 'C:/Users/dshgu/AppData/Local/hermes/cache/scratch/stop_out.csv';
end
T = readtable(configfile);

E=210e9; RHO=7850; NU=0.3;
pd = cantilever_beam(1.0,1.0);
md = struct('degree',[3 3],'regularity',[2 2],'nsub',[64 32],'nquad',[4 4]);
[gm,mh,sp] = buildSpaces(pd,md);
H = harmonic_kd_setup(mh,sp,gm,E,NU,RHO,'cpu','double','plate');
F = zeros(H.N,1); [~,ic] = min(abs((1:H.N).'-round(H.N/2))); F(ic)=1;
fprintf('setup done, N=%d\n', H.N);

NMAX = 150;
fprintf('%-4s %5s %6s %10s %7s %10s %12s %12s %8s\n', ...
    'idx','kmax','growth','raa0eps','iters','stopIt','tailspread','bestJ','reason');
for r = 1:height(T)
    km = T.kmax(r); gr = T.growth(r); re = T.raa0eps(r);
    opt = struct('tmin',0.4,'tmax',1.0,'volfrac',0.5,'nmax',NMAX,'stop',true, ...
                 'objscale','log','optimizer','gcmma', ...
                 'kmax',km,'growth',gr,'raa0eps',re,'raaeps',re,'gcmma_impl','fast', ...
                 'ftol',1e-4,'xchgtol',1e-4,'nhold',3,'window',10);
    [~,hist] = harmonic_topopt(H,F,2*pi*12.09,25,0,opt);
    nit = numel(hist.J);
    fprintf('%-4d %5d %6.2f %10.1e %7d %10d %12.4e %12.4e %8s\n', ...
        T.idx(r), km, gr, re, nit, nit, hist.tailspread, min(hist.J), hist.stopReason);
    fid = fopen(outfile,'a');
    if fid > 0
        fprintf(fid, '%d,%d,%.2f,%.1e,%d,%.6e,%.6e,%s\n', ...
            T.idx(r), km, gr, re, nit, hist.tailspread, min(hist.J), hist.stopReason);
        fclose(fid);
    end
end
fprintf('STOP_TEST_DONE\n');
end
