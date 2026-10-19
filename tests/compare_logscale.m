function compare_logscale(km, gr, raa0eps, outfile)
% COMPARE_LOGSCALE Fair MMA vs GCMMA at the SAME objective scaling (objscale=
% 'log') and the SAME budget (nmax=60, stop=false). Reports tail spread /
% best J / wall time for both. GCMMA uses the supplied (kmax,growth,raa0eps).
%
% Usage: compare_logscale(8, 1.3, 1e-5, out.csv)

addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
if nargin < 4 || isempty(outfile)
    outfile = 'C:/Users/dshgu/AppData/Local/hermes/cache/scratch/mma_vs_gcmma.csv';
end

E=210e9; RHO=7850; NU=0.3;
pd = cantilever_beam(1.0,1.0);
md = struct('degree',[3 3],'regularity',[2 2],'nsub',[64 32],'nquad',[4 4]);
[gm,mh,sp] = buildSpaces(pd,md);
H = harmonic_kd_setup(mh,sp,gm,E,NU,RHO,'cpu','double','plate');
F = zeros(H.N,1); [~,ic] = min(abs((1:H.N).'-round(H.N/2))); F(ic)=1;
fprintf('setup done, N=%d\n', H.N);
NMAX = 60;

fprintf('%-6s %5s %6s %10s %7s %9s %12s %12s %12s %12s %8s\n', ...
    'who','kmax','growth','raa0eps','iters','evals','wall_s','tailspread','bestJ','lastJ','reason');
for who = {'gcmma','mma'}
    opt = struct('tmin',0.4,'tmax',1.0,'volfrac',0.5,'nmax',NMAX,'stop',false, ...
                 'objscale','log','optimizer',who{1});
    if strcmp(who{1},'gcmma')
        opt.kmax = km; opt.growth = gr; opt.raa0eps = raa0eps; opt.raaeps = raa0eps;
        opt.gcmma_impl = 'fast';
    end
    t0 = tic;
    [~,hist] = harmonic_topopt(H,F,2*pi*12.09,25,0,opt);
    wall = toc(t0);
    nit = numel(hist.J);
    if isempty(hist.innerit), evals = nit; epo = 1; else
        evals = nit + sum(hist.innerit); epo = evals/max(nit,1);
    end
    fprintf('%-6s %5d %6.2f %10.1e %7d %9d %12.2f %12.4e %12.4e %12.4e %8s\n', ...
        who{1}, km, gr, raa0eps, nit, evals, wall, hist.tailspread, min(hist.J), hist.J(end), hist.stopReason);
    fid = fopen(outfile,'a');
    if fid > 0
        fprintf(fid, '%s,%d,%.2f,%.1e,%d,%d,%.2f,%.6e,%.6e,%.6e,%s\n', ...
            who{1}, km, gr, raa0eps, nit, evals, wall, hist.tailspread, min(hist.J), hist.J(end), hist.stopReason);
        fclose(fid);
    end
    save(sprintf('C:/Users/dshgu/AppData/Local/hermes/cache/scratch/mv_%s.mat', who{1}), 'hist');
end
fprintf('COMPARE_LOGSCALE_DONE\n');
end
