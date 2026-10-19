function sweep_gcmma_conservatism(rowlo, rowhi, configfile, outfile)
% SWEEP_GCMMA_CONSERVATISM Re-sweep GCMMA conservatism knobs now that the inner
% loop is LIVE (opt.eval is supplied by the driver). Runs config rows
% [rowlo,rowhi] of a CSV (idx,kmax,growth,raa0eps) at 64x32, nmax=60, stop=false
% and APPENDS one result line per config to outfile (CSV). Idempotent per row so
% the sweep can be split across MATLAB sessions.
%
% Usage: sweep_gcmma_conservatism(1, 6, cfg.csv, out.csv)

addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');

if nargin < 1 || isempty(rowlo), rowlo = 1; end
if nargin < 2 || isempty(rowhi), rowhi = rowlo; end
if nargin < 3 || isempty(configfile)
    configfile = 'C:/Users/dshgu/AppData/Local/hermes/cache/scratch/gcmma_cfg.csv';
end
if nargin < 4 || isempty(outfile)
    outfile = 'C:/Users/dshgu/AppData/Local/hermes/cache/scratch/gcmma_sweep_out.csv';
end

T = readtable(configfile);
if rowhi > height(T), rowhi = height(T); end

% ---- one-time mesh setup (shared across all rows) ----
E=210e9; RHO=7850; NU=0.3;
pd = cantilever_beam(1.0,1.0);
md = struct('degree',[3 3],'regularity',[2 2],'nsub',[64 32],'nquad',[4 4]);
[gm,mh,sp] = buildSpaces(pd,md);
H = harmonic_kd_setup(mh,sp,gm,E,NU,RHO,'cpu','double','plate');
F = zeros(H.N,1); [~,ic] = min(abs((1:H.N).'-round(H.N/2))); F(ic)=1;
tsetup = tic; fprintf('setup done (%.2f s), N=%d\n', toc(tsetup), H.N);

NMAX = 60;
fprintf('%-4s %5s %6s %10s %6s %8s %10s %9s %12s %12s %8s\n', ...
    'idx','kmax','growth','raa0eps','iters','evals','eval/out','wall_s','tailspread','bestJ','conserv');
for r = rowlo:rowhi
    km = T.kmax(r); gr = T.growth(r); re = T.raa0eps(r);
    opt = struct('tmin',0.4,'tmax',1.0,'volfrac',0.5,'nmax',NMAX,'stop',false, ...
                 'objscale','log','optimizer','gcmma', ...
                 'kmax',km,'growth',gr,'raa0eps',re,'raaeps',re,'gcmma_impl','fast');
    t0 = tic;
    [~,hist] = harmonic_topopt(H,F,2*pi*12.09,25,0,opt);
    wall = toc(t0);
    nit   = numel(hist.J);
    evals = nit + sum(hist.innerit);
    epo   = evals/max(nit,1);
    bJ    = min(hist.J);
    % last-config conservatism fraction (fraction of outer iters where the
    % accepted point was already conservative with innerit==0 is NOT stored;
    % report mean innerit instead via epo). tailspread from hist.
    fprintf('%-4d %5d %6.2f %10.1e %6d %8d %10.2f %9.2f %12.4e %12.4e %8s\n', ...
        r, km, gr, re, nit, evals, epo, wall, hist.tailspread, bJ, hist.stopReason);
    % append CSV
    fid = fopen(outfile,'a');
    if fid > 0
        fprintf(fid, '%d,%d,%.2f,%.1e,%d,%d,%.3f,%.2f,%.6e,%.6e,%s\n', ...
            r, km, gr, re, nit, evals, epo, wall, hist.tailspread, bJ, hist.stopReason);
        fclose(fid);
    end
end
fprintf('SWEEP_DONE rows %d..%d\n', rowlo, rowhi);
end
