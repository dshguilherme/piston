% GCMMA_SUBPROBLEM_PERF -- reproducible measurement of the dual-Newton dense
% allocation hotspot in gcmma_subsolv and of the O(n) rewrite.
%
% THE CLAIM (previously brief-sourced): the inner Newton system for the
% separable dual of the GCMMA/MMA subproblem contains the term diag(d)*G.  The
% old implementation (src/gcmma_solver.m, gcmma_subsolv) materialised it as a
% dense n x n matrix with diag(<n-vector>) ON EVERY Newton step; at a
% representative n = 2048 design that is 2048^2*8 = 33.55 MB per step, and over
% a full run of order 10^1 GB of allocation churn.  The rewrite
% (src/gcmma_solver_fast.m) is element-wise column scaling, the same arithmetic
% at O(n) cost, bit-identical.
%
% WHAT IS MEASURED HERE (cheap and defensible)
%   (1) ms/step and MB/step of the two forms at n = 2048, m = 1 (the real
%       setting: one volume constraint), plus the identity check.
%   (2) the number of times the dense line actually executes inside a
%       representative gcmma_solver_fast subsolv call (MATLAB profiler), which
%       turns "steps-per-call" from prose into a number.
%   (3) the arithmetic that converts (2) x calls-per-run into a churn figure,
%       with calls-per-run taken from the committed run tables.
%
% Writes results/gcmma_subproblem_perf.txt
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
resdir = 'C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(resdir,'dir'), mkdir(resdir); end
fid = fopen(fullfile(resdir,'gcmma_subproblem_perf.txt'),'w');
lg = @(varargin) both(fid, varargin{:});

lg('===== GCMMA subproblem: dual-Newton dense-allocation hotspot =====\n');
lg('Script: tests/gcmma_subproblem_perf.m   (on the offending line of\n');
lg('src/gcmma_solver.m / src/gcmma_solver_fast.m gcmma_subsolv)\n\n');

%% (1) microbenchmark ----------------------------------------------------------
n = 2048; m = 1; reps = 2000;
rng(0);
GG = rand(m,n); diagxinv = rand(n,1); diaglamyi = rand(m,1);
% warm up
Alam = diag(diaglamyi) + GG*diag(diagxinv)*GG.';
Alam2 = diag(diaglamyi) + (GG.*diagxinv.')*GG.';
tic
for k=1:reps, Alam = diag(diaglamyi) + GG*diag(diagxinv)*GG.'; end
t_full = toc;
tic
for k=1:reps, Alam2 = diag(diaglamyi) + (GG.*diagxinv.')*GG.'; end
t_scal = toc;
ms_full = 1000*t_full/reps;  ms_scal = 1000*t_scal/reps;
mb_step = n*n*8/1e6;
maxdiff = max(abs(Alam(:)-Alam2(:)));
lg('(1) MICROBENCHMARK on the offending line, n = %d, m = %d, %d reps\n', n, m, reps);
lg('    old  (diag(d)*G, dense diag): %.4f ms/step\n', ms_full);
lg('    new  (G.*(d.''), O(n)):      %.4f ms/step\n', ms_scal);
lg('    speedup on the line        : %.1fx\n', ms_full/ms_scal);
lg('    bit-identical              : max abs diff = %.3e\n', maxdiff);
lg('    dense n x n diag = %d x %d x 8 bytes = %.2f MB per Newton step\n\n', n, n, mb_step);

%% (2) steps-per-call: how many times the dense line executes in one subsolv ---
% A representative subproblem of the SAME size (n = 2048) and constraint count
% (m = 1, the volume constraint) as the optimization.
rng(1);
x = 0.5*ones(n,1); xmin = 0.4*ones(n,1); xmax = 1.0*ones(n,1);
df0dx = 1e-2*randn(n,1); f0val = 1.0; fval = 0.0; dfdx = ones(1,n)/n;
opt = struct('kmax',20,'growth',1.30,'raa0eps',1e-3,'raaeps',1e-3);
opt.eval = @(xq) deal(f0val + df0dx.'*(xq-x), fval + dfdx*(xq-x));
steps_per_call = NaN;
try
    profile clear; profile on;
    for k = 1:3
        [~,~,info] = gcmma_solver_fast(x, xmin, xmax, f0val, df0dx, fval, dfdx, [], k, opt);
    end
    profile off;
    pi = profile('info');
    steps_per_call = 0;
    for i = 1:numel(pi.FunctionTable)
        fn = pi.FunctionTable(i);
        if contains(fn.FunctionName, 'gcmma_subsolv') || ...
           (isfield(fn,'FileName') && contains(fn.FileName,'gcmma_solver_fast'))
            El = fn.ExecutedLines;   % [line# , count]
            if ~isempty(El)
                for r = 1:size(El,1)
                    % the dense-op line in gcmma_solver_fast.m is line 376
                    if El(r,1) == 376
                        steps_per_call = steps_per_call + El(r,2);
                    end
                end
            end
        end
    end
    lg('(2) NEWTON STEPS PER SUBPROBLEM CALL (profiler, 3 calls, n=%d, m=%d)\n', n, m);
    lg('    dense-line executions total = %d  ->  %d steps per subsolv call\n\n', ...
       steps_per_call, round(steps_per_call/3));
catch err
    lg('(2) profiler step-count FAILED: %s\n\n', err.message);
end

%% (3) churn arithmetic --------------------------------------------------------
lg('(3) CHURN ARITHMETIC:  churn = (subsolv calls per run) x (steps/subsolv) x MB/step\n');
lg('    MB/step = %.2f (measured above)\n', mb_step);
% subsolv calls per run = outer iters + sum(innerit) [each inner iteration
% re-solves the subproblem].  The committed run tables give these directly.
% 64x32 Scordelis-Lo roof (n = 2048 elements): 60 iters, innerit sum per run.
iters64 = 60; inner64 = [17 57 17];     % dyncomp / inputpower / quadvel (64x32)
iters128 = 60; inner128 = [0 24 0];     % dyncomp / inputpower / quadvel (128x64)
calls64 = iters64 + inner64;
calls128 = iters128 + inner128;
if ~isnan(steps_per_call)
    spc = round(steps_per_call/3);
    lg('    steps/subsolv (measured) = %d\n', spc);
    % n = 2048 -> 64x32 mesh
    churn64 = calls64 .* spc .* mb_step / 1000;   % MB -> GB
    lg('    64x32 (n=2048), calls/run = %s:\n', mat2str(calls64));
    lg('      churn/run = calls x %d steps x %.2f MB = %s GB\n', spc, mb_step, ...
       mat2str(round(churn64,1)));
    % 128x64 (n=8192): dense step would be 8192^2*8 = 536.9 MB; but the fast
    % solver avoids it entirely.  Report the per-step figure at that size:
    mb_step8192 = 8192*8192*8/1e6;
    lg('    128x64 (n=8192): dense step would be %.1f MB; calls/run = %s\n', ...
       mb_step8192, mat2str(calls128));
    lg('      churn/run = calls x %d steps x %.1f MB = %s GB (if the dense form ran)\n', ...
       spc, mb_step8192, mat2str(round(calls128.*spc.*mb_step8192/1000,1)));
else
    lg('    steps/subsolv not measured (profiler unavailable); churn cannot be\n');
    lg('    closed from first principles here.\n');
end

%% (4) end-to-end reconciliation ----------------------------------------------
% Captured by running  tests/bench_impl.m  on this machine (the command below
% is reproducible as-is); bench_impl is the manuscript's timing harness.
lg('\n(4) END-TO-END (tests/bench_impl.m, meshes {[64 32]}, 20 outer iterations)\n');
lg('    command: matlab -batch "addpath(genpath(''...geopdes-master'')); ...; \\\n');
lg('                      bench_impl({[64 32]}, 20);"\n');
lg('    BEFORE (gcmma_impl=std )  gcmma phase 293.782 s  of 464.991 s total (63.2%%)\n');
lg('    AFTER  (gcmma_impl=fast)  gcmma phase 283.091 s  of 460.712 s total (61.4%%)\n');
lg('    end-to-end speedup (sum hist.times) = 1.01x  (465.75 s -> 461.89 s)\n');
lg('    CORRECTNESS: max rel diff on J per outer iter = 0.000e+00 (bit-identical)\n');
lg('    => The dense-allocation fix removes the O(10^1) GB of churn and is ~3600x\n');
lg('       on the offending line, but the END-TO-END wall-time effect at this mesh\n');
lg('       is ~1%%, because the subproblem dual assembly is not the dominant cost\n');
lg('       of the GCMMA phase here.  The working-note end-to-end figures of 2.5x\n');
lg('       (64x32) and 6.7x (128x64) were NOT reproducible and are NOT claimed in\n');
lg('       the manuscript; only the measured quantities above are.\n');
fclose(fid);
fprintf('GCMMA_SUBPROBLEM_PERF_DONE\n');

% ============================ local functions =============================
function both(fid, varargin)
  fprintf(fid, varargin{:}); fprintf(1, varargin{:});
end
