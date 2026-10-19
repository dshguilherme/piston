function [topt, hist, meta] = run_metric_opt(H, F, w, alpha, beta, metric, opt)
% RUN_METRIC_OPT Thickness topology optimization for ONE harmonic objective,
% reproducing the outer loop of src/harmonic_topopt.m (which is hard-wired to
% dynamic compliance) so that the OTHER two metrics ('inputpower', 'quadvel')
% can be optimised with the SAME building blocks and the SAME dB scaling.
%
%   [topt,hist,meta] = run_metric_opt(H,F,w,alpha,beta,metric,opt)
%
% Building blocks (per outer iteration), exactly as the driver:
%   KS = harmonic_kd_form(H,t,w,alpha,beta,'cpu_sparse');  U = KS\F;
%   J  = obj_raw(...);   dJ = harmonic_kd_sens(H,t,F,U,w,alpha,beta,metric,R0);
%   dB scaling {P}_dB = Abar + fac*log10(J), f0 = 100*{P}_dB/p0,
%   df0 = 100*(fac/ln10) * dJ/J / p0   (Olavo thesis eqs. 6.113-6.115),
%   then one gcmma_solver call with opt.eval wired to the SAME scaling.
%
% metric: 'dyncomp' | 'inputpower' | 'quadvel'   (case-insensitive)
% opt.R0: quadvel weighting spec ([] identity, 'mass' -> unit-thickness M0).

if nargin < 7, opt = struct(); end
if ~isfield(opt,'tmin'),     opt.tmin     = 0.4;  end
if ~isfield(opt,'tmax'),     opt.tmax     = 1.0;  end
if ~isfield(opt,'volfrac'),  opt.volfrac  = 0.5;  end
if ~isfield(opt,'nmax'),     opt.nmax     = 100;  end
if ~isfield(opt,'tol'),      opt.tol      = 0;    end
if ~isfield(opt,'stop'),     opt.stop     = true; end
if ~isfield(opt,'ftol'),     opt.ftol     = 1e-4; end
if ~isfield(opt,'xchgtol'),  opt.xchgtol  = 1e-4; end
if ~isfield(opt,'nhold'),    opt.nhold    = 3;    end
if ~isfield(opt,'feastol'),  opt.feastol  = 1e-6; end
if ~isfield(opt,'window'),   opt.window   = 10;   end
if ~isfield(opt,'Abar'),     opt.Abar     = 100;  end
if ~isfield(opt,'dbfac'),    opt.dbfac    = 10;   end
if ~isfield(opt,'R0'),       opt.R0       = [];   end

nel1 = H.nel_dir(1); nel2 = H.nel_dir(2);
t      = opt.volfrac * opt.tmax * ones(nel1, nel2);
tmin_e = opt.tmin * ones(nel1, nel2);
tmax_e = opt.tmax * ones(nel1, nel2);
V0 = opt.volfrac * sum(tmax_e(:));

st = []; Jbest = inf; tbest = t; p0 = []; nholdcnt = 0;
hist = struct('J',[],'Jbest',[],'f0',[],'vol',[],'chg',[],'times',[],'innerit',[]);

% GCMMA tuning identical to the driver (stiff conservatism).
gopt = opt;
gopt.kmax = 20; gopt.growth = 1.30; gopt.raa0eps = 1e-3; gopt.raaeps = 1e-3;

stopReason = 'nmax';
t_start = tic;
for it = 1:opt.nmax
    t0 = tic;
    KS = harmonic_kd_form(H, t, w, alpha, beta, 'cpu_sparse');
    U  = KS \ F;
    J  = obj_raw(U, F, w, metric, opt.R0);
    dJ = harmonic_kd_sens(H, t, F, U, w, alpha, beta, metric, opt.R0);

    gv  = sum(t(:))/V0 - 1;
    dgv = (1/V0) * ones(nel1, nel2);

    Jdb = opt.Abar + opt.dbfac * log10(max(J, eps));
    if isempty(p0), p0 = Jdb; end
    f0s  = 100 * Jdb / p0;
    df0s = 100 * (opt.dbfac/log(10)) * dJ / max(J, eps) / p0;
    hist.f0(end+1) = f0s;

    % Conservative-check callback: SAME scaling as f0s.
    gopt.eval = @(tq) deal( ...
        obj_scaled(tq, H, F, w, alpha, beta, metric, opt, V0, nel1, nel2, p0), ...
        sum(tq(:))/V0 - 1 );

    [tnew, st, info] = gcmma_solver(t, tmin_e, tmax_e, f0s, df0s, gv, dgv, st, it, gopt);
    hist.innerit(it) = info.innerit;

    chg = max(abs(tnew(:) - t(:)));
    if J < Jbest, Jbest = J; tbest = t; end          % track the BEST iterate
    hist.J(end+1)     = J;
    hist.Jbest(end+1) = Jbest;
    hist.vol(end+1)   = sum(t(:))/sum(tmax_e(:));
    hist.chg(end+1)   = chg;
    hist.times(end+1) = toc(t0);

    fprintf('  [%-10s] it %3d: J=%.6e  f0(dB)=%+.4f  Jbest=%.6e  vol=%.4f  chg=%.2e  inner=%d\n', ...
            metric, it, J, f0s, Jbest, hist.vol(end), chg, info.innerit);

    t = reshape(tnew, nel1, nel2);

    if opt.stop && numel(hist.J) >= 2*opt.window + 1
        W = opt.window; fo = hist.f0;
        fstat = abs(mean(fo(end-W+1:end)) - mean(fo(end-2*W+1:end-W))) / ...
                max(abs(mean(fo(end-2*W+1:end-W))), eps);
        xstat = chg / max(opt.tmax - opt.tmin, eps);
        if xstat < opt.xchgtol && fstat < opt.ftol && gv <= opt.feastol
            nholdcnt = nholdcnt + 1;
        else
            nholdcnt = 0;
        end
        if nholdcnt >= opt.nhold
            stopReason = 'converged';
            fprintf('  [%s] STOP at it %d: converged (xstat %.2e, fstat %.2e)\n', ...
                    metric, it, xstat, fstat);
            break;
        end
    end
end

topt = tbest;
hist.stopReason = stopReason;
[~, ib] = min(hist.J); hist.iters_to_best = ib;
hist.upjump = max(0, max(diff(hist.J)));
if numel(hist.J) > 5
    tl = hist.J(max(1,end-49):end);
    hist.tailspread = max(tl)/max(min(tl), eps);
else
    hist.tailspread = NaN;
end
hist.walltime    = toc(t_start);
hist.numeval     = sum(hist.innerit);
hist.metric      = metric;
hist.p0          = p0;
hist.tmin        = opt.tmin;  hist.tmax = opt.tmax;
hist.final_vol_best = sum(tbest(:))/sum(tmax_e(:));
hist.final_vol_last = sum(t(:))/sum(tmax_e(:));

meta = struct('metric',metric, 'niter',numel(hist.J), 'numeval',sum(hist.innerit), ...
              'walltime',hist.walltime, 'final_vol',hist.final_vol_best, ...
              'tailspread',hist.tailspread, 'stopReason',stopReason, 'Jbest',Jbest);
end

% ---------------------------------------------------------------------------
function J = obj_raw(U, F, w, metric, R0)
% Raw (unscaled) objective for the chosen metric -- identical to the value
% harmonic_kd_sens assumes when it returns dJ for that metric.
U = U(:); F = F(:);
switch lower(metric)
    case 'dyncomp'
        J = abs(F.' * U);
    case 'inputpower'
        J = -0.5 * w * imag(conj(F).' * U);     % = 0.5*Re(F^H V), V=j*w*U
    case 'quadvel'
        if isempty(R0)
            J = w^2 * real(U' * U);             % R0 = I
        elseif ischar(R0)
            error('run_metric_opt:obj_raw','char R0 must be resolved by caller.');
        else
            J = w^2 * real(U' * (R0*U));
        end
    otherwise
        error('run_metric_opt:obj_raw','unknown metric ''%s''.', metric);
end
end

% ---------------------------------------------------------------------------
function f0q = obj_scaled(tq, H, F, w, alpha, beta, metric, opt, V0, nel1, nel2, p0)
% Scaled objective at a trial thickness field for GCMMA's conservatism check.
tq = reshape(tq, nel1, nel2);
KS = harmonic_kd_form(H, tq, w, alpha, beta, 'cpu_sparse');
Uq = KS \ F;
Jq = obj_raw(Uq, F, w, metric, opt.R0);
f0q = 100 * (opt.Abar + opt.dbfac * log10(max(Jq, eps))) / p0;
end
