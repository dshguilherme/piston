function [t_opt, hist] = harmonic_topopt(H, F, w, alpha, beta, opt)
% HARMONIC_TOPOPT Thickness optimization for the harmonic objectives, using the
% WQ matrix-free operator, the strain-energy-trick sensitivities and MMA.
%
% opt fields: tmin, tmax, volfrac, nmax (default 100), tol (default 1e-4),
%             penal (SIMP-like exponent on the thickness, default 1)
%
% Objective: dynamic compliance J = |F.'U| (thesis Eq. 3.31), minimised.
% Constraint: total thickness (volume) <= volfrac * sum(tmax).
%
% Each iteration: one sparse harmonic solve (the operator is verified against the
% matrix-free matvec to 1e-16, so either path gives the same state), the
% gradient from harmonic_kd_sens, one MMA step.

% NOTE: this function takes SIX inputs, so the guard must be `nargin < 6`.
% (An earlier `nargin < 7` discarded the caller's opt struct entirely and
% silently ran the defaults -- the same off-by-one as in harmonic_kd_sens.)
if nargin < 6, opt = struct(); end
if ~isfield(opt,'tmin'),    opt.tmin = 0.4; end
if ~isfield(opt,'tmax'),    opt.tmax = 1.0; end
if ~isfield(opt,'volfrac'), opt.volfrac = 0.5; end
if ~isfield(opt,'nmax'),    opt.nmax = 100; end
if ~isfield(opt,'tol'),     opt.tol = 0;    end   % legacy raw design-change stop
% STOPPING CRITERION (opt.stop, default on). One instantaneous test cannot work
% here: the harmonic objective oscillates, so |J_k - J_{k-1}| stays O(1) even
% after the design has settled. Require THREE conditions at once, for opt.nhold
% consecutive iterations:
%   (1) windowed objective stationarity: mean of the last W iterates matches the
%       mean of the previous W to within opt.ftol;
%   (2) design stationarity: max|dx|/(xmax-xmin) < opt.xchgtol;
%   (3) feasibility: the volume constraint holds to opt.feastol.
% hist.stopReason records 'converged' vs 'nmax', so a reported iteration count
% always says WHY it stopped. The KKT residual is reported (hist.kkt) but NOT used
% to stop: it saturates when most elements sit at a bound, so it is not yet
% trustworthy.
% OBJECTIVE SCALING (opt.objscale). For vibration/acoustics the objective spans
% many orders of magnitude (here J goes from ~1e-8 to ~1e-15), which wrecks the
% moving-asymptote models: the conservatism raa0 is in objective units, so an
% 8-order dynamic range makes it either dominant or irrelevant, and the stopping
% tests see O(1) oscillations. Standard practice in vibro-acoustics topology
% optimization is a LOGARITHMIC (dB) objective:
%       f = 10*log10(J / Jref)        (Jref = the first-iterate value)
%       df/dx = (10/ln(10)) * (dJ/dx) / J
% NOTE: this is a MONOTONE transform, so the minimiser is IDENTICAL -- it changes
% the conditioning, the model behaviour and the stopping behaviour, not the
% problem. Use objscale='linear' for the plain |F.'U| and 'none' for raw values.
if ~isfield(opt,'objscale'), opt.objscale = 'log'; end
% SOLVER BACKEND (opt.solver). The article's enabling technology is the
% matrix-free operator + COCG in the complex-symmetric BILINEAR inner product
% (no conjugation). 'sparse' keeps the assembled direct path as a reference for
% verification and for the timing table.
% DEFAULT IS 'sparse' ON PURPOSE: measured at 64x32 near resonance the iterative
% path needs ~thousands of COCG iterations (cond(K_D) ~ 3.4e5) where the direct
% path does ONE factorisation -- 0.37 s/outer-iter vs not finishing one iteration
% in 4 min. Matrix-free pays only where factorisation is infeasible (large DOF)
% or the conditioning is tamed (damping / better preconditioner). Do not flip
% this default without re-measuring at the target mesh.
if ~isfield(opt,'solver'), opt.solver = 'sparse'; end
if ~isfield(opt,'stop'),    opt.stop = true; end
if ~isfield(opt,'ftol'),    opt.ftol = 1e-4; end   % thesis Sec. 4.5 tolerance
if ~isfield(opt,'xchgtol'), opt.xchgtol = 1e-4; end
if ~isfield(opt,'nhold'),   opt.nhold = 3; end
if ~isfield(opt,'feastol'), opt.feastol = 1e-6; end
if ~isfield(opt,'window'),  opt.window = 10; end
% GCMMA IMPLEMENTATION (opt.gcmma_impl). 'std' (default) calls the shared
% gcmma_solver.m; 'fast' calls the drop-in gcmma_solver_fast.m, which is
% numerically identical but does not form an n x n dense diagonal on every
% Newton step of the dual solve (the measured hot spot -- see
% tests/microbench_diag.m). Default is 'std' so behaviour is unchanged unless
% explicitly requested.
if ~isfield(opt,'gcmma_impl'), opt.gcmma_impl = 'std'; end

nel1 = H.nel_dir(1); nel2 = H.nel_dir(2);
t = opt.volfrac * opt.tmax * ones(nel1, nel2);
tmin_e = opt.tmin * ones(nel1, nel2); tmax_e = opt.tmax * ones(nel1, nel2);
st = [];
hist = struct('J', [], 'Jbest', [], 'f0', [], 'vol', [], 'chg', [], 'times', [], ...
              'kkt', [], 'optimizer', '', 'objscale', '', 'itsolve', [], 'solver', '', 'innerit', []);
% PER-PHASE TIMER ACCUMULATOR (instrumentation). Sums wall time (tic/toc) over
% the whole run around each independent phase of the outer iteration so the
% phase split can be read off hist.tphase after >= 50 iterations.
Pacc = struct('form',0,'solve',0,'sens',0,'kkt',0,'gcmma',0,'attach',0);
hist.gcmma_impl = opt.gcmma_impl;
% Choose the GCMMA implementation.
gsolver = @gcmma_solver;
if strcmpi(opt.gcmma_impl,'fast')
    if exist('gcmma_solver_fast','file') == 2
        gsolver = @gcmma_solver_fast;
        hist.gcmma_impl = 'fast';
    else
        warning('harmonic_topopt:noFastGcmma', ...
            'gcmma_impl=''fast'' requested but gcmma_solver_fast.m not on the path -- using std.');
        hist.gcmma_impl = 'std';
    end
end
hist.optimizer = 'gcmma'; if isfield(opt,'optimizer'), hist.optimizer = opt.optimizer; end
hist.stopReason = 'nmax';
hist.objscale = opt.objscale;
Jbest = inf; tbest = t;   % track the BEST iterate, not just the last
nholdcnt = 0;

for it = 1:opt.nmax
    t0 = tic;
    tp = tic;
    switch lower(opt.solver)
        case 'sparse'
            KD = harmonic_kd_form(H, t, w, alpha, beta, 'cpu_sparse');
            Pacc.form = Pacc.form + toc(tp); tp = tic;
            U  = KD \ F;
            Pacc.solve = Pacc.solve + toc(tp);
            itsolve = 0;
        otherwise                       % matrix-free + COCG (bilinear inner product)
            Afun = @(x) harmonic_kd_matvec(H, t, x, w, alpha, beta);
            dd   = harmonic_kd_diag(H, t, w, alpha, beta);
            [U, fl, rr, itsolve] = cocg(Afun, F, @(x) x ./ dd, 1e-10, 4*H.N, []);
            Pacc.solve = Pacc.solve + toc(tp);
            if fl ~= 0
                fprintf('  [warn] COCG flag %d at outer it %d (relres %.2e, %d its)\n', ...
                        fl, it, rr, itsolve);
            end
    end
    hist.itsolve(end+1) = itsolve;
    hist.solver = opt.solver;
    g = F.'*U; J = abs(g);
    tp = tic;
    dJ = harmonic_kd_sens(H, t, F, U, w, alpha, beta, 'dyncomp');
    Pacc.sens = Pacc.sens + toc(tp);
    % volume constraint  g = sum(t)/V0 - 1 <= 0
    tp = tic;
    V0 = opt.volfrac * sum(tmax_e(:));
    gv = sum(t(:))/V0 - 1;
    dgv = (1/V0) * ones(nel1, nel2);
    % NORMALISATION (why): the moving-asymptote conservatism parameter raa0 is
    % initialised in OBJECTIVE UNITS (Svanberg uses 1e-5). Our dynamic compliance
    % is ~1e-8, so an un-normalised objective lets the conservatism term dominate
    % the gradients -- the models then barely move. Scaling f0 and df0 by the
    % first-iterate value makes the model gradient-driven again and leaves the
    % optimum unchanged.
    switch lower(opt.objscale)
        case 'log'
            % THESIS-EXACT dB SCALING (Olavo's thesis, eqs. 6.113-6.115):
            %   {P}_dB = Abar + 10*log10(P)     (Abar keeps it positive for P < 1)
            %   p0     = {P}_dB at the FIRST iteration
            %   min    100*{P}_dB/p0,   d/dx {P}_dB = 10/(ln10 * P) * dP/dx
            % (for velocity-type metrics the thesis uses {R}_dB = 100 + 100*log10(R);
            %  the same pattern with a larger factor -- set opt.dbfac = 100 for that.)
            Abar = 100; if isfield(opt,'Abar'), Abar = opt.Abar; end
            fac  = 10;  if isfield(opt,'dbfac'), fac = opt.dbfac; end
            Jdb = Abar + fac*log10(max(J, eps));
            if ~isfield(opt,'p0') || isempty(opt.p0), opt.p0 = Jdb; end
            f0s  = 100 * Jdb / opt.p0;
            df0s = 100 * (fac/log(10)) * dJ / max(J, eps) / opt.p0;
        case 'linear'
            if ~isfield(opt,'Jscale') || isempty(opt.Jscale), opt.Jscale = max(J, eps); end
            f0s = J / opt.Jscale;  df0s = dJ / opt.Jscale;
        otherwise
            f0s = J;  df0s = dJ;
    end
    hist.f0(end+1) = f0s;

    % stationarity measure: reduced gradient wrt the active set (box + volume).
    % For components at a bound the KKT condition allows an outward gradient, so
    % those are excluded; the rest must satisfy grad f0 + mu grad g = 0.
    dL = df0s(:) + max(0, -(df0s(:).'*dgv(:))/(dgv(:).'*dgv(:))) * dgv(:);
    tflat = t(:); atlo = tflat <= tmin_e(:) + 1e-12; athi = tflat >= tmax_e(:) - 1e-12;
    keep = ~((atlo & dL > 0) | (athi & dL < 0));
    kkt = max(abs(dL(keep))) / max(max(abs(df0s(:))), eps);
    Pacc.kkt = Pacc.kkt + toc(tp);

    % GCMMA TUNING. The conservatism sweep showed the stiff setting cuts the
    % objective thrashing ~9x at equal objective, so it is the default:
    %   kmax 20, growth 1.30, raa0eps/raaeps 1e-3.
    % GHINIT -- Svanberg's asymptote-initialisation constant (the 0.5 in his
    % eq. for L_jk/U_jk, tuned in the thesis' Fig. 5.45 between 0.5 and 0.01) --
    % is opt.asyinit; a SMALLER value tightens the initial asymptotes (stiffer).
    gopt = opt;
    if ~isfield(gopt,'kmax'),    gopt.kmax    = 20;   end
    if ~isfield(gopt,'growth'),  gopt.growth  = 1.30; end
    if ~isfield(gopt,'raa0eps'), gopt.raa0eps = 1e-3; end
    if ~isfield(gopt,'raaeps'),  gopt.raaeps  = 1e-3; end

    % TRUE-EVALUATION CALLBACK -- REQUIRED. Without it gcmma_solver cannot check
    % whether its model is conservative, and it used to fall back silently to a
    % tangent model whose conservatism test is vacuously true: the inner loop
    % never ran and kmax/growth were inert (root cause found by investigation).
    % The callback MUST return the SAME scaling the driver passes in f0s --
    % feeding raw |F.'U| while f0s is dB-scaled compares incompatible quantities.
    gopt.eval = @(tq) deal( ...
        obj_scaled(tq, H, F, w, alpha, beta, opt, V0, nel1, nel2), ...
        sum(tq(:))/V0 - 1 );
    use_mma = isfield(opt,'optimizer') && strcmpi(opt.optimizer,'mma');
    if ~use_mma && exist('gcmma_solver', 'file') == 2
        tp = tic;
        [tnew, st, info] = gsolver(t, tmin_e, tmax_e, f0s, df0s, gv, dgv, st, it, gopt);
        Pacc.gcmma = Pacc.gcmma + toc(tp);
        hist.innerit(it) = info.innerit;   % true evaluations this outer iteration
    else
        warning('harmonic_topopt:gcmmaMissing', ...
            'gcmma_solver.m not found on the path -- falling back to MMA, which is NOT recommended for vibration objectives.');
        [tnew, st] = mma_solver(t, tmin_e, tmax_e, f0s, df0s, gv, dgv, st, it);
    end
    chg = max(abs(tnew(:) - t(:)));
    if J < Jbest, Jbest = J; tbest = t; end
    hist.J(end+1) = J; hist.Jbest(end+1) = Jbest; hist.kkt(end+1) = kkt;
    hist.vol(end+1) = sum(t(:))/sum(tmax_e(:));
    hist.chg(end+1) = chg; hist.times(end+1) = toc(t0);
    t = reshape(tnew, nel1, nel2);   % solvers return a flattened n-vector; the
                                     % operator (harmonic_kd_form) wants the grid.
    if opt.stop && numel(hist.J) >= 2*opt.window + 1
        W = opt.window;
        % test on the SCALED objective actually being optimised: with dB scaling
        % the raw J oscillates by orders of magnitude and would never satisfy it
        fo = hist.f0;
        if numel(fo) < 2*W+1, fo = hist.J; end
        fstat = abs(mean(fo(end-W+1:end)) - mean(fo(end-2*W+1:end-W))) / ...
                max(abs(mean(fo(end-2*W+1:end-W))), eps);
        xstat = chg / max(opt.tmax - opt.tmin, eps);
        if xstat < opt.xchgtol && fstat < opt.ftol && gv <= opt.feastol
            nholdcnt = nholdcnt + 1;
        else
            nholdcnt = 0;
        end
        if nholdcnt >= opt.nhold
            hist.stopReason = 'converged';
            fprintf(['  STOP at iteration %d: converged (xstat %.2e < %.1e, ' ...
                     'fstat %.2e < %.1e, feasible)\n'], it, xstat, opt.xchgtol, fstat, opt.ftol);
            break;
        end
    end
    if chg < opt.tol
        hist.stopReason = 'converged';
        fprintf('  converged at iteration %d (legacy design change %.2e < %.1e)\n', it, chg, opt.tol);
        break;
    end
    if mod(it, 10) == 0 || it == 1
        fprintf('  it %3d: J = %.6e  f0(dB) = %+.4f  Jbest = %.6e  volfrac = %.4f  chg = %.2e\n', ...
                it, J, f0s, Jbest, hist.vol(end), chg);
    end
end
if strcmp(hist.stopReason,'nmax')
    fprintf('  STOP: budget exhausted (nmax = %d) WITHOUT meeting the criterion\n', opt.nmax);
end
t_opt = tbest;   % return the BEST design found (the last one is often mid-oscillation)
% NOTE: no monotonicity metric. GCMMA's approximations are conservative rather
% than descent-based, so objective increases between outer iterations are BY
% DESIGN (they are how it leaves a local basin). Report best-so-far J and the
% stationarity measure instead.
[~, ib] = min(hist.J); hist.iters_to_best = ib;
% TRAJECTORY STABILITY diagnostics. upjump = largest single upward step.
% tailspread = max/min over the last 50 iterations -- but that window is only
% meaningful for runs of AT LEAST 50 iterations. With the conservatism loop live
% the driver now stops early (typically ~29), and a 'last 50' window then spans
% the whole startup descent: it reported 1e6-1e9 on a perfectly settled design.
% Two independent investigations hit this, so both windows are reported:
%   hist.tailspread    long window; NaN unless the run reached 50 iterations
%   hist.tailspread10  max/min over the last 10; valid for short runs
hist.upjump = max(0, max(diff(hist.J)));
tl10 = hist.J(max(1,end-9):end);
hist.tailspread10 = max(tl10)/max(min(tl10), eps);
if numel(hist.J) >= 50
    tl = hist.J(end-49:end);
    hist.tailspread = max(tl)/max(min(tl), eps);
else
    hist.tailspread = NaN;   % <50 iters: the long-window number is an artifact
end
hist.tbest = tbest;
% PER-PHASE SPLIT (instrumentation): totals, per-iteration means and percent of
% the summed phase time. Read hist.tphase_mean / hist.tphase_pct after a run.
n_it = numel(hist.times);
tot_ph = Pacc.form + Pacc.solve + Pacc.sens + Pacc.kkt + Pacc.gcmma;
hist.tphase = Pacc; hist.tphase_tot = tot_ph;
if n_it > 0
    denom = max(tot_ph, eps);
    hist.tphase_mean = struct('form',Pacc.form/n_it, 'solve',Pacc.solve/n_it, ...
        'sens',Pacc.sens/n_it, 'kkt',Pacc.kkt/n_it, 'gcmma',Pacc.gcmma/n_it, ...
        'total',tot_ph/n_it);
    hist.tphase_pct = struct('form',100*Pacc.form/denom, 'solve',100*Pacc.solve/denom, ...
        'sens',100*Pacc.sens/denom, 'kkt',100*Pacc.kkt/denom, 'gcmma',100*Pacc.gcmma/denom);
end
end


function f0q = obj_scaled(tq, H, F, w, alpha, beta, opt, V0, nel1, nel2)
% Scaled objective at a trial thickness field for GCMMA's conservatism check.
% Uses the SAME transform as the main loop, so the model is compared against
% like-with-like. Costs one weighted-quadrature formation + solve per call,
% which is the irreducible price of a correct conservatism check.
    tq  = reshape(tq, nel1, nel2);
    KS  = harmonic_kd_form(H, tq, w, alpha, beta, 'cpu_sparse');
    Jq  = abs(F.' * (KS \ F));
    switch lower(opt.objscale)
        case 'log'
            Abar = 100; if isfield(opt,'Abar'),  Abar = opt.Abar;  end
            fac  = 10;  if isfield(opt,'dbfac'), fac  = opt.dbfac; end
            f0q  = 100 * (Abar + fac*log10(max(Jq, eps))) / opt.p0;
        case 'linear'
            f0q  = Jq / opt.Jscale;
        otherwise
            f0q  = Jq;
    end
end
