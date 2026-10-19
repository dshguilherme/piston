function [xnew, st, info] = gcmma_solver(x, xmin, xmax, f0val, df0dx, fval, dfdx, st, it, opt)
%GCMMA_SOLVER One OUTER iteration of Svanberg's Globally Convergent Method of
% Moving Asymptotes (GCMMA, 2002), specialised to the harmonic/thickness
% topology-optimization driver of this project.
%
% WHY GCMMA AND NOT MMA FOR VIBRATION PROBLEMS
% ---------------------------------------------------------------------------
% The plain MMA subproblem is only a LOCAL model: it is built from the value and
% gradient at x^k and, away from x^k, it need not lie above the true function.
% For resonance-sensitive objectives (J = |F.'U| near a natural frequency) the
% true function can be far above its MMA model between iterations, so MMA steps
% can be accepted although they make J worse: the iteration then stalls or
% oscillates. GCMMA closes that gap with an INNER loop that keeps the model
% accountable: it solves the subproblem, evaluates the TRUE objective and
% constraints at the trial point, and only ACCEPTS the step once the separable
% model is provably CONSERVATIVE there (model >= true). When it is not, the
% curvature of the model is raised (raa0, raa are increased) and the subproblem
% is re-solved -- a monotonically growing trust region around the current point.
% This is what makes GCMMA globally convergent and what makes it monotone on the
% harmonic objective where MMA is not. (The user's PhD thesis used GCMMA for all
% vibration problems for exactly this reason.)
%
% FAITHFUL IMPLEMENTATION
% ---------------------------------------------------------------------------
% This file is a direct MATLAB port of Krister Svanberg's reference code
%   mmasub.m / gcmmasub.m / subsolv.m / asymp.m / raaupdate.m / concheck.m
% (KTH, 2007/2008; equivalently the line-by-line Python translation by Arjen
% Deetman, GCMMA-MMA-Python). The outer-iteration flowchart is
%
%   (1) asymp():  update asymptotes low/upp with the oscillation heuristic and
%                 RESET the conservatism parameters raa0, raa from the gradient
%                 (this reset is the "reduction of conservatism" the task asks
%                 for -- it happens every outer iteration, not inside the inner
%                 loop).
%   (2) gcmmasub(): build the separable convex model
%                     fhat_i(x) = r_i + sum_j [ p_ij/(U_j-x_j) + q_ij/(x_j-L_j) ]
%                 with r_i pinned so that fhat_i(x^k) = f_i(x^k), and solve its
%                 dual through subsolv() (a primal-dual Newton method on the
%                 multipliers). The artificial variables y_i >= 0 and z >= 0
%                 make the subproblem ALWAYS FEASIBLE (y_i absorbs any
%                 constraint violation; z relaxes all constraints when a_i>0),
%                 so the subproblem never fails no matter how aggressive the
%                 model is.
%   (3) concheck(): evaluate the TRUE f0, f_i at the trial point x^mma and test
%                 fhat_0(x^mma) >= f_0(x^mma) and fhat_i >= f_i for all i.
%   (4) inner loop: while the model is NOT conservative and innerit < kmax,
%                 raaupdate() raises raa0/raa (curvature up => model up) and
%                 gcmmasub() is re-solved from the (frozen) current point x^k.
%                 The first conservative trial point is the accepted x^{k+1}.
%
% SIGNATURE
%   [xnew, st, info] = gcmma_solver(x, xmin, xmax, f0val, df0dx, fval, dfdx, st, it, opt)
%     x       (n x 1)        current design (flattened thickness field)
%     xmin,xmax (n x 1)      box bounds
%     f0val   scalar         objective J at x
%     df0dx   (n x 1)        dJ/dx at x   (from harmonic_kd_sens)
%     fval    (m x 1)        constraint values at x  (volume: sum(t)/V0 - 1)
%     dfdx    (m x n)        constraint Jacobian at x (volume: 1/V0)
%     st      persistent state from the previous call (asymptotes low/upp,
%             raa0/raa, xold1/xold2); pass [] on the first call.
%     it      1-based OUTER iteration number
%     opt     optional knobs (see below)
%   xnew, st, info as described below.
%
% REQUIRED opt.eval FOR THE TRUE INNER LOOP
%   The signature carries no objective/constraint callback, so the inner loop
%   gets it from opt:
%       opt.eval = @(xq) deal(f0q, fq)   % true objective and constraints at xq
%   e.g. for the driver:  opt.eval = @(tq) deal(Jfun(tq), volfun(tq)).
%   opt.eval is REQUIRED. If it is omitted the inner loop has no true values to
%   test against, so gcmma_solver raises 'gcmma_solver:missingEval' rather than
%   silently degenerating: a "linear" fallback compared the model against its own
%   tangent, which is provably always "conservative", making the inner loop dead
%   (info.innerit == 0 forever) and voiding GCMMA's global-convergence guarantee.
%   NOTE: the values returned by opt.eval MUST use the SAME objective/constraint
%   scaling as f0val/fval (e.g. the dB-scaled objective the harmonic driver
%   passes), otherwise the conservatism test compares incompatible quantities.
%   harmonic_topopt.m does NOT construct opt.eval itself; the caller must set it.
%
% opt knobs (all optional):
%   opt.eval      function handle @(xq)[f0q,fq]  (see above)
%   opt.kmax      max inner iterations per outer iteration      (default 10)
%   opt.raa0eps   minimum of the objective conservatism raa0      (default 1e-6)
%   opt.raaeps    minimum of the constraint conservatism raa       (default 1e-6)
%   opt.growth    conservatism growth factor in raaupdate          (default 1.1)
%   opt.epsimin   inner-loop tolerance                             (default 1e-7)
%   opt.albefa    bound shrink factor (alfa/beta)                  (default 0.1)
%   opt.move,opt.asyinit,opt.asydecr,opt.asyincr,opt.asymin,opt.asymax
%                 asymptote/move-limit knobs (Svanberg defaults)
%   opt.a0,opt.a,opt.c,opt.d  artificial-variable costs. Recommended values:
%                 a0=1, a=zeros(m,1), c=1000*ones(m,1), d=ones(m,1).
%                 a_i>0 lets z relax constraints i; large c_i penalises the
%                 feasibility slacks y_i so the true (unrelaxed) constraint is
%                 effectively enforced -- important for the volume constraint.
%
% info fields: innerit (inner iterations actually used), raa0, raa (final
%   conservatism parameters), conserv (1 if the accepted point is conservative),
%   f0app/fapp (model values at the accepted point), f0new/fnew (true values),
%   mode ('gcmma'; opt.eval is mandatory, see above).
%
% This file also contains gcmma_asymp, gcmma_gcmmasub, gcmma_subsolv,
% gcmma_concheck, gcmma_raaupdate -- the faithful port of Svanberg's routines.

%% ---------------------------------------------------------------- opt defaults
if nargin < 10 || isempty(opt), opt = struct(); end
if ~isfield(opt,'kmax'),    opt.kmax    = 10;   end
if ~isfield(opt,'raa0eps'), opt.raa0eps = 1e-6; end
if ~isfield(opt,'raaeps'),  opt.raaeps  = 1e-6; end
if ~isfield(opt,'growth'),  opt.growth  = 1.1;  end   % Svanberg's raaupdate factor
if ~isfield(opt,'epsimin'), opt.epsimin = 1e-7; end
if ~isfield(opt,'albefa'),  opt.albefa  = 0.1;  end
if ~isfield(opt,'move'),    opt.move    = 0.5;  end
if ~isfield(opt,'asyinit'), opt.asyinit = 0.5;  end
if ~isfield(opt,'asydecr'), opt.asydecr = 0.7;  end
if ~isfield(opt,'asyincr'), opt.asyincr = 1.2;  end
if ~isfield(opt,'asymin'),  opt.asymin  = 0.01; end
if ~isfield(opt,'asymax'),  opt.asymax  = 10;   end
if ~isfield(opt,'a0'),      opt.a0      = 1;    end

%% ---------------------------------------------------------------- shapes
n = numel(x);
x     = x(:);
xmin  = xmin(:);
xmax  = xmax(:);
df0dx = df0dx(:);
m     = numel(fval);
fval  = fval(:);
dfdx  = reshape(dfdx, m, n);              % m x n even when the caller passes n x 1

if ~isfield(opt,'a'), opt.a = zeros(m,1);       end
if ~isfield(opt,'c'), opt.c = 1000*ones(m,1);   end
if ~isfield(opt,'d'), opt.d = ones(m,1);        end
opt.a = opt.a(:); opt.c = opt.c(:); opt.d = opt.d(:);

% The inner conservatism loop can only work with the TRUE objective/constraints at
% the trial point, and the signature cannot supply them: they must come from
% opt.eval. A missing callback used to fall back SILENTLY to a first-order (linear)
% model of the true functions, and gcmma_concheck() then compared the Svanberg
% model against its OWN tangent at x. That comparison is provably ALWAYS true:
%   model - tangent = (xmma-x)^2 .* [ df0+/(upp-xmma) + df0-/(xmma-low)
%       + (0.001*(df0+ + df0-) + raa0/xmami) * (upp-low)/((upp-xmma).*(xmma-low)) ]
% is elementwise >= 0 (every factor is non-negative on the alfa/beta box), and
% identically for each constraint. So `conserv` was always 1, the body of
% `while conserv == 0` never ran (info.innerit == 0 forever), kmax/growth became
% inert, and GCMMA's global-convergence guarantee was silently voided -- it
% degenerated to a plain MMA step. Refuse to run without the true-eval callback.
if ~isfield(opt,'eval') || isempty(opt.eval)
    error('gcmma_solver:missingEval', ...
        ['opt.eval is required: the GCMMA inner loop must evaluate the TRUE ', ...
         'objective/constraints at the trial point or it cannot verify ', ...
         'conservatism and would silently degenerate to MMA. Pass ', ...
         'opt.eval = @(xq) deal(f0q, fq) (same units as f0val/fval).']);
end
mode = 'gcmma';

%% ---------------------------------------------------------------- persistent state
% On the first call initialise the asymptotes around x and seed raa0/raa; the
% oscillation heuristic needs two previous points, so xold1=xold2=x here and
% asymp() uses the "outeriter<=2" branch for the first two outer iterations.
if isempty(st) || ~isstruct(st) || ~isfield(st,'low') || isempty(st.low)
    st.low   = x - opt.asyinit*(xmax - xmin);
    st.upp   = x + opt.asyinit*(xmax - xmin);
    st.xold1 = x;
    st.xold2 = x;
    st.raa0  = opt.raa0eps;
    st.raa   = opt.raaeps*ones(m,1);
end

%% ---------------------------------------------------------------- (1) asymptotes + raa reset
[low, upp, raa0, raa] = gcmma_asymp(it, n, x, st.xold1, st.xold2, xmin, xmax, ...
                                    st.low, st.upp, st.raa0, st.raa, ...
                                    opt.raa0eps, opt.raaeps, df0dx, dfdx, opt);

%% ---------------------------------------------------------------- (2..4) inner conservatism loop
innerit = 0;
[xmma, f0app, fapp] = gcmma_gcmmasub(m, n, opt.epsimin, x, xmin, xmax, low, upp, ...
                                     raa0, raa, f0val, df0dx, fval, dfdx, opt);
[f0new, fnew] = eval_true(opt, xmma, x, f0val, df0dx, fval, dfdx);
conserv = gcmma_concheck(m, opt.epsimin, f0app, f0new, fapp, fnew);

while conserv == 0 && innerit < opt.kmax
    innerit = innerit + 1;
    % Raise the model curvature around the (frozen) current point x, re-solve,
    % re-evaluate the TRUE functions, re-test conservativeness.
    [raa0, raa] = gcmma_raaupdate(xmma, x, xmin, xmax, low, upp, ...
                                  f0new, fnew, f0app, fapp, raa0, raa, opt);
    [xmma, f0app, fapp] = gcmma_gcmmasub(m, n, opt.epsimin, x, xmin, xmax, low, upp, ...
                                         raa0, raa, f0val, df0dx, fval, dfdx, opt);
    [f0new, fnew] = eval_true(opt, xmma, x, f0val, df0dx, fval, dfdx);
    conserv = gcmma_concheck(m, opt.epsimin, f0app, f0new, fapp, fnew);
end

%% ---------------------------------------------------------------- state + info
xnew = xmma;
st.xold2 = st.xold1;
st.xold1 = x;
st.low   = low;
st.upp   = upp;
st.raa0  = raa0;
st.raa   = raa;

info = struct('innerit', innerit, 'raa0', raa0, 'raa', raa, 'conserv', conserv, ...
              'f0app', f0app, 'fapp', fapp, 'f0new', f0new, 'fnew', fnew, 'mode', mode);
end

% =========================================================================
% Local functions: faithful port of Svanberg's GCMMA-MMA reference routines
% =========================================================================

function [f0new, fnew] = eval_true(opt, xmma, x, f0val, df0dx, fval, dfdx)
% True objective/constraint values at the trial point. Uses opt.eval when given
% (the real GCMMA check); otherwise a first-order model about x (documented
% degradation -- see header).
if ~isempty(opt.eval)
    [f0new, fnew] = opt.eval(xmma);
    fnew = fnew(:);
else
    f0new = f0val + df0dx.'*(xmma - x);
    fnew  = fval  + dfdx  *(xmma - x);
end
end

function [low, upp, raa0, raa] = gcmma_asymp(outeriter, n, xval, xold1, xold2, ...
        xmin, xmax, low, upp, raa0, raa, raa0eps, raaeps, df0dx, dfdx, opt)
% Svanberg asymp(): asymptote update (oscillation heuristic) + raa reset.
% WHY the raa reset: at the start of every OUTER iteration the conservatism
% parameters are recomputed from the current gradient magnitude. This is the
% "reduction" step: a fresh outer iteration starts with a SMALL (gradient-scaled)
% curvature, and raaupdate() only ever raises it inside the inner loop.
xmami = max(xmax - xmin, 1e-5);            % guard against zero-width boxes
raa0  = max(raa0eps, (0.1/n) * (abs(df0dx).' * xmami));
raa   = max(raaeps, (0.1/n) * (abs(dfdx)  * xmami));   % m x 1

if outeriter <= 2
    % No history: symmetric asymptotes at asyinit times the box width.
    low = xval - opt.asyinit * xmami;
    upp = xval + opt.asyinit * xmami;
else
    % Oscillation heuristic: if x_j moved in the same direction twice, widen the
    % asymptotes (faster, less conservative); if it oscillated, shrink them.
    xxx = (xval - xold1) .* (xold1 - xold2);
    factor = ones(n,1);
    factor(xxx > 0) = opt.asyincr;
    factor(xxx < 0) = opt.asydecr;
    low = xval - factor .* (xold1 - low);
    upp = xval + factor .* (upp - xold1);
    % Clamp the asymptote distances to [asymin, asymax] box widths.
    low = max(low, xval - opt.asymax * xmami);
    low = min(low, xval - opt.asymin * xmami);
    upp = min(upp, xval + opt.asymax * xmami);
    upp = max(upp, xval + opt.asymin * xmami);
end
end

function [xmma, f0app, fapp] = gcmma_gcmmasub(m, n, epsimin, xval, xmin, xmax, ...
        low, upp, raa0, raa, f0val, df0dx, fval, dfdx, opt)
% Svanberg gcmmasub(): build the separable convex model and solve its dual.
%
% Model (i = 0 objective, i = 1..m constraints):
%   fhat_i(x) = r_i + sum_j [ p_ij/(U_j - x_j) + q_ij/(x_j - L_j) ]
% with Svanberg's coefficient rules p_ij = (U_j-x_j^k)^2 (1.001 d^+ + 0.001 d^- + raa_i/dx_j)
% (and 0.001/1.001 swapped for q_ij), r_i pinned so fhat_i(x^k)=f_i(x^k).
% The 1.001/0.001 split biases the model slightly so its gradient at x^k matches
% the true gradient to first order; raa_i/dx_j adds the conservatism curvature.

xval = xval(:); xmin = xmin(:); xmax = xmax(:); low = low(:); upp = upp(:);

xmami     = max(xmax - xmin, 1e-5);
xmami_inv = 1 ./ xmami;

ux1 = upp - xval;  ux2 = ux1.^2;  ux_inv = 1 ./ ux1;
xl1 = xval - low;  xl2 = xl1.^2;  xl_inv = 1 ./ xl1;

% --- objective coefficients p0, q0 and the offset r0 ------------------------
p0 = max(df0dx, 0);
q0 = max(-df0dx, 0);
pq0 = p0 + q0;
p0 = p0 + 0.001*pq0 + raa0*xmami_inv;
q0 = q0 + 0.001*pq0 + raa0*xmami_inv;
p0 = p0 .* ux2;
q0 = q0 .* xl2;
r0 = f0val - p0.'*ux_inv - q0.'*xl_inv;

% --- constraint coefficients P, Q (m x n) and the offsets r, b --------------
P = max(dfdx, 0);
Q = max(-dfdx, 0);
PQ = P + Q;
raa_col = raa * xmami_inv.';                 % m x n outer product: raa_i / dx_j
P = P + 0.001*PQ + raa_col;
Q = Q + 0.001*PQ + raa_col;
P = P .* (ux2.');                            % scale column j by (U_j-x_j)^2
Q = Q .* (xl2.');
r = fval - P*ux_inv - Q*xl_inv;              % m x 1 model offset
b = -r;                                      % RHS handed to subsolv

% --- inner box [alfa,beta]: keep the trial point strictly inside (low,upp) ----
alfa = max(low + opt.albefa*(xval - low), xmin);
beta = min(upp - opt.albefa*(upp - xval), xmax);
beta = max(beta, alfa + 1e-12);              % guard a degenerate box

% --- solve the dual (Newton/bisection on the multipliers) -------------------
xmma = gcmma_subsolv(m, n, epsimin, low, upp, alfa, beta, p0, q0, P, Q, ...
                     opt.a0, opt.a, b, opt.c, opt.d);

% --- model values at the trial point (for concheck) -------------------------
ux1 = upp - xmma;  xl1 = xmma - low;
f0app = r0 + p0.'*(1./ux1) + q0.'*(1./xl1);
fapp  = r  + P*(1./ux1)  + Q*(1./xl1);
end

function [x, y, z, lam, xsi, eta, mu, zet, s] = gcmma_subsolv(m, n, epsimin, ...
        low, upp, alfa, beta, p0, q0, P, Q, a0, a, b, c, d)
% Svanberg subsolv(): primal-dual Newton method for the separable subproblem
%   min  sum_j[p0j/(Uj-xj)+q0j/(xj-Lj)] + a0*z + sum_i[ci*yi+0.5*di*(yi)^2]
%   s.t. sum_j[pij/(Uj-xj)+qij/(xj-Lj)] - ai*z - yi <= bi, alfa<=x<=beta, y,z>=0.
% It iterates x, y, z, lam, xsi, eta, mu, zet, s to the KKT point, shrinking the
% barrier parameter epsi by 0.1 each pass and taking a damped Newton step with a
% residual-decreasing line search. Two linear-system variants are used: a small
% m x m dual system when m<n (the common case, e.g. one volume constraint) and
% an n x n primal system otherwise. This is a faithful port of the reference.

een = ones(n,1);  eem = ones(m,1);
epsi = 1;
epsvecn = epsi*een;  epsvecm = epsi*eem;
x   = 0.5*(alfa + beta);
y   = eem;
z   = 1;
lam = eem;
xsi = max(een./(x - alfa), een);
eta = max(een./(beta - x), een);
mu  = max(0.5*c, eem);
zet = 1;
s   = eem;
itera = 0;                                                   %#ok<NASGU>

while epsi > epsimin
    epsvecn = epsi*een;  epsvecm = epsi*eem;
    ux1 = upp - x;  xl1 = x - low;  ux2 = ux1.^2;  xl2 = xl1.^2;
    uxinv1 = 1./ux1;  xlinv1 = 1./xl1;
    plam = p0 + P.'*lam;
    qlam = q0 + Q.'*lam;
    gvec = P*uxinv1 + Q*xlinv1;
    dpsidx = plam./ux2 - qlam./xl2;
    rex   = dpsidx - xsi + eta;
    rey   = c + d.*y - mu - lam;
    rez   = a0 - zet - a.'*lam;
    relam = gvec - a*z - y + s - b;
    rexsi = xsi.*(x - alfa) - epsvecn;
    reeta = eta.*(beta - x) - epsvecn;
    remu  = mu.*y - epsvecm;
    rezet = zet*z - epsi;
    res   = lam.*s - epsvecm;
    residu    = [rex; rey; rez; relam; rexsi; reeta; remu; rezet; res];
    residunorm = norm(residu);
    residumax  = max(abs(residu));
    ittt = 0;

    while residumax > 0.9*epsi && ittt < 200
        ittt = ittt + 1;  itera = itera + 1;
        ux1 = upp - x;  xl1 = x - low;
        ux2 = ux1.^2;   xl2 = xl1.^2;
        ux3 = ux1.*ux2; xl3 = xl1.*xl2;
        uxinv1 = 1./ux1; xlinv1 = 1./xl1; uxinv2 = 1./ux2; xlinv2 = 1./xl2;
        plam = p0 + P.'*lam;
        qlam = q0 + Q.'*lam;
        gvec = P*uxinv1 + Q*xlinv1;
        GG   = P.*(uxinv2.') - Q.*(xlinv2.');        % d(gvec)/dx, m x n
        dpsidx = plam./ux2 - qlam./xl2;
        delx   = dpsidx - epsvecn./(x - alfa) + epsvecn./(beta - x);
        dely   = c + d.*y - lam - epsvecm./y;
        delz   = a0 - a.'*lam - epsi/z;
        dellam = gvec - a*z - y - b + epsvecm./lam;
        diagx  = plam./ux3 + qlam./xl3;
        diagx  = 2*diagx + xsi./(x - alfa) + eta./(beta - x);
        diagxinv = 1./diagx;
        diagy  = d + mu./y;
        diagyinv = 1./diagy;
        diaglam = s./lam;
        diaglamyi = diaglam + diagyinv;

        % ---- Newton system: dual (m<n) or primal (m>=n) variant ------------
        if m < n
            blam = dellam + dely./diagy - GG*(delx./diagx);
            bb   = [blam; delz];
            Alam = diag(diaglamyi) + (GG.*diagxinv.')*GG.';   % O(n) column scaling: was materialising a dense n x n diag (4110x slower)
            AA   = [Alam, a; a.', -zet/z];
            solut = AA\bb;
            dlam = solut(1:m);
            dz   = solut(m+1);
            dx   = -delx./diagx - (GG.'*dlam)./diagx;
        else
            diaglamyiinv = 1./diaglamyi;
            dellamyi = dellam + dely./diagy;
            Axx = diag(diagx) + (GG.'.*diaglamyiinv.')*GG;      % O(n) column scaling, same identity
            azz = zet/z + a.'*(a./diaglamyi);
            axz = -GG.'*(a./diaglamyi);
            bx  = delx + GG.'*(dellamyi./diaglamyi);
            bz  = delz - a.'*(dellamyi./diaglamyi);
            AA  = [Axx, axz; axz.', azz];
            bb  = [-bx; -bz];
            solut = AA\bb;
            dx   = solut(1:n);
            dz   = solut(n+1);
            dlam = (GG*dx)./diaglamyi - dz*(a./diaglamyi) + dellamyi./diaglamyi;
        end

        dy   = -dely./diagy + dlam./diagy;
        dxsi = -xsi + epsvecn./(x - alfa) - xsi.*dx./(x - alfa);
        deta = -eta + epsvecn./(beta - x) + eta.*dx./(beta - x);
        dmu  = -mu + epsvecm./y - mu.*dy./y;
        dzet = -zet + epsi/z - zet.*dz./z;
        ds   = -s + epsvecm./lam - s.*dlam./lam;

        xx  = [y; z; lam; xsi; eta; mu; zet; s];
        dxx = [dy; dz; dlam; dxsi; deta; dmu; dzet; ds];

        % ---- fraction-to-the-boundary step --------------------------------
        stepxx = -1.01*dxx./xx;          stmxx   = max(stepxx);
        stepalfa = -1.01*dx./(x - alfa); stmalfa = max(stepalfa);
        stepbeta =  1.01*dx./(beta - x); stmbeta = max(stepbeta);
        stmalbe   = max(stmalfa, stmbeta);
        stmalbexx = max(stmalbe, stmxx);
        stminv    = max(stmalbexx, 1.0);
        steg      = 1.0/stminv;

        % ---- damped update + residual-decreasing line search --------------
        xold=x; yold=y; zold=z; lamold=lam; xsiold=xsi; etaold=eta;
        muold=mu; zetold=zet; sold=s;
        itto = 0;  resinew = 2*residunorm;
        while resinew > residunorm && itto < 50
            itto = itto + 1;
            x   = xold   + steg*dx;
            y   = yold   + steg*dy;
            z   = zold   + steg*dz;
            lam = lamold + steg*dlam;
            xsi = xsiold + steg*dxsi;
            eta = etaold + steg*deta;
            mu  = muold  + steg*dmu;
            zet = zetold + steg*dzet;
            s   = sold   + steg*ds;
            ux1 = upp - x;  xl1 = x - low;  ux2 = ux1.^2;  xl2 = xl1.^2;
            uxinv1 = 1./ux1;  xlinv1 = 1./xl1;
            plam = p0 + P.'*lam;
            qlam = q0 + Q.'*lam;
            gvec = P*uxinv1 + Q*xlinv1;
            dpsidx = plam./ux2 - qlam./xl2;
            rex   = dpsidx - xsi + eta;
            rey   = c + d.*y - mu - lam;
            rez   = a0 - zet - a.'*lam;
            relam = gvec - a*z - y + s - b;
            rexsi = xsi.*(x - alfa) - epsvecn;
            reeta = eta.*(beta - x) - epsvecn;
            remu  = mu.*y - epsvecm;
            rezet = zet*z - epsi;
            res   = lam.*s - epsvecm;
            residu = [rex; rey; rez; relam; rexsi; reeta; remu; rezet; res];
            resinew = norm(residu);
            steg = steg/2;
        end
        residunorm = resinew;
        residumax  = max(abs(residu));
        steg = 2*steg;
    end
    epsi = 0.1*epsi;
end
end

function conserv = gcmma_concheck(m, epsimin, f0app, f0new, fapp, fnew)
% Svanberg concheck(): is the model conservative at the trial point?
% Conservative iff model >= true for the objective AND every constraint, each
% with a half-epsimin cushion so near-ties count as conservative.
%   f0app + epsimin >= f0new   AND   fapp + epsimin >= fnew (elementwise)
conserv = all([f0app + epsimin; fapp + epsimin*ones(m,1)] >= [f0new; fnew]);
end

function [raa0, raa] = gcmma_raaupdate(xmma, xval, xmin, xmax, low, upp, ...
        f0new, fnew, f0app, fapp, raa0, raa, opt)
% Svanberg raaupdate(): raise the model curvature when it was not conservative.
% raacof measures how far x^mma moved inside (low,upp): the further the trial
% point, the less curvature must be added to lift the model over the true value.
% Only the entries with a genuine shortfall (f > fapp) get increased; increases
% are capped at 10x to avoid a runaway.

raacofmin = 1e-12;
xmami = max(xmax - xmin, 1e-5);
xxux  = (xmma - xval)./(upp - xmma);
xxxl  = (xmma - xval)./(xmma - low);
xxul  = xxux.*xxxl;
ulxx  = (upp - low)./xmami;
raacof = max(sum(xxul.*ulxx), raacofmin);

f0appe = f0app + 0.5*opt.epsimin;
if all(f0new > f0appe)
    deltaraa0 = (1.0/raacof)*(f0new - f0app);
    raa0 = min(opt.growth*(raa0 + deltaraa0), 10*raa0);
end

fappe  = fapp + 0.5*opt.epsimin;
fdelta = fnew - fappe;
deltaraa = (1.0/raacof)*(fnew - fapp);
zzz = min(opt.growth*(raa + deltaraa), 10*raa);
raa(fdelta > 0) = zzz(fdelta > 0);
end
