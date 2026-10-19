function [xnew, st] = mma_solver(x, xmin, xmax, f0, df0, fval, dfdx, st, it)
% MMA_SOLVER One outer iteration of Svanberg's Method of Moving Asymptotes
%
% NOTE: kept for COMPARISON and for smooth problems only. Do NOT use MMA for the
% vibration/harmonic objectives of this project: its subproblem is not
% conservative away from the current point, so it stalls or oscillates when the
% objective is resonance-sensitive. Use GCMMA (gcmma_solver.m) there -- that is
% also what the underlying PhD thesis used for all vibration problems.
% (mma_solver is the mmasub+subsolv pair specialised to ONE constraint, which
% is all the thickness formulation needs: the volume constraint).
%
%   x     (n x 1) current design      xmin, xmax  box bounds
%   f0    scalar objective            df0  (n x 1) its gradient
%   fval  scalar constraint value     dfdx (n x 1) its gradient
%   st    solver state (asymptotes) returned by the previous call, or []
%   it    outer iteration number (1-based)
%
% The separable convex subproblem is
%   min  sum_j [ p0j/(Uj-xj) + q0j/(xj-Lj) ]
%   s.t. sum_j [ p1j/(Uj-xj) + q1j/(xj-Lj) ] <= b1 ,  alpha <= x <= beta
% and is solved through its dual: for a Lagrange multiplier y >= 0 the inner
% minimiser has the closed form
%   xj(y) = (Lj*sqrt(Bj) + Uj*sqrt(Aj)) / (sqrt(Aj) + sqrt(Bj)),  Aj = p0j+y*p1j,
%                                                                Bj = q0j+y*q1j
% clamped to the box; y is found by bisection so that the constraint is met.
% The asymptote update uses Svanberg's standard oscillation heuristic.

n = numel(x);
x = x(:); xmin = xmin(:); xmax = xmax(:); df0 = df0(:); dfdx = dfdx(:);

% --- asymptotes ---
if it <= 2 || isempty(st)
    st.low = x - 0.5*(xmax - xmin);
    st.upp = x + 0.5*(xmax - xmin);
    st.xold1 = x; st.xold2 = x;
else
    zzz = (x - st.xold1) .* (st.xold1 - st.xold2);
    st.low = x - 0.7*(st.xold1 - st.low);
    st.upp = x + 0.7*(st.upp - st.xold1);
    st.low(zzz < 0) = x(zzz < 0) - 0.5*(xmax(zzz < 0) - xmin(zzz < 0));
    st.upp(zzz > 0) = x(zzz > 0) + 0.5*(xmax(zzz > 0) - xmin(zzz > 0));
    % Asymptote-distance clamp (Svanberg asymin/asymax): without it the
    % asymptotes can collapse onto x, beta -> x becomes active and MMA stalls
    % ~1e-4 short of the optimum (observed on case (a)). Keeping each asymptote
    % within [asymin, asymax] box-widths of x restores full convergence.
    dw = xmax - xmin;
    st.low = max(st.low, x - 10*dw);   st.low = min(st.low, x - 0.01*dw);
    st.upp = min(st.upp, x + 10*dw);   st.upp = max(st.upp, x + 0.01*dw);
end
L = st.low; U = st.upp;
raa0 = 1e-5;
dx = xmax - xmin;
p0 = (U-x).^2 .* (1.001*max(df0,0) + 0.001*max(-df0,0) + raa0./dx);
q0 = (x-L).^2 .* (0.001*max(df0,0) + 1.001*max(-df0,0) + raa0./dx);
p1 = (U-x).^2 .* (1.001*max(dfdx,0) + 0.001*max(-dfdx,0) + raa0./dx);
q1 = (x-L).^2 .* (0.001*max(dfdx,0) + 1.001*max(-dfdx,0) + raa0./dx);

% box for the subproblem (Svanberg's inward shrink)
alpha = max(xmin, L + 0.1*(x - L));
beta  = min(xmax, U - 0.1*(U - x));

% constraint bound: chosen so that x = current design satisfies it with equality
b1 = sum(p1./(U-x) + q1./(x-L)) - fval;

con = @(y) sum(p1./(U - xsub(y,p0,q0,p1,q1,L,U,alpha,beta)) + ...
               q1./(xsub(y,p0,q0,p1,q1,L,U,alpha,beta) - L)) - b1;
if con(0) <= 0
    y = 0;
else
    hi = 1;
    for k = 1:60
        if con(hi) <= 0 || hi > 1e12, break; end
        hi = hi * 2;
    end
    lo = 0;
    for k = 1:80
        mid = 0.5*(lo + hi);
        if con(mid) > 0, lo = mid; else, hi = mid; end
    end
    y = 0.5*(lo + hi);
end
xnew = xsub(y, p0, q0, p1, q1, L, U, alpha, beta);
st.xold2 = st.xold1; st.xold1 = x;
end

function xv = xsub(y, p0, q0, p1, q1, L, U, alpha, beta)
% Closed-form inner minimiser of  p/(U-x) + q/(x-L):
%   dM/dx = p/(U-x)^2 - q/(x-L)^2 = 0  =>  sqrt(p)/(U-x) = sqrt(q)/(x-L)
%   =>  x = (L*sqrt(p) + U*sqrt(q)) / (sqrt(p) + sqrt(q)).
% NOTE (bug fix): the original code had this as (L*sqrt(B)+U*sqrt(A)) with
% A=p, B=q, i.e. p and q SWAPPED -- it returned the FAR end of the box, so
% MMA stepped the wrong way and stalled (e.g. case (a) drove to (0,0)).
A = p0 + y*p1;  B = q0 + y*q1;
xv = (L.*sqrt(A) + U.*sqrt(B)) ./ (sqrt(A) + sqrt(B));
xv = min(max(xv, alpha), beta);
end
