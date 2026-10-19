function [x, flag, relres, iter, resvec] = cocg(Afun, b, Mfun, tol, maxit, x0)
% COCG  Conjugate Orthogonal Conjugate Gradient for COMPLEX-SYMMETRIC systems.
%
%   A x = b  with  A.' = A  (transpose, NOT conjugate transpose)
%
% K_D = (t^3 Kb + t Km) + j*w*(alpha*M + beta*(t^3 Kb + t Km)) - w^2 M is
% complex symmetric, so the Krylov process must use the BILINEAR inner product
%     <u,v> = u.' * v        (no conjugation)
% Using the Hermitian inner product instead (as standard PCG does) is wrong for
% this operator and destroys the short recurrence.
%
% The preconditioner Mfun (optional) must also be complex symmetric; the Jacobi
% choice is Mfun = @(r) r ./ d with d = diag(A) from harmonic_kd_diag.
%
% Inputs:  Afun   function handle y = Afun(x)
%          b      right-hand side (N x 1)
%          Mfun   preconditioner handle, or [] / omitted for none
%          tol    relative residual tolerance (default 1e-8)
%          maxit  iteration cap (default 500)
%          x0     initial guess (default zero)
% Outputs: x, flag (0 converged, 1 maxit, 2 breakdown), relres, iter, resvec
%
% Breakdown: for complex symmetric A the bilinear denominator <p,Ap> can vanish
% even when the system is solvable; this is detected and reported rather than
% producing NaNs.

if nargin < 3, Mfun = []; end
if nargin < 4 || isempty(tol), tol = 1e-8; end
if nargin < 5 || isempty(maxit), maxit = 500; end
if nargin < 6 || isempty(x0), x0 = zeros(size(b)); end

b = b(:); x = x0(:);
nb = norm(b);
if nb == 0, flag = 0; relres = 0; iter = 0; resvec = 0; return; end

r = b - Afun(x);
resvec = zeros(maxit+1, 1);
resvec(1) = norm(r)/nb;
if resvec(1) <= tol, flag = 0; relres = resvec(1); iter = 0; return; end

if isempty(Mfun), z = r; else, z = Mfun(r); end
p = z;
rho = r.' * z;                       % BILINEAR (no conjugation)

flag = 1; relres = resvec(1); iter = maxit;
small = 1e-300;
for k = 1:maxit
    Ap = Afun(p);
    den = p.' * Ap;                  % bilinear denominator
    if abs(den) < small * norm(p) * norm(Ap)
        flag = 2; iter = k-1; relres = norm(r)/nb; break;
    end
    alp = rho / den;
    x = x + alp * p;
    r = r - alp * Ap;
    rr = norm(r)/nb;
    resvec(k+1) = rr;
    if rr <= tol
        flag = 0; iter = k; relres = rr; break;
    end
    if isempty(Mfun), z = r; else, z = Mfun(r); end
    rho_new = r.' * z;
    if abs(rho) < small * abs(rho_new), flag = 2; iter = k; relres = rr; break; end
    bet = rho_new / rho;
    p = z + bet * p;
    rho = rho_new;
    relres = rr;
end
resvec = resvec(1:min(iter+1, numel(resvec)));
end
