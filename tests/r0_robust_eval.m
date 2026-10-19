function [J, dJ, info] = r0_robust_eval(ctx, t, op, R)
% R0_ROBUST_EVAL  Robust (sample-average) radiated-power objective under random relative
% thickness errors, with its adjoint gradient.  Experiment 10.
%
%   h_i = t .* (1 + sigma * xi_i),   xi_i = R.Xi(:,i) ~ N(0, I)  (independent per element)
%   P_i = w^2 || F u(h_i) ||^2,      u(h) solves (K_D(h) [+ iwZ]) u = F on the free dofs
%   J   = mean_i P_i                              (R.mode = 'mean', sample-average approximation)
%   J   = mean_i P_i + k std_i P_i                (R.mode = 'meanstd', R.k = k; std with 1/(N-1))
%   dJ/dt_e = sum_i wt_i (1 + sigma xi_ie) dP_i/dh_e,  wt_i = 1/N [+ k (P_i - mean)/((N-1) std)]
%
% Per sample the state solve, J = w^2||factor U||^2, the adjoint K_D z = a and the element
% sensitivities are those of r0_pilot_eval (= articleB_r0_eval(...,'quadvel',op,'consistent')),
% evaluated at the perturbed thickness h_i; the chain rule through h_i = t (1 + sigma xi_i) gives
% the factor (1 + sigma xi_i).  With sigma = 0 (or an empty Xi) and solver 'direct' it is one
% nominal sample and reproduces r0_pilot_eval to round-off (checked in r0_robust_unit).
%
% op    : struct('factor', F) [, 'iwZ', dense i w s Z_plus on the free dofs: fluid loading in the
%         state solve, as tests/r0_fsi_eval.m; Z is thickness-independent and symmetric, so the
%         adjoint still uses the same matrix]
% R     : struct with fields sigma, Xi (ne x N), mode ('mean'), k (0), solver ('gmres' | 'direct'),
%         tol (1e-10, left-preconditioned relative residual), extra (cell of additional factors
%         evaluated on the same sample states, returned in info.Pext, N x numel(extra); no gradient)
% Solver 'gmres' (default): one sparse LU of the NOMINAL K_D(t) (no fluid loading) preconditions
% every perturbed system (relative perturbation O(3 sigma) of the bending term, O(1e-3) for the
% fluid term); the initial guess is the nominal solution.  A GMRES whose residual stays above
% 10*tol after 60 iterations falls back to a direct solve (counted in info.nfallback).
% 'direct' does one LU per sample (sparse, or dense with iwZ).
% Cache: one entry keyed on (t, sigma, Xi, solver, fluid loading): the sample matrices and states of
% the last call are kept, so the gradient call at GCMMA's accepted point needs only the adjoints.
% Factor products are batched over the samples (one matrix-matrix product per factor).
persistent C PC
t = t(:);
if ~isfield(R,'mode') || isempty(R.mode), R.mode = 'mean'; end
if ~isfield(R,'k') || isempty(R.k), R.k = 0; end
if ~isfield(R,'solver') || isempty(R.solver), R.solver = 'gmres'; end
if ~isfield(R,'tol') || isempty(R.tol), R.tol = 1e-10; end
if ~isfield(R,'extra'), R.extra = {}; end
fsi = isfield(op,'iwZ') && ~isempty(op.iwZ);
if R.sigma == 0 || isempty(R.Xi), Xi = zeros(numel(t),1); sig = 0; else, Xi = R.Xi; sig = R.sigma; end
N = size(Xi,2);
ndof = ctx.ndof; fr = ctx.free; nf = numel(fr); ne = ctx.ne; lm = ctx.lm; n = size(lm,2);
cK = 1 + 1i*ctx.w*ctx.beta;
cM = 1i*ctx.w*ctx.alpha - ctx.w^2;
w2 = ctx.w^2;
pkey = sprintf('%d_%d_%.17g_%.17g', ndof, nf, ctx.w, ctx.beta);
if isempty(PC) || ~strcmp(PC.key, pkey)
    map = zeros(ndof,1); map(fr) = 1:nf;
    keep = map(ctx.Rows(:)) > 0 & map(ctx.Cols(:)) > 0;
    PC = struct('key',pkey,'keep',keep,'rf',map(ctx.Rows(keep)),'cf',map(ctx.Cols(keep)), ...
        'C0',cK*ctx.Ske + cM*ctx.Me);
end
useg = strcmp(R.solver,'gmres') && (sig > 0 || fsi);
H = t .* (1 + sig*Xi);                       % ne x N perturbed thickness fields
key = {t, sig, size(Xi), R.solver, fsi, nf, R.tol};
hit = ~isempty(C) && isequal(C.key, key) && isequal(C.Xi, Xi);
nfb = 0;
if ~hit
    C = struct('key',{key},'Xi',Xi,'U',[],'A',{cell(N,1)},'L0',[],'nit',zeros(N,1));
    if useg
        A0 = assemble_free(ctx, PC, t, cK, cM);
        [L0,U0,P0,Q0] = lu(A0);
        C.L0 = struct('L',L0,'U',U0,'P',P0,'Q',Q0);
        b = ctx.F(fr); u0 = Q0*(U0\(L0\(P0*b)));
    end
    Us = complex(zeros(nf,N));
    for i = 1:N
        A = assemble_free(ctx, PC, H(:,i), cK, cM);
        if useg
            [Us(:,i), C.nit(i), fb] = gsolve(A, op, fsi, ctx.F(fr), C.L0, u0, R.tol); nfb = nfb + fb;
        else
            Us(:,i) = dsolve(A, op, fsi, ctx.F(fr));
        end
        C.A{i} = A;
    end
    C.U = Us;
end
Ufull = complex(zeros(ndof,N)); Ufull(fr,:) = C.U;
FU = op.factor*Ufull;                        % batched over samples
P = w2*real(sum(conj(FU).*FU,1)).';
m = mean(P);
if N > 1, s = std(P); else, s = 0; end
switch R.mode
    case 'mean',    J = m;
    case 'meanstd', J = m + R.k*s;
    otherwise, error('unknown mode %s', R.mode);
end
info = struct('P',P,'mean',m,'std',s,'N',N,'gmres_it',C.nit,'nfallback',nfb,'cachehit',hit);
if ~isempty(R.extra)
    info.Pext = zeros(N,numel(R.extra));
    for j = 1:numel(R.extra)
        V = R.extra{j}*Ufull; info.Pext(:,j) = w2*real(sum(conj(V).*V,1)).';
    end
end
if nargout < 2, return; end
wt = ones(N,1)/N;
if strcmp(R.mode,'meanstd') && N > 1 && s > 0
    wt = wt + R.k*(P - m)/((N-1)*s);
end
Aall = 2*w2*conj(op.factor'*FU);             % adjoint right-hand sides, batched
dJ = zeros(ne,1);
if useg && isempty(C.L0)
    A0 = assemble_free(ctx, PC, t, cK, cM);
    [L0,U0,P0,Q0] = lu(A0); C.L0 = struct('L',L0,'U',U0,'P',P0,'Q',Q0);
end
for i = 1:N
    h = H(:,i);
    a = Aall(fr,i);
    if useg
        L0 = C.L0; z0 = L0.Q*(L0.U\(L0.L\(L0.P*a)));
        [zf, ~, fb] = gsolve(C.A{i}, op, fsi, a, L0, z0, R.tol); nfb = nfb + fb;
    else
        zf = dsolve(C.A{i}, op, fsi, a);
    end
    z = zeros(ndof,1); z(fr) = zf;
    ui = Ufull(:,i);
    ZU = reshape(ui(lm),ne,n,1) .* reshape(z(lm),ne,1,n);
    ZU = reshape(ZU,ne,n*n);
    q = 3*cK*(h.^2).*sum(ZU.*ctx.Bke,2) + sum(ZU.*PC.C0,2);
    dJ = dJ + wt(i) * (1 + sig*Xi(:,i)) .* (-real(q));
end
info.nfallback = nfb;
end

function A = assemble_free(ctx, PC, h, cK, cM)
EmatK = ctx.Bke.*(h.^3) + ctx.Ske.*h;
EmatM = ctx.Me.*h;
v = cK*EmatK(:) + cM*EmatM(:);
nf = numel(ctx.free);
A = sparse(PC.rf, PC.cf, v(PC.keep), nf, nf);
end

function x = dsolve(A, op, fsi, b)
if fsi
    [L,U,P] = lu(full(A) + op.iwZ); x = U\(L\(P*b));
else
    [L,U,P,Q] = lu(A); x = Q*(U\(L\(P*b)));
end
end

function [x, it, fb] = gsolve(A, op, fsi, b, L0, x0, tol)
if fsi, afun = @(v) A*v + op.iwZ*v; else, afun = @(v) A*v; end
mfun = @(v) L0.Q*(L0.U\(L0.L\(L0.P*v)));
[x, it, relres] = pgmres(afun, mfun, b, x0, tol, 60);
fb = 0;
if relres > 10*tol
    x = dsolve(A, op, fsi, b); fb = 1;
end
end

function [x, k, relres] = pgmres(afun, mfun, b, x0, tol, maxit)
% left-preconditioned GMRES (modified Gram-Schmidt Arnoldi, no restart); relres is
% ||M\(b - A x)|| / ||M\b||, the quantity MATLAB's gmres reports for a left preconditioner
r = mfun(b - afun(x0));
nb = norm(mfun(b));
beta = norm(r); relres = beta/nb; k = 0;
if relres <= tol, x = x0; return; end
nn = numel(b);
V = complex(zeros(nn,maxit+1)); Hm = zeros(maxit+1,maxit);
V(:,1) = r/beta;
e1 = zeros(maxit+1,1); e1(1) = beta;
y = [];
for k = 1:maxit
    w = mfun(afun(V(:,k)));
    for j = 1:k
        hj = V(:,j)'*w; w = w - hj*V(:,j); Hm(j,k) = hj;
    end
    Hm(k+1,k) = norm(w);
    if Hm(k+1,k) > 0, V(:,k+1) = w/Hm(k+1,k); end
    Hk = Hm(1:k+1,1:k);
    y = Hk\e1(1:k+1);
    relres = norm(Hk*y - e1(1:k+1))/nb;
    if relres <= tol || Hm(k+1,k) == 0, break; end
end
x = x0 + V(:,1:k)*y;
end
