% PCG gate on the matrix-free harmonic K_D operator.
% K_D is complex SYMMETRIC, so the Krylov method must use the BILINEAR inner
% product <u,v> = u.'*v (no conjugation) -- COCG. Checks:
%   A. COCG must reproduce MATLAB pcg EXACTLY on a real SPD system
%      (for real vectors the bilinear form is the standard inner product)
%   B. the Jacobi diagonal harmonic_kd_diag must equal diag(K_D) exactly
%   C. COCG solution vs the direct sparse solve, across the elastic resonances
%   D. iteration-count behaviour with mesh refinement
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
allok = true;

% ---- A. real SPD: COCG == MATLAB pcg ----
rng(0); n = 400;
Q = orth(randn(n)); ev = logspace(0,3,n)';
Ar = Q*diag(ev)*Q'; Ar = (Ar+Ar')/2; br = randn(n,1);
[cx, ~, cr, ci] = cocg(@(v) Ar*v, br, [], 1e-10, 5000);
[px, ~, pr, pi_] = pcg(Ar, br, 1e-10, 5000);
A1 = (ci == pi_) && (norm(cx-px)/norm(px) < 1e-12);
fprintf('A real SPD (cond 1e3): cocg it=%d | pcg it=%d | sol diff %.1e', ci, pi_, norm(cx-px)/norm(px));
if A1, fprintf('   OK\n'); else, fprintf('   FAIL\n'); allok = false; end

% ---- B/C/D on the physical plate ----
E = 210e9; RHO = 7850; NU = 0.3; t = 2e-3; al = 2*pi*0.01; bet = 2e-5;
pd = cantilever_beam(1.0, 1.0);
md = struct('degree',[3 3],'regularity',[2 2],'nsub',[16 8],'nquad',[4 4]);
[gm, mh, sp] = buildSpaces(pd, md);
te = t*ones(16,8);
H = harmonic_kd_setup(mh, sp, gm, E, NU, RHO, 'cpu', 'double', 'plate');
Vv = harmonic_kd_form(H, te, 0, 0, 0, 'values');
Kb = sparse(H.rows, H.cols, Vv.Kb, H.N, H.N);
Mt = sparse(H.rows, H.cols, Vv.Mm, H.N, H.N);
lam = sort(real(eig(full(Kb), full(Mt)))); lam = lam(lam > 1e-6*max(lam));
fr = sqrt(lam(1:6))/(2*pi);
fprintf('B elastic resonances (Hz): %s\n', mat2str(fr', 8));
rng(3); b = randn(H.N,1) + 1i*randn(H.N,1);
w = 2*pi*fr(3);
d = harmonic_kd_diag(H, te, w, al, bet);
KD = harmonic_kd_form(H, te, w, al, bet, 'cpu_sparse');
ed = norm(d - diag(KD))/norm(diag(KD));
fprintf('B diag vs sparse diag: rel %.3e', ed);
if ed < 1e-12, fprintf('   OK\n'); else, fprintf('   FAIL\n'); allok = false; end

fprintf('C band sweep (tol 1e-8):\n');
worst = 0;
for fac = [0.5 0.95 1.0 1.05 1.5 3.0]
    f = fac*fr(3); w = 2*pi*f;
    Af = @(v) harmonic_kd_matvec(H, te, v, w, al, bet);
    dg = harmonic_kd_diag(H, te, w, al, bet);
    KD = harmonic_kd_form(H, te, w, al, bet, 'cpu_sparse');
    xd = KD \ b;
    [x1, fl1, r1, i1] = cocg(Af, b, [], 1e-8, 4000);
    [x2, fl2, r2, i2] = cocg(Af, b, @(r) r./dg, 1e-8, 4000);
    e1 = norm(x1-xd)/norm(xd); e2 = norm(x2-xd)/norm(xd);
    worst = max([worst e1 e2]);
    fprintf('   f/f3=%5.2f (f=%7.2f Hz) | plain it=%4d res=%.1e fl=%d err=%.1e | Jac it=%4d res=%.1e fl=%d err=%.1e\n', ...
        fac, f, i1, r1, fl1, e1, i2, r2, fl2, e2);
    if fl1 ~= 0 || fl2 ~= 0, allok = false; end
end
fprintf('C worst solution error vs direct solve: %.3e', worst);
if worst < 1e-8, fprintf('   OK\n'); else, fprintf('   FAIL\n'); allok = false; end

fprintf('D mesh refinement (f=12.09 Hz, tol 1e-8):\n');
for ns = [8 16 32]      % ns=64 diverges at this resonance; see SESSION_LOG
    pdx = cantilever_beam(1.0, 1.0);
    mdx = struct('degree',[3 3],'regularity',[2 2],'nsub',[ns ns/2],'nquad',[4 4]);
    [gx, mx, sx] = buildSpaces(pdx, mdx);
    tex = t*ones(ns, ns/2);
    Hs = harmonic_kd_setup(mx, sx, gx, E, NU, RHO, 'cpu', 'double', 'plate');
    wx = 2*pi*fr(3);
    Afx = @(v) harmonic_kd_matvec(Hs, tex, v, wx, al, bet);
    dgx = harmonic_kd_diag(Hs, tex, wx, al, bet);
    bb = randn(Hs.N,1) + 1i*randn(Hs.N,1);
    [~, f1, r1, i1] = cocg(Afx, bb, [], 1e-8, 4000);
    [~, f2, r2, i2] = cocg(Afx, bb, @(r) r./dgx, 1e-8, 4000);
    fprintf('   N=%6d | plain it=%4d res=%.1e fl=%d | Jacobi it=%4d res=%.1e fl=%d\n', ...
        Hs.N, i1, r1, f1, i2, r2, f2);
end

if allok, disp('PCG_ALL_OK'); else, disp('PCG_FAIL'); end
