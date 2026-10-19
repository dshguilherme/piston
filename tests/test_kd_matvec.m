% Matrix-free matvec gate: y = K_D*x (pages) vs the sparse form.
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
problem_data = cantilever_beam(1.0, 1.0);
method_data.degree = [3 3]; method_data.regularity = [2 2];
method_data.nsub = [8 4]; method_data.nquad = [4 4];
[geometry, msh, space] = buildSpaces(problem_data, method_data);
te = 1e-3 * ones(8,4); w = 2*pi*100; alpha = 1e-5; beta = 1e-7;
rng(0);
for model = {'T1','plate'}
  H = harmonic_kd_setup(msh, space, geometry, 1.0, 0.3, 1.0, 'cpu', 'double', model{1});
  KD = harmonic_kd_form(H, te, w, alpha, beta, 'cpu_sparse');
  x = randn(H.N,1) + 1i*randn(H.N,1);
  y = harmonic_kd_matvec(H, te, x, w, alpha, beta);
  yref = KD * x;
  e = norm(y - yref)/norm(yref);
  fprintf('matvec %-5s : rel err %.3e', model{1}, e);
  if e < 1e-12, fprintf('   OK\n'); else, fprintf('   FAIL\n'); end
  % complex-SYMMETRIC (not Hermitian) bilinear pairing: x2.'*K*x == (K*x2).'*x
  x2 = randn(H.N,1) + 1i*randn(H.N,1);
  y2 = harmonic_kd_matvec(H, te, x2, w, alpha, beta);
  r1 = x2.'*y; r2 = (y2).'*x;
  fprintf('  bilinear pairing rel diff = %.3e\n', abs(r1-r2)/abs(r1));
end
disp('MATVEC_DONE');
