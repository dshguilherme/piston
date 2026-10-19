% MICROBENCH Confirm the diag() full-matrix allocation is the subsolv hotspot.
n = 2048; m = 1;
rng(0);
GG = rand(m,n); diagxinv = rand(n,1); diaglamyi = rand(m,1); a = 1;
reps = 2000;
tic
for k=1:reps
    Alam = diag(diaglamyi) + GG*diag(diagxinv)*GG.';
end
t_full = toc;
tic
for k=1:reps
    Alam2 = diag(diaglamyi) + (GG.*diagxinv.')*GG.';
end
t_scal = toc;
fprintf('diag() form   : %.4f ms/iter  (x%d)\n', 1000*t_full/reps, reps);
fprintf('scaling form  : %.4f ms/iter\n', 1000*t_scal/reps);
fprintf('identical? max abs diff = %.3e\n', max(abs(Alam(:)-Alam2(:))));
fprintf('full 2048x2048 diag bytes = %.1f MB\n', n*n*8/1e6);