function run_all_tests()
% RUN_ALL_TESTS Regression suite for the harmonic K_D module.
%   tests/test_wq2_exactness.m   - exactness of all new WQ rule types
%   tests/test_kd_vs_assembled.m - K_D (bending + mass) vs assembled reference
%   tests/test_kd_matvec.m       - matrix-free matvec vs the sparse form
%   tests/test_kd_pcg.m          - COCG correctness, Jacobi diag, band sweep
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..', 'src'));
fprintf('\n=== WQ rule exactness ===\n');
run(fullfile(here, 'test_wq2_exactness.m'));
fprintf('\n=== K_D vs assembled reference ===\n');
run(fullfile(here, 'test_kd_vs_assembled.m'));
fprintf('\n=== matrix-free matvec ===\n');
run(fullfile(here, 'test_kd_matvec.m'));
fprintf('\n=== PCG/COCG ===\n');
run(fullfile(here, 'test_kd_pcg.m'));
end
