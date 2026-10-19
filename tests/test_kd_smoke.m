addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
problem_data = cantilever_beam(1.0, 0.5);
method_data.degree = [3 3]; method_data.regularity = [2 2];
method_data.nsub = [8 4]; method_data.nquad = [4 4];
[geometry, msh, space] = buildSpaces(problem_data, method_data);
H = harmonic_kd_setup(msh, space, geometry, 1.0, 0.3, 1.0, 'cpu', 'double');
te = 1e-3 * ones(8,4);   % thin plate
w = 2*pi*100;
KD = harmonic_kd_form(H, te, w, 1e-5, 1e-7, 'cpu_sparse');
fprintf('KD size %d nnz %d symerr %.2e\n', size(KD,1), nnz(KD), norm(KD-KD.',1)/max(norm(KD,1),eps));
v = harmonic_kd_form(H, te, w, 1e-5, 1e-7, 'values');
Kv = sparse(H.rows, H.cols, v, H.N, H.N);
fprintf('values-vs-sparse rel diff %.2e\n', norm(full(KD(:)-Kv(:)),inf)/max(norm(full(KD(:)),inf),eps));
disp('FORM_OK');
