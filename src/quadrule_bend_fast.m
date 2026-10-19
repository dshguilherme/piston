function Q2 = quadrule_bend_fast(space)
% QUADRULE_BEND_FAST 1D weighted-quadrature rules of all types (k1,k2) with
% k1,k2 in {0,1,2}, built as exact clones of the official
% quadrule_stiff_fast / quadrule_gen construction (source of truth).
%
% OFFICIAL CONVENTION (quadrule_gen, lines 121-147 of quadrule_stiff_fast.m):
%   the rule of type (k1,k2) for test function a satisfies
%       sum_q w_q * phi_c(x_q) = int N_a^{(k1)}(x) * phi_c(x) dx
%   for every trial function phi_c of the (p-k2) space whose support touches
%   supp(N_a); the assembly PAGE stores  w_a(q) * N_b^{(k2)}(x_q)  with N_b the
%   ORIGINAL-space basis, so the contraction reproduces
%       sum_q w_a(q) N_b^{(k2)}(x_q) C(x_q) = int N_a^{(k1)} N_b^{(k2)} C dx.
%
% TYPES
%   official (re-exported from quadrule_stiff_fast):
%     (0,0) values x values          (0,1) values x first derivatives
%     (1,0) first derivs x values    (1,1) first derivs x first derivs
%   new (this file):
%     (2,2) N_a'' x N_b''            (2,0) N_a'' x N_b
%     (0,2) N_a   x N_b''            (2,1) N_a'' x N_b'
%     (1,2) N_a'  x N_b''
%   (2,1) and (1,2) are what the metric expansion of the physical Hessian needs
%   (cross terms w_xi xi * w_xi eta and w_xi eta * w_eta eta).
%
% All new rules are built on the ORIGINAL-space neighbour set, so every page
% shares the J / len / mask sparsity structure.
%
% Collocations (original space):
%   Gauss grid, p+1 points per element (the official derived-space trick):
%     B_gauss  values, B1_gauss first derivatives, B2_gauss second derivatives
%   Reduced points (all_points of the stiffness rule):
%     B_global values, B1_global first derivatives, B2_global second derivatives
% Derivatives come from the exact difference chains
%   N' = (p-1) values * D1 ,  N'' = (p-2) values * (D2*D1)
% (verified against basisfunder to 1e-13).
%
% Output Q2 fields: all_points, ind_points, nquad_points, quad_points, plus
% quad_weights_{00,01,10,11,22,20,02,21,12} and the collocations above.

knots  = space.knots;
ndof   = space.ndof;
degree = space.degree;
assert(degree >= 3, 'quadrule_bend_fast: requires degree >= 3 (p-2 space must be >= 1).');
distinct_knots = unique(knots);
nel = numel(distinct_knots) - 1;
nqn_g = degree + 1;                 % Gauss points per element (matches space)

%% 0. Official stiffness rule: points, ind_points, weights (source of truth)
Qs = quadrule_stiff_fast(space);
all_points = Qs.all_points;
npts = numel(all_points);

%% 1. Derived spaces on the SAME Gauss grid
% (p-1) space
knots_d1 = knots(2:end-1);
p1 = degree - 1;
nd1 = ndof - 1;
qn_d1 = msh_set_quad_nodes(distinct_knots, msh_gauss_nodes(nqn_g));
sp1 = sp_bspline_1d_param(knots_d1, p1, qn_d1);
% (p-2) space
knots_d2 = knots(3:end-2);
p2 = degree - 2;
nd2 = ndof - 2;
qn_d2 = msh_set_quad_nodes(distinct_knots, msh_gauss_nodes(nqn_g));
sp2 = sp_bspline_1d_param(knots_d2, p2, qn_d2);
assert(sp1.ndof == nd1 && sp2.ndof == nd2, 'quadrule_bend_fast: derived space dof mismatch.');

rows_g1 = repmat(reshape(1:nqn_g*nel, [nqn_g, 1, nel]), [1, sp1.nsh_max, 1]);
cols_g1 = repmat(reshape(sp1.connectivity, [1, sp1.nsh_max, nel]), [nqn_g, 1, 1]);
B_gauss_der1 = sparse(rows_g1(:), cols_g1(:), sp1.shape_functions(:), nqn_g*nel, nd1);

rows_g2 = repmat(reshape(1:nqn_g*nel, [nqn_g, 1, nel]), [1, sp2.nsh_max, 1]);
cols_g2 = repmat(reshape(sp2.connectivity, [1, sp2.nsh_max, nel]), [nqn_g, 1, 1]);
B_gauss_der2 = sparse(rows_g2(:), cols_g2(:), sp2.shape_functions(:), nqn_g*nel, nd2);

%% 2. Exact difference operators
c1 = degree     ./ (knots(degree+2 : ndof+degree)   - knots(2:ndof))';      % [ndof-1]
c2 = (degree-1) ./ (knots(degree+2 : ndof+degree-1) - knots(3:ndof))';      % [ndof-2]
D1 = spdiags([-c1 c1], [0 1], ndof-1, ndof);
D2 = spdiags([-c2 c2], [0 1], ndof-2, ndof-1);

%% 3. Collocations of the ORIGINAL space
% Gauss grid
ivg = findspan(ndof-1, degree, qn_d2(:)', knots);
Bgv = basisfun(ivg, qn_d2(:)', degree, knots);
nbg = numbasisfun(ivg, qn_d2(:)', degree, knots) + 1;
rowsg = repmat((1:nqn_g*nel).', 1, degree+1);
B_gauss  = sparse(rowsg(:), nbg(:), Bgv(:), nqn_g*nel, ndof);
B1_gauss = B_gauss_der1 * D1;                 % first derivatives
B2_gauss = B_gauss_der2 * (D2 * D1);          % second derivatives

% Reduced points (all_points)
iv = findspan(ndof-1, degree, all_points, knots);
Bv = basisfun(iv, all_points, degree, knots);
nbv = numbasisfun(iv, all_points, degree, knots) + 1;
rowsv = repmat((1:npts).', 1, degree+1);
B_global = sparse(rowsv(:), nbv(:), Bv(:), npts, ndof);

% derived-space values at the reduced points, then the chains
iv1 = findspan(nd1-1, p1, all_points, knots_d1);
Bv1 = basisfun(iv1, all_points, p1, knots_d1);
nbv1 = numbasisfun(iv1, all_points, p1, knots_d1) + 1;
rowsv1 = repmat((1:npts).', 1, p1+1);
B_global_der1 = sparse(rowsv1(:), nbv1(:), Bv1(:), npts, nd1);

iv2 = findspan(nd2-1, p2, all_points, knots_d2);
Bv2 = basisfun(iv2, all_points, p2, knots_d2);
nbv2 = numbasisfun(iv2, all_points, p2, knots_d2) + 1;
rowsv2 = repmat((1:npts).', 1, p2+1);
B_global_der2 = sparse(rowsv2(:), nbv2(:), Bv2(:), npts, nd2);

B1_global = B_global_der1 * D1;
B2_global = B_global_der2 * (D2 * D1);

%% 4. Per-test-function weights, all new types (quadrule_gen clone)
gauss_rule = msh_gauss_nodes(nqn_g);
supp = space.supp;
conn = space.connectivity;

Q2 = struct('all_points', all_points, 'ind_points', {Qs.ind_points}, ...
            'nquad_points', Qs.nquad_points, 'quad_points', {Qs.quad_points});
newtypes = {[2 2],[2 0],[0 2],[2 1],[1 2]};
for it = 1:numel(newtypes)
    Q2.(sprintf('quad_weights_%d%d', newtypes{it}(1), newtypes{it}(2))) = cell(ndof,1);
end

for ii = 1:ndof
    neighbors = unique(conn(:, supp{ii}));          % ORIGINAL-space neighbours
    iel = supp{ii}';
    [~, qw] = msh_set_quad_nodes(distinct_knots(iel(1):iel(end)+1), gauss_rule);
    w_gauss = qw(:).';
    local_points = ((iel(1)-1)*nqn_g + 1) : iel(end)*nqn_g;
    ip = Qs.ind_points{ii};

    % collocation matrices per trial derivative order k2
    Coll = {B_global, B1_global, B2_global};
    CollG = {B_gauss, B1_gauss, B2_gauss};
    % test arrays per test derivative order k1
    Test = {B_gauss, B1_gauss, B2_gauss};

    for it = 1:numel(newtypes)
        k1 = newtypes{it}(1); k2 = newtypes{it}(2);
        local_B   = full(Coll{k2+1}(ip, neighbors).');                % [nc x nq]
        local_rhs = full(CollG{k2+1}(local_points, neighbors)).' * ...
                    (full(Test{k1+1}(local_points, ii)).' .* w_gauss).';
        % Minimum-norm least-squares solution. The collocation matrix can be
        % rank-deficient for the derivative types (the set {N_b'} or {N_b''}
        % over the support need not be linearly independent), so a plain
        % QR backslash is not well posed: it returns an arbitrary member of the
        % solution set and warns about singularity. lsqminnorm returns the
        % minimum-norm solution of the consistent system, which keeps the
        % reduced weights bounded and the page values stable.
        % system: local_B * w = local_rhs, local_B is [nc x nq]
        Q2.(sprintf('quad_weights_%d%d', k1, k2)){ii} = ...
            lsqminnorm(local_B, local_rhs).';
    end
end

%% 5. Re-export official stiffness weights + expose collocations
Q2.quad_weights_00 = Qs.quad_weights_00;
Q2.quad_weights_01 = Qs.quad_weights_01;
Q2.quad_weights_10 = Qs.quad_weights_10;
Q2.quad_weights_11 = Qs.quad_weights_11;
Q2.B_global = B_global;   Q2.B1_global = B1_global;  Q2.B2_global = B2_global;
Q2.B_gauss  = B_gauss;    Q2.B1_gauss  = B1_gauss;   Q2.B2_gauss  = B2_gauss;
Q2.D1 = D1; Q2.D2 = D2;
end
