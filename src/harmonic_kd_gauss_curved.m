function [Kb, info] = harmonic_kd_gauss_curved(msh, space, geometry, YOUNG, POISSON, bend_model, nqg)
% HARMONIC_KD_GAUSS_CURVED  Per-element GAUSS assembly of the curved
% (metric-carrying, 5-component) Kirchhoff-Love bending operator
%
%     Kb = int D det(J) Phi' G Phi dx ,
%     Phi = [w_xi_xi ; w_xi_eta ; w_eta_eta ; w_xi ; w_eta]  (parametric derivs)
%
% with the EXACT second-order chain rule of quadrule_bend_curved.m (Christoffel
% / inverse-map second-derivative term S included).  It is the CONSISTENT
% REMEDY for the weighted-quadrature failure: instead of the reduced rule of
% quadrule_bend_fast (exact only for a CONSTANT coefficient, see the moment
% theorem), every element is integrated with a genuine tensor-product Gauss
% rule and the spatially VARYING coefficient C = D det(J) B' G B is evaluated
% at the Gauss points.  The resulting quadrature is consistent: for degree-p
% splines and nqg Gauss points per direction the polynomial part of the
% integrand is integrated exactly (affine maps) and the curved-map error
% converges under h-refinement.
%
% INPUTS
%   msh, space, geometry : geopdes 2D mesh / scalar space / geometry.
%   YOUNG, POISSON       : material constants.
%   bend_model           : 'T1' (default) or 'plate'.
%   nqg                  : Gauss points per element per direction (default p+1).
%
% OUTPUT
%   Kb   : sparse NxN real bending operator (bending ONLY), symmetrized.
%   info : .nqg, .D, .bend_model, .nel, .N
%
% NOTE: no existing src/ file is modified; this is a new file.

if nargin < 6 || isempty(bend_model), bend_model = 'T1'; end

sp_s = space.scalar_spaces{1};
n    = sp_s.ndof_dir;              % [n1 n2]
N    = prod(n);

% degree per direction
deg = cell(1,2);
for d = 1:2, deg{d} = sp_s.sp_univ(d).degree; end
if nargin < 7 || isempty(nqg), nqg = max(deg{1},deg{2}) + 1; end

D = YOUNG / (12 * (1 - POISSON^2));
switch lower(bend_model)
    case 't1',    G = diag([1 1 2]);
    case 'plate', G = [1 POISSON 0; POISSON 1 0; 0 0 2*(1-POISSON)];
    otherwise, error('harmonic_kd_gauss_curved: unknown bend_model ''%s''.', bend_model);
end

%% 1. Per-direction Gauss grid and 1D collocations (values / 1st / 2nd derivs)
B0 = cell(1,2); B1 = cell(1,2); B2 = cell(1,2);
qd = cell(1,2); qw = cell(1,2); conn = cell(1,2);
for d = 1:2
    kv    = sp_s.sp_univ(d).knots;
    knots = unique(kv);
    rule  = msh_gauss_nodes(nqg);
    [qd{d}, qw{d}] = msh_set_quad_nodes(knots, rule);
    s0 = sp_bspline_1d_param(kv, deg{d}, qd{d});
    s1 = sp_bspline_1d_param(kv, deg{d}, qd{d}, 'gradient', true);
    s2 = sp_bspline_1d_param(kv, deg{d}, qd{d}, 'hessian',  true);
    B0{d} = colloc(s0.shape_functions(:),         s0.connectivity, nqg, s0.nsh_max, s0.ndof);
    B1{d} = colloc(s1.shape_function_gradients(:),s1.connectivity, nqg, s1.nsh_max, s1.ndof);
    B2{d} = colloc(s2.shape_function_hessians(:), s2.connectivity, nqg, s2.nsh_max, s2.ndof);
    conn{d} = s0.connectivity;      % [nsh x nel]
end

comps = [2 0; 1 1; 0 2; 1 0; 0 1];   % the five Phi components

%% 2. Element assembly
nel1 = size(conn{1},2); nel2 = size(conn{2},2);
rows = []; cols = []; vals = [];
for el2 = 1:nel2
    loc2 = conn{2}(:,el2);      % nsh2 dofs in direction 2
    g2idx = (el2-1)*nqg + (1:nqg);
    for el1 = 1:nel1
        loc1 = conn{1}(:,el1);  % nsh1 dofs in direction 1
        g1idx = (el1-1)*nqg + (1:nqg);
        gdof = loc1(:) + (loc2(:).' - 1) * n(1);   % nsh1*nsh2 (dir1 fast)
        nloc = numel(gdof);
        Ke = zeros(nloc, nloc);
        for g2 = 1:nqg
            eta = qd{2}(g2,el2); w2 = qw{2}(g2,el2);
            v0b = full(B0{2}(g2idx(g2), loc2)); v0b = v0b(:);
            v1b = full(B1{2}(g2idx(g2), loc2)); v1b = v1b(:);
            v2b = full(B2{2}(g2idx(g2), loc2)); v2b = v2b(:);
            for g1 = 1:nqg
                xi = qd{1}(g1,el1); w1 = qw{1}(g1,el1);
                P = [xi; eta];
                [~, J] = geometry.map_der(P);   J = squeeze(J);
                [~, ~, H] = geometry.map_der2(P); H = squeeze(H);
                dJ = J(1,1)*J(2,2) - J(1,2)*J(2,1);
                Ji = inv(J);  S = inv_map_der2(Ji, H);
                b_xx = [Ji(1,1)^2; 2*Ji(1,1)*Ji(2,1); Ji(2,1)^2; S(1,1,1); S(2,1,1)];
                b_yy = [Ji(1,2)^2; 2*Ji(1,2)*Ji(2,2); Ji(2,2)^2; S(1,2,2); S(2,2,2)];
                b_xy = [Ji(1,1)*Ji(1,2); Ji(1,1)*Ji(2,2)+Ji(2,1)*Ji(1,2); Ji(2,1)*Ji(2,2); S(1,1,2); S(2,1,2)];
                Bm = [b_xx.'; b_yy.'; b_xy.'];      % 3x5
                M = D * dJ * (Bm.' * G * Bm);       % 5x5
                v0a = full(B0{1}(g1idx(g1), loc1)); v0a = v0a(:);
                v1a = full(B1{1}(g1idx(g1), loc1)); v1a = v1a(:);
                v2a = full(B2{1}(g1idx(g1), loc1)); v2a = v2a(:);
                % Phi(:,k) : component k = (i derivs in dir1, j in dir2).
                % gdof is indexed dir1-fast, so build Phi with dir1-fast kron.
                Phi = [ kron(v0b, v2a), kron(v1b, v1a), kron(v2b, v0a), ...
                        kron(v0b, v1a), kron(v1b, v0a) ];   % nloc x 5
                Ke = Ke + (w1*w2) * (Phi * M * Phi.');
            end
        end
        [ii, jj] = ndgrid(1:nloc, 1:nloc);
        rows = [rows; gdof(ii(:))];   %#ok<AGROW>
        cols = [cols; gdof(jj(:))];   %#ok<AGROW>
        vals = [vals; reshape(Ke, [], 1)];  %#ok<AGROW>
    end
end
Kb = sparse(rows, cols, vals, N, N);
Kb = 0.5 * (Kb + Kb.');

info = struct('nqg', nqg, 'D', D, 'bend_model', lower(bend_model), ...
              'nel', [nel1 nel2], 'N', N);
end

% ---------------------------------------------------------------------------
function Bs = colloc(vals, connectivity, nqn, nsh, ndof)
% Sparse [nqn*nel x ndof] collocation from a [nqn x nsh x nel] array.
nel = size(connectivity, 2);
rows = repmat(reshape(1:nqn*nel, [nqn, 1, nel]), [1, nsh, 1]);
cols = repmat(reshape(connectivity, [1, nsh, nel]),  [nqn, 1, 1]);
Bs = sparse(rows(:), cols(:), vals(:), nqn*nel, ndof);
end

% ---------------------------------------------------------------------------
function S = inv_map_der2(Ji, H)
S = zeros(2,2,2);
for alpha = 1:2
    for i = 1:2
        for j = 1:2
            s = 0;
            for m = 1:2
                for b = 1:2
                    for c = 1:2
                        s = s - Ji(alpha,m) * H(m,b,c) * Ji(b,i) * Ji(c,j);
                    end
                end
            end
            S(alpha,i,j) = s;
        end
    end
end
end
