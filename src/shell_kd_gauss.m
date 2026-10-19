function [Kb, info] = shell_kd_gauss(msh, space, geometry, YOUNG, POISSON, bend_model, nqg)
% SHELL_KD_GAUSS  Per-element GAUSS assembly of the metric-carrying bending
% operator for a THIN SHELL: a 2D manifold embedded in R^3 (rdim = 3).
%
% This is the genuine shell counterpart of harmonic_kd_gauss_curved.m (the
% rdim=2 curved-chart Gauss remedy) and the CONSISTENT scheme for the shell:
% instead of the reduced weighted-quadrature rule of quadrule_bend_fast /
% quadrule_bend_shell (exact only for a CONSTANT coefficient, see the moment
% theorem), every element is integrated with a tensor-product Gauss rule and the
% spatially VARYING surface coefficient
%
%     M(xi,eta) = D * sqrt(a) [ nu t' t + (1-nu) Q ]        ('plate')
%     M(xi,eta) = D * sqrt(a) [ Q ]                          ('T1')
%
% is evaluated at every Gauss point of every element from the genuine surface
% kinematics (surface metric a, sqrt(a), second fundamental form b, Christoffel
% symbols) exactly as in quadrule_bend_shell.m / shell_curvature.m. The resulting
% quadrature is consistent: the polynomial part of the integrand is integrated
% exactly and the curved-geometry error converges under h-refinement.
%
% OPERATOR
%   Kb = int D sqrt(a) Phi' G_shell Phi dxi deta ,
%   Phi = [w_xixi; w_xieta; w_etaeta; w_xi; w_eta]   (parametric components)
% with the 5x5 surface coefficient G_shell built exactly as in
% quadrule_bend_shell.m, so Kb here is the SAME operator the WQ path assembles
% (only the quadrature differs) -- the three schemes are directly comparable.
%
% The 1D shape functions / gradients / hessians are taken from the ANALYSIS
% SPACE itself (sp_scalar.sp_univ) evaluated on a Gauss mesh, so rational NURBS
% geometry and rational spaces are handled correctly (unlike the B-spline-only
% sp_bspline_1d_param path of the rdim=2 helper).
%
% INPUTS
%   msh, space, geometry : geopdes 2D mesh / sp_vector space / rdim=3 geometry.
%   YOUNG, POISSON       : material constants.
%   bend_model           : 'plate' (default) or 'T1'.
%   nqg                  : Gauss points per element per direction (default p+1).
%
% OUTPUT
%   Kb   : sparse Ns x Ns real bending operator on the SCALAR space (bending
%          ONLY), symmetrized. Ns = prod(space.scalar_spaces{1}.ndof_dir).
%   info : .nqg, .D, .bend_model, .nel, .N, .npts_shell
%
% NOTE: new file; no existing src/ file is modified.

if nargin < 6 || isempty(bend_model), bend_model = 'plate'; end

sp_s = space.scalar_spaces{1};
n    = sp_s.ndof_dir;              % [n1 n2]
N    = prod(n);

deg = cell(1,2);
for d = 1:2, deg{d} = sp_s.sp_univ(d).degree; end
if nargin < 7 || isempty(nqg), nqg = max(deg{1},deg{2}) + 1; end

D = YOUNG / (12 * (1 - POISSON^2));
switch lower(bend_model)
    case 't1',    nu = 0;
    case 'plate', nu = POISSON;
    otherwise, error('shell_kd_gauss: unknown bend_model ''%s''.', bend_model);
end

%% 1. Gauss mesh and 1D collocations from the analysis space
knots_u = geometry.nurbs.knots;
rule = msh_gauss_nodes([nqg nqg]);
[qn_g, qw_g] = msh_set_quad_nodes(knots_u, rule);
mshg = msh_cartesian(knots_u, qn_g, qw_g, geometry);
spg  = sp_nurbs(geometry.nurbs, mshg);

B0 = cell(1,2); B1 = cell(1,2); B2 = cell(1,2); conn = cell(1,2);
for d = 1:2
    sug = spg.sp_univ(d);
    B0{d}  = sug.shape_functions;            % [nqg x nsh x nel]
    B1{d}  = sug.shape_function_gradients;
    B2{d}  = sug.shape_function_hessians;
    conn{d} = sug.connectivity;              % [nsh x nel]
end

%% 2. Shell coefficient M (5x5) at every Gauss point, one shell_curvature call
nel1 = size(conn{1},2); nel2 = size(conn{2},2);
qn1 = reshape(qn_g{1}, nqg, nel1);   % [nqg x nel1]
qn2 = reshape(qn_g{2}, nqg, nel2);
[G1g, G2g] = ndgrid(1:nqg, 1:nqg);
XiE  = reshape(qn1, [], 1);                  % dir1-fast ordering of the first column
% full point list: index i1 (dir1, over nqg*nel1) fastest, then i2
xiall = qn1(:);  etall = qn2(:);             % each already dir-fast over its own elements
[I1all, I2all] = ndgrid(1:numel(xiall), 1:numel(etall));
P_all = [xiall(I1all(:)).'; etall(I2all(:)).'];
npt1 = numel(xiall);
M_all = shell_M(geometry, P_all, D, nu);     % (npt1*npt2) x 25

%% 3. Element assembly
comps = [2 0; 1 1; 0 2; 1 0; 0 1];   % Phi components
rows = []; cols = []; vals = [];
for el2 = 1:nel2
    loc2 = conn{2}(:,el2);
    for el1 = 1:nel1
        loc1 = conn{1}(:,el1);
        gdof = loc1(:) + (loc2(:).' - 1) * n(1);   % nsh1*nsh2 (dir1 fast)
        nloc = numel(gdof);
        Ke = zeros(nloc, nloc);
        for g2 = 1:nqg
            eta = qn2(g2,el2); w2 = qw_g{2}(g2,el2);
            i2g = (el2-1)*nqg + g2;
            v0b = full(B0{2}(g2, :, el2)); v1b = full(B1{2}(g2, :, el2)); v2b = full(B2{2}(g2, :, el2));
            v0b = v0b(:); v1b = v1b(:); v2b = v2b(:);
            for g1 = 1:nqg
                w1 = qw_g{1}(g1,el1);
                i1g = (el1-1)*nqg + g1;
                pid = i1g + (i2g-1)*npt1;
                Mm  = reshape(M_all(pid,:), 5, 5);
                v0a = full(B0{1}(g1, :, el1)); v1a = full(B1{1}(g1, :, el1)); v2a = full(B2{1}(g1, :, el1));
                v0a = v0a(:); v1a = v1a(:); v2a = v2a(:);
                Phi = [ kron(v0b, v2a), kron(v1b, v1a), kron(v2b, v0a), ...
                        kron(v0b, v1a), kron(v1b, v0a) ];   % nloc x 5
                Ke = Ke + (w1*w2) * (Phi * Mm * Phi.');
            end
        end
        [ii, jj] = ndgrid(1:nloc, 1:nloc);
        rows = [rows; gdof(ii(:))]; %#ok<AGROW>
        cols = [cols; gdof(jj(:))]; %#ok<AGROW>
        vals = [vals; reshape(Ke, [], 1)]; %#ok<AGROW>
    end
end
Kb = sparse(rows, cols, vals, N, N);
Kb = 0.5 * (Kb + Kb.');

info = struct('nqg', nqg, 'D', D, 'bend_model', lower(bend_model), ...
              'nel', [nel1 nel2], 'N', N, 'npts_shell', numel(P_all));
end

% ---------------------------------------------------------------------------
function M = shell_M(geometry, P, D, nu)
% Shell bending coefficient page [npts x 25] at parametric points P (2 x npts),
% identical in form to quadrule_bend_shell.m (surface metric / curvature).
Sk = shell_curvature(geometry, P, 'analytic');
npts = Sk.npts;
comps = [2 0; 1 1; 0 2; 1 0; 0 1]; %#ok<NASGU>
M = zeros(npts, 25);
for p = 1:npts
    ai = Sk.ainv(:,:,p);  Sq = Sk.sqrta(p);  Gm = Sk.Gamma(:,:,:,p);
    A11 = [1 0 0  -Gm(1,1,1)  -Gm(1,1,2)];
    A12 = [0 1 0  -Gm(1,2,1)  -Gm(1,2,2)];
    A22 = [0 0 1  -Gm(2,2,1)  -Gm(2,2,2)];
    Aof = cell(2,2); Aof{1,1}=A11; Aof{1,2}=A12; Aof{2,1}=A12; Aof{2,2}=A22;
    t = ai(1,1)*A11 + ai(1,2)*A12 + ai(2,1)*A12 + ai(2,2)*A22;
    Q = zeros(5);
    for al = 1:2
        for be = 1:2
            for g = 1:2
                for de = 1:2
                    Q = Q + ai(al,g) * ai(be,de) * (Aof{al,be}.' * Aof{g,de});
                end
            end
        end
    end
    Mm = D * Sq * (nu * (t.'*t) + (1-nu) * Q);
    M(p,:) = reshape(Mm, 1, 25);
end
end
