function Q2 = quadrule_bend_cmoment(space, nqg, nr)
% QUADRULE_BEND_CMOMENT  C-weighted / moment-enriched weighted-quadrature rule.
%
% The reduced rule of quadrule_bend_fast is exact against the (p-k2)-spline
% space only, which is why a spatially VARYING coefficient C(x) breaks it (see
% the moment-condition theorem). This rule ENRICHES the moment conditions so the
% rule is exact for the leading (polynomial) variation of C: on the per-element
% GAUSS point set it enforces, for every test DOF a and every neighbour c of the
% (p-k2) space,
%     sum_q w_a(q) phi_c(x_q) x_q^m = int N_a^{(k1)}(x) phi_c(x) x^m dx ,
%     m = 0,1,...,nr-1 ,
% i.e. exactness against the coefficient-multiplied multiplier space
% phi_c * {1,x,...,x^{nr-1}}.  With nr=2 the rule is exact for C affine in the
% direction (a + b x), hence for the leading order of any smooth C; the residual
% defect for a smooth C then vanishes with h (instead of pinning at O(1)).
%
% Points: nqg Gauss nodes per element (union over the mesh, deduplicated). The
% number of points in the support must be >= nc*nr (else lsqminnorm returns the
% best fit and exactness of the extra moments is only partial near boundaries).
%
% Output fields mirror quadrule_bend_fast: all_points, ind_points, nquad_points,
% quad_points, quad_weights_{00,01,10,11,22,20,02,21,12}, B_global, B1_global,
% B2_global, D1, D2.

if nargin < 2 || isempty(nqg), nqg = space.degree + 5; end
if nargin < 3 || isempty(nr),  nr  = 2; end

knots = space.knots; ndof = space.ndof; degree = space.degree;
assert(degree >= 3, 'quadrule_bend_cmoment: requires degree >= 3.');
distinct = unique(knots); nel = numel(distinct) - 1;

% Gauss point set (union over elements)
rule = msh_gauss_nodes(nqg);
[qd, qw] = msh_set_quad_nodes(distinct, rule);        %#ok<ASGLU>
all_points = unique(qd(:).');
npts = numel(all_points);

% derived spaces / difference operators
kd1 = knots(2:end-1); p1 = degree-1; nd1 = ndof-1;
kd2 = knots(3:end-2); p2 = degree-2; nd2 = ndof-2;
c1 = degree    ./ (knots(degree+2 : ndof+degree)   - knots(2:ndof))';
c2 = (degree-1)./ (knots(degree+2 : ndof+degree-1) - knots(3:ndof))';
D1 = spdiags([-c1 c1], [0 1], ndof-1, ndof);
D2 = spdiags([-c2 c2], [0 1], ndof-2, ndof-1);

% collocations of the (test/trial) spaces at all_points (via difference chains)
iv = findspan(ndof-1, degree, all_points, knots);
Bv = basisfun(iv, all_points, degree, knots); nbg = numbasisfun(iv, all_points, degree, knots)+1;
B_global = spcol(npts, degree, Bv, nbg, ndof);
[iv1, Bv1, nbg1] = coll1d(nd1, p1, all_points, kd1);
[iv2, Bv2, nbg2] = coll1d(nd2, p2, all_points, kd2);
Bgd1 = spcol(npts, p1, Bv1, nbg1, nd1);
Bgd2 = spcol(npts, p2, Bv2, nbg2, nd2);
B1_global = Bgd1 * D1;
B2_global = Bgd2 * (D2 * D1);
Coll = {B_global, B1_global, B2_global};      % indexed by k2+1

% fine Gauss grid over ALL elements for the exact RHS moments
ngf = 2*degree + 2;
[xf, wf] = msh_set_quad_nodes(distinct, msh_gauss_nodes(ngf));
xf = xf(:).'; wf = wf(:).';
nxf = numel(xf);
ivf = findspan(ndof-1, degree, xf, knots); Bf = basisfun(ivf, xf, degree, knots); nbf = numbasisfun(ivf, xf, degree, knots)+1;
Bf0 = spcol(nxf, degree, Bf, nbf, ndof);
[~, Bf1, nbf1] = coll1d(nd1, p1, xf, kd1); Bfd1 = spcol(nxf, p1, Bf1, nbf1, nd1);
[~, Bf2, nbf2] = coll1d(nd2, p2, xf, kd2); Bfd2 = spcol(nxf, p2, Bf2, nbf2, nd2);
Btf  = {Bf0, Bfd1*D1, Bfd2*(D2*D1)};          % test side at fine points
Bcol = {Bf0, Bfd1*D1, Bfd2*(D2*D1)};          % trial side at fine points

conn = space.connectivity; supp = space.supp;
Q2 = struct('all_points', all_points, 'ind_points', {cell(ndof,1)}, ...
            'nquad_points', zeros(ndof,1), 'quad_points', {cell(ndof,1)});
types = [0 0;0 1;1 0;1 1;2 0;0 2;2 1;1 2;2 2];
for it = 1:size(types,1)
    Q2.(sprintf('quad_weights_%d%d', types(it,1), types(it,2))) = cell(ndof,1);
end

for a = 1:ndof
    els = supp{a};
    lo = distinct(els(1)); hi = distinct(els(end)+1);
    ip = find(all_points >= lo - 1e-12 & all_points <= hi + 1e-12);
    xq = all_points(ip);                      % row
    nbrs = unique(conn(:, els)).';
    nc = numel(nbrs);
    Q2.ind_points{a} = ip(:).'; Q2.nquad_points(a) = numel(ip);
    Q2.quad_points{a} = xq;
    for it = 1:size(types,1)
        k1 = types(it,1); k2 = types(it,2);
        colB = full(Coll{k2+1}(ip, nbrs));        % npoints x nc
        LB = []; rh = [];
        for m = 0:nr-1
            LB = [LB; (colB .* (xq(:).^m)).'];                       %#ok<AGROW>  nc x npoints
            rhs_m = full(Bcol{k2+1}(:, nbrs)).' * ...
                    (full(Btf{k1+1}(:, a)) .* (xf(:).^m) .* wf(:));    % nc x 1
            rh = [rh; rhs_m];                                          %#ok<AGROW>
        end
        w = lsqminnorm(LB, rh);
        Q2.(sprintf('quad_weights_%d%d', k1, k2)){a} = w(:).';
    end
end

Q2.B_global = B_global; Q2.B1_global = B1_global; Q2.B2_global = B2_global;
Q2.D1 = D1; Q2.D2 = D2;
end

% ---------------------------------------------------------------------------
function Bs = spcol(nrow, p, vals, nbidx, ndof)
rows = repmat((1:nrow).', 1, p+1);
Bs = sparse(rows(:), nbidx(:), vals(:), nrow, ndof);
end

function [iv, Bv, nbidx] = coll1d(nd, p, pts, kn)
iv = findspan(nd-1, p, pts, kn);
Bv = basisfun(iv, pts, p, kn);
nbidx = numbasisfun(iv, pts, p, kn) + 1;
end
