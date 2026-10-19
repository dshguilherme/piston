function [bendC, info] = quadrule_bend_curved(geometry, qn, YOUNG, POISSON, bend_model, mapmode, fdh)
% QUADRULE_BEND_CURVED Metric-carrying bending coefficient page for a GENERAL
% (non-affine / curved) map x(xi,eta), including the map's own second
% derivatives (the Christoffel / inverse-map second-derivative terms that the
% affine page in harmonic_kd_setup.m omits).
%
% THE PHYSICS. The Kirchhoff-Love plate (or shell, in a 2D chart) bending
% density is, in PHYSICAL Cartesian coordinates (x,y) = x(xi,eta),
%     e_bend = D * ( w_xx^2 + w_yy^2 + 2 nu w_xx w_yy + 2(1-nu) w_xy^2 )   ('plate')
%     e_bend = D * ( w_xx^2 + w_yy^2 + 2 w_xy^2 )                          ('T1')
%     D = YOUNG/(12*(1-POISSON^2))
% and the energy is int e_bend * det(J) dxi deta.
%
% THE DERIVATION (full second-order chain rule). Let xi^alpha = xi^alpha(x,y)
% be the inverse map, Ji(alpha,i) = d xi^alpha / d x^i (i=1 -> x, i=2 -> y),
% h = [w_xixi; w_xieta; w_etaeta], g = [w_xi; w_eta]. Then
%     w_{,ij} = w_{,alpha beta} Ji(alpha,i) Ji(beta,j) + w_{,alpha} S(alpha,i,j)
% where the SECOND term is the missing piece and
%     S(alpha,i,j) = d^2 xi^alpha / d x^i d x^j = - Ji(alpha,m) H(m,beta,gamma) Ji(beta,i) Ji(gamma,j)
% with H(m,beta,gamma) = d^2 x^m / d xi^beta d xi^gamma the map Hessian
% (verified against geopdes_inv_der2__, the construction geo_load uses;
%  max abs deviation 5.7e-14 over random J,H).
%
% Collecting the five parametric components Phi = [w_xixi; w_xieta; w_etaeta;
% w_xi; w_eta], the physical Hessian components are the 3x5 matrix B with rows
%     beta^{xx} = [Ji11^2 ; 2 Ji11 Ji21 ; Ji21^2 ; S(1,1,1) ; S(2,1,1)]
%     beta^{yy} = [Ji12^2 ; 2 Ji12 Ji22 ; Ji22^2 ; S(1,2,2) ; S(2,2,2)]
%     beta^{xy} = [Ji11 Ji12 ; Ji11 Ji22+Ji21 Ji12 ; Ji21 Ji22 ; S(1,1,2) ; S(2,1,2)]
% (the first three entries reproduce the affine vectors a,b,c of
%  harmonic_kd_setup.m exactly), and the 5x5 bending coefficient is
%     M = D * det(J) * B' * G * B ,
%     G = diag(1,1,2) ('T1')  or  [1 nu 0; nu 1 0; 0 0 2(1-nu)] ('plate').
% For an AFFINE map S=0 and the last two columns of B vanish: M collapses to
% the affine 3x3 block (regression-exact).
%
% COMPONENT MAP (used by the contraction): component 1 = xi xi = (2,0),
% 2 = xi eta = (1,1), 3 = eta eta = (0,2), 4 = xi = (1,0), 5 = eta = (0,1).
% For test component alpha=(i,j) and trial component beta=(k,l), direction 1
% uses WQ rule (i,k) and direction 2 uses rule (j,l) -- all (k1,k2) in {0,1,2}
% are already produced by quadrule_bend_fast.
%
% POINT ORDERING. The page is returned FLAT, index p = i1 + (i2-1)*nqtot(1),
% varying the direction-1 point index fastest. This is the ordering required by
% the LIN gather of harmonic_kd_setup(_curved). NOTE: the existing
% harmonic_kd_setup.m uses meshgrid(qn{1},qn{2}) + reshape, which orders the
% direction-2 index fastest; for a NON-affine map that mixes the two directions
% (harmless for the affine maps it was written for, since the Jacobian is then
% constant). This function therefore builds P_all with ndgrid (dir-1 fastest).
%
% INPUTS
%   geometry : geopdes geometry (needs .map, .map_der, and .map_der2 for
%              mapmode='analytic'; only .map for mapmode='fd').
%   qn       : {q1, q2}, the WQ point coordinates per direction (cell).
%   YOUNG, POISSON : material constants.
%   bend_model : 'T1' (default) or 'plate'.
%   mapmode  : 'analytic' (default) uses geometry.map_der/map_der2;
%              'fd' uses central finite differences of geometry.map.
%   fdh      : FD step for mapmode='fd' (default 1e-4).
%
% OUTPUTS
%   bendC (nqtot(1)*nqtot(2) x 25) flat coefficient page.
%   info  : .comps (5x2), .detJ (flat), .D, .bend_model.

if nargin < 5 || isempty(bend_model), bend_model = 'T1'; end
if nargin < 6 || isempty(mapmode),    mapmode    = 'analytic'; end
if nargin < 7 || isempty(fdh),        fdh        = 1e-4; end

nq1 = numel(qn{1}); nq2 = numel(qn{2});
npts = nq1 * nq2;

% P_all with the direction-1 point index varying fastest (ndgrid, NOT meshgrid)
[I1g, I2g] = ndgrid(1:nq1, 1:nq2);
P1v = qn{1}(I1g(:)); P2v = qn{2}(I2g(:));
P_all = [P1v(:).'; P2v(:).'];

switch lower(mapmode)
    case 'analytic'
        [~, Jall] = geometry.map_der(P_all);            % (2,2,npts)
        Jall = reshape(Jall, 2, 2, nq1, nq2);
        Hall = geometry.map_der2(P_all);                % (2,2,2,npts)
        Hall = reshape(Hall, 2, 2, 2, nq1, nq2);
    case 'fd'
        Jall = zeros(2,2,nq1,nq2);
        Hall = zeros(2,2,2,nq1,nq2);
        for i2 = 1:nq2
            for i1 = 1:nq1
                [Jk, Hk] = fd_map_derivatives(geometry, qn{1}(i1), qn{2}(i2), fdh);
                Jall(:,:,i1,i2) = Jk;
                Hall(:,:,:,i1,i2) = Hk;
            end
        end
    otherwise
        error('quadrule_bend_curved: unknown mapmode ''%s''.', mapmode);
end

D = YOUNG / (12 * (1 - POISSON^2));
switch lower(bend_model)
    case 't1',    G = diag([1 1 2]);
    case 'plate', G = [1 POISSON 0; POISSON 1 0; 0 0 2*(1-POISSON)];
    otherwise, error('quadrule_bend_curved: unknown bend_model ''%s''.', bend_model);
end

comps = [2 0; 1 1; 0 2; 1 0; 0 1];   % the five Phi components
bendC = zeros(npts, 25);
detJ  = zeros(npts, 1);

for i2 = 1:nq2
    for i1 = 1:nq1
        Jm = Jall(:,:,i1,i2);
        Hm = Hall(:,:,:,i1,i2);
        dJ = Jm(1,1)*Jm(2,2) - Jm(1,2)*Jm(2,1);
        Ji = inv(Jm);
        S  = inv_map_der2(Ji, Hm);
        b_xx = [Ji(1,1)^2; ...
                Ji(1,1)*Ji(2,1) + Ji(2,1)*Ji(1,1); ...
                Ji(2,1)^2; S(1,1,1); S(2,1,1)];
        b_yy = [Ji(1,2)^2; ...
                2*Ji(1,2)*Ji(2,2); ...
                Ji(2,2)^2; S(1,2,2); S(2,2,2)];
        b_xy = [Ji(1,1)*Ji(1,2); ...
                Ji(1,1)*Ji(2,2) + Ji(2,1)*Ji(1,2); ...
                Ji(2,1)*Ji(2,2); S(1,1,2); S(2,1,2)];
        Bm = [b_xx.'; b_yy.'; b_xy.'];          % 3x5
        Mm = D * dJ * (Bm.' * G * Bm);          % 5x5
        p_flat = i1 + (i2-1)*nq1;
        bendC(p_flat, :) = reshape(Mm, 1, 25);
        detJ(p_flat) = dJ;
    end
end

info = struct('comps', comps, 'detJ', detJ, 'D', D, ...
              'bend_model', lower(bend_model), 'mapmode', lower(mapmode));
end

% ---------------------------------------------------------------------------
function S = inv_map_der2(Ji, H)
% S(alpha,i,j) = - sum_{m,b,c} Ji(alpha,m) H(m,b,c) Ji(b,i) Ji(c,j)
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

% ---------------------------------------------------------------------------
function [J, H] = fd_map_derivatives(geometry, xi, eta, h)
% Central finite differences of the map, from x(xi,eta) alone.
x0  = geometry.map([xi; eta]);
xp  = geometry.map([xi+h; eta]);  xm = geometry.map([xi-h; eta]);
yp  = geometry.map([xi; eta+h]);  ym = geometry.map([xi; eta-h]);
xpp = geometry.map([xi+h; eta+h]); xpm = geometry.map([xi+h; eta-h]);
xmp = geometry.map([xi-h; eta+h]); xmm = geometry.map([xi-h; eta-h]);
J = [(xp - xm)/(2*h), (yp - ym)/(2*h)];              % columns: d/dxi, d/deta
H = zeros(2,2,2);
H(:,1,1) = (xp - 2*x0 + xm)/h^2;                     % d2/dxi2
H(:,2,2) = (yp - 2*x0 + ym)/h^2;                     % d2/deta2
H(:,1,2) = (xpp - xpm - xmp + xmm)/(4*h^2);          % d2/dxideta
H(:,2,1) = H(:,1,2);
end
