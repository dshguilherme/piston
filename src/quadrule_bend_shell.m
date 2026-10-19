function [bendC, info] = quadrule_bend_shell(geometry, qn, YOUNG, POISSON, bend_model, mapmode, fdh)
% QUADRULE_BEND_SHELL  Metric-carrying bending coefficient page for a THIN
% SHELL: a 2D manifold embedded in R^3 (rdim = 3) with surface metric and
% curvature. This is the genuine shell counterpart of quadrule_bend_curved.m
% (which is a planar chart, rdim = 2, physical Cartesian Hessian only).
%
% PHYSICS (isotropic Kirchhoff-Love bending, physical components on the surface).
% With the covariant second derivative of the scalar field w,
%   kappa_{alpha beta} = w_{;alpha beta} = w_{,alpha beta} - Gamma^g_{alpha beta} w_{,g} ,
% the bending energy density per unit AREA is
%   e_bend = D [ nu (a^{alpha beta} kappa_{alpha beta})^2
%               + (1-nu) a^{alpha gamma} a^{beta delta} kappa_{alpha beta} kappa_{gamma delta} ]
%          = D [ nu I1^2 + (1-nu) I2 ] ,   D = YOUNG/(12 (1-POISSON^2)) ,
% and the energy is  int e_bend dA = int e_bend sqrt(a) dxi deta .
%
% IDENTIFICATION WITH THE PLATE. On a flat surface with a Cartesian chart
% (a_{alpha beta}=delta, b=0, Gamma=0): I1 = w_xx + w_yy and
% I2 = w_xx^2 + w_yy^2 + 2 w_xy^2, so
%   e_bend = D [ w_xx^2 + w_yy^2 + 2 nu w_xx w_yy + 2(1-nu) w_xy^2 ]  ('plate'),
% i.e. EXACTLY the rdim=2 density used by quadrule_bend_curved.m. The bend_model
% 'T1' of that file is the nu = 0 member e_bend = D I2 (D unchanged).
%
% PARAMETRIC EXPANSION (the five components the WQ pages carry).
%   Phi = [w_xi_xi; w_xi_eta; w_eta_eta; w_xi; w_eta] ,  components
%   1: xi xi = (2,0)  2: xi eta = (1,1)  3: eta eta = (0,2)  4: xi = (1,0)  5: eta = (0,1).
% Since kappa_{alpha beta} is linear in Phi,  kappa_{alpha beta} = A_{alpha beta} Phi with
%   A11 = [1 0 0  -Gamma^1_{11}  -Gamma^2_{11}]
%   A12 = [0 1 0  -Gamma^1_{12}  -Gamma^2_{12}]
%   A22 = [0 0 1  -Gamma^1_{22}  -Gamma^2_{22}]
% we get I1 = t' Phi  with  t = a^{alpha beta} A_{alpha beta}  (1 x 5)  and
%        I2 = Phi' Q Phi  with  Q = sum_{al,be,g,d} a^{al g} a^{be d} A_{al be}' A_{g d}  (5 x 5).
% Hence the 5 x 5 coefficient page
%   M = D sqrt(a) [ nu t' t + (1-nu) Q ]   ('plate')
%   M = D sqrt(a) [ Q ]                    ('T1', nu=0 member)
% exactly mirroring the 5 x 5 page of quadrule_bend_curved.m, so the SAME
% assembly (harmonic_kd_form_curved.m) consumes it.
%
% POINT ORDERING. Returned FLAT, index p = i1 + (i2-1)*nqtot(1), direction-1
% point index fastest (ndgrid), matching the LIN gather of shell_kd_setup.m.
%
% INPUTS
%   geometry : geopdes geometry (needs .map_der/.map_der2 for 'analytic',
%              only .map for 'fd'). rdim should be 3.
%   qn       : {q1,q2}, WQ point coordinates per direction.
%   YOUNG, POISSON : material constants.
%   bend_model : 'plate' (default) or 'T1'.
%   mapmode  : 'analytic' (default) or 'fd'.
%   fdh      : FD step for mapmode='fd' (default 1e-5).
%
% OUTPUTS
%   bendC (nq1*nq2 x 25) flat coefficient page.
%   info  : .comps (5x2), .detJ (=sqrt(a) flat), .D, .bend_model, .mapmode.

if nargin < 5 || isempty(bend_model), bend_model = 'plate'; end
if nargin < 6 || isempty(mapmode),    mapmode    = 'analytic'; end
if nargin < 7 || isempty(fdh),        fdh        = 1e-5; end

nq1 = numel(qn{1}); nq2 = numel(qn{2});
npts = nq1 * nq2;

[I1g, I2g] = ndgrid(1:nq1, 1:nq2);
P_all = [qn{1}(I1g(:)).'; qn{2}(I2g(:)).'];

Sk = shell_curvature(geometry, P_all, mapmode, fdh);

D = YOUNG / (12 * (1 - POISSON^2));
comps = [2 0; 1 1; 0 2; 1 0; 0 1];
bendC = zeros(npts, 25);
detJ  = zeros(npts, 1);


for p = 1:npts
    ai = Sk.ainv(:,:,p);
    Sq = Sk.sqrta(p);
    Gm = Sk.Gamma(:,:,:,p);     % Gamma(alpha,beta,gamma)

    A11 = [1 0 0  -Gm(1,1,1)  -Gm(1,1,2)];
    A12 = [0 1 0  -Gm(1,2,1)  -Gm(1,2,2)];
    A22 = [0 0 1  -Gm(2,2,1)  -Gm(2,2,2)];
    Aof = cell(2,2);                       % Aof{alpha,beta}, with Aof{2,1} = Aof{1,2}
    Aof{1,1} = A11; Aof{1,2} = A12; Aof{2,1} = A12; Aof{2,2} = A22;

    % t = a^{alpha beta} A_{alpha beta}   (1 x 5)
    t = ai(1,1)*A11 + ai(1,2)*A12 + ai(2,1)*A12 + ai(2,2)*A22;

    % Q = sum_{al,be,g,d = 1..2} a^{al g} a^{be d} A_{al be}' A_{g d}   (5 x 5)
    % NOTE: the FULL sum over all four indices is required; the off-diagonal
    % (1,2)/(2,1) pair appears with multiplicity two (A is symmetric).
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

    switch lower(bend_model)
        case 'plate'
            Mm = D * Sq * (POISSON * (t.'*t) + (1-POISSON) * Q);
        case 't1'
            Mm = D * Sq * Q;
        otherwise
            error('quadrule_bend_shell: unknown bend_model ''%s''.', bend_model);
    end

    p_flat = I1g(p) + (I2g(p)-1)*nq1;   % == p (ndgrid order); kept explicit
    bendC(p_flat, :) = reshape(Mm, 1, 25);
    detJ(p_flat) = Sq;
end

info = struct('comps', comps, 'detJ', detJ, 'D', D, ...
              'bend_model', lower(bend_model), 'mapmode', lower(mapmode));
end
