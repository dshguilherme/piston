function S = shell_curvature(geometry, P, mapmode, fdh)
% SHELL_CURVATURE  Kirchhoff-Love surface kinematics of a 2D manifold embedded
% in R^3 (rdim = 3), evaluated at parametric points P = [xi; eta].
%
% This is the genuine SHELL counterpart of the rdim=2 "curved chart" kinematics
% of quadrule_bend_curved.m. Where that file computes the physical CARTESIAN
% Hessian w_{,ij} of a planar chart, a thin-walled shell is a curved surface and
% its bending strain is built from the SURFACE metric, the SECOND FUNDAMENTAL
% FORM and the Christoffel symbols of the surface.
%
% GEOMETRY. Let x(xi,eta) in R^3 be the surface map (our geometry.map), with
%   a_alpha    = x_{,alpha}                 (covariant / tangent basis, 3-vectors)
%   J          = [x_{,xi} x_{,eta}]         (3 x 2 Jacobian; columns tangent vectors)
%   H(:,al,be) = x_{,alpha beta}            (3 x 2 x 2 map Hessian)
% Then
%   surface metric            a_{albe}  = a_al . a_be = J' J           (2 x 2)
%   area element              sqrt(a)   = sqrt(det(a_{albe})) = |x_{,xi} x x_{,eta}|
%   unit normal               n         = (x_{,xi} x x_{,eta}) / sqrt(a)
%   second fundamental form   b_{albe}  = n . x_{,alpha beta}          (2 x 2)
%   Christoffel (2nd kind)    Gamma^g_{albe} = a^{g d} (x_{,d} . x_{,alpha beta})
% with a^{gd} = (a^{-1})^{gd}. (First kind Gamma_{d,albe} = a_d . x_{,alpha beta}.)
%
% PHYSICAL SECOND DERIVATIVES. For a scalar field w(xi,eta) on the surface the
% covariant (physical) second derivative is
%   w_{;alpha beta} = w_{,alpha beta} - Gamma^g_{alpha beta} w_{,g} .
% The unit-normal sign follows the parameterization orientation (x_{,xi} x x_{,eta}).
%
% INDEPENDENT CHECKS (used by tests/test_shell_curvature.m):
%   (i)  b via the Weingarten route -n_{,alpha} . x_{,beta} (FD of the normal),
%   (ii) Gamma via FD of the metric, Gamma_{d,albe} = (a_{d al,be}+a_{d be,al}-a_{al be,d})/2,
%   (iii) the parameterization-invariant curvatures K = det(b)/det(a) and
%         Hm = (a^{albe} b_{albe})/2 against the analytic values of the surface.
%
% INPUTS
%   geometry : geopdes geometry struct (.map, .map_der, .map_der2, .rdim=3).
%   P        : 2 x npts matrix of parametric coordinates (xi; eta).
%   mapmode  : 'analytic' (default) uses geometry.map_der / geometry.map_der2;
%              'fd' uses central finite differences of geometry.map only
%              (independent of the analytic derivatives).
%   fdh      : FD step for mapmode='fd' (default 1e-5; clamped to stay in domain).
%
% OUTPUT S : struct with fields, each over the npts points:
%   a (2,2,npts), ainv (2,2,npts), sqrta (1,npts), normal (3,npts),
%   b (2,2,npts), Gamma (2,2,2,npts) stored as Gamma(alpha,beta,gamma),
%   J (3,2,npts), H (3,2,2,npts), rdim, npts, mapmode.

if nargin < 3 || isempty(mapmode), mapmode = 'analytic'; end
if nargin < 4 || isempty(fdh),     fdh     = 1e-5;       end

P = double(P);
assert(size(P,1) == 2, 'shell_curvature: P must be 2 x npts (xi; eta).');
npts = size(P,2);

% physical space dimension (rdim): prefer geometry.rdim, else probe the map
if isfield(geometry,'rdim') && ~isempty(geometry.rdim)
    rdim = geometry.rdim;
else
    rdim = size(geometry.map(P(:,1)), 1);
end
assert(rdim >= 2, 'shell_curvature: physical dimension must be >= 2.');

switch lower(mapmode)
    case 'analytic'
        Jr = geometry.map_der(P);   J = reshape(Jr, rdim, 2, []);
        Hr = geometry.map_der2(P);  H = reshape(Hr, rdim, 2, 2, []);
        assert(size(J,3)==npts && size(H,4)==npts, ...
            'shell_curvature: map_der/map_der2 returned unexpected sizes.');
    case 'fd'
        [J, H] = fd_derivs(geometry, P, fdh, rdim);
    otherwise
        error('shell_curvature: unknown mapmode ''%s''.', mapmode);
end

a    = zeros(2,2,npts);
ainv = zeros(2,2,npts);
sqrta = zeros(1,npts);
nrm  = zeros(rdim,npts);
b    = zeros(2,2,npts);
Gam  = zeros(2,2,2,npts);

for p = 1:npts
    Jp = J(:,:,p);
    Hp = H(:,:,:,p);
    ap = Jp.' * Jp;                       % a_{alpha beta}
    a(:,:,p) = ap;
    dt = ap(1,1)*ap(2,2) - ap(1,2)*ap(2,1);
    if dt <= 0
        error('shell_curvature: non-positive metric determinant at point %d.', p);
    end
    sq = sqrt(dt);
    sqrta(p) = sq;
    ai = ap \ eye(2);
    ainv(:,:,p) = ai;

    c1 = Jp(:,1); c2 = Jp(:,2);           % x_{,xi}, x_{,eta}
    if rdim == 3
        nn = cross(c1, c2) / sq;
    else                                   % rdim==2: flat plane, +z normal
        nn = [0; 0; 1];
    end
    nrm(:,p) = nn;

    % b_{alpha beta} = n . x_{,alpha beta}
    for al = 1:2
        for be = 1:2
            b(al,be,p) = nn.' * Hp(:,al,be);
        end
    end

    % Gamma^gamma_{alpha beta} = a^{gamma delta} (x_{,delta} . x_{,alpha beta})
    for al = 1:2
        for be = 1:2
            for g = 1:2
                s = 0;
                for d = 1:2
                    s = s + ai(g,d) * (Jp(:,d).' * Hp(:,al,be));
                end
                Gam(al,be,g,p) = s;
            end
        end
    end
end

S = struct('a',a, 'ainv',ainv, 'sqrta',sqrta, 'normal',nrm, 'b',b, ...
           'Gamma',Gam, 'J',J, 'H',H, 'rdim',rdim, 'npts',npts, ...
           'mapmode',lower(mapmode));
end

% ---------------------------------------------------------------------------
function [J, H] = fd_derivs(geometry, P, h0, rdim)
% Central finite differences of geometry.map at every point, with the step
% clamped so the 5-point stencil stays inside [0,1]^2 (the parametric domain
% of the surfaces used here). Independent of map_der / map_der2.
npts = size(P,2);
J = zeros(rdim,2,npts);
H = zeros(rdim,2,2,npts);
for p = 1:npts
    xi = P(1,p); eta = P(2,p);
    h = min(h0, 0.4*min([xi, 1-xi, eta, 1-eta]));
    h = max(h, 1e-5);
    x0  = geometry.map([xi; eta]);
    xp  = geometry.map([xi+h; eta]);   xm = geometry.map([xi-h; eta]);
    yp  = geometry.map([xi; eta+h]);   ym = geometry.map([xi; eta-h]);
    xpp = geometry.map([xi+h; eta+h]); xpm = geometry.map([xi+h; eta-h]);
    xmp = geometry.map([xi-h; eta+h]); xmm = geometry.map([xi-h; eta-h]);
    J(:,1,p) = (xp - xm) / (2*h);
    J(:,2,p) = (yp - ym) / (2*h);
    H(:,1,1,p) = (xp - 2*x0 + xm) / h^2;
    H(:,2,2,p) = (yp - 2*x0 + ym) / h^2;
    H(:,1,2,p) = (xpp - xpm - xmp + xmm) / (4*h^2);
    H(:,2,1,p) = H(:,1,2,p);
end
end
