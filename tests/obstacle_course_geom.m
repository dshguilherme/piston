function problems = obstacle_course_geom(which)
% OBSTACLE_COURSE_GEOM  Geometry / BC / load / target definitions for the
% classical thin-shell OBSTACLE COURSE (Belytschko, Wong & Stolarski 1989),
% as analysed by Chapelle & Bathe (2003) and in the IGA literature
% (Hughes, Cottrell & Bazilevs 2005 ; Kiendl et al. 2009 ; Benson et al. 2010).
%
%   problems = obstacle_course_geom()        -> all three
%   problems = obstacle_course_geom('cylinder'|'hemisphere'|'roof')
%
% Each entry is a struct with the fields
%   .name        short identifier
%   .nurbs       NURBS *mid-surface* (NURBS-toolbox struct, rdim = 3), FULL
%                domain in the canonical orientation (reduce by symmetry when
%                assembling -- symmetry info in .symmetry)
%   .rdim        3 (a shell is a surface embedded in R^3)
%   .param_domain {[0 1],[0 1]} ; .deg, .order  surface degrees / orders
%   .thickness   t
%   .E, .nu       isotropic material
%   .body_load   @(x,y,z) -> [3 x N] distributed load per unit MID-SURFACE AREA
%   .point_loads struct array: .x (3x1 physical point), .F (3x1 force), .P mag
%   .dirichlet   struct: .sides (GeoPDEs side ids), .components (1..3),
%                .value (0 = homogeneous), .physical  (description)
%   .target      struct: .x (3x1), .component (1..3), .ref (benchmark value),
%                .normalisation (string), .quantity (string)
%
% -------------------------------------------------------------------------
% THE THREE BENCHMARKS  (dimensions, BCs, load, target -- see also the
% citation block at the bottom and tests/obstacle_course_references.md)
% -------------------------------------------------------------------------
% (1) PINCHED CYLINDER WITH END DIAPHRAGMS
%     R = 300, L = 600, t = 3 ; E = 3.0e6, nu = 0.3.
%     Closed circular cylinder, clamped-free by RIGID DIAPHRAGMS on the two
%     end circles (in-plane, i.e. transverse, displacements = 0; axial free).
%     Two radially OPPOSING unit point loads P = 1 at mid-length (z = L/2).
%     Severity: inextensional bending + a complex membrane state.
%     TARGET: radial displacement under a load point.
%     REF  u = 1.8248e-5  (exact  w = 164.24 * P/(E*t) = 0.18248e-4 ;
%          Belytschko/Wong/Stolarski 1989, p.404 ; Bouclier/Elguedj/
%          Combescure 2013 CMAME 267:86-110).
%
% (2) PINCHED HEMISPHERE WITH AN 18-DEGREE CIRCULAR HOLE
%     R = 10, t = 0.04 ; E = 6.825e7, nu = 0.3.
%     Hemisphere whose pole is cut out to a polar angle of 18 deg
%     (hole edge radius R*sin18 = 3.0902). Equator is otherwise FREE.
%     Four equal radial point loads at the four cardinal equator points,
%     alternating inward/outward (two in, two out), |P| = 1.
%     Severity: inextensional bending and rigid-body-motion capability.
%     TARGET: radial displacement under a load point (sign along the load).
%     REF  |u| = 0.0938  (MacNeal & Harder 1985 ; the analytic value of
%          Simo, Fox & Rifai 1989 is 0.093 ; some authors quote 0.094).
%
% (3) SCORDELIS-LO ROOF
%     R = 25, L = 50, half-angle 40 deg (total opening 80 deg), t = 0.25 ;
%     E = 4.32e8, nu = 0.0.
%     Cylindrical section; the two end arcs rest on RIGID DIAPHRAGMS
%     (in-plane displacement zero), the two straight (longitudinal) edges are
%     FREE. Loaded by its own weight, a uniform vertical load g = 90 per unit
%     area (downward, -z). Membrane dominated.
%     TARGET: vertical displacement (u_z) at the MIDPOINT of a free edge.
%     REF  u_z = 0.3024  (Belytschko/Wong/Stolarski 1989 ; Lee & Bathe 2002
%          Comput Struct 80:1441-1455 report 0.3024 at t = 0.25 ;
%          Ko, Lee & Bathe 2017 Comput Struct 193:187-206 also use 0.3024.)
% -------------------------------------------------------------------------

if nargin < 1 || isempty(which), which = 'all'; end

problems = struct('name', {}, 'nurbs', {}, 'rdim', {}, 'param_domain', {}, ...
                  'deg', {}, 'order', {}, 'thickness', {}, 'E', {}, 'nu', {}, ...
                  'body_load', {}, 'point_loads', {}, 'dirichlet', {}, ...
                  'target', {}, 'symmetry', {});

want = @(s) any(strcmpi(which, {'all', s}));

%% (1) PINCHED CYLINDER  --------------------------------------------------
if want('cylinder')
  R = 300; L = 600;
  % Full circular cylinder, axis = +z, from z=0 to z=L, centred on the z-axis.
  nurbs = nrbcylind(L, R, [0 0 0], 0, 2*pi);
  % load points (mid-length, diametrically opposed)
  Pc = [R; 0; L/2];   Pm = [-R; 0; L/2];
  pl = struct('x', {}, 'F', {}, 'P', {});
  pl(1) = struct('x', Pc, 'F', [-1;0;0], 'P', 1);   % inward at +x
  pl(2) = struct('x', Pm, 'F', [ 1;0;0], 'P', 1);   % inward at -x
  prob = struct( ...
    'name', 'pinched_cylinder', 'nurbs', nurbs, 'rdim', 3, ...
    'param_domain', {{[0 1],[0 1]}}, 'deg', [2 1], 'order', [3 2], ...
    'thickness', 3, 'E', 3.0e6, 'nu', 0.3, ...
    'body_load', @(x,y,z) zeros(3, numel(x)), ...
    'point_loads', pl, ...
    'dirichlet', struct('sides', [3 4], 'components', [1 2], 'value', 0, ...
        'physical', ['rigid end diaphragms at z=0 and z=L (eta=0,1): ' ...
                     'transverse in-plane displacements u_x = u_y = 0, ' ...
                     'axial u_z free']), ...
    'target', struct('x', Pc, 'component', 1, 'ref', 1.8248e-5, ...
        'quantity', 'radial (x) displacement under the point load', ...
        'normalisation', 'u = 164.24 * P/(E*t); P=1, E=3e6, t=3'), ...
    'symmetry', struct('planes', 'x=0, y=0, z=L/2 (eighth model)', ...
        'note', 'GeoPDEs sides: 3,4 = eta=0,1 (the two end circles)'));
  problems(end+1) = prob; %#ok<AGROW>
end

%% (2) PINCHED HEMISPHERE WITH 18-DEGREE HOLE  -----------------------------
if want('hemisphere')
  R = 10;
  hole = 18*pi/180;                     % polar angle removed at the pole
  % Meridian arc (in the x-y plane, angles measured from +x axis) from the
  % hole edge (polar angle 18 deg) to the equator (90 deg); revolve a full
  % turn about the x axis -> hemisphere with an 18 deg cap removed.
  crv  = nrbcirc(R, [0 0 0], hole, pi/2);
  nurbs = nrbrevolve(crv, [0 0 0], [1 0 0]);
  % Four cardinal equator points: (0, R cos(phi), R sin(phi)), phi = 0:90:270
  php  = [0 pi/2 pi -pi/2];
  eq   = R*[zeros(1,4); cos(php); sin(php)];
  sgn  = [-1 +1 -1 +1];                 % inward, outward, inward, outward
  pl = struct('x', {}, 'F', {}, 'P', {});
  for k = 1:4
    nrm = [0; cos(php(k)); sin(php(k))];   % radial direction at the equator
    pl(k) = struct('x', eq(:,k), 'F', sgn(k)*nrm, 'P', 1);
  end
  prob = struct( ...
    'name', 'pinched_hemisphere_hole', 'nurbs', nurbs, 'rdim', 3, ...
    'param_domain', {{[0 1],[0 1]}}, 'deg', [2 2], 'order', [3 3], ...
    'thickness', 0.04, 'E', 6.825e7, 'nu', 0.3, ...
    'body_load', @(x,y,z) zeros(3, numel(x)), ...
    'point_loads', pl, ...
    'dirichlet', struct('sides', [], 'components', [], 'value', 0, ...
        'physical', ['equator is FREE (no essential BC); the pole hole ' ...
                     'edge is traction-free. Symmetry (two planes through ' ...
                     'the loads) is used to retain one quadrant']), ...
    'target', struct('x', eq(:,1), 'component', 2, 'ref', 0.0938, ...
        'quantity', 'radial displacement under a load point', ...
        'normalisation', 'unit load P = 1, E = 6.825e7, t = 0.04'), ...
    'symmetry', struct('planes', 'y=0 and z=0 (quarter model)', ...
        'note', 'pole hole edge radius R*sin(18deg) = 3.0902'));
  problems(end+1) = prob; %#ok<AGROW>
end

%% (3) SCORDELIS-LO ROOF  -------------------------------------------------
if want('roof')
  R = 25; L = 50; th = 2.0/9.0*pi;      % half opening angle = 40 deg
  % GeoPDEs construction (ex_proof): circular arc 50deg..130deg in the x-y
  % plane, rotated by 90deg about x so the cylinder axis runs along +y, then
  % extruded a length L=50 along +y.  Cross-section is in the (x,z) plane.
  nrb   = nrbreverse(nrbcirc(R, [0 0 0], pi/2-th, pi/2+th));
  nurbs = nrbextrude(nrbtform(nrb, vecrotx(pi/2)), [0 L 0]);
  % Free-edge midpoint (u=0 -> arc end at 50 deg; v=0.5 -> y=L/2).
  xmid = [R*cos(pi/2-th); L/2; R*sin(pi/2-th)];   % = (R cos50, 25, R sin50)
  prob = struct( ...
    'name', 'scordelis_lo_roof', 'nurbs', nurbs, 'rdim', 3, ...
    'param_domain', {{[0 1],[0 1]}}, 'deg', [2 1], 'order', [3 2], ...
    'thickness', 0.25, 'E', 4.32e8, 'nu', 0.0, ...
    'body_load', @(x,y,z) -90*ones(3, numel(x)) .* [0;0;1], ...   % g_z = -90
    'point_loads', struct('x', {}, 'F', {}, 'P', {}), ...
    'dirichlet', struct('sides', [3 4], 'components', [1 3], 'value', 0, ...
        'physical', ['rigid diaphragms on the two end arcs (eta=0,1): ' ...
                     'u_x = u_z = 0; the two straight longitudinal edges ' ...
                     '(xi=0,1) are free']), ...
    'target', struct('x', xmid, 'component', 3, 'ref', 0.3024, ...
        'quantity', 'vertical (z) displacement at the midpoint of a free edge', ...
        'normalisation', 'E=4.32e8, nu=0, t=0.25, gravity g=90 per unit area'), ...
    'symmetry', struct('planes', 'x=0 (crown) and y=L/2 (quarter model)', ...
        'note', 'GeoPDEs sides: 3,4 = eta=0,1 (the two diaphragm end arcs)'));
  problems(end+1) = prob; %#ok<AGROW>
end
end

% =========================================================================
%  CITED REQUIREMENTS / REFERENCE DATA
% =========================================================================
%  T. Belytschko, B. L. Wong, H. Stolarski,
%    "Assumed strain stabilization procedure for the 9-node Lagrange shell
%     element", Int. J. Numer. Methods Engrg. 28 (1989) 385-414.
%    -> defines the three-problem shell obstacle course and the target values.
%
%  D. Chapelle, K.-J. Bathe,
%    "Fundamental considerations for the finite element analysis of shell
%     structures", Computers & Structures 66(1) (1998) 19-36;
%    "The Finite Element Analysis of Shells - Fundamentals", Springer, 2003
%     (2nd ed. 2011).  Extensive analysis of these three benchmarks, the
%     asymptotic (membrane vs bending) classification, and the reference
%     solutions used.
%
%  P.-S. Lee, K.-J. Bathe,
%    "On the asymptotic behavior of shell structures and the evaluation in
%     finite element solutions", Computers & Structures 80 (2002) 1441-1455.
%    -> Scordelis-Lo target 0.3024 at t=0.25 (p. ~1449).
%
%  Y. Ko, P.-S. Lee, K.-J. Bathe,
%    "Performance of the MITC3+ and MITC4+ shell elements in widely-used
%     benchmark problems", Computers & Structures 193 (2017) 187-206.
%    -> Scordelis-Lo reference displacement -0.3024.
%
%  R. H. MacNeal, R. L. Harder,
%    "A proposed standard set of problems to test finite element accuracy",
%     Finite Elements in Analysis and Design 1 (1985) 3-20.
%    -> pinched hemisphere with 18 deg hole, |u| = 0.0938.
%  J. C. Simo, D. D. Fox, M. S. Rifai,
%    "On a stress resultant geometrically exact shell model ...",
%     Comput. Methods Appl. Mech. Engrg. 73 (1989) 53-92.
%    -> analytic value 0.093 for the same problem.
%
%  IGA treatments of the same obstacle course:
%   - T. J. R. Hughes, J. A. Cottrell, Y. Bazilevs, "Isogeometric analysis:
%     CAD, finite elements, NURBS, exact geometry and mesh refinement",
%     Comput. Methods Appl. Mech. Engrg. 194 (2005) 4135-4195, Sec. 4.5
%     "Shell obstacle course" (all three problems; convergence to the
%     Belytschko et al. benchmark solutions).
%   - J. Kiendl, K.-U. Bletzinger, J. Linhard, R. Wuechner, "Isogeometric
%     shell analysis with Kirchhoff-Love elements", CMAME 198 (2009)
%     3902-3914 (Scordelis-Lo, pinched cylinder, pinched hemisphere).
%   - J. Kiendl, Y. Bazilevs, M.-C. Hsu, R. Wuechner, K.-U. Bletzinger, "The
%     bending strip method for IGA of Kirchhoff-Love shells ...", CMAME 199
%     (2010) 2403-2416 (pinched cylinder 1.8248e-5; hemispherical shell).
%   - D. J. Benson, Y. Bazilevs, M.-C. Hsu, T. J. R. Hughes, "Isogeometric
%     shell analysis: the Reissner-Mindlin shell", CMAME 199 (2010) 276-289
%     (obstacle course with the RM shell).
%   - R. Bouclier, T. Elguedj, A. Combescure, "Efficient isogeometric
%     NURBS-based solid-shell elements: mixed formulation and B-bar method",
%     CMAME 267 (2013) 86-110 (pinched cylinder 1.8248e-5, roof 0.3024).
%   - Geometrically exact/mixed KL: M. Hosseini et al. arXiv:1803.02721
%     (pinched cylinder 1.8248e-5, citing Kiendl 2010).
%
%  Origin/critical review of the roof value:
%   - M. A. Sanz, et al., "Benchmarking Computational Shell Models",
%     Archives of Computational Methods in Engineering (2022),
%     doi:10.1007/s11831-022-09798-5.  Traces 0.3024 to Belytschko et al.
%     (1985) and notes it superseded Carr's 0.3080/0.3086.
% =========================================================================
