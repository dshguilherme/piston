% VERIFICATION of the curved (non-affine) bending kinematics.
%
% PART A -- REGRESSION: on an AFFINE map the 5-component curved construction
%   must agree with the existing 3-component affine one (harmonic_kd_setup /
%   harmonic_kd_form) to machine precision (~1e-14), both at the coefficient
%   page level and at the operator level.
%
% PART B -- MANUFACTURED SOLUTION on a CURVED NURBS chart. A cubic manufactured
%   displacement w(xi,eta) (in the degree-3 spline space on every mesh) is
%   interpolated to nodal coefficients c; the WQ bending energy c'Kb c is
%   compared with a HIGH-ORDER DENSE GAUSS quadrature of the exact integrand
%   D det(J) (w_xx^2+w_yy^2+2nu w_xx w_yy+2(1-nu) w_xy^2) whose map derivatives
%   J,H are obtained by CENTRAL FINITE DIFFERENCES of x(xi,eta) (i.e. from the
%   geometry map alone, independently of geometry.map_der/map_der2). Refining
%   the mesh drives the relative error down: the WQ rule is exact only for a
%   CONSTANT coefficient, so a curved coefficient C(xi,eta) costs a quadrature
%   error that vanishes with h.
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');

E = 1.0; NU = 0.3; RHO = 1.0;
D = E/(12*(1-NU^2));
fprintf('==== PART A: affine regression vs harmonic_kd_setup (must be ~1e-14) ====\n');
geoms = {'rect 1.0x0.5 [5 4]', 1.0, 0.5, [5 4]; ...
         'square 1x1 [6 6]',   1.0, 1.0, [6 6]};
Aok = true;
for ig = 1:size(geoms,1)
 for model = {'plate','T1'}
  problem_data = cantilever_beam(geoms{ig,2}, geoms{ig,3});
  method_data.degree = [3 3]; method_data.regularity = [2 2];
  method_data.nsub = geoms{ig,4}; method_data.nquad = [4 4];
  [geometry, msh, space] = buildSpaces(problem_data, method_data);
  te = ones(method_data.nsub);
  He = harmonic_kd_setup(msh, space, geometry, E, NU, RHO, 'cpu', 'double', model{1});
  Vv = harmonic_kd_form(He, te, 1, 0, 0, 'values');
  Kb_aff = sparse(He.rows, He.cols, Vv.Kb, He.N, He.N); Kb_aff = 0.5*(Kb_aff+Kb_aff.');
  Hc = harmonic_kd_setup_curved(msh, space, geometry, E, NU, RHO, 'cpu', 'double', model{1});
  Kb_cur = harmonic_kd_form_curved(Hc, te);
  e_op = norm(Kb_cur - Kb_aff, 'fro') / norm(Kb_aff, 'fro');
  epage = 0;
  for ia = 1:3
    for ib = 1:3
      Ce = He.bendKp{(ia-1)*3 + ib}; Cc = Hc.bendKp{(ia-1)*5 + ib};
      if numel(Ce) == 1 || numel(Cc) == 1
        epage = max(epage, abs(double(Ce(1)) - double(Cc(1))) / max(abs(double(Ce(1))), 1));
      else
        epage = max(epage, max(abs(double(Ce(:)) - double(Cc(:)))) / max(max(abs(double(Ce(:)))), 1));
      end
    end
  end
  ex = 0;
  for ia = 1:5
    for ib = 1:5
      if ia<=3 && ib<=3, continue; end
      Cc = Hc.bendKp{(ia-1)*5+ib};
      ex = max(ex, max(abs(double(Cc(:)))));
    end
  end
  fprintf('  %-18s %-6s : op %.2e | block %.2e | new16 %.2e\n', geoms{ig,1}, model{1}, e_op, epage, ex);
  Aok = Aok && (max([e_op epage]) < 1e-12) && (ex < 1e-12);
 end
end
if Aok, fprintf('  PART A: OK\n'); else, fprintf('  PART A: FAIL\n'); end

fprintf('\n==== PART B: manufactured solution on a curved NURBS chart ====\n');
model = 'plate';
% curved degree-3 chart: bilinear -> elevate -> perturb interior control points
n0 = nrb4surf([0 0],[1 0],[0 1],[1 1]);
n0 = nrbdegelev(n0,[2 2]);           % order 4 (degree 3), 4x4 net
n0.coefs(1,2,2) = n0.coefs(1,2,2) + 0.15;
n0.coefs(2,2,2) = n0.coefs(2,2,2) - 0.15;
n0.coefs(1,3,3) = n0.coefs(1,3,3) - 0.12;
n0.coefs(2,3,3) = n0.coefs(2,3,3) + 0.12;
n0.coefs(1,2,3) = n0.coefs(1,2,3) + 0.06;
n0.coefs(2,3,2) = n0.coefs(2,3,2) - 0.06;

grid = linspace(0,1,21);
[GG1,GG2] = meshgrid(grid,grid);
detg = zeros(size(GG1));
for a=1:numel(GG1)
  Jp = n0_mapder([GG1(a);GG2(a)]);
  detg(a) = Jp(1,1)*Jp(2,2)-Jp(1,2)*Jp(2,1);
end
fprintf('  curved chart det(J) range on [0,1]^2 grid: [%.4f, %.4f]\n', min(detg(:)), max(detg(:)));

nsubs = [2 2; 4 4; 8 8; 16 16; 32 32];
ng = 14;   % dense Gauss points per element per direction (reference)
fprintf('\n  %-8s %6s %14s %14s %10s %12s %12s %7s\n', ...
        'nsub','N','E_WQ','E_ref(FD)','rel.err','rel.noS','interp.res','order');
fprintf('  %s\n', repmat('-',1,96));
prev = NaN; prevh = NaN;
for is = 1:size(nsubs,1)
  pd = struct('geo_name', n0);
  md = struct('degree',[3 3],'regularity',[2 2],'nsub',nsubs(is,:),'nquad',[5 5]);
  [geometry, msh, space] = buildSpaces(pd, md);
  te = ones(nsubs(is,:));

  Hc = harmonic_kd_setup_curved(msh, space, geometry, E, NU, RHO, 'cpu','double',model);
  Kb = harmonic_kd_form_curved(Hc, te);
  Kb_nos = form_noS(Hc, te);   % affine (3x3) block only: ignores the S (map Hessian) term

  % nodal coefficients of the manufactured cubic w (exact in the degree-3 space)
  sp1 = space.scalar_spaces{1}.sp_univ(1);
  sp2 = space.scalar_spaces{1}.sp_univ(2);
  B1g = sp_colloc(sp1); B2g = sp_colloc(sp2);
  X1 = repmat(msh.qn{1}(:), 1, numel(msh.qn{2}));
  X2 = repmat(msh.qn{2}(:).', numel(msh.qn{1}), 1);
  W = wfun(X1, X2);
  C = B1g \ W / B2g.';
  interp_res = norm(B1g*C*B2g.' - W) / norm(W);
  c = C(:);
  E_wq = c.' * Kb * c;
  E_nos = c.' * Kb_nos * c;

  % independent dense-Gauss reference with FINITE-DIFFERENCED map derivatives
  E_ref = dense_gauss_ref(msh, geometry, ng, D, NU, model);
  rel = abs(E_wq - E_ref)/abs(E_ref);
  rel_nos = abs(E_nos - E_ref)/abs(E_ref);
  h = 1/nsubs(is,1);
  if isnan(prev), order = NaN; else, order = log(prev/rel)/log(prevh/h); end
  fprintf('  %-8s %6d %14.8e %14.8e %10.3e %12.3e %12.2e %7.2f\n', ...
          mat2str(nsubs(is,:)), Hc.N, E_wq, E_ref, rel, rel_nos, interp_res, order);
  prev = rel; prevh = h;
end

fprintf('\n==== PART C: pointwise physical Hessian of the curved chart ====\n');
fprintf('  chain rule (Ji, S=inv-map 2nd deriv from analytic map) vs\n');
fprintf('  central FD of w(x,y) in PHYSICAL coords (Newton inversion of the map)\n');
pd = struct('geo_name',n0); md = struct('degree',[3 3],'regularity',[2 2],'nsub',[6 6],'nquad',[5 5]);
[geometry, msh, space] = buildSpaces(pd, md);
rng(7); nx = 200; e = zeros(nx,3);
for k = 1:nx
  xi = 0.08 + 0.84*rand; eta = 0.08 + 0.84*rand;
  [~, J]  = geometry.map_der([xi;eta]);
  [~, ~, H] = geometry.map_der2([xi;eta]);
  Ji = inv(J); S = inv_map_der2(Ji, H);
  Phi = [3*xi-2*eta; -2*xi+1.4*eta+1; 1.4*xi+1.8*eta; ...
         1.5*xi^2-2*xi*eta+0.7*eta^2+eta+0.2; -xi^2+1.4*xi*eta+0.9*eta^2+xi+0.1];
  bxx=[Ji(1,1)^2;2*Ji(1,1)*Ji(2,1);Ji(2,1)^2;S(1,1,1);S(2,1,1)];
  byy=[Ji(1,2)^2;2*Ji(1,2)*Ji(2,2);Ji(2,2)^2;S(1,2,2);S(2,2,2)];
  bxy=[Ji(1,1)*Ji(1,2);Ji(1,1)*Ji(2,2)+Ji(2,1)*Ji(1,2);Ji(2,1)*Ji(2,2);S(1,1,2);S(2,1,2)];
  wxx_a=bxx.'*Phi; wyy_a=byy.'*Phi; wxy_a=bxy.'*Phi;
  P0 = geometry.map([xi;eta]); d = 1e-3;
  W = @(xx,yy) wphys(geometry, xx, yy, xi, eta);
  w0 = wfun(xi,eta);
  wxx_f = (W(P0(1)+d,P0(2)) - 2*w0 + W(P0(1)-d,P0(2)))/d^2;
  wyy_f = (W(P0(1),P0(2)+d) - 2*w0 + W(P0(1),P0(2)-d))/d^2;
  wxy_f = (W(P0(1)+d,P0(2)+d)-W(P0(1)+d,P0(2)-d)-W(P0(1)-d,P0(2)+d)+W(P0(1)-d,P0(2)-d))/(4*d^2);
  sc = max(abs([wxx_a wyy_a wxy_a]));
  e(k,:) = abs([wxx_a-wxx_f, wyy_a-wyy_f, wxy_a-wxy_f])/sc;
end
fprintf('  max rel diff (w_xx, w_yy, w_xy) = %.3e, %.3e, %.3e\n', max(e(:,1)), max(e(:,2)), max(e(:,3)));
if max(e(:)) < 1e-5, fprintf('  PART C: OK (kinematic transformation validated)\n'); else, fprintf('  PART C: FAIL\n'); end

fprintf('\n==== PART D: is the w_xx term needed? (affine a,b,c vectors ONLY) ====\n');
% Single sample point: physical w_xx computed with and without the S term.
xi=0.5; eta=0.5;
[~,J]=geometry.map_der([xi;eta]); [~,~,H]=geometry.map_der2([xi;eta]); Ji=inv(J); S=inv_map_der2(Ji,H);
Phi=[3*xi-2*eta; -2*xi+1.4*eta+1; 1.4*xi+1.8*eta; 1.5*xi^2-2*xi*eta+0.7*eta^2+eta+0.2; -xi^2+1.4*xi*eta+0.9*eta^2+xi+0.1];
bxx_nos=[Ji(1,1)^2;2*Ji(1,1)*Ji(2,1);Ji(2,1)^2;0;0];
bxx_full=[Ji(1,1)^2;2*Ji(1,1)*Ji(2,1);Ji(2,1)^2;S(1,1,1);S(2,1,1)];
wxx_nos=bxx_nos.'*Phi; wxx_full=bxx_full.'*Phi;
fprintf('  at (xi,eta)=(0.5,0.5):  w_xx without S = %.6e , with S = %.6e , |S-part| = %.3e\n', ...
        wxx_nos, wxx_full, abs(wxx_full-wxx_nos));
disp('CURVED_TEST_DONE');

% ---------------- local functions ----------------
function B1g = sp_colloc(sp1)
  nqn = size(sp1.shape_functions,1); nel = size(sp1.shape_functions,3);
  rows = repmat(reshape(1:nqn*nel,[nqn,1,nel]),[1,sp1.nsh_max,1]);
  cols = repmat(reshape(sp1.connectivity,[1,sp1.nsh_max,nel]),[nqn,1,1]);
  B1g = sparse(rows(:), cols(:), sp1.shape_functions(:), nqn*nel, sp1.ndof);
end

function Knos = form_noS(H, t_e)
  % The AFFINE operator a curved map would get if the map Hessian (S) term is
  % dropped -- i.e. what harmonic_kd_setup.m computes. Uses only the 3x3 block.
  N = H.N;
  tq = t_e(H.e1, H.e2); s3 = tq(H.LIN);
  comps = H.comps; V = [];
  for ia = 1:3
    for ib = 1:3
      C = H.bendKp{(ia-1)*5 + ib};
      i=comps(ia,1); j=comps(ia,2); k=comps(ib,1); l=comps(ib,2);
      P1 = H.WBe{1}{i+1,k+1}; P2 = H.WBe{2}{j+1,l+1};
      T = s3 .* C;
      R = pagemtimes(T,'none',P2,'transpose'); R = pagemtimes(P1,R);
      if isempty(V), V = R; else, V = V + R; end
    end
  end
  Vb = V(H.mask);
  Knos = sparse(H.rows, H.cols, double(Vb), N, N);
  Knos = 0.5*(Knos + Knos.');
end

function [Jp] = n0_mapder(ps)
  persistent g
  if isempty(g)
    n = nrb4surf([0 0],[1 0],[0 1],[1 1]);
    n = nrbdegelev(n,[2 2]);
    n.coefs(1,2,2) = n.coefs(1,2,2) + 0.15;
    n.coefs(2,2,2) = n.coefs(2,2,2) - 0.15;
    n.coefs(1,3,3) = n.coefs(1,3,3) - 0.12;
    n.coefs(2,3,3) = n.coefs(2,3,3) + 0.12;
    n.coefs(1,2,3) = n.coefs(1,2,3) + 0.06;
    n.coefs(2,3,2) = n.coefs(2,3,2) - 0.06;
    g = geo_load(n);
  end
  [~, Jp] = g.map_der(ps);
end

function w = wfun(xi, eta)
  w = 0.5*xi.^3 - xi.^2.*eta + 0.7*xi.*eta.^2 + 0.3*eta.^3 + xi.*eta + 0.2*xi + 0.1*eta;
end

function E = dense_gauss_ref(msh, geometry, ng, D, NU, model)
  rule = msh_gauss_nodes(ng);
  [qd1, qw1] = msh_set_quad_nodes(unique(msh.breaks{1}), rule);
  [qd2, qw2] = msh_set_quad_nodes(unique(msh.breaks{2}), rule);
  nel1 = size(qd1,2); nel2 = size(qd2,2);
  E = 0;
  for el2 = 1:nel2
    for el1 = 1:nel1
      for g2 = 1:ng
        for g1 = 1:ng
          xi = qd1(g1,el1); eta = qd2(g2,el2);
          wq = qw1(g1,el1)*qw2(g2,el2);
          [J,H] = fd_map_derivatives(geometry, xi, eta, 2e-3);
          dJ = J(1,1)*J(2,2)-J(1,2)*J(2,1);
          Ji = inv(J);
          S = inv_map_der2(Ji, H);
          % analytic parametric derivatives of the manufactured field
          wxi   = 1.5*xi^2 - 2*xi*eta + 0.7*eta^2 + eta + 0.2;
          weta  = -xi^2 + 1.4*xi*eta + 0.9*eta^2 + xi + 0.1;
          wxixi = 3*xi - 2*eta;
          wxieta= -2*xi + 1.4*eta + 1;
          wetaeta=1.4*xi + 1.8*eta;
          bxx = [Ji(1,1)^2; 2*Ji(1,1)*Ji(2,1); Ji(2,1)^2; S(1,1,1); S(2,1,1)];
          byy = [Ji(1,2)^2; 2*Ji(1,2)*Ji(2,2); Ji(2,2)^2; S(1,2,2); S(2,2,2)];
          bxy = [Ji(1,1)*Ji(1,2); Ji(1,1)*Ji(2,2)+Ji(2,1)*Ji(1,2); Ji(2,1)*Ji(2,2); S(1,1,2); S(2,1,2)];
          Phi = [wxixi; wxieta; wetaeta; wxi; weta];
          wxx = bxx.'*Phi; wyy = byy.'*Phi; wxy = bxy.'*Phi;
          if strcmpi(model,'plate')
            dens = wxx^2 + wyy^2 + 2*NU*wxx*wyy + 2*(1-NU)*wxy^2;
          else
            dens = wxx^2 + wyy^2 + 2*wxy^2;
          end
          E = E + wq * D * dJ * dens;
        end
      end
    end
  end
end

function [J, H] = fd_map_derivatives(geometry, xi, eta, h0)
  % adaptive step: keep the 5-point stencil inside the parametric domain [0,1]^2
  h = min(h0, 0.45*min([xi, 1-xi, eta, 1-eta]));
  h = max(h, 5e-6);
  x0  = geometry.map([xi; eta]);
  xp  = geometry.map([xi+h; eta]);  xm = geometry.map([xi-h; eta]);
  yp  = geometry.map([xi; eta+h]);  ym = geometry.map([xi; eta-h]);
  xpp = geometry.map([xi+h; eta+h]); xpm = geometry.map([xi+h; eta-h]);
  xmp = geometry.map([xi-h; eta+h]); xmm = geometry.map([xi-h; eta-h]);
  J = [(xp-xm)/(2*h), (yp-ym)/(2*h)];
  H = zeros(2,2,2);
  H(:,1,1) = (xp - 2*x0 + xm)/h^2;
  H(:,2,2) = (yp - 2*x0 + ym)/h^2;
  H(:,1,2) = (xpp - xpm - xmp + xmm)/(4*h^2);
  H(:,2,1) = H(:,1,2);
end

function S = inv_map_der2(Ji, H)
  S = zeros(2,2,2);
  for alpha=1:2
    for i=1:2
      for j=1:2
        s=0;
        for m=1:2, for b=1:2, for c=1:2
          s = s - Ji(alpha,m)*H(m,b,c)*Ji(b,i)*Ji(c,j);
        end, end, end
        S(alpha,i,j)=s;
      end
    end
  end
end

function w = wphys(geometry, xx, yy, xi0, eta0)
  q = newton_inv(geometry, xx, yy, xi0, eta0);
  w = wfun(q(1), q(2));
end

function q = newton_inv(geometry, xx, yy, xi0, eta0)
  q = [xi0; eta0];
  for it = 1:30
    F = geometry.map(q) - [xx; yy];
    [~, Jp] = geometry.map_der(q);
    J = reshape(Jp, 2, 2);
    dq = J \ F;
    q = q - dq;
    if norm(dq) < 1e-14, break; end
  end
end
