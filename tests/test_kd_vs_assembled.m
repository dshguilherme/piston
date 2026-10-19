% T1 gate: harmonic K_D (WQ pages, metric-carrying bending) vs an assembled
% PHYSICAL Gauss reference, on BOTH an identity map and a non-identity affine
% map. The reference expands the parametric Hessian into the physical Hessian
% through the inverse map Jacobian, exactly as the setup does:
%   w_xx = cxx'h , w_yy = cyy'h , w_xy = cxy'h ,  h = [w_xi xi; w_xi eta; w_eta eta]
% with cxx = [Ji11^2; 2 Ji11 Ji21; Ji21^2] etc.
% Models: 'T1'    D(w_xx^2 + w_yy^2 + 2 w_xy^2)
%         'plate' D(w_xx^2 + w_yy^2 + 2 nu w_xx w_yy + 2(1-nu) w_xy^2)
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');

geoms = { 'unit square (identity)', 1.0, 1.0, [8 4];
          'rectangle 1.0 x 0.5',   1.0, 0.5, [8 4] };
NU = 0.3; E = 1.0; RHO = 1.0;
w = 2*pi*100; alpha = 1e-5; beta = 1e-7;
D = E/(12*(1-NU^2));

for ig = 1:size(geoms,1)
  gname = geoms{ig,1}; Lx = geoms{ig,2}; Ly = geoms{ig,3}; nsub = geoms{ig,4};
  problem_data = cantilever_beam(Lx, Ly);
  method_data.degree = [3 3]; method_data.regularity = [2 2];
  method_data.nsub = nsub; method_data.nquad = [4 4];
  [geometry, msh, space] = buildSpaces(problem_data, method_data);
  te = 1e-3 * ones(nsub);

  sp1 = space.scalar_spaces{1}.sp_univ(1);
  sp2 = space.scalar_spaces{1}.sp_univ(2);
  sp1h = sp_bspline_1d_param(sp1.knots, sp1.degree, msh.qn{1}, 'hessian', true);
  sp2h = sp_bspline_1d_param(sp2.knots, sp2.degree, msh.qn{2}, 'hessian', true);
  B1  = sp1.shape_functions;           B2  = sp2.shape_functions;
  Bx1 = sp1.shape_function_gradients;  Bx2 = sp2.shape_function_gradients;
  Bxx1= sp1h.shape_function_hessians;  Bxx2= sp2h.shape_function_hessians;
  qw1 = msh.qw{1}; qw2 = msh.qw{2};
  conn1 = sp1.connectivity; conn2 = sp2.connectivity;
  n1 = sp1.ndof;
  m1 = msh.nqn_dir(1); m2 = msh.nqn_dir(2);

  % metric coefficients at every Gauss point
  CXX = zeros(m1,m2,msh.nel_dir(1),msh.nel_dir(2),3);
  CYY = CXX; CXY = CXX; DJ = zeros(m1,m2,msh.nel_dir(1),msh.nel_dir(2));
  for el1 = 1:msh.nel_dir(1)
    for el2 = 1:msh.nel_dir(2)
      for g1 = 1:m1
        for g2 = 1:m2
          P = [msh.qn{1}(g1,el1); msh.qn{2}(g2,el2)];
          [~, jp] = geometry.map_der(P);
          dJ = jp(1,1)*jp(2,2) - jp(1,2)*jp(2,1);
          Ji = inv(jp);
          DJ(g1,g2,el1,el2) = dJ;
          CXX(g1,g2,el1,el2,:) = [Ji(1,1)^2, 2*Ji(1,1)*Ji(2,1), Ji(2,1)^2];
          CYY(g1,g2,el1,el2,:) = [Ji(1,2)^2, 2*Ji(1,2)*Ji(2,2), Ji(2,2)^2];
          CXY(g1,g2,el1,el2,:) = [Ji(1,1)*Ji(1,2), Ji(1,1)*Ji(2,2)+Ji(1,2)*Ji(2,1), Ji(2,1)*Ji(2,2)];
        end
      end
    end
  end
  fprintf('\n--- %s (detJ in [%.4f, %.4f]) ---\n', gname, min(DJ(:)), max(DJ(:)));

  for model = {'T1','plate'}
    H = harmonic_kd_setup(msh, space, geometry, E, NU, RHO, 'cpu', 'double', model{1});
    KD = harmonic_kd_form(H, te, w, alpha, beta, 'cpu_sparse');
    [Kb_ref, Mm_ref] = assemble_ref(msh, te, D, NU, B1,B2,Bx1,Bx2,Bxx1,Bxx2, ...
                                    qw1,qw2, conn1,conn2, n1, CXX,CYY,CXY,DJ, model{1});
    KDref = (1 + 1i*w*beta)*Kb_ref + (1i*w*alpha - w^2)*Mm_ref;
    KDref = 0.5*(KDref + KDref.');
    Vv = harmonic_kd_form(H, te, w, alpha, beta, 'values');
    Kb_ours = sparse(H.rows, H.cols, Vv.Kb, H.N, H.N);
    Mm_ours = sparse(H.rows, H.cols, Vv.Mm, H.N, H.N);
    y  = harmonic_kd_matvec(H, te, ones(H.N,1), w, alpha, beta);
    ym = KD * ones(H.N,1);
    eb = norm(Kb_ours - Kb_ref, 'fro')/norm(Kb_ref,'fro');
    em = norm(Mm_ours - Mm_ref, 'fro')/norm(Mm_ref,'fro');
    ed = norm(KD - KDref, 'fro')/norm(KDref,'fro');
    ev = norm(y - ym)/norm(ym);
    fprintf('  %-5s : bending %.2e | mass %.2e | K_D %.2e | matvec %.2e', model{1}, eb, em, ed, ev);
    if max([eb em ed ev]) < 1e-9, fprintf('   OK\n'); else, fprintf('   FAIL\n'); end
  end
end
disp('VERIFY_DONE');

function [Kb, Mm] = assemble_ref(msh, te, D, NU, B1,B2,Bx1,Bx2,Bxx1,Bxx2, ...
                                 qw1,qw2, conn1,conn2, n1, CXX,CYY,CXY,DJ, model)
  N = n1 * max(conn2(:));
  Kb = zeros(N); Mm = zeros(N);
  m1 = size(B1,1); m2 = size(B2,1);
  for el1 = 1:msh.nel_dir(1)
    for el2 = 1:msh.nel_dir(2)
      te_el = te(el1, el2);
      idx1 = conn1(:,el1); idx2 = conn2(:,el2);
      nA = numel(idx1); nB = numel(idx2);
      I = (idx1 - 1) + (idx2 - 1)' * n1 + 1; I = I(:);
      A = zeros(numel(I)); Am = zeros(numel(I));
      for g2 = 1:m2
        for g1 = 1:m1
          dJ  = DJ(g1,g2,el1,el2);
          cxx = reshape(CXX(g1,g2,el1,el2,:), 1, 3);
          cyy = reshape(CYY(g1,g2,el1,el2,:), 1, 3);
          cxy = reshape(CXY(g1,g2,el1,el2,:), 1, 3);
          wgt = qw1(g1,el1)*qw2(g2,el2)*dJ*D*te_el^3;
          % physical Hessian components for each local basis function
          wxr = zeros(nA,nB); wyr = wxr; wxyr = wxr;
          for a = 1:nA
            for b = 1:nB
              h = [Bxx1(g1,a,el1)*B2(g2,b,el2);
                   Bx1(g1,a,el1)*Bx2(g2,b,el2);
                   B1(g1,a,el1)*Bxx2(g2,b,el2)];
              wxr(a,b)  = cxx*h;
              wyr(a,b)  = cyy*h;
              wxyr(a,b) = cxy*h;
            end
          end
          for a = 1:nA
            for b = 1:nB
              r = a + (b-1)*nA;
              for c = 1:nA
                for d = 1:nB
                  s = c + (d-1)*nA;
                  if strcmp(model,'plate')
                    bend = wxr(a,b)*wxr(c,d) + wyr(a,b)*wyr(c,d) + ...
                           NU*(wxr(a,b)*wyr(c,d) + wyr(a,b)*wxr(c,d)) + ...
                           2*(1-NU)*wxyr(a,b)*wxyr(c,d);
                  else
                    bend = wxr(a,b)*wxr(c,d) + wyr(a,b)*wyr(c,d) + ...
                           2*wxyr(a,b)*wxyr(c,d);
                  end
                  A(r,s) = A(r,s) + wgt*bend;
                  Am(r,s) = Am(r,s) + te_el*(B1(g1,a,el1)*B2(g2,b,el2))* ...
                            (B1(g1,c,el1)*B2(g2,d,el2))*qw1(g1,el1)*qw2(g2,el2)*dJ;
                end
              end
            end
          end
        end
      end
      Kb(I,I) = Kb(I,I) + A; Mm(I,I) = Mm(I,I) + Am;
    end
  end
end
