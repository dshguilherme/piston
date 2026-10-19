% RULE-ACCURACY TEST (1D) -- WHY the curved bending pages do not converge.
%
% The reduced weighted-quadrature rules of quadrule_bend_fast are exact for a
% CONSTANT coefficient C(x). This test measures the error of the assembly
% contraction
%     sum_q w_a(q) C(x_q) N_b^{(k2)}(x_q)   vs   int N_a^{(k1)} C N_b^{(k2)} dx
% (integrated exactly by a 20-point Gauss rule per element) as the mesh refines,
% for C == const, C linear, and C quadratic. Result: for a VARYING C the
% derivative-bearing rule types (2,2),(2,1),(1,0),(0,1) do NOT converge -- their
% error is O(1) (and (2,2) grows). This is the mechanism behind the curved-map
% non-convergence of test_kd_curved.m; it is a property of the reduced rule
% (constant-coefficient exactness only), not of the curved coefficient page.
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');

p = 3; nqn = p+1;
types = [2 2;2 1;1 1;2 0;1 0;0 2;0 1;0 0];
nsubs = [4 8 16 32];
cfuns = {@(x) 1+0*x, @(x) 1+0.5*x, @(x) 1+0.5*x+0.3*x.^2};
cnames = {'const','linear','quadratic'};
summary = struct();
for ic = 1:3
  Cfun = cfuns{ic};
  fprintf('\nC = %s\n', cnames{ic});
  fprintf('%6s','nsub'); for it=1:size(types,1), fprintf('  (%d,%d)',types(it,1),types(it,2)); end; fprintf('\n');
  E = zeros(numel(nsubs), size(types,1));
  for is = 1:numel(nsubs)
    nsub = nsubs(is);
    kn = [zeros(1,p) linspace(0,1,nsub+1) ones(1,p)];
    qn = msh_set_quad_nodes(unique(kn), msh_gauss_nodes(nqn));
    sp = sp_bspline_1d_param(kn,p,qn); Q2 = quadrule_bend_fast(sp); ndof=sp.ndof;
    sg=findspan(ndof-1,p,Q2.all_points,kn); B0=basisfun(sg,Q2.all_points,p,kn); nbg=numbasisfun(sg,Q2.all_points,p,kn)+1;
    rows=repmat((1:numel(Q2.all_points)).',1,p+1); Bg=sparse(rows(:),nbg(:),B0(:),numel(Q2.all_points),ndof);
    Coll={Bg,Q2.B1_global,Q2.B2_global};
    ng=20; dk=unique(kn); [qall,wall]=msh_set_quad_nodes(dk,msh_gauss_nodes(ng)); qall=qall(:).'; wall=wall(:);
    sgg=findspan(ndof-1,p,qall,kn); Bgg=basisfun(sgg,qall,p,kn); nbgg=numbasisfun(sgg,qall,p,kn)+1;
    rowg=repmat((1:numel(qall)).',1,p+1); Bgg=sparse(rowg(:),nbgg(:),Bgg(:),numel(qall),ndof);
    kd1=kn(2:end-1); p1=p-1; nd1=ndof-1; v1=findspan(nd1-1,p1,qall,kd1); bv1=basisfun(v1,qall,p1,kd1); n1=numbasisfun(v1,qall,p1,kd1)+1;
    r1=repmat((1:numel(qall)).',1,p1+1); Bv1=sparse(r1(:),n1(:),bv1(:),numel(qall),nd1);
    kd2=kn(3:end-2); p2=p-2; nd2=ndof-2; v2=findspan(nd2-1,p2,qall,kd2); bv2=basisfun(v2,qall,p2,kd2); n2=numbasisfun(v2,qall,p2,kd2)+1;
    r2=repmat((1:numel(qall)).',1,p2+1); Bv2=sparse(r2(:),n2(:),bv2(:),numel(qall),nd2);
    CollG={Bgg, Bv1*Q2.D1, Bv2*(Q2.D2*Q2.D1)};
    Cq=Cfun(qall(:)); Cw=Cfun(Q2.all_points(:));
    for it=1:size(types,1)
      k1=types(it,1); k2=types(it,2); W=Q2.(sprintf('quad_weights_%d%d',k1,k2)); mx=0;
      for a=1:ndof
        w=W{a}; if isempty(w), continue; end
        nbrs=unique(sp.connectivity(:,sp.supp{a})); ip=Q2.ind_points{a};
        for b=nbrs(:)'
          lhs=w*(Cw(ip).*full(Coll{k2+1}(ip,b)));
          rhs=sum(wall.*full(CollG{k1+1}(:,a)).*Cq.*full(CollG{k2+1}(:,b)));
          if abs(rhs)>1e-12, mx=max(mx,abs(lhs-rhs)/abs(rhs)); end
        end
      end
      E(is,it)=mx;
    end
    fprintf('%6d',nsub); for it=1:size(types,1), fprintf('  %.1e',E(is,it)); end; fprintf('\n');
  end
  summary.(cnames{ic}) = E;
end
fprintf('\nVERDICT:\n');
fprintf('  const C : rules exact to %.1e (design property; (2,2) degrades at nsub>=16)\n', max(max(summary.const(1:2,:))));
fprintf('  varying C: rule (2,2) error %.2e -> %.2e ; rule (1,0) error %.2e -> %.2e\n', ...
        summary.linear(1,types(:,1)==2&types(:,2)==2), summary.quadratic(end,types(:,1)==2&types(:,2)==2), ...
        summary.linear(1,types(:,1)==1&types(:,2)==0), summary.quadratic(end,types(:,1)==1&types(:,2)==0));
disp('RULE_TEST_DONE');
