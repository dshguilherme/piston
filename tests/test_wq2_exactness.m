% EXACTNESS TEST for the new WQ rules, all types (k1,k2) with k1,k2 in {0,1,2}.
% Validates each rule exactly as the assembly contraction uses it:
%   sum_q w_a(q) * N_b^{(k2)}(x_q) == int N_a^{(k1)} N_b^{(k2)} dx
% for every original-space trial b sharing support with a, with the RHS
% integrated per element (the integrand is piecewise polynomial with breaks at
% interior knots, so a single Gauss rule over the whole support is NOT exact).
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
problem_data = cantilever_beam(1.0, 1.0);
method_data.degree = [4 4]; method_data.regularity = [3 3];
method_data.nsub = [6 3]; method_data.nquad = [5 5];
[~, msh, space] = buildSpaces(problem_data, method_data);
sp1 = space.scalar_spaces{1}.sp_univ(1);
knots = sp1.knots; p = sp1.degree; ndof = sp1.ndof;
dk = unique(knots);

Q2 = quadrule_bend_fast(sp1);
nqn = p+1;

% per-element Gauss grid over the whole domain (official layout)
[qn_all, qw_all] = msh_set_quad_nodes(dk, msh_gauss_nodes(nqn));
qn_all = qn_all(:); qw_all = qw_all(:);
sg = findspan(ndof-1, p, qn_all', knots);
Bg = basisfun(sg, qn_all', p, knots);
nbg = numbasisfun(sg, qn_all', p, knots) + 1;
rowsg = repmat((1:numel(qn_all)).', 1, p+1);
B_gauss = sparse(rowsg(:), nbg(:), Bg(:), numel(qn_all), ndof);
B2_gauss = full(Q2.B2_gauss);
B1_gauss = full(Q2.B1_gauss);

% cross-check the chains against basisfunder on the Gauss grid
tdg = basisfunder(sg, p, qn_all', knots, 2);
chk1 = 0; chk2 = 0;
for ig = 1:numel(qn_all)
    for k = 1:p+1
        chk2 = max(chk2, abs(B2_gauss(ig, nbg(ig,k)) - tdg(ig,3,k)));
        chk1 = max(chk1, abs(B1_gauss(ig, nbg(ig,k)) - tdg(ig,2,k)));
    end
end
fprintf('chain vs basisfunder (Gauss grid): 1st %.3e  2nd %.3e\n', chk1, chk2);

types = [2 2; 2 0; 0 2; 2 1; 1 2];
Coll   = {full(Q2.B_global),  full(Q2.B1_global),  full(Q2.B2_global)};
CollG  = {full(B_gauss),      full(B1_gauss),      full(B2_gauss)};
TestG  = {full(B_gauss),      full(B1_gauss),      full(B2_gauss)};
allok = true;
for it = 1:size(types,1)
    k1 = types(it,1); k2 = types(it,2);
    W = Q2.(sprintf('quad_weights_%d%d', k1, k2));
    maxerr = 0; npairs = 0; worst = [];
    for a = 1:ndof
        w = W{a}; if isempty(w), continue; end
        nbrs = unique(sp1.connectivity(:, sp1.supp{a}));
        iel = sp1.supp{a}';
        lp = ((iel(1)-1)*nqn + 1) : (iel(end)*nqn);
        Ta = TestG{k1+1}(lp, a);
        for b = nbrs(:)'
            Tb = CollG{k2+1}(lp, b);
            lhs = w * Coll{k2+1}(Q2.ind_points{a}, b);
            rhs = sum(qw_all(lp) .* Ta .* Tb);
            err = abs(lhs - rhs);
            if abs(rhs) > 1e-13 || abs(lhs) > 1e-13
                npairs = npairs + 1;
                if err > maxerr, maxerr = err; worst = [a b lhs rhs]; end
            end
        end
    end
    fprintf('rule (%d,%d): %3d pairs, max |lhs-rhs| = %.3e', k1, k2, npairs, maxerr);
    if maxerr < 1e-9
        fprintf('   OK\n');
    else
        allok = false;
        fprintf('   FAIL  worst (a,b)=(%d,%d) lhs=%.6e rhs=%.6e\n', worst(1), worst(2), worst(3), worst(4));
    end
end
if allok, disp('WQ_ALL_TYPES_OK'); else, disp('WQ_ALL_TYPES_FAIL'); end
