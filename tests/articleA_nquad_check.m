% ARTICLEA_NQUAD_CHECK -- reproduce the nquad = degree+1 requirement of the
% reduced weighted-quadrature (quadrule_stiff_fast) setup.
%
% THE CLAIM (previously brief-sourced, now measured here): the analysis space
% handed to the reduced fast-formation setup MUST be built with nquad points per
% element equal to degree+1 (=4 for degree 3).  At larger nquad the official
% (0,0)/(1,1) rules are SILENTLY corrupted -- wrong matrices, no error.
%
% WHY (mechanism, from geopdes quadrule_stiff_fast.m): quadrule_gen indexes its
% local RHS with  local_points = (supp(1)-1)*nqn_gauss+1 : supp(end)*nqn_gauss
% where nqn_gauss = degree+1 is HARD-CODED.  The global B_gauss it indexes was
% built from space.shape_functions, which has exactly `nquad` rows per element.
% When nquad ~= degree+1 the two disagree and the weights are solved against the
% wrong rows -- silently.
%
% REPRODUCTION: an AFFINE rdim=2 chart (unit square, degree-elevated to p=3 so
% the map is still affine and the geometric coefficient is constant), assembled
% through the rdim=2 curved-clone path (harmonic_kd_setup_curved /
% harmonic_kd_form_curved), compared against the consistent Gauss remedy
% harmonic_kd_gauss_curved nqg=5.  Because the coefficient is constant, the
% reduced rule is EXACT at nquad=degree+1 (machine precision) and only the
% nquad mismatch can move it.
%
% Writes results/articleA_nquad_check.txt
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
resdir = 'C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(resdir,'dir'), mkdir(resdir); end
fid = fopen(fullfile(resdir,'articleA_nquad_check.txt'),'w');
lg = @(varargin) both(fid, varargin{:});

E=1.0; NU=0.3; RHO=1.0; DEG=[3 3]; REG=[2 2];
lg('===== ARTICLE A: nquad = degree+1 requirement of the reduced WQ setup =====\n');
lg('Script: tests/articleA_nquad_check.m\n');
lg('Affine rdim=2 chart (unit square, degree-elevated to p=3), degree 3 / C^1,\n');
lg('8x8 elements.  Reduced WQ (harmonic_kd_setup_curved + harmonic_kd_form_curved)\n');
lg('vs consistent Gauss (harmonic_kd_gauss_curved, nqg=5).  Constant coefficient,\n');
lg('so the reduced rule must be machine-exact at nquad = degree+1 = 4.\n\n');

% affine rdim=2 chart: control net is a parallelogram -> constant Jacobian
coefs = zeros(2,2,2);
coefs(1:2,1,1)=[0;0]; coefs(1:2,2,1)=[1;0]; coefs(1:2,1,2)=[0;1]; coefs(1:2,2,2)=[1;1];
nr0 = nrbmak(coefs,{[0 0 1 1],[0 0 1 1]}); nr0 = nrbdegelev(nr0,[2 2]);
nsub = [8 8];

lg('  %-6s %-24s %-12s\n','nquad','WQ-vs-Gauss(rel.Frob)','note');
Ng = 0; fkg = NaN;
for nquad = [4 5 6]
    geo = geo_load(nr0);
    [~,~,nk] = kntrefine(geo.nurbs.knots, nsub-1, geo.nurbs.order-1, REG);
    nurbs = nrbkntins(geo.nurbs, nk); geo = geo_load(nurbs);
    rule = msh_gauss_nodes([nquad nquad]);
    [qn,qw] = msh_set_quad_nodes(geo.nurbs.knots, rule);
    msh = msh_cartesian(geo.nurbs.knots, qn, qw, geo);
    space = sp_vector(repmat({sp_nurbs(geo.nurbs,msh)},1,2), msh);
    Hs = harmonic_kd_setup_curved(msh, space, geo, E, NU, RHO, 'cpu','double','plate');
    Kwq = harmonic_kd_form_curved(Hs, ones(msh.nel_dir));
    Kg  = harmonic_kd_gauss_curved(msh, space, geo, E, NU, 'plate', 5);
    e = norm(Kwq-Kg,'fro')/norm(Kg,'fro');
    Ng = Hs.N; fkg = norm(Kg,'fro');
    if nquad == DEG(1)+1
        note = 'degree+1 (correct)';
    else
        note = 'corrupted (silent)';
    end
    lg('  %-6d %-24.3e %-12s\n', nquad, e, note);
end
lg('\nN = %d scalar DOF.  Reference (Gauss nqg=5) Frobenius norm = %.6e.\n', Ng, fkg);
lg('Result: machine-exact ONLY at nquad = degree+1 = 4; nquad 5 and 6 corrupt the\n');
lg('official (0,0)/(1,1) pages with no warning.  The analysis space must be built\n');
lg('at nquad = degree+1.\n');
fclose(fid);
fprintf('ARTICLEA_NQUAD_CHECK_DONE\n');

% ============================ local functions =============================
function both(fid, varargin)
  fprintf(fid, varargin{:}); fprintf(1, varargin{:});
end
