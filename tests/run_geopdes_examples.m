% RUN_GEOPDES_EXAMPLES -- item 5 of the fourth directive.
%
% Locate and RUN, UNMODIFIED, the shell examples that GeoPDEs itself ships, on
% the corresponding obstacle problems.  Their printed result is the ground
% truth for the geometry/BC/load/normalisation setup.
%
% WHAT GEOpdes SHIPS (surveyed under geopdes-master):
%   * geopdes/inst/examples/kirchhoff_love_shell/ex_KL_shell_scordelis_lo_roof.m
%   * geopdes/inst/examples/kirchhoff_love_shell/ex_KL_shell_square.m
%   * geopdes/inst/examples/geometry_files/sphere_octant.m
%   * topopt/Functions/pinched_hemisphere.m   (a topopt SCRIPT, not a clean example)
%   * topopt/Functions/scordelis_problem.m
%   * NO pinched-cylinder example.
%   * NO solid-shell / Kirchhoff-solid obstacle example (only standard 3D
%     elasticity: op_lin_elast_* / ex_lin_elast_horseshoe / ex_plane_strain_*).
%
% It also performs the item-3(a) LIKE-FOR-LIKE operator comparison: assemble the
% stiffness matrix on the SAME space/geometry/quadrature through geopdes' native
% buildSpaces() and through our own build, and compare element-by-element.
%
% Saves results/geopdes_examples.txt.
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
resdir = 'C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(resdir,'dir'), mkdir(resdir); end
scratch = 'C:/Users/dshgu/AppData/Local/hermes/cache/scratch';
fid=fopen(fullfile(resdir,'geopdes_examples.txt'),'w');
lg=@(varargin) both(fid,varargin{:});
cwd0 = pwd; cd(scratch);   % examples write VTK into cwd

lg('==== GeoPDEs OWN EXAMPLES (run unmodified) ====\n\n');

%% ---- (5A) Scordelis-Lo roof example -------------------------------------
lg('---- ex_KL_shell_scordelis_lo_roof.m (unmodified) ----\n');
lg('(it prints its own numerical vs reference value; reference_displacement=-0.3005925)\n');
try
  clear u space geometry msh
  run('C:/Users/dshgu/OneDrive/Documents/geopdes-master/geopdes/inst/examples/kirchhoff_love_shell/ex_KL_shell_scordelis_lo_roof.m');
  pts{1}=0.0; pts{2}=0.5; [eu,~]=sp_eval(u,space,geometry,pts);
  lg('  [captured] u_z(edge mid) = %.8g   (vs obstacle ref -0.3024)\n', eu(3));
catch ME
  lg('  FAILED: %s\n', ME.message);
end
lg('\n');

%% ---- (5B) pinched hemisphere script -------------------------------------
lg('---- topopt/Functions/pinched_hemisphere.m (unmodified) ----\n');
lg('NOTE: geometry = sphere_octant(10) = nrbrevolve(nrbcirc(10,0,0,pi/2),...,pi/2):\n');
lg('      the meridian arc runs 0..90deg, so the POLAR HOLE IS NOT CUT. This is a\n');
lg('      full-octant hemisphere, NOT the 18deg-hole obstacle problem. It pins ALL\n');
lg('      dofs of side 3 and loads two single dofs.  It prints NO displacement\n');
lg('      value (only writes VTK); we capture u at its load dofs.\n');
try
  clear u sp geometry msh F_y F_z K F
  run('C:/Users/dshgu/OneDrive/Documents/geopdes-master/topopt/Functions/pinched_hemisphere.m');
  lg('  [captured] u(F_y) = %.6g   u(F_z) = %.6g\n', full(u(F_y)), full(u(F_z)));
  lg('  ndof=%d  #dofs_3(pinned)=%d  F_y=%d F_z=%d\n', sp.ndof, numel(dofs_3), F_y, F_z);
catch ME
  lg('  FAILED: %s\n', ME.message);
end
lg('\n');

%% ---- (5C) ex_KL_shell_square (sanity, unmodified) ------------------------
lg('---- ex_KL_shell_square.m (unmodified, sanity only) ----\n');
try
  clear u space geometry msh
  run('C:/Users/dshgu/OneDrive/Documents/geopdes-master/geopdes/inst/examples/kirchhoff_love_shell/ex_KL_shell_square.m');
  lg('  ran OK (see its own printed output above)\n');
catch ME
  lg('  FAILED: %s\n', ME.message);
end
lg('\n');

cd(cwd0);

%% ---- (3a) LIKE-FOR-LIKE operator comparison on the hemisphere -----------
lg('==== item 3(a): like-for-like K comparison (our build vs geopdes buildSpaces) ====\n');
lg('Geometry = sphere_octant(10); p=3, C^2, nsub=[16 16], nquad=[4 4].\n');
lg('The obstacle harness calls op_KL_shells_tp DIRECTLY, so the operator is the\n');
lg('library one; what can differ is our space/geometry/quadrature construction.\n');
E=6.825e7; nu=0.3; t=0.04; degree=3; regularity=2; nsub=[16 16]; nquad=[4 4];

% (i) geopdes native path
clear problem_data method_data
problem_data.geo_name = sphere_octant(10);
problem_data.E_coeff = @(x,y,z) E*ones(size(x));
problem_data.nu_coeff = @(x,y,z) nu*ones(size(x));
problem_data.thickness = t;
method_data.degree=[degree degree]; method_data.regularity=[regularity regularity];
method_data.nsub=nsub; method_data.nquad=nquad;
[gN,mN,spN]=buildSpaces(problem_data,method_data);
KN=op_KL_shells_tp(spN,spN,mN,@(x,y,z)E*ones(size(x)),@(x,y,z)nu*ones(size(x)),t);

% (ii) our path
[gO,mO,spO]=build_shell_rdim3(sphere_octant(10),[degree degree],[regularity regularity],nsub,nquad,3);
KO=op_KL_shells_tp(spO,spO,mO,@(x,y,z)E*ones(size(x)),@(x,y,z)nu*ones(size(x)),t);

lg('  geometry same points?      max|gN.map-gO.map| = %.3e\n', max(abs(gN.map([0.3;0.7])-gO.map([0.3;0.7]))));
lg('  ndof native=%d  ours=%d\n', spN.ndof, spO.ndof);
lg('  ||KN-KO||_F / ||KN||_F = %.3e\n', norm(KN-KO,'fro')/norm(KN,'fro'));
% per-element Frobenius norms
ne1=mN.nel_dir(1); ne2=mN.nel_dir(2);
dE=zeros(ne1,ne2);
% element dof blocks via connectivity
for e2=1:ne2
  for e1=1:ne1
    idxN = elem_dofs(spN,e1,e2); idxO = elem_dofs(spO,e1,e2);
    Kn = KN(idxN,idxN); Ko = KO(idxO,idxO);
    dE(e1,e2) = norm(Kn-Ko,'fro')/max(norm(Kn,'fro'),eps);
  end
end
lg('  per-element rel diff: max=%.3e  mean=%.3e\n', max(dE(:)), mean(dE(:)));
lg('  interior vs boundary (max rel diff):\n');
interior = dE(2:end-1,2:end-1);
lg('     interior(1,1)=%.3e  all-interior max=%.3e  boundary max=%.3e\n', ...
   dE(1,1), max(interior(:)), max([dE(1,:),dE(end,:),dE(:,1).',dE(:,end).']));

lg('\nGEOPDES_EXAMPLES_DONE\n');
fclose(fid);
fprintf('GEOPDES_EXAMPLES_DONE\n');

%% ============================== helpers ==================================
function idx=elem_dofs(sp,e1,e2)
  s1=sp.scalar_spaces{1}; n1=s1.ndof_dir(1); Ns=n1*s1.ndof_dir(2);
  c1=s1.sp_univ(1).connectivity(:,e1); c2=s1.sp_univ(2).connectivity(:,e2);
  idx=[];
  for c=1:3
    for b=1:numel(c2)
      idx=union(idx, c1+(c2(b)-1)*n1+(c-1)*Ns);
    end
  end
end
function [geometry,msh,space]=build_shell_rdim3(nr,degree,regularity,nsub,nquad,rdim)
  geometry=geo_load(nr);
  degelev=max(degree-(geometry.nurbs.order-1),0);
  nurbs=nrbdegelev(geometry.nurbs,degelev);
  [rk,zeta,nk]=kntrefine(nurbs.knots,nsub-1,nurbs.order-1,regularity);
  nurbs=nrbkntins(nurbs,nk); geometry=geo_load(nurbs);
  rule=msh_gauss_nodes(nquad); [qn,qw]=msh_set_quad_nodes(geometry.nurbs.knots,rule);
  msh=msh_cartesian(geometry.nurbs.knots,qn,qw,geometry);
  sp_scalar=sp_nurbs(geometry.nurbs,msh); space=sp_vector(repmat({sp_scalar},1,rdim),msh);
end
function both(fid,varargin)
  fprintf(fid,varargin{:}); fprintf(1,varargin{:});
end
