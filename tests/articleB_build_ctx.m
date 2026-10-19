function ctx = articleB_build_ctx(nsub, E, NU, RHO, t0, freq, zeta)
% ARTICLEB_BUILD_CTX  Curved-shell harmonic context for Article B.
%
% Builds the SCORDELIS-LO ROOF (a genuinely CURVED cylindrical shell, R=25,
% L=50, opening 80 deg) as a KIRCHHOFF-LOVE shell: geometry, space, and the
% precomputed element matrices
%     Bke  bending   (carries E)      -> scales as t^3
%     Ske  membrane  (carries E)      -> scales as t
%     Me   mass      (carries RHO)    -> scales as t
% plus the vectorised assembly indices, the Dirichlet set (rigid end diaphragms
% u_x=u_z=0 on the two end arcs, plus one axial u_y dof pinned to remove the
% rigid axial translation), and the harmonic body load (self-weight direction).
%
% WHY A CURVED SHELL AND NOT THE FLAT PLATE: on a flat plate all three harmonic
% objectives are DEGENERATE (the optimiser exploits a local antiresonance at the
% excitation point). The curvature couples membrane and bending so the design
% must redistribute material globally. The roof is an OPEN shell (no closed
% parametric seam), so the C0-seam hinge that afflicts closed rational circles
% in a Kirchhoff-Love model cannot arise here; the operator is assembled by
% Gauss quadrature (op_KL_shells_elements), which handles rational NURBS exactly.
%
% freq: excitation frequency (Hz), or NaN to build without a frequency.

addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');

pd = scordelis_problem;
pd.E_coeff  = @(x,y,z) E*ones(size(x));
pd.nu_coeff = @(x,y,z) NU*ones(size(x));
md = struct('degree',[2 2],'regularity',[1 1],'nsub',nsub,'nquad',[3 3]);
[geometry,msh,space]=buildSpaces(pd,md);
msh1 = msh_precompute(msh);
sp1 = sp_precompute_param(space,msh1,'value',true,'gradient',true,'hessian',true);
for idim=1:msh1.rdim, x{idim}=reshape(msh1.geo_map(idim,:,:),msh1.nqn,msh1.nel); end
[Bke,Ske,~]=op_KL_shells_elements(sp1,sp1,msh1,pd.E_coeff(x{:}),pd.nu_coeff(x{:}));
Me = op_u_v_elements(sp1,sp1,msh1,RHO);
lm = sp1.connectivity.';
ne = msh1.nel; n = size(lm,2); ndof = space.ndof;

Rows = zeros(ne, n*n); Cols = zeros(ne, n*n);
for a=1:n
  for b=1:n
    Rows(:,(a-1)*n+b) = lm(:,a);
    Cols(:,(a-1)*n+b) = lm(:,b);
  end
end

dr=[]; for s=pd.drchlt_sides
  comps = pd.drchlt_components{find(pd.drchlt_sides==s,1)};
  for ic=comps, dr=union(dr,space.boundary(s).dofs(space.boundary(s).comp_dofs{ic})); end
end
dy3 = space.boundary(3).dofs(space.boundary(3).comp_dofs{2});
dr = union(dr, dy3(1));
free = setdiff(1:ndof, dr);
F = op_f_v_tp(space,msh,pd.f);

ctx = struct();
ctx.Bke=Bke; ctx.Ske=Ske; ctx.Me=Me; ctx.lm=lm; ctx.ne=ne; ctx.ndof=ndof;
ctx.Rows=Rows; ctx.Cols=Cols; ctx.free=free; ctx.dr=dr; ctx.F=F;
ctx.E=E; ctx.RHO=RHO; ctx.t0=t0; ctx.nsub=nsub;
ctx.space=space; ctx.msh=msh1; ctx.geometry=geometry;
if ~isnan(freq)
  ctx.w = 2*pi*freq; ctx.zeta = zeta;
  ctx.alpha = 2*zeta*ctx.w; ctx.beta = 0;
else
  ctx.w=NaN; ctx.alpha=0; ctx.beta=0;
end
ctx.R0 = [];
end
