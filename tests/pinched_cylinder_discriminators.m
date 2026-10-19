% PINCHED_CYLINDER_DISCRIMINATORS -- targeted probes for the pinched-cylinder
% discrepancy (the full-model KL shell converges at ~2.7x the reference
% 1.8248e-5; see results/obstacle_course.txt).  Raw evidence, no forcing.
%
%  PART B  MODEL PROBE (plate vs T1).  The pipeline has two bending-coefficient
%          models: the full Kirchhoff-Love plate coefficient
%              M = a a' + b b' + nu (a b' + b a') + 2(1-nu) c c'   ('plate')
%          (this is what op_KL_shells_tp uses, via op_KL_bending_stress with the
%          Poisson term) and the shear-deformable/Reissner-Mindlin-type member
%              M = a a' + b b' + 2 c c'                            ('T1')
%          obtained by dropping the Poisson coupling in the BENDING term only
%          (membrane unchanged).  The IGA book attributes its "very slightly
%          softer than the benchmark" cylinder to TRANSVERSE SHEAR -- which only
%          acts in a shear-deformable formulation.  We run BOTH models on the
%          pinched cylinder and compare.
%
%  PART C  BC SIDE IDENTIFICATION.  Print the GeoPDEs side table (midpoints /
%          coordinate ranges / dof counts) for the cylindrical patch and VERIFY
%          that the diaphragm constraints u_x=u_y=0 land on the two END CIRCLES,
%          not on a longitudinal (seam) edge.
%
%  PART D  FREE-END DISCRIMINATOR.  Solve with the end diaphragms REMOVED
%          (all ends free).  If the diaphragm BC is really landing and really
%          stiffening, removing it must make the cylinder dramatically softer.
%
% Saves results/pinched_cylinder_discriminators.txt.
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
resdir = 'C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(resdir,'dir'), mkdir(resdir); end
fid=fopen(fullfile(resdir,'pinched_cylinder_discriminators.txt'),'w');
lg=@(varargin) both(fid,varargin{:});
REF=1.8248e-5; R=300; L=600; t=3; E=3.0e6; nu=0.3;

lg('==== PINCHED CYLINDER DISCRIMINATORS  (reference 1.8248e-5) ====\n\n');

%% ============================ PART C: SIDE TABLE ==========================
lg('---- PART C: GeoPDEs side table for the cylindrical patch ----\n');
lg('nrbcylind(L=600, R=300, 0..2pi): dir1 = theta (revolution), dir2 = z (axial).\n');
nru = nrbcylind(L, R, [0 0 0], 0, 2*pi);
nrc = nrbtform(nru, vecrotz(pi/2));
gu  = geo_load(nru); gc = geo_load(nrc);
lg('%-6s %-10s %-26s %-26s %-26s %-8s\n','patch','side','x-range','y-range','z-range','#bdof');
for ip=1:2
  if ip==1, g=gu; tag='raw'; else, g=gc; tag='rot90z'; end
  for s=1:numel(g.boundary)
    tt=linspace(0,1,21); P=g.boundary(s).map(tt); P=reshape(P,3,[]);
    ex=g.boundary(s).map(0); ex=ex(:);
    lg('%-6s %-10d [%8.2f,%8.2f] %15s [%8.2f,%8.2f] %15s [%8.2f,%8.2f] %15s\n', ...
       tag, s, min(P(1,:)),max(P(1,:)),'',min(P(2,:)),max(P(2,:)),'',min(P(3,:)),max(P(3,:)),'');
  end
end
lg('Interpretation: the two END CIRCLES are constant-z planes (z=0 and z=L=600);\n');
lg('the two SEAM/longitudinal edges lie on the theta=0 / theta=2pi generator.\n');
% identify which sides the diaphragm DOFs of the rotated patch fall on
[g,m,sp]=build_shell_rdim3(nrc,[3 3],[2 2],[16 16],[5 5],3);
lg('\nDiaphragm dofs (u_x=u_y) are taken from sides 3,4.  Per-side dof counts:\n');
for s=1:numel(sp.boundary)
  dx=numel(side_comp(sp,s,1)); dy=numel(side_comp(sp,s,2)); dz=numel(side_comp(sp,s,3));
  % geometry check: midpoint of side s
  P=g.boundary(s).map(0.5); P=P(:);
  lg('  side %d: midpoint=[%8.3f %8.3f %8.3f]  #ux=%5d #uy=%5d #uz=%5d\n', s, P(1),P(2),P(3),dx,dy,dz);
end
lg('CHECK: sides 3,4 must be the constant-z end circles (midpoint z=0 and z=600).\n\n');

%% ========================== PART B: PLATE vs T1 ==========================
lg('---- PART B: bending-coefficient model probe (plate vs T1) ----\n');
lg('%-10s %-8s %-9s %-14s %-10s\n','model','nsub','ndof','u_x(load)','ratio');
% self-check: custom assembler with 'plate' must reproduce op_KL_shells_tp
[g,m,sp]=build_shell_rdim3(nrc,[3 3],[2 2],[8 8],[5 5],3);
Kref=op_KL_shells_tp(sp,sp,m,@(x,y,z)E*ones(size(x)),@(x,y,z)nu*ones(size(x)),t);
Kcus=op_KL_shells_tp_model(sp,sp,m,@(x,y,z)E*ones(size(x)),@(x,y,z)nu*ones(size(x)),t,'plate');
lg('  [self-check] custom plate assembler vs op_KL_shells_tp: rel err = %.3e\n', ...
   norm(Kcus-Kref,'fro')/norm(Kref,'fro'));
for model={'plate','T1'}
  for nsub=[16 32 60]
    [g,m,sp]=build_shell_rdim3(nrc,[3 3],[2 2],[nsub nsub],[5 5],3);
    if strcmp(model{1},'plate')
      K=op_KL_shells_tp(sp,sp,m,@(x,y,z)E*ones(size(x)),@(x,y,z)nu*ones(size(x)),t);
    else
      K=op_KL_shells_tp_model(sp,sp,m,@(x,y,z)E*ones(size(x)),@(x,y,z)nu*ones(size(x)),t,'T1');
    end
    rhs=zeros(sp.ndof,1);
    rhs=rhs+point_load_rhs(sp,newton_inv(g,[R;0;L/2]),[-1;0;0]);
    rhs=rhs+point_load_rhs(sp,newton_inv(g,[-R;0;L/2]),[ 1;0;0]);
    dd=union([side_comp(sp,3,1);side_comp(sp,3,2)],[side_comp(sp,4,1);side_comp(sp,4,2)]);
    int=setdiff(1:sp.ndof,dd); u=zeros(sp.ndof,1); u(int)=K(int,int)\rhs(int);
    q=newton_inv(g,[R;0;L/2]); val=sp_eval(u,sp,g,{q(1),q(2)}); val=val(1);
    lg('%-10s %-8s %-9d %-14.6g %-10.4f\n', model{1}, mat2str([nsub nsub]), sp.ndof, val, val/REF);
  end
end
lg('\n');

%% ========================== PART D: FREE ENDS ============================
lg('---- PART D: free-end discriminator (diaphragms removed) ----\n');
lg('%-16s %-10s %-14s %-10s\n','BC','nsub','u_x(load)','ratio');
for nsub=[16 32]
  [g,m,sp]=build_shell_rdim3(nrc,[3 3],[2 2],[nsub nsub],[5 5],3);
  K=op_KL_shells_tp(sp,sp,m,@(x,y,z)E*ones(size(x)),@(x,y,z)nu*ones(size(x)),t);
  rhs=zeros(sp.ndof,1);
  rhs=rhs+point_load_rhs(sp,newton_inv(g,[R;0;L/2]),[-1;0;0]);
  rhs=rhs+point_load_rhs(sp,newton_inv(g,[-R;0;L/2]),[ 1;0;0]);
  % (i) diaphragm (baseline), (ii) free ends (no essential BC but pin the axial
  %     rigid translation at one end node so the system is solvable)
  cases = {'diaphragm', true; 'free(pin u_z)', false};
  for ic=1:2
    if cases{ic,2}
      dd=union([side_comp(sp,3,1);side_comp(sp,3,2)],[side_comp(sp,4,1);side_comp(sp,4,2)]);
      tagn='diaphragm';
    else
      dd=dof_at(sp,newton_inv(g,[R;0;0]),3);   % only kill axial rigid translation
      tagn='free(pin u_z)';
    end
    int=setdiff(1:sp.ndof,dd); u=zeros(sp.ndof,1); u(int)=K(int,int)\rhs(int);
    q=newton_inv(g,[R;0;L/2]); val=sp_eval(u,sp,g,{q(1),q(2)}); val=val(1);
    lg('%-16s %-10s %-14.6g %-10.4f\n', tagn, mat2str([nsub nsub]), val, val/REF);
  end
end
lg('\nNOTE: if the diaphragm is landing correctly, the free-end value must be\n');
lg('much softer (larger |u_x|).  A ~unchanged value would mean the diaphragm is\n');
lg('NOT landing or is not contributing as assumed.\n');

lg('\nPINCHED_CYL_DISCRIM_DONE\n');
fclose(fid);
fprintf('PINCHED_CYL_DISCRIM_DONE\n');

%% ============================== helpers ==================================
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
function A=op_KL_shells_tp_model(space_u,space_v,msh,E_coeff,nu_coeff,t_coeff,bend_model)
  % KL shell stiffness, tensor-product, with the BENDING Poisson coupling
  % optionally dropped ('T1'): bending nu=0, membrane nu=nu_coeff.
  A = spalloc(space_v.ndof, space_u.ndof, 3*space_u.ndof);
  for iel=1:msh.nel_dir(1)
    msh_col=msh_evaluate_col(msh,iel);
    sp_col_u=sp_evaluate_col_param(space_u,msh_col,'value',true,'gradient',true,'hessian',true);
    sp_col_v=sp_evaluate_col_param(space_v,msh_col,'value',true,'gradient',true,'hessian',true);
    for idim=1:msh.rdim
      x{idim}=reshape(msh_col.geo_map(idim,:,:),msh_col.nqn,msh_col.nel); %#ok<AGROW>
    end
    Ec=E_coeff(x{:}); nuc=nu_coeff(x{:});
    if strcmpi(bend_model,'T1'), nub=zeros(size(nuc)); else, nub=nuc; end
    A = A + assemble_KL(sp_col_u,sp_col_v,msh_col,Ec,nuc,nub,t_coeff);
  end
end
function mat=assemble_KL(ref_sp_u,ref_sp_v,msh,E_coeff,nu_mem,nu_bend,t_coeff)
  % copy of geopdes op_KL_shells with separate membrane/bending Poisson ratios.
  wq = msh.jacdet .* msh.quad_weights;
  [bending_stress,Kappa]   = op_KL_bending_stress (ref_sp_u,ref_sp_v,msh,E_coeff.*wq,nu_bend,t_coeff);
  [membrane_stress,Epsilon]= op_KL_membrane_stress(ref_sp_u,ref_sp_v,msh,E_coeff.*wq,nu_mem ,t_coeff);
  bending_stress  = reshape(bending_stress ,[],msh.nqn,ref_sp_u.nsh_max,msh.nel);
  membrane_stress = reshape(membrane_stress,[],msh.nqn,ref_sp_u.nsh_max,msh.nel);
  Kappa   = reshape(Kappa  ,[],msh.nqn,ref_sp_v.nsh_max,msh.nel);
  Epsilon = reshape(Epsilon,[],msh.nqn,ref_sp_v.nsh_max,msh.nel);
  rows=zeros(msh.nel*ref_sp_u.nsh_max*ref_sp_v.nsh_max,1);
  cols=rows; values=rows; ncounter=0;
  for iel=1:msh.nel
    if all(msh.jacdet(:,iel))
      ndof_u=ref_sp_u.nsh_max; ndof_v=ref_sp_v.nsh_max;
      bs=reshape(bending_stress (:, :, :, iel),[],msh.nqn,1,ndof_u);
      ms=reshape(membrane_stress(:, :, :, iel),[],msh.nqn,1,ndof_u);
      Ki=reshape(Kappa  (:, :, :, iel),[],msh.nqn,ndof_v,1);
      Ei=reshape(Epsilon(:, :, :, iel),[],msh.nqn,ndof_v,1);
      t1=sum(bsxfun(@times,bs,Ki),1); t2=sum(bsxfun(@times,ms,Ei),1);
      ev=reshape(sum(t1,2),ndof_v,ndof_u)+reshape(sum(t2,2),ndof_v,ndof_u);
      [rl,cl]=ndgrid(ref_sp_v.connectivity(:,iel),ref_sp_u.connectivity(:,iel));
      idx=rl & cl;
      rows(ncounter+(1:ref_sp_u.nsh(iel)*ref_sp_v.nsh(iel)))=rl(idx);
      cols(ncounter+(1:ref_sp_u.nsh(iel)*ref_sp_v.nsh(iel)))=cl(idx);
      values(ncounter+(1:ref_sp_u.nsh(iel)*ref_sp_v.nsh(iel)))=ev(idx);
      ncounter=ncounter+ref_sp_u.nsh(iel)*ref_sp_v.nsh(iel);
    end
  end
  mat=sparse(rows(1:ncounter),cols(1:ncounter),values(1:ncounter),ref_sp_v.ndof,ref_sp_u.ndof);
end
function d=side_comp(space,s,c)
  d=space.boundary(s).dofs(space.boundary(s).comp_dofs{c});
end
function gd=dof_at(space,q,c)
  sp=space.scalar_spaces{1}; Ns=prod(sp.ndof_dir); idx=cell(1,2); val=cell(1,2);
  for d=1:2
    su=sp.sp_univ(d); kn=su.knots(:).'; p=su.degree; nd=su.ndof;
    spn=findspan(nd-1,p,q(d),kn); b=basisfun(spn,q(d),p,kn); nb=numbasisfun(spn,q(d),p,kn)+1;
    idx{d}=nb(:); val{d}=b(:);
  end
  best=0; sd=1;
  for a=1:numel(idx{1})
    for b=1:numel(idx{2})
      w=val{1}(a)*val{2}(b); if w>best, best=w; sd=idx{1}(a)+(idx{2}(b)-1)*sp.ndof_dir(1); end
    end
  end
  gd=sd+(c-1)*Ns;
end
function q=newton_inv(g,xt)
  starts=[0.5 0.5;0.25 0.5;0.75 0.5;0.5 0.25;0.5 0.75;0.15 0.5;0.85 0.5;0.5 0.15;0.5 0.85];
  best=inf; qbest=[0.5;0.5];
  for s=1:size(starts,1)
    q=starts(s,:).';
    for it=1:80
      qc=min(max(q,1e-9),1-1e-9);
      F=g.map(qc)-xt; Jp=g.map_der(qc); J=reshape(Jp,3,2);
      dq=(J.'*J)\(J.'*F); q=min(max(qc-dq,1e-9),1-1e-9);
      if norm(dq)<1e-14,break;end
    end
    r=norm(g.map(q)-xt);
    if r<best, best=r; qbest=q; end
  end
  q=qbest;
end
function rhs=point_load_rhs(space,q,F)
  sp=space.scalar_spaces{1}; rhs=zeros(space.ndof,1);
  idx=cell(1,2); val=cell(1,2);
  for d=1:2
    su=sp.sp_univ(d); kn=su.knots(:).'; p=su.degree; nd=su.ndof;
    spn=findspan(nd-1,p,q(d),kn); b=basisfun(spn,q(d),p,kn); nb=numbasisfun(spn,q(d),p,kn)+1;
    idx{d}=nb(:); val{d}=b(:);
  end
  Ns=prod(sp.ndof_dir);
  for a=1:numel(idx{1})
    for b=1:numel(idx{2})
      sd=idx{1}(a)+(idx{2}(b)-1)*sp.ndof_dir(1); w=val{1}(a)*val{2}(b);
      for c=1:3, rhs(sd+(c-1)*Ns)=rhs(sd+(c-1)*Ns)+w*F(c); end
    end
  end
end
function both(fid,varargin)
  fprintf(fid,varargin{:}); fprintf(1,varargin{:});
end
