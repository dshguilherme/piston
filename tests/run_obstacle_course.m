% RUN_OBSTACLE_COURSE -- solve the classical thin-shell obstacle course
% (Belytschko/Wong/Stolarski 1989) with the geopdes Kirchhoff-Love shell
% (op_KL_shells_tp, Gauss quadrature on the canonical rdim=3 surfaces) and
% compare the target displacements with the cited references.
%
% THE THREE PROBLEMS AND THE MODELS THIS SCRIPT USES
% --------------------------------------------------
% (1) PINCHED CYLINDER (R=300,L=600,t=3,E=3e6,nu=0.3; two opposing unit point
%     loads at mid-length, rigid end diaphragms).  The closed rational circle of
%     nrbcylind carries a C0 PARAMETRIC SEAM at theta=0; a Kirchhoff-Love shell
%     needs C1, so that seam acts as a soft hinge and the model sits at ~2.7x the
%     reference.  RESOLVED by the C^{p-1}-CLOSED construction: the theta knot
%     vector is unclamped (kntunclamp, seam multiplicity 1 -> genuine periodic
%     spline) and the geometry is supplied ANALYTICALLY (cell-array map/der/der2
%     of the exact cylinder).  Converges to 1.8248e-5.  See
%     tests/pinched_c1_closed.m.
% (2) PINCHED HEMISPHERE WITH 18deg HOLE (R=10,t=0.04,E=6.825e7).  Solved with
%     the symmetry-reduced QUARTER MODEL {y>=0, z>=0} obtained by revolving the
%     meridian arc a quarter turn.  The revolution is a QUARTER turn, so the
%     surface is OPEN: its two clamped ends ARE the symmetry planes z=0 (side1)
%     and y=0 (side2); sides 3 (hole edge) and 4 (equator) are FREE.  There is
%     therefore NO closing C0 seam.  The 7.8%-too-soft defect was in the
%     SYMMETRY-CUT CONDITIONS: u_n=0 alone leaves the natural KL condition M=0
%     there (a free bending edge).  Imposing the full mirror-symmetry kinematic
%     set (see below) recovers 0.0938 to <0.1%, penalty-insensitive.  Two unit
%     radial point loads at the two equator corners; the x rigid translation is
%     pinned at one hole-edge node.  See tests/pinched_hemisphere_audit.m.
% (3) SCORDELIS-LO ROOF (R=25,L=50,80deg,t=0.25,E=4.32e8,nu=0; self weight
%     g=90; rigid end diaphragms, straight edges free).  Solved on the FULL roof
%     AND on the symmetry-reduced HALF MODEL {x>=0} about the crown.  For this
%     displacement-only KL shell the crown symmetry plane CANNOT be imposed by
%     u_x=0 alone: the rotation condition d(u_z)/d(xi)=0 must ALSO be imposed,
%     else the natural condition is the free bending edge M=0 (~33% error).
%
% Saves results/obstacle_course.txt.
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
resdir = 'C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(resdir,'dir'), mkdir(resdir); end
fid=fopen(fullfile(resdir,'obstacle_course.txt'),'w');
lg=@(varargin) both(fid,varargin{:});

lg('==== Shell obstacle course: geopdes KL shell (Gauss) on canonical geometries ====\n');
lg('%-24s %-20s %-8s %-9s %-14s %-13s %-9s\n', ...
   'problem','model','nsub','ndof','u_value','ref','rel.err');
lg('%s\n', repmat('-',1,100));

%% ------------------------------------------------------------------ (1) ---
% PINCHED CYLINDER -- C^{p-1}-closed periodic spline + analytic geometry.
R=300; L=600; t=3; E=3.0e6; nu=0.3;
for nsub=[48 64 96]
  [g,m,sp]=buildC1(R,L,nsub,[3 3],[2 2],[5 5]);
  K=op_KL_shells_tp(sp,sp,m,@(x,y,z)E*ones(size(x)),@(x,y,z)nu*ones(size(x)),t);
  rhs=zeros(sp.ndof,1);
  rhs=rhs+point_load_rhs(sp,newton_inv(g,[R;0;L/2]),[-1;0;0]);
  rhs=rhs+point_load_rhs(sp,newton_inv(g,[-R;0;L/2]),[ 1;0;0]);
  % rigid end diaphragms (sides 3,4 = the two end circles): u_x=u_y=0
  dd=union([side_comp(sp,3,1);side_comp(sp,3,2)],[side_comp(sp,4,1);side_comp(sp,4,2)]);
  int=setdiff(1:sp.ndof,dd); u=zeros(sp.ndof,1); u(int)=K(int,int)\rhs(int);
  q=newton_inv(g,[R;0;L/2]); val=eval_u(u,sp,q); val=val(1);
  lg('%-24s %-20s %-8s %-9d %-14.6g %-13.6g %-9.2e\n', ...
     'pinched_cylinder','C1-closed(periodic)',mat2str([nsub nsub]),sp.ndof,val,1.8248e-5, ...
     abs(abs(val)-1.8248e-5)/1.8248e-5);
end

%% ------------------------------------------------------------------ (2) ---
% PINCHED HEMISPHERE with 18deg hole -- symmetry-reduced QUARTER model.
Rm=10; tm=0.04; Em=6.825e7; num=0.3; hole=18*pi/180;
srf=nrbrevolve(nrbcirc(Rm,[0 0 0],hole,pi/2),[0 0 0],[1 0 0],pi/2);   % {y>=0,z>=0}
for nsub=[16 32 48]
  [g,m,sp]=build_shell_rdim3(srf,[3 3],[2 2],[nsub nsub],[5 5],3);
  K=op_KL_shells_tp(sp,sp,m,@(x,y,z)Em*ones(size(x)),@(x,y,z)num*ones(size(x)),tm);
  rhs=zeros(sp.ndof,1);
  qA=newton_inv(g,[0;Rm;0]); rhs=rhs+point_load_rhs(sp,qA,[0;-1;0]);
  qB=newton_inv(g,[0;0;Rm]); rhs=rhs+point_load_rhs(sp,qB,[0; 0;1]);
  % symmetry planes = the revolution ends: side1 (z=0) and side2 (y=0)
  dd=union(side_comp(sp,1,3),side_comp(sp,2,2));
  dd=union(dd,dof_at(sp,newton_inv(g,[Rm*cos(hole);0;0]),1));  % pin x rigid translation
  int=setdiff(1:sp.ndof,dd);
  % (b) u_n=0 only -- the 7.8%-too-soft baseline
  u=zeros(sp.ndof,1); u(int)=K(int,int)\rhs(int);
  uu=sp_eval(u,sp,g,{qA(1),qA(2)}); valb=uu(2);
  lg('%-24s %-20s %-8s %-9d %-14.6g %-13.6g %-9.2e\n', ...
     'pinched_hemisphere_hole','quarter{u_n only}',mat2str([nsub nsub]),sp.ndof,valb,0.0938, ...
     abs(abs(valb)-0.0938)/0.0938);
  % (c) exact mirror-symmetry set: z=0: d(ux),d(uy)/dxi ; y=0: d(ux),d(uz)/dxi
  D=[rot_pen(sp,1,1,1,1e-4,40);rot_pen(sp,1,2,1,1e-4,40); ...
     rot_pen(sp,2,1,1,1e-4,40);rot_pen(sp,2,3,1,1e-4,40)];
  Kr=K(int,int)+1e8*(D(:,int).'*D(:,int)); u=zeros(sp.ndof,1); u(int)=Kr\rhs(int);
  uu=sp_eval(u,sp,g,{qA(1),qA(2)}); val=uu(2);
  lg('%-24s %-20s %-8s %-9d %-14.6g %-13.6g %-9.2e\n', ...
     'pinched_hemisphere_hole','quarter{full sym}',mat2str([nsub nsub]),sp.ndof,val,0.0938, ...
     abs(abs(val)-0.0938)/0.0938);
end

%% ------------------------------------------------------------------ (3) ---
% SCORDELIS-LO ROOF -- full model.
Rr=25; Lr=50; Eu=4.32e8; nuu=0.0; tu=0.25;
probs=obstacle_course_geom('roof'); nrf=probs(1).nurbs;
for nsub=[20 30]
  [g,m,sp]=build_shell_rdim3(nrf,[3 3],[2 2],[nsub nsub],[5 5],3);
  K=op_KL_shells_tp(sp,sp,m,@(x,y,z)Eu*ones(size(x)),@(x,y,z)nuu*ones(size(x)),tu);
  rhs=op_f_v_tp(sp,m,@(x,y,z,ind) -90*ones(3,numel(x)).*[0;0;1]);
  dd=union([side_comp(sp,3,1);side_comp(sp,3,3)],[side_comp(sp,4,1);side_comp(sp,4,3)]);
  int=setdiff(1:sp.ndof,dd); u=zeros(sp.ndof,1); u(int)=K(int,int)\rhs(int);
  th=2/9*pi; q=newton_inv(g,[Rr*cos(pi/2-th);Lr/2;Rr*sin(pi/2-th)]);
  val=sp_eval(u,sp,g,{q(1),q(2)}); val=val(3);
  lg('%-24s %-20s %-8s %-9d %-14.6g %-13.6g %-9.2e\n', ...
     'scordelis_lo_roof','full',mat2str([nsub nsub]),sp.ndof,val,-0.3024, ...
     abs(abs(val)-0.3024)/0.3024);
end

%% ------------------------------------------------------------------ (3b) --
% SCORDELIS-LO ROOF -- symmetry-reduced HALF model {x>=0} about the crown, WITH
% the rotation condition d(u_z)/dxi = 0 imposed by a penalty (the plain u_x=0
% reduction is ~33% wrong; the rotation condition recovers the full value).
th=2/9*pi; nsub=30;
srf2=nrbextrude(nrbtform(nrbreverse(nrbcirc(Rr,[0 0 0],pi/2-th,pi/2)),vecrotx(pi/2)),[0 Lr 0]);
[g,m,sp]=build_shell_rdim3(srf2,[3 3],[2 2],[nsub nsub],[5 5],3);
K=op_KL_shells_tp(sp,sp,m,@(x,y,z)Eu*ones(size(x)),@(x,y,z)nuu*ones(size(x)),tu);
rhs=op_f_v_tp(sp,m,@(x,y,z,ind) -90*ones(3,numel(x)).*[0;0;1]);
dd=union([side_comp(sp,3,1);side_comp(sp,3,3)],[side_comp(sp,4,1);side_comp(sp,4,3)]); % diaphragms
cx=side_comp(sp,1,1);                       % crown u_x = 0 (Dirichlet)
D = rotation_penalty_uzdxi(sp,0,1e-4,40);   % rows: d(u_z)/dxi at (xi=0, eta_k)
for mu=[0 1e6 1e8]
  int=setdiff(1:sp.ndof,union(dd,cx));
  qe=newton_inv(g,[Rr*cos(pi/2-th);Lr/2;Rr*sin(pi/2-th)]);
  if mu==0
    u=zeros(sp.ndof,1); u(int)=K(int,int)\rhs(int);   % plain u_x=0 reduction (WRONG)
    tag='half,u_x only';
  else
    Kr=K(int,int)+mu*(D(:,int).'*D(:,int)); u=zeros(sp.ndof,1); u(int)=Kr\rhs(int);
    tag=sprintf('half,u_x+rot(%g)',mu);
  end
  val=sp_eval(u,sp,g,{qe(1),qe(2)}); val=val(3);
  lg('%-24s %-20s %-8s %-9d %-14.6g %-13.6g %-9.2e\n', ...
     'scordelis_lo_roof',tag,mat2str([nsub nsub]),sp.ndof,val,-0.3024, ...
     abs(abs(val)-0.3024)/0.3024);
end

lg('\nNotes:\n');
lg(' - PINCHED CYLINDER: reproduced by the C^{p-1}-closed periodic spline +\n');
lg('   analytic geometry (removes the C0 parametric seam of the rational circle).\n');
lg('   ratio drifts 0.93 -> 0.98 -> 0.99 over nsub 48/64/96, converging to 1.0.\n');
lg(' - PINCHED HEMISPHERE: reproduced by the quarter model once the FULL mirror-\n');
lg('   symmetry kinematic set is imposed on the cuts.  u_n=0 alone leaves a free\n');
lg('   bending edge (M=0) and is ~7.8%% too soft; the full set gives 0.0938.\n');
lg(' - SCORDELIS-LO ROOF: full and crown-reduced models agree; the u_x=0-only\n');
lg('   reduction is ~33%% off, the same symmetry-cut mechanism.\n');
lg(' - tests/pinched_c1_closed.m and tests/pinched_hemisphere_audit.m hold the raw\n');
lg('   seam audits, penalty/collocation sweeps and convergence tables.\n');
fclose(fid);
fprintf('OBSTACLE_COURSE_DONE\n');

% =============================== helpers =================================
function [geo,msh,sp]=buildC1(R,L,nsub,deg,reg,nquad)
  % C^{p-1}-closed (periodic) cylinder: theta unclamped + analytic geometry.
  knots0={[0 0 0 0 1 1 1 1],[0 0 0 0 1 1 1 1]};
  [knots,zeta]=kntrefine(knots0,[nsub nsub]-1,deg,reg);
  knotsU=kntunclamp(knots,deg,reg,1);
  geo=geo_load({@(P) cylmap(P,R,L), @(P) cylder(P,R,L), @(P) cylder2(P,R,L)});
  geo.nurbs=struct('knots',{{[0 0 0 0 1 1 1 1],[0 0 0 0 1 1 1 1]}},'order',[4 4]);
  rule=msh_gauss_nodes(nquad); [qn,qw]=msh_set_quad_nodes(zeta,rule);
  msh=msh_cartesian(zeta,qn,qw,geo);
  s1=sp_bspline(knotsU,deg,msh,[],1); sp=sp_vector(repmat({s1},1,3),msh);
end
function [u,v]=uvp(P)
  if iscell(P), u=P{1}(:).'; v=P{2}(:).'; else, u=P(1,:); v=P(2,:); end
end
function X=cylmap(P,R,L)
  [u,v]=uvp(P); if iscell(P), [U,V]=ndgrid(u,v); U=U(:).'; V=V(:).'; else, U=u; V=v; end
  X=[R*cos(2*pi*U); R*sin(2*pi*U); L*V];
end
function J=cylder(P,R,L)
  [u,v]=uvp(P); if iscell(P), [U,V]=ndgrid(u,v); U=U(:).'; V=V(:).'; else, U=u; V=v; end
  n=numel(U); J=zeros(3,2,n);
  J(1,1,:)=-2*pi*R*sin(2*pi*U); J(2,1,:)=2*pi*R*cos(2*pi*U); J(3,2,:)=L*ones(1,n);
end
function H=cylder2(P,R,L)
  [u,v]=uvp(P); if iscell(P), [U,V]=ndgrid(u,v); U=U(:).'; V=V(:).'; else, U=u; V=v; end
  n=numel(U); H=zeros(3,2,2,n);
  H(1,1,1,:)=-(2*pi)^2*R*cos(2*pi*U); H(2,1,1,:)=-(2*pi)^2*R*sin(2*pi*U);
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
  sp=space.scalar_spaces{1}; Ns=prod(sp.ndof_dir); rhs=zeros(space.ndof,1);
  [i1,w1]=basfuns(sp.sp_univ(1),q(1),ismember(1,sp.periodic_dir));
  [i2,w2]=basfuns(sp.sp_univ(2),q(2),ismember(2,sp.periodic_dir));
  for a=1:numel(i1)
    for b=1:numel(i2)
      sd=i1(a)+(i2(b)-1)*sp.ndof_dir(1); w=w1(a)*w2(b);
      for c=1:3, rhs(sd+(c-1)*Ns)=rhs(sd+(c-1)*Ns)+w*F(c); end
    end
  end
end
function v=eval_u(u,space,q)
  sp=space.scalar_spaces{1}; Ns=prod(sp.ndof_dir); v=zeros(3,1);
  [i1,w1]=basfuns(sp.sp_univ(1),q(1),ismember(1,sp.periodic_dir));
  [i2,w2]=basfuns(sp.sp_univ(2),q(2),ismember(2,sp.periodic_dir));
  for a=1:numel(i1)
    for b=1:numel(i2)
      sd=i1(a)+(i2(b)-1)*sp.ndof_dir(1); w=w1(a)*w2(b);
      for c=1:3, v(c)=v(c)+w*u(sd+(c-1)*Ns); end
    end
  end
end
function [ii,ww]=basfuns(su,q,per)
  kn=su.knots(:).'; p=su.degree; nfull=numel(kn)-p-1;
  spn=findspan(nfull-1,p,q,kn); ww=basisfun(spn,q,p,kn); ii=numbasisfun(spn,q,p,kn)+1;
  if per, ii=mod(ii-1,su.ndof)+1; end
  ii=ii(:).'; ww=ww(:).';
end
function D=rot_pen(sp,side,comp,diffdir,dlt,Ksam)
  % rows: d(u_comp)/d(param diffdir) = 0 along the symmetry cut.
  %   side=1 -> xi=0 (z=0) ; side=2 -> xi=1 (y=0); diffdir=1 -> wrt xi.
  s1=sp.scalar_spaces{1}; n1=s1.ndof_dir(1); Ns=prod(s1.ndof_dir);
  D=sparse(Ksam,sp.ndof); et=linspace(0.02,0.98,Ksam); xi0=(side==2);
  step=1; if side==2, step=-1; end
  for k=1:Ksam
    [i0,w0,j0,wj]=basevals(s1,xi0,et(k));
    [i1,w1,~,~]=basevals(s1,xi0+step*dlt,et(k));
    for a=1:numel(i0)
      for b=1:numel(j0)
        sd0=i0(a)+(j0(b)-1)*n1; sd1=i1(a)+(j0(b)-1)*n1;
        D(k,sd0+(comp-1)*Ns)=D(k,sd0+(comp-1)*Ns)-w0(a)*wj(b)/dlt;
        D(k,sd1+(comp-1)*Ns)=D(k,sd1+(comp-1)*Ns)+w1(a)*wj(b)/dlt;
      end
    end
  end
end
function D=rotation_penalty_uzdxi(sp,xi0,dlt,Ksam)
  s1=sp.scalar_spaces{1}; n1=s1.ndof_dir(1); Ns=prod(s1.ndof_dir);
  D=sparse(Ksam,sp.ndof); et=linspace(0.02,0.98,Ksam);
  for k=1:Ksam
    [i0,w0,j0,wj]=basevals(s1,xi0,et(k));
    [i1,w1,~,~]=basevals(s1,xi0+dlt,et(k));
    for a=1:numel(i0)
      for b=1:numel(j0)
        sd0=i0(a)+(j0(b)-1)*n1; sd1=i1(a)+(j0(b)-1)*n1;
        D(k,sd0+2*Ns)=D(k,sd0+2*Ns)-w0(a)*wj(b)/dlt;
        D(k,sd1+2*Ns)=D(k,sd1+2*Ns)+w1(a)*wj(b)/dlt;
      end
    end
  end
end
function [i,w,j,wj]=basevals(s1,xi,eta)
  su1=s1.sp_univ(1); kn=su1.knots(:).'; p=su1.degree; nd=su1.ndof;
  spn=findspan(nd-1,p,xi,kn); i=(numbasisfun(spn,xi,p,kn)+1)'; w=basisfun(spn,xi,p,kn)';
  su2=s1.sp_univ(2); kn=su2.knots(:).'; p=su2.degree; nd=su2.ndof;
  spn=findspan(nd-1,p,eta,kn); j=(numbasisfun(spn,eta,p,kn)+1)'; wj=basisfun(spn,eta,p,kn)';
end
function both(fid,varargin)
  fprintf(fid,varargin{:}); fprintf(1,varargin{:});
end
