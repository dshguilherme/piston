% PINCHED_SEAM_MEASURE -- direct measurement of the computed solution across the
% parametric seam, plus the closed-ring benchmark.
%
% For a genuinely C^{p-1} closed representation the one-sided tangential
% derivative of the solution must match across the seam (jump -> 0 as the
% sampling step dth -> 0).  For a C0 (clamped, seam-bearing) representation the
% slope jump, if present, is O(1) independent of dth.  We sweep dth and report
% the value / slope / curvature jumps, decomposed in the radial and tangential
% directions, for the C1-closed (periodic) and C0 (nrbcirc) models.
%
% Saves results/pinched_seam_measure.txt.
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
resdir='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
fid=fopen(fullfile(resdir,'pinched_seam_measure.txt'),'w');
lg=@(varargin) both(fid,varargin{:});
DEG=[3 3]; REG=[2 2]; NQ=[5 5];

%% C1-closed periodic cylinder, nsub=64 (converged value ~0.98 ref)
[geo,msh,sp]=buildC1(300,600,64,DEG,REG,NQ);
K=op_KL_shells_tp(sp,sp,msh,@(x,y,z)3e6*ones(size(x)),@(x,y,z)0.3*ones(size(x)),3);
rhs=zeros(sp.ndof,1);
rhs=rhs+point_load_rhs(sp,newton_inv(geo,[300;0;300]),[-1;0;0]);
rhs=rhs+point_load_rhs(sp,newton_inv(geo,[-300;0;300]),[1;0;0]);
dd=union([side_comp(sp,3,1);side_comp(sp,3,2)],[side_comp(sp,4,1);side_comp(sp,4,2)]);
int=setdiff(1:sp.ndof,dd); u1=zeros(sp.ndof,1); u1(int)=K(int,int)\rhs(int);

%% C0 nrbcirc cylinder, nsub=32
nrc=nrbtform(nrbcylind(600,300,[0 0 0],0,2*pi),vecrotz(pi/2));
[g,m,spC]=buildC0(nrc,DEG,REG,[32 32],NQ);
KC=op_KL_shells_tp(spC,spC,m,@(x,y,z)3e6*ones(size(x)),@(x,y,z)0.3*ones(size(x)),3);
rhsC=zeros(spC.ndof,1);
rhsC=rhsC+point_load_rhs(spC,newton_inv(g,[300;0;300]),[-1;0;0]);
rhsC=rhsC+point_load_rhs(spC,newton_inv(g,[-300;0;300]),[1;0;0]);
ddC=union([side_comp(spC,3,1);side_comp(spC,3,2)],[side_comp(spC,4,1);side_comp(spC,4,2)]);
intC=setdiff(1:spC.ndof,ddC); u0=zeros(spC.ndof,1); u0(intC)=KC(intC,intC)\rhsC(intC);

lg('==== SEAM JUMP SWEEP (mid-length circle) ====\n');
lg('C1-closed: seam at theta=0 (periodic, no seam); C0: seam at theta=90deg.\n');
lg('%-12s %-9s %-13s %-13s %-13s %-13s\n','model','dth','val-jump-r','val-jump-t','slope-jump-r','slope-jump-t');
for dth=[1e-2 1e-3 1e-4]
  [vr,vt,sr,st]=seam_sweep(geo,[],sp,u1,300,0,dth);
  lg('%-12s %-9.0e %-13.3e %-13.3e %-13.3e %-13.3e\n','C1-closed',dth,vr,vt,sr,st);
  [vr,vt,sr,st]=seam_sweep([],g,spC,u0,300,pi/2,dth);
  lg('%-12s %-9.0e %-13.3e %-13.3e %-13.3e %-13.3e\n','C0',dth,vr,vt,sr,st);
end
lg('(slope-jump* is the one-sided d(u)/dtheta difference; if it is O(1) as dth->0\n');
lg(' it is a genuine slope discontinuity; if it shrinks with dth it is smooth.)\n\n');

%% CLOSED RING (per-point ref 0.4464).  The closed ring is free-free, so the
%% stiffness is singular (6 rigid modes); the load is self-equilibrated, so we
%% solve the minimum-norm solution.  Note the CONDITIONING: the min-norm solve
%% loses all accuracy once rcond < ~1e-20 (nsub>=48).
lg('==== CLOSED RING (R=100 b=2 t=1 E=1e6 nu=0.3 P=1; per-point ref 0.4464) ====\n');
Rr=100; br=2; tr=1; Er=1e6; nur=0.3; Pr=1; Ir=br*tr^3/12; dref=(0.1488*Pr*Rr^3/(Er*Ir))/2;
lg('%-24s %-8s %-9s %-12s %-14s %-9s\n','model','nsub','ndof','rcond','u_radial(load)','ratio');
for nsub=[16 24 32]
  [geo,msh,sp]=buildC1(Rr,br,nsub,DEG,REG,NQ);
  K=op_KL_shells_tp(sp,sp,msh,@(x,y,z)Er*ones(size(x)),@(x,y,z)nur*ones(size(x)),tr);
  rhs=zeros(sp.ndof,1);
  qA=newton_inv(geo,[Rr;0;br/2]); rhs=rhs+point_load_rhs(sp,qA,[-Pr;0;0]);
  qB=newton_inv(geo,[-Rr;0;br/2]); rhs=rhs+point_load_rhs(sp,qB,[Pr;0;0]);
  rc=1/condest(K); u=lsqminnorm(K,rhs); val=eval_u(u,sp,qA); val=val(1);
  lg('%-24s %-8s %-9d %-12.2e %-14.6g %-9.4f\n','C1-closed (min-norm)',mat2str([nsub nsub]),sp.ndof,rc,val,val/dref);
end
lg('\nPINCHED_SEAM_MEASURE_DONE\n');
fclose(fid); fprintf('PINCHED_SEAM_MEASURE_DONE\n');

%% ============================== helpers ==================================
function [vr,vt,sr,st]=seam_sweep(geo,g,sp,u,R,th0,dth)
  if ~isempty(geo), gg=geo; else, gg=g; end
  U=@(j) uv(gg,sp,u,R,th0+j*dth);
  up1=U(1); up2=U(2); um1=U(-1); um2=U(-2);
  dir=@(th,val) [cos(th) sin(th)]*[val(1);val(2)];   % radial
  tan=@(th,val) [-sin(th) cos(th)]*[val(1);val(2)];  % tangential
  vr=abs(dir(th0,up1)-dir(th0,um1));  vt=abs(tan(th0,up1)-tan(th0,um1));
  sr=abs(dir(th0,up2)-dir(th0,up1)-dir(th0,um1)+dir(th0,um2))/dth;
  st=abs(tan(th0,up2)-tan(th0,up1)-tan(th0,um1)+tan(th0,um2))/dth;
end
function val=uv(geo,sp,u,R,th)
  xt=[R*cos(th);R*sin(th);300]; q=newton_inv(geo,xt); val=eval_u(u,sp,q);
end
function [geo,msh,sp]=buildC1(R,L,nsub,deg,reg,nquad)
  knots0={[0 0 0 0 1 1 1 1],[0 0 0 0 1 1 1 1]};
  [knots,zeta]=kntrefine(knots0,[nsub nsub]-1,deg,reg); knotsU=kntunclamp(knots,deg,reg,1);
  geo=geo_load({@(P) cmap(P,R,L), @(P) cder(P,R,L), @(P) cder2(P,R,L)});
  geo.nurbs=struct('knots',{{[0 0 0 0 1 1 1 1],[0 0 0 0 1 1 1 1]}},'order',[4 4]);
  rule=msh_gauss_nodes(nquad); [qn,qw]=msh_set_quad_nodes(zeta,rule);
  msh=msh_cartesian(zeta,qn,qw,geo); s1=sp_bspline(knotsU,deg,msh,[],1); sp=sp_vector(repmat({s1},1,3),msh);
end
function [geometry,msh,space]=buildC0(nr,degree,regularity,nsub,nquad)
  geometry=geo_load(nr); degelev=max(degree-(geometry.nurbs.order-1),0);
  nurbs=nrbdegelev(geometry.nurbs,degelev);
  [rk,zeta,nk]=kntrefine(nurbs.knots,nsub-1,nurbs.order-1,regularity);
  nurbs=nrbkntins(nurbs,nk); geometry=geo_load(nurbs);
  rule=msh_gauss_nodes(nquad); [qn,qw]=msh_set_quad_nodes(geometry.nurbs.knots,rule);
  msh=msh_cartesian(geometry.nurbs.knots,qn,qw,geometry);
  sp_scalar=sp_nurbs(geometry.nurbs,msh); space=sp_vector(repmat({sp_scalar},1,3),msh);
end
function rhs=point_load_rhs(space,q,F)
  sp=space.scalar_spaces{1}; Ns=prod(sp.ndof_dir); rhs=zeros(space.ndof,1);
  [i1,w1]=basfuns(sp.sp_univ(1),q(1),ismember(1,sp.periodic_dir));
  [i2,w2]=basfuns(sp.sp_univ(2),q(2),ismember(2,sp.periodic_dir));
  for a=1:numel(i1); for b=1:numel(i2)
      sd=i1(a)+(i2(b)-1)*sp.ndof_dir(1); w=w1(a)*w2(b);
      for c=1:3, rhs(sd+(c-1)*Ns)=rhs(sd+(c-1)*Ns)+w*F(c); end
  end; end
end
function v=eval_u(u,space,q)
  sp=space.scalar_spaces{1}; Ns=prod(sp.ndof_dir); v=zeros(3,1);
  [i1,w1]=basfuns(sp.sp_univ(1),q(1),ismember(1,sp.periodic_dir));
  [i2,w2]=basfuns(sp.sp_univ(2),q(2),ismember(2,sp.periodic_dir));
  for a=1:numel(i1); for b=1:numel(i2)
      sd=i1(a)+(i2(b)-1)*sp.ndof_dir(1); w=w1(a)*w2(b);
      for c=1:3, v(c)=v(c)+w*u(sd+(c-1)*Ns); end
  end; end
end
function [ii,ww]=basfuns(su,q,per)
  kn=su.knots(:).'; p=su.degree; nfull=numel(kn)-p-1;
  spn=findspan(nfull-1,p,q,kn); ww=basisfun(spn,q,p,kn); ii=numbasisfun(spn,q,p,kn)+1;
  if per, ii=mod(ii-1,su.ndof)+1; end
  ii=ii(:).'; ww=ww(:).';
end
function d=side_comp(space,s,c); d=space.boundary(s).dofs(space.boundary(s).comp_dofs{c}); end
function [u,v]=uvp(P); if iscell(P), u=P{1}(:).'; v=P{2}(:).'; else, u=P(1,:); v=P(2,:); end; end
function X=cmap(P,R,L)
  [u,v]=uvp(P); if iscell(P), [U,V]=ndgrid(u,v); U=U(:).'; V=V(:).'; else, U=u; V=v; end
  X=[R*cos(2*pi*U); R*sin(2*pi*U); L*V];
end
function J=cder(P,R,L)
  c=iscell(P); [u,v]=uvp(P); if c, [U,V]=ndgrid(u,v); U=U(:).'; V=V(:).'; else, U=u; V=v; end %#ok<NASGU>
  n=numel(U); J=zeros(3,2,n);
  J(1,1,:)=-2*pi*R*sin(2*pi*U); J(2,1,:)=2*pi*R*cos(2*pi*U); J(3,2,:)=L*ones(1,n);
end
function H=cder2(P,R,L)
  c=iscell(P); [u,v]=uvp(P); if c, [U,V]=ndgrid(u,v); U=U(:).'; V=V(:).'; else, U=u; V=v; end %#ok<NASGU>
  n=numel(U); H=zeros(3,2,2,n);
  H(1,1,1,:)=-(2*pi)^2*R*cos(2*pi*U); H(2,1,1,:)=-(2*pi)^2*R*sin(2*pi*U);
end
function q=newton_inv(g,xt)
  starts=[0.5 0.5;0.25 0.5;0.75 0.5;0.5 0.25;0.5 0.75;0.15 0.5;0.85 0.5;0.5 0.15;0.5 0.85];
  best=inf; qbest=[0.5;0.5];
  for s=1:size(starts,1)
    q=starts(s,:).';
    for it=1:80
      qc=min(max(q,1e-9),1-1e-9); F=g.map(qc)-xt; Jp=g.map_der(qc); J=reshape(Jp,3,2);
      dq=(J.'*J)\(J.'*F); q=min(max(qc-dq,1e-9),1-1e-9);
      if norm(dq)<1e-14,break;end
    end
    r=norm(g.map(q)-xt); if r<best, best=r; qbest=q; end
  end
  q=qbest;
end
function both(fid,varargin); fprintf(fid,varargin{:}); fprintf(1,varargin{:}); end
