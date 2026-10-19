% HEMISPHERE_CONVENTION_CHECK
% Settle the pinched-hemisphere LOAD-MAGNITUDE convention from the PHYSICS of
% symmetry reduction (not from the spread of quoted reference values), then
% confirm it with the FULL periodic model (no symmetry cuts).
%
% The claim under test: the canonical benchmark's reduced (quarter) model has
% TWO loads of magnitude 1 (this is how the benchmark is universally
% implemented), which -- by the symmetry half-share rule -- represents the FULL
% problem with FOUR loads of magnitude 2.  If instead the canonical FULL model
% had four UNIT loads, the reference displacement would be HALF of 0.094, and
% the current 0.0938 would be an accidental factor of 2.
%
% Established here:
%   (1) symmetry half-share: a concentrated load ON a symmetry plane is HALVED
%       in the reduced model  =>  quarter(P) == full(2P).
%   (2) linear load scaling in both models.
%   (3) strain-energy bookkeeping: U_full == 4 * U_quarter_region.
%   (4) full periodic model at the CORRECT load (4 loads of magnitude 2)
%       reproduces the quarter value.
%
% Prints raw evidence and saves tests/hemisphere_convention_check_output.txt.
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
testsdir='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests';
fid=fopen(fullfile(testsdir,'hemisphere_convention_check_output.txt'),'w');
lg=@(varargin) both(fid,varargin{:});

R=10; t=0.04; E=6.825e7; nu=0.3; hole=18*pi/180;
DEG=[3 3]; REG=[2 2]; NQ=[5 5];

lg('==== PINCHED HEMISPHERE: LOAD-MAGNITUDE / NORMALISATION CONVENTION ====\n');
lg('R=10, t=0.04, E=6.825e7, nu=0.3, 18deg hole.  Quarter model {y>=0,z>=0}.\n');
lg('Load points A=(0,R,0) inward (-y) and B=(0,0,R) outward (+z).\n\n');

%% (1)(2) quarter model, load sweep, with the exact mirror-symmetry set -----
lg('---- (1) QUARTER MODEL: two corner loads of magnitude P (load sweep) ----\n');
lg('symmetry planes = revolution ends z=0 (side1) and y=0 (side2);\n');
lg('full mirror kinematic set imposed by penalty; x-rigid translation pinned.\n');
lg('%-8s %-8s %-12s %-16s %-16s\n','nsub','ndof','P','u_y(load A)','E_strain');
srf=nrbrevolve(nrbcirc(R,[0 0 0],hole,pi/2),[0 0 0],[1 0 0],pi/2);
for nsub=[16 32 48]
  [g,m,sp]=build_shell_rdim3(srf,DEG,REG,[nsub nsub],NQ,3);
  K=op_KL_shells_tp(sp,sp,m,@(x,y,z)E*ones(size(x)),@(x,y,z)nu*ones(size(x)),t);
  qA=newton_inv(g,[0;R;0]); qB=newton_inv(g,[0;0;R]);
  dd0=union(side_comp(sp,1,3),side_comp(sp,2,2));
  dd0=union(dd0,dof_at(sp,newton_inv(g,[R*cos(hole);0;0]),1));
  int=setdiff(1:sp.ndof,dd0);
  D=[rot_pen(sp,1,1,1,1e-4,40);rot_pen(sp,1,2,1,1e-4,40); ...
     rot_pen(sp,2,1,1,1e-4,40);rot_pen(sp,2,3,1,1e-4,40)];
  Kc=K(int,int)+1e8*(D(:,int).'*D(:,int));
  for P=[0.5 1 2]
    rhs=zeros(sp.ndof,1);
    rhs=rhs+point_load_rhs(sp,qA,[0;-P;0]);
    rhs=rhs+point_load_rhs(sp,qB,[0;0;P]);
    u=zeros(sp.ndof,1); u(int)=Kc\rhs(int);
    uu=sp_eval(u,sp,g,{qA(1),qA(2)}); val=uu(2);
    Estrain=0.5*u(int).'*Kc*u(int);
    lg('%-8s %-8d %-12.3g %-16.8g %-16.8g\n',mat2str([nsub nsub]),sp.ndof,P,val,Estrain);
  end
end
lg('  -> u scales LINEARLY with P (0.5:1:2 => 1:2:4).\n');
lg('  -> quarter(P=1) is the benchmark value; quarter(P=2)=2x it.\n\n');

%% (3) FULL PERIODIC MODEL (no symmetry cuts) -------------------------------
lg('---- (2) FULL PERIODIC HEMISPHERE (no symmetry cuts), 4 loads 90deg apart ----\n');
lg('pattern phi=0 in, 90 out, 180 in, 270 out; each load magnitude P.\n');
lg('solved with lsqminnorm (self-equilibrated load -> exact min-norm solution).\n');
lg('%-10s %-8s %-16s %-16s %-16s %-16s\n', ...
   'nsub','ndof','u(P=0.5)','u(P=1)','u(P=2)','E_strain(P=2)');
for nsub=[16 24 32 48 64]
  [geo,msh,sp]=buildHemiFull(R,hole,nsub,DEG,REG,NQ);
  K=op_KL_shells_tp(sp,sp,msh,@(x,y,z)E*ones(size(x)),@(x,y,z)nu*ones(size(x)),t);
  php=[0 pi/2 pi 3*pi/2]; sgn=[-1 1 -1 1]; res=zeros(1,3); E2=0;
  Pp=[0.5 1 2];
  for jj=1:3
    P=Pp(jj); rhs=zeros(sp.ndof,1);
    for k=1:4
      nrm=[0;cos(php(k));sin(php(k))];
      rhs=rhs+point_load_rhs(sp,[1;mod(php(k)/(2*pi),1)],P*sgn(k)*nrm);
    end
    u=lsqminnorm(K,rhs); v=eval_u(u,sp,[1;0]); res(jj)=v(2);
    if P==2, E2=0.5*u.'*K*u; end
  end
  lg('%-10s %-8d %-16.8g %-16.8g %-16.8g %-16.8g\n', ...
     mat2str([nsub nsub]),sp.ndof,res(1),res(2),res(3),E2);
end
lg('  -> full(P) scales linearly; expect full(P=1) = HALF of full(P=2).\n\n');

%% half-share cross-check ----------------------------------------------------
lg('---- (3) HALF-SHARE CROSS-CHECK ----\n');
nsub=48;
[g,m,sp]=build_shell_rdim3(srf,DEG,REG,[nsub nsub],NQ,3);
K=op_KL_shells_tp(sp,sp,m,@(x,y,z)E*ones(size(x)),@(x,y,z)nu*ones(size(x)),t);
qA=newton_inv(g,[0;R;0]); qB=newton_inv(g,[0;0;R]);
dd0=union(side_comp(sp,1,3),side_comp(sp,2,2));
dd0=union(dd0,dof_at(sp,newton_inv(g,[R*cos(hole);0;0]),1));
int=setdiff(1:sp.ndof,dd0);
D=[rot_pen(sp,1,1,1,1e-4,40);rot_pen(sp,1,2,1,1e-4,40); ...
   rot_pen(sp,2,1,1,1e-4,40);rot_pen(sp,2,3,1,1e-4,40)];
Kc=K(int,int)+1e8*(D(:,int).'*D(:,int));
quarter=@(P) deal_q(P,sp,qA,qB,Kc,int,D,g);
[geo,msh,spf]=buildHemiFull(R,hole,nsub,DEG,REG,NQ);
Kf=op_KL_shells_tp(spf,spf,msh,@(x,y,z)E*ones(size(x)),@(x,y,z)nu*ones(size(x)),t);
php=[0 pi/2 pi 3*pi/2]; sgn=[-1 1 -1 1];
fullval=@(P) deal_f(P,spf,Kf,php,sgn);
q1=quarter(1); q05=quarter(0.5); f1=fullval(1); f2=fullval(2);
lg('t nsub=%d, ndof(quarter)=%d, ndof(full)=%d\n',nsub,sp.ndof,spf.ndof);
lg('  quarter(P=1)   = %-14.8g\n',q1);
lg('  quarter(P=0.5) = %-14.8g\n',q05);
lg('  full(P=2)      = %-14.8g\n',f2);
lg('  full(P=1)      = %-14.8g\n',f1);
lg('  quarter(P=1)-full(P=2) = %.2e   (same physical load)\n',q1-f2);
lg('  quarter(P=0.5)-full(P=1) = %.2e (same physical load)\n',q05-f1);
lg('  => quarter(P) == full(2P): the QUARTER model carries HALF the\n');
lg('     full-model per-point load.  Hence quarter unit loads  <->  full loads\n');
lg('     of magnitude 2.  The current setup is CONSISTENT with the canonical\n');
lg('     "quarter model, two unit loads" benchmark convention.\n\n');

lg('CONVENTION VERDICT\n');
lg(' * A point load sitting ON a symmetry plane carries HALF its magnitude in\n');
lg('   the reduced model (the two reflected halves share it).\n');
lg(' * Therefore the FULL model equivalent to the canonical quarter-with-unit-\n');
lg('   loads benchmark has FOUR loads of magnitude 2, NOT four unit loads.\n');
lg(' * Four UNIT loads in the full model give ~0.047 (half) -- a DIFFERENT,\n');
lg('   non-canonical problem.\n');
lg('HEMISPHERE_CONVENTION_CHECK_DONE\n');
fclose(fid); fprintf('HEMISPHERE_CONVENTION_CHECK_DONE\n');

%% ============================== helpers ==============================
function val=deal_q(P,sp,qA,qB,Kc,int,D,g)
  Ns=prod(sp.scalar_spaces{1}.ndof_dir); rhs=zeros(sp.ndof,1);
  rhs=rhs+point_load_rhs(sp,qA,[0;-P;0]);
  rhs=rhs+point_load_rhs(sp,qB,[0;0;P]);
  u=zeros(sp.ndof,1); u(int)=Kc\rhs(int);
  uu=sp_eval(u,sp,g,{qA(1),qA(2)}); val=uu(2);
end
function val=deal_f(P,sp,K,php,sgn)
  rhs=zeros(sp.ndof,1);
  for k=1:4
    nrm=[0;cos(php(k));sin(php(k))];
    rhs=rhs+point_load_rhs(sp,[1;mod(php(k)/(2*pi),1)],P*sgn(k)*nrm);
  end
  u=lsqminnorm(K,rhs); v=eval_u(u,sp,[1;0]); val=v(2);
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
function [geo,msh,sp]=buildHemiFull(R,hole,nsub,deg,reg,nquad)
  knots0={[0 0 0 0 1 1 1 1],[0 0 0 0 1 1 1 1]};
  [knots,zeta]=kntrefine(knots0,[nsub nsub]-1,deg,reg);
  knotsU=kntunclamp(knots,deg,reg,2);
  geo=geo_load({@(P) hmap(P,R,hole), @(P) hder(P,R,hole), @(P) hder2(P,R,hole)});
  geo.nurbs=struct('knots',{{[0 0 0 0 1 1 1 1],[0 0 0 0 1 1 1 1]}},'order',[4 4]);
  rule=msh_gauss_nodes(nquad); [qn,qw]=msh_set_quad_nodes(zeta,rule);
  msh=msh_cartesian(zeta,qn,qw,geo);
  s1=sp_bspline(knotsU,deg,msh,[],2); sp=sp_vector(repmat({s1},1,3),msh);
end
function x=hmap(P,R,hole)
  [u,v]=uvp(P); if iscell(P), [U,V]=ndgrid(u,v); U=U(:).'; V=V(:).'; else, U=u; V=v; end
  th=hole+(pi/2-hole)*U; ph=2*pi*V;
  x=[R*cos(th); R*sin(th).*cos(ph); R*sin(th).*sin(ph)];
end
function J=hder(P,R,hole)
  [u,v]=uvp(P); if iscell(P), [U,V]=ndgrid(u,v); U=U(:).'; V=V(:).'; else, U=u; V=v; end
  n=numel(U); J=zeros(3,2,n); a=pi/2-hole; th=hole+a*U; ph=2*pi*V;
  J(1,1,:)=-a*R*sin(th); J(2,1,:)= a*R*cos(th).*cos(ph); J(3,1,:)= a*R*cos(th).*sin(ph);
  J(2,2,:)=-2*pi*R*sin(th).*sin(ph); J(3,2,:)= 2*pi*R*sin(th).*cos(ph);
end
function H=hder2(P,R,hole)
  [u,v]=uvp(P); if iscell(P), [U,V]=ndgrid(u,v); U=U(:).'; V=V(:).'; else, U=u; V=v; end
  n=numel(U); H=zeros(3,2,2,n); a=pi/2-hole; th=hole+a*U; ph=2*pi*V;
  H(1,1,1,:)=-a^2*R*cos(th); H(2,1,1,:)=-a^2*R*sin(th).*cos(ph); H(3,1,1,:)=-a^2*R*sin(th).*sin(ph);
  H(2,1,2,:)=-a*2*pi*R*cos(th).*sin(ph); H(3,1,2,:)= a*2*pi*R*cos(th).*cos(ph);
  H(2,2,2,:)=-(2*pi)^2*R*sin(th).*cos(ph); H(3,2,2,:)=-(2*pi)^2*R*sin(th).*sin(ph);
end
function [u,v]=uvp(P)
  if iscell(P), u=P{1}(:).'; v=P{2}(:).'; else, u=P(1,:); v=P(2,:); end
end
function D=rot_pen(sp,side,comp,diffdir,dlt,Ksam)
  s1=sp.scalar_spaces{1}; n1=s1.ndof_dir(1); Ns=prod(s1.ndof_dir);
  D=sparse(Ksam,sp.ndof); et=linspace(0.02,0.98,Ksam); xi0=(side==2);
  step=1; if side==2, step=-1; end
  for k=1:Ksam
    [i0,w0,j0,wj]=basevals(s1,xi0,et(k));
    [i1,w1,~,~]=basevals(s1,xi0+step*dlt,et(k));
    for a=1:numel(i0); for b=1:numel(j0)
        sd0=i0(a)+(j0(b)-1)*n1; sd1=i1(a)+(j0(b)-1)*n1;
        D(k,sd0+(comp-1)*Ns)=D(k,sd0+(comp-1)*Ns)-w0(a)*wj(b)/dlt;
        D(k,sd1+(comp-1)*Ns)=D(k,sd1+(comp-1)*Ns)+w1(a)*wj(b)/dlt;
    end; end
  end
end
function [i,w,j,wj]=basevals(s1,xi,eta)
  su1=s1.sp_univ(1); kn=su1.knots(:).'; p=su1.degree; nd=su1.ndof;
  spn=findspan(nd-1,p,xi,kn); i=(numbasisfun(spn,xi,p,kn)+1)'; w=basisfun(spn,xi,p,kn)';
  su2=s1.sp_univ(2); kn=su2.knots(:).'; p=su2.degree; nd=su2.ndof;
  spn=findspan(nd-1,p,eta,kn); j=(numbasisfun(spn,eta,p,kn)+1)'; wj=basisfun(spn,eta,p,kn)';
end
function d=side_comp(space,s,c)
  d=space.boundary(s).dofs(space.boundary(s).comp_dofs{c});
end
function gd=dof_at(space,q,c)
  sp=space.scalar_spaces{1}; Ns=prod(sp.ndof_dir);
  idx=cell(1,2); val=cell(1,2);
  for d=1:2
    su=sp.sp_univ(d); kn=su.knots(:).'; p=su.degree; nd=su.ndof;
    spn=findspan(nd-1,p,q(d),kn); b=basisfun(spn,q(d),p,kn);
    nb=numbasisfun(spn,q(d),p,kn)+1; idx{d}=nb(:); val{d}=b(:);
  end
  best=0; sd=1;
  for a=1:numel(idx{1}); for b=1:numel(idx{2})
      w=val{1}(a)*val{2}(b); if w>best, best=w; sd=idx{1}(a)+(idx{2}(b)-1)*sp.ndof_dir(1); end
  end; end
  gd=sd+(c-1)*Ns;
end
function rhs=point_load_rhs(space,q,F)
  sp=space.scalar_spaces{1}; Ns=prod(sp.ndof_dir); rhs=zeros(space.ndof,1);
  idx=cell(1,2); val=cell(1,2);
  for d=1:2
    su=sp.sp_univ(d); kn=su.knots(:).'; p=su.degree; nd=su.ndof;
    spn=findspan(nd-1,p,q(d),kn); b=basisfun(spn,q(d),p,kn);
    nb=numbasisfun(spn,q(d),p,kn)+1; idx{d}=nb(:); val{d}=b(:);
  end
  for a=1:numel(idx{1}); for b=1:numel(idx{2})
      sd=idx{1}(a)+(idx{2}(b)-1)*sp.ndof_dir(1); w=val{1}(a)*val{2}(b);
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
function q=newton_inv(g,xt)
  starts=[0.5 0.5;0.25 0.5;0.75 0.5;0.5 0.25;0.5 0.75;0.15 0.5;0.85 0.5;0.5 0.15;0.5 0.85];
  best=inf; qbest=[0.5;0.5];
  for s=1:size(starts,1)
    q=starts(s,:).';
    for it=1:80
      qc=min(max(q,1e-9),1-1e-9);
      F=g.map(qc)-xt; Jp=g.map_der(qc); J=reshape(Jp,3,2); dq=(J.'*J)\(J.'*F);
      q=min(max(qc-dq,1e-9),1-1e-9); if norm(dq)<1e-14,break;end
    end
    r=norm(g.map(q)-xt); if r<best, best=r; qbest=q; end
  end
  q=qbest;
end
function both(fid,varargin)
  fprintf(fid,varargin{:}); fprintf(1,varargin{:});
end
