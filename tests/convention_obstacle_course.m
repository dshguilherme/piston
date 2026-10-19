% CONVENTION_OBSTACLE_COURSE
% Re-run the classical thin-shell obstacle course (Belytschko/Wong/Stolarski
% 1989; problems from MacNeal & Harder 1985) with CORRECT, self-consistent load
% conventions and normalisations, each carrying its citation.  This is a
% convention-corrected companion to tests/run_obstacle_course.m; it exists
% because the pinched-hemisphere reference in that script (0.0938) was the
% computed value itself (self-referential: 0.0162% is not a real error against
% any published benchmark).
%
% LOAD-CONVENTION STATEMENT (pinched hemisphere) -- see
% tests/hemisphere_convention_check.m for the raw proof:
%   * Canonical benchmark = QUARTER model {y>=0,z>=0} with TWO unit point loads
%     at the equator corners (MacNeal-Harder 1985 Fig.9; Bentley STAAD.Pro V20
%     verification model applies "two unit forces on the quadrant", theory
%     0.094; DTIC AD-A369167 applies +/-1.0 on the quarter, theory 0.0940).
%   * A point load sitting ON a symmetry plane is HALVED in the reduced model
%     (the reflected halves share it); equivalently the FULL model carries FOUR
%     loads of magnitude 2.  (COMSOL 6.4 verification doc: full = 4 x 200 N ->
%     quarter = 2 x 100 N.)  Four UNIT loads in the FULL model would give ~0.047,
%     a DIFFERENT problem.  The quarter-with-unit-loads setup below is correct.
%
% Saves results/articleA_convention.txt (and refreshes results/obstacle_course.txt).
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
resdir = 'C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(resdir,'dir'), mkdir(resdir); end
fid=fopen(fullfile(resdir,'articleA_convention.txt'),'w');
lg=@(varargin) both(fid,varargin{:});

lg('==== Shell obstacle course: geopdes KL shell (Gauss) -- CORRECTED CONVENTIONS ====\n');
lg('problems after MacNeal & Harder, Finite Elem. Anal. Des. 1 (1985) 3-20, as\n');
lg('restated in Belytschko, Wong & Stolarski, Int. J. Numer. Meth. Engng 28 (1989) 385-414.\n\n');

%% ------------------------------------------------------------------ (1) ---
% PINCHED CYLINDER: R=300, L=600, t=3, E=3e6, nu=0.3; two diametrically
% opposite UNIT radial point loads at mid-length; rigid end diaphragms.
% Reference |u| = 1.8248e-5 (MacNeal & Harder 1985; Belytschko et al. 1989).
R=300; L=600; t=3; E=3.0e6; nu=0.3;
lg('-- (1) PINCHED CYLINDER (unit loads, R=300,L=600,t=3,E=3e6,nu=0.3) --\n');
lg('   ref |u| = 1.8248e-5  [MacNeal & Harder 1985; Belytschko/Wong/Stolarski 1989]\n');
lg('%-26s %-8s %-9s %-14s %-13s %-9s\n','model','nsub','ndof','u_value','ref','rel.err');
for nsub=[48 64 96]
  [g,m,sp]=buildC1(R,L,nsub,[3 3],[2 2],[5 5]);
  K=op_KL_shells_tp(sp,sp,m,@(x,y,z)E*ones(size(x)),@(x,y,z)nu*ones(size(x)),t);
  rhs=zeros(sp.ndof,1);
  rhs=rhs+point_load_rhs(sp,newton_inv(g,[R;0;L/2]),[-1;0;0]);
  rhs=rhs+point_load_rhs(sp,newton_inv(g,[-R;0;L/2]),[ 1;0;0]);
  dd=union([side_comp(sp,3,1);side_comp(sp,3,2)],[side_comp(sp,4,1);side_comp(sp,4,2)]);
  int=setdiff(1:sp.ndof,dd); u=zeros(sp.ndof,1); u(int)=K(int,int)\rhs(int);
  q=newton_inv(g,[R;0;L/2]); val=eval_u(u,sp,q); val=val(1);
  lg('%-26s %-8s %-9d %-14.6g %-13.6g %-9.2e\n', ...
     'C1-closed(periodic)',mat2str([nsub nsub]),sp.ndof,val,1.8248e-5, ...
     abs(abs(val)-1.8248e-5)/1.8248e-5);
end
lg('\n');

%% ------------------------------------------------------------------ (2) ---
% PINCHED HEMISPHERE (18 deg hole): QUARTER model, two UNIT corner loads.
Rm=10; tm=0.04; Em=6.825e7; num=0.3; hole=18*pi/180;
srf=nrbrevolve(nrbcirc(Rm,[0 0 0],hole,pi/2),[0 0 0],[1 0 0],pi/2);   % {y>=0,z>=0}
lg('-- (2) PINCHED HEMISPHERE, 18deg hole (quarter model, TWO unit corner loads) --\n');
lg('   primary ref 0.0940  [MacNeal & Harder 1985, Table 5(b), open-top hemisphere]\n');
lg('   alt refs: 0.093 [Simo/Fox/Rifai, CMAME 73 (1989) 53-92, sec 7.5.1, 18deg hole];\n');
lg('             0.093597 [Krysl & Chen, Arch. Comput. Methods Eng. 29 (2022)]\n');
lg('   (Note: 0.0924 [Morley & Morris 1978] is the CLOSED sphere, NOT the 18deg-hole model.)\n');
lg('%-26s %-8s %-9s %-14s %-13s %-9s\n','model','nsub','ndof','u_value','ref(0.0940)','rel.err');
for nsub=[16 32 48]
  [g,m,sp]=build_shell_rdim3(srf,[3 3],[2 2],[nsub nsub],[5 5],3);
  K=op_KL_shells_tp(sp,sp,m,@(x,y,z)Em*ones(size(x)),@(x,y,z)num*ones(size(x)),tm);
  rhs=zeros(sp.ndof,1);
  qA=newton_inv(g,[0;Rm;0]); rhs=rhs+point_load_rhs(sp,qA,[0;-1;0]);
  qB=newton_inv(g,[0;0;Rm]); rhs=rhs+point_load_rhs(sp,qB,[0; 0;1]);
  dd=union(side_comp(sp,1,3),side_comp(sp,2,2));
  dd=union(dd,dof_at(sp,newton_inv(g,[Rm*cos(hole);0;0]),1));
  int=setdiff(1:sp.ndof,dd);
  D=[rot_pen(sp,1,1,1,1e-4,40);rot_pen(sp,1,2,1,1e-4,40); ...
     rot_pen(sp,2,1,1,1e-4,40);rot_pen(sp,2,3,1,1e-4,40)];
  Kr=K(int,int)+1e8*(D(:,int).'*D(:,int)); u=zeros(sp.ndof,1); u(int)=Kr\rhs(int);
  uu=sp_eval(u,sp,g,{qA(1),qA(2)}); val=uu(2);
  lg('%-26s %-8s %-9d %-14.6g %-13.6g %-9.2e\n', ...
     'quarter{full sym}',mat2str([nsub nsub]),sp.ndof,val,0.0940, ...
     abs(abs(val)-0.0940)/0.0940);
end
% independent FULL periodic model, no symmetry cuts, FOUR loads of magnitude 2
% (the full-model equivalent of the canonical quarter-with-unit-loads setup)
lg('%-26s %-8s %-9s %-14s %-13s %-9s\n','full(4 loads, mag 2)','nsub','ndof','u_value','ref(0.0940)','rel.err');
for nsub=[32 48 64]
  [geo,msh,sp]=buildHemiFull(Rm,hole,nsub,[3 3],[2 2],[5 5]);
  K=op_KL_shells_tp(sp,sp,msh,@(x,y,z)Em*ones(size(x)),@(x,y,z)num*ones(size(x)),tm);
  rhs=zeros(sp.ndof,1); php=[0 pi/2 pi 3*pi/2]; sgn=[-1 1 -1 1];
  for k=1:4
    nrm=[0;cos(php(k));sin(php(k))];
    rhs=rhs+point_load_rhs(sp,[1;mod(php(k)/(2*pi),1)],2*sgn(k)*nrm);
  end
  u=lsqminnorm(K,rhs); v=eval_u(u,sp,[1;0]); val=v(2);
  lg('%-26s %-8s %-9d %-14.6g %-13.6g %-9.2e\n', ...
     'full(no cuts)',mat2str([nsub nsub]),sp.ndof,val,0.0940,abs(abs(val)-0.0940)/0.0940);
end
lg('\n');

%% ------------------------------------------------------------------ (3) ---
% SCORDELIS-LO ROOF: R=25, L=50, 80deg, t=0.25, E=4.32e8, nu=0; self weight
% g=90 per unit area; rigid end diaphragms.  Ref -0.3024 [MacNeal & Harder 1985].
% (Krysl & Chen 2022 argue the shear-rigid/KL limit is ~0.3006.)
Rr=25; Lr=50; Eu=4.32e8; nuu=0.0; tu=0.25;
probs=obstacle_course_geom('roof'); nrf=probs(1).nurbs;
lg('-- (3) SCORDELIS-LO ROOF (self-weight g=90, R=25,L=50,t=0.25,E=4.32e8,nu=0) --\n');
lg('   ref -0.3024  [MacNeal & Harder 1985];  KL-limit ~-0.3006 [Krysl & Chen 2022]\n');
lg('%-26s %-8s %-9s %-14s %-13s %-9s\n','model','nsub','ndof','u_value','ref(-0.3024)','rel.err');
for nsub=[20 30]
  [g,m,sp]=build_shell_rdim3(nrf,[3 3],[2 2],[nsub nsub],[5 5],3);
  K=op_KL_shells_tp(sp,sp,m,@(x,y,z)Eu*ones(size(x)),@(x,y,z)nuu*ones(size(x)),tu);
  rhs=op_f_v_tp(sp,m,@(x,y,z,ind) -90*ones(3,numel(x)).*[0;0;1]);
  dd=union([side_comp(sp,3,1);side_comp(sp,3,3)],[side_comp(sp,4,1);side_comp(sp,4,3)]);
  int=setdiff(1:sp.ndof,dd); u=zeros(sp.ndof,1); u(int)=K(int,int)\rhs(int);
  th=2/9*pi; q=newton_inv(g,[Rr*cos(pi/2-th);Lr/2;Rr*sin(pi/2-th)]);
  val=sp_eval(u,sp,g,{q(1),q(2)}); val=val(3);
  lg('%-26s %-8s %-9d %-14.6g %-13.6g %-9.2e\n', ...
     'full',mat2str([nsub nsub]),sp.ndof,val,-0.3024,abs(abs(val)-0.3024)/0.3024);
end
lg('\n');

lg('NORMALISATION / LOAD-CONVENTION NOTES (self-consistent across all three):\n');
lg(' - CYLINDER: two opposing UNIT point loads at mid-length; u = radial\n');
lg('   component at the load point; ref 1.8248e-5 [MH85; BWS89].\n');
lg(' - HEMISPHERE: canonical QUARTER model with TWO unit corner loads (MH85\n');
lg('   Fig.9; STAAD.Pro V20 verification, theory 0.094).  Half-share: a load on\n');
lg('   a symmetry plane is halved, so this equals the FULL model with FOUR\n');
lg('   loads of magnitude 2 (COMSOL 6.4, 200 N -> 100 N).  ref 0.0940 (open-top)\n');
lg('   / 0.093 [SFR89] / 0.0936 [KC22].  The full row above reproduces the\n');
lg('   quarter value with NO symmetry reduction, at the equivalent load.\n');
lg(' - ROOF: distributed self-weight 90/unit area; u_z at the mid free edge;\n');
lg('   ref -0.3024 [MH85].\n');
lg('CONVENTION_OBSTACLE_COURSE_DONE\n');
fclose(fid);

% refresh the primary obstacle-course results file with the corrected table
% (read the articleA file we just wrote and copy it verbatim)
src=fopen(fullfile(resdir,'articleA_convention.txt'),'r'); txt=fread(src,'*char'); fclose(src);
dst=fopen(fullfile(resdir,'obstacle_course.txt'),'w'); fwrite(dst,txt); fclose(dst);
fprintf('CONVENTION_OBSTACLE_COURSE_DONE\n');

% =============================== helpers =================================
function [geo,msh,sp]=buildC1(R,L,nsub,deg,reg,nquad)
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
  sp=space.scalar_spaces{1}; Ns=prod(sp.ndof_dir); idx=cell(1,2); val=cell(1,2);
  for d=1:2
    su=sp.sp_univ(d); kn=su.knots(:).'; p=su.degree; nd=su.ndof;
    spn=findspan(nd-1,p,q(d),kn); b=basisfun(spn,q(d),p,kn); nb=numbasisfun(spn,q(d),p,kn)+1;
    idx{d}=nb(:); val{d}=b(:);
  end
  best=0; sd=1;
  for a=1:numel(idx{1}); for b=1:numel(idx{2})
      w=val{1}(a)*val{2}(b); if w>best, best=w; sd=idx{1}(a)+(idx{2}(b)-1)*sp.ndof_dir(1); end
  end; end
  gd=sd+(c-1)*Ns;
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
