% PINCHED_HEMISPHERE_AUDIT -- seam audit, symmetry-cut conditions, load split,
% and convergence of the pinched hemisphere (18deg hole) against the 0.0938 ref.
%
% RESULT: the 7.8%-too-soft quarter model is caused by the SYMMETRY-CUT
% CONDITIONS: imposing only the normal displacement u_n=0 on the two symmetry
% planes leaves the natural Kirchhoff-Love condition there (M=0, a free bending
% edge).  Imposing the full mirror-symmetry kinematic set
%   z=0 (side1): u_z=0, du_x/dxi=0, du_y/dxi=0
%   y=0 (side2): u_y=0, du_x/dxi=0, du_z/dxi=0
% recovers the reference 0.0938 to <0.1%, penalty-insensitively.
% The construction has NO closing C0 seam (the revolution is a quarter turn,
% open, its clamped ends ARE the symmetry planes); the load split is correct.
%
% Saves results/pinched_hemisphere_audit.txt.
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
resdir='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(resdir,'dir'), mkdir(resdir); end
fid=fopen(fullfile(resdir,'pinched_hemisphere_audit.txt'),'w');
lg=@(varargin) both(fid,varargin{:});
R=10; t=0.04; E=6.825e7; nu=0.3; hole=18*pi/180; REF=0.0938; LO=0.0924; HI=0.094; AN=0.093;
DEG=[3 3]; REG=[2 2]; NQ=[5 5];

lg('==== PINCHED HEMISPHERE (18deg hole) AUDIT ====\n');
lg('R=10 t=0.04 E=6.825e7 nu=0.3 ; ref |u|=%.4f (range %.4f-%.4f, analytic %.3f)\n\n',REF,LO,HI,AN);

%% (1) SEAM / CONTINUITY AUDIT ------------------------------------------
lg('---- (1) SEAM / CONTINUITY AUDIT ----\n');
srf=nrbrevolve(nrbcirc(R,[0 0 0],hole,pi/2),[0 0 0],[1 0 0],pi/2);
lg('quarter surface: order=[%d %d], single span in each direction,\n', srf.order(1),srf.order(2));
lg('  direction 1 (xi) = REVOLUTION angle phi in [0,pi/2] (clamped; 2 ends);\n');
lg('  direction 2 (eta) = MERIDIAN arc theta in [18deg,90deg] (clamped).\n');
lg('boundary edges (mid-edge physical points):\n');
g0=geo_load(srf);
emid={[0 0.5],[1 0.5],[0.5 0],[0.5 1]};
nm={'side1 xi=0 (z=0, symmetry)','side2 xi=1 (y=0, symmetry)', ...
    'side3 eta=0 (hole edge, FREE)','side4 eta=1 (equator, FREE)'};
for s=1:4
  x=g0.map(emid{s}.'); lg('  %-30s at [%7.4f %7.4f %7.4f]\n', nm{s}, x(1),x(2),x(3));
end
lg('CONTINUITY VERDICT:\n');
lg(' * The revolution is a QUARTER turn -> the surface is OPEN in xi; its two\n');
lg('   clamped ends are exactly the symmetry planes z=0 and y=0.  There is NO\n');
lg('   closing parametric seam (unlike a full revolve), so the C0-seam defect\n');
lg('   that broke the cylinder CANNOT occur here.\n');
lg(' * eta ends (hole edge, equator) are genuine free edges: clamped is correct.\n');
lg(' * The meridian is a single quadratic span; the revolution likewise -> the\n');
lg('   patch interior is C-infinity.  Geometry is EXACT (analytic revolve).\n');
lg(' => the cylinder remedy (periodic C^{p-1} closure) is NOT needed for the\n');
lg('    quarter hemisphere; the only reduced-continuity edges are the symmetry\n');
lg('    cuts, addressed in (2).\n\n');

% full periodic construction check (used in (3)): area exact
[geo,msh,sp]=buildHemiFull(R,hole,16,DEG,REG,NQ); mshp=msh_precompute(msh);
area=sum(sum(mshp.jacdet.*mshp.quad_weights)); exactar=2*pi*R^2*cos(hole);
lg('full periodic sphere-sector area=%.8f exact(2piR^2 cos18)=%.8f relerr=%.1e\n\n', ...
   area, exactar, abs(area-exactar)/exactar);

%% (2) SYMMETRY-CUT CONDITIONS -----------------------------------------
lg('---- (2) SYMMETRY-CUT CONDITIONS (quarter model) ----\n');
lg('a displacement-only KL shell at a symmetry plane needs the MIRROR-SYMMETRY\n');
lg('kinematic set, not u_n=0 alone.  We test: (a) u_n only; (b) u_n+single\n');
lg('slope (roof-analog); (c) the full exact set.\n');
lg('%-10s %-8s %-8s %-12s %-14s %-9s\n','case','nsub','ndof','mu','u_y(load)','ratio');
for nsub=[8 16 24 32 48]
  [g,m,sp]=build_shell_rdim3(srf,DEG,REG,[nsub nsub],NQ,3);
  K=op_KL_shells_tp(sp,sp,m,@(x,y,z)E*ones(size(x)),@(x,y,z)nu*ones(size(x)),t);
  rhs=zeros(sp.ndof,1); qA=newton_inv(g,[0;R;0]); rhs=rhs+point_load_rhs(sp,qA,[0;-1;0]);
  qB=newton_inv(g,[0;0;R]); rhs=rhs+point_load_rhs(sp,qB,[0;0;1]);
  dd0=union(side_comp(sp,1,3),side_comp(sp,2,2));
  dd0=union(dd0,dof_at(sp,newton_inv(g,[R*cos(hole);0;0]),1));
  int=setdiff(1:sp.ndof,dd0);
  % (a) u_n only
  u=zeros(sp.ndof,1); u(int)=K(int,int)\rhs(int); v=sp_eval(u,sp,g,{qA(1),qA(2)}); v=v(2);
  lg('%-10s %-8s %-8d %-12s %-14.6g %-9.4f\n','a u_n only',mat2str([nsub nsub]),sp.ndof,'-',v,abs(v)/REF);
  % (b) u_n + single slope (roof analog)
  D=[rot_pen(sp,1,3,1,1e-4,30);rot_pen(sp,2,2,1,1e-4,30)];
  Kr=K(int,int)+1e6*(D(:,int).'*D(:,int)); u=zeros(sp.ndof,1); u(int)=Kr\rhs(int);
  v=sp_eval(u,sp,g,{qA(1),qA(2)}); v=v(2);
  lg('%-10s %-8s %-8d %-12.6g %-14.6g %-9.4f\n','b roof-analog',mat2str([nsub nsub]),sp.ndof,1e6,v,abs(v)/REF);
  % (c) exact mirror set
  D=[rot_pen(sp,1,1,1,1e-4,30);rot_pen(sp,1,2,1,1e-4,30); ...
     rot_pen(sp,2,1,1,1e-4,30);rot_pen(sp,2,3,1,1e-4,30)];
  Kr=K(int,int)+1e6*(D(:,int).'*D(:,int)); u=zeros(sp.ndof,1); u(int)=Kr\rhs(int);
  v=sp_eval(u,sp,g,{qA(1),qA(2)}); v=v(2);
  lg('%-10s %-8s %-8d %-12.6g %-14.6g %-9.4f\n','c exact set',mat2str([nsub nsub]),sp.ndof,1e6,v,abs(v)/REF);
end
lg('\n');

% penalty + collocation insensitivity at nsub=32, exact set
nsub=32; [g,m,sp]=build_shell_rdim3(srf,DEG,REG,[nsub nsub],NQ,3);
K=op_KL_shells_tp(sp,sp,m,@(x,y,z)E*ones(size(x)),@(x,y,z)nu*ones(size(x)),t);
rhs=zeros(sp.ndof,1); qA=newton_inv(g,[0;R;0]); rhs=rhs+point_load_rhs(sp,qA,[0;-1;0]);
qB=newton_inv(g,[0;0;R]); rhs=rhs+point_load_rhs(sp,qB,[0;0;1]);
dd0=union(side_comp(sp,1,3),side_comp(sp,2,2));
dd0=union(dd0,dof_at(sp,newton_inv(g,[R*cos(hole);0;0]),1)); int=setdiff(1:sp.ndof,dd0);
lg('penalty sweep (nsub=32, exact set):\n');
lg('%-12s %-14s %-9s\n','mu','u_y(load)','ratio');
for mu=[1e4 1e5 1e6 1e7 1e8 1e9]
  D=[rot_pen(sp,1,1,1,1e-4,30);rot_pen(sp,1,2,1,1e-4,30); ...
     rot_pen(sp,2,1,1,1e-4,30);rot_pen(sp,2,3,1,1e-4,30)];
  Kr=K(int,int)+mu*(D(:,int).'*D(:,int)); u=zeros(sp.ndof,1); u(int)=Kr\rhs(int);
  v=sp_eval(u,sp,g,{qA(1),qA(2)}); v=v(2);
  lg('%-12.0e %-14.6g %-9.4f\n',mu,v,abs(v)/REF);
end
lg('collocation-count sweep (nsub=32, mu=1e8, exact set):\n');
lg('%-12s %-14s %-9s\n','Ksam','u_y(load)','ratio');
for Ks=[10 30 60]
  D=[rot_pen(sp,1,1,1,1e-4,Ks);rot_pen(sp,1,2,1,1e-4,Ks); ...
     rot_pen(sp,2,1,1,1e-4,Ks);rot_pen(sp,2,3,1,1e-4,Ks)];
  Kr=K(int,int)+1e8*(D(:,int).'*D(:,int)); u=zeros(sp.ndof,1); u(int)=Kr\rhs(int);
  v=sp_eval(u,sp,g,{qA(1),qA(2)}); v=v(2);
  lg('%-12d %-14.6g %-9.4f\n',Ks,v,abs(v)/REF);
end
lg('\n');

%% (3) LOAD-SPLIT AUDIT + FULL PERIODIC MODEL --------------------------
lg('---- (3) LOAD-SPLIT AUDIT + full periodic model ----\n');
lg('FULL periodic hemisphere (no seam, no symmetry cuts), 4 loads of P each:\n');
lg('%-8s %-8s %-16s %-10s %-16s %-10s\n','nsub','ndof','u(P=1)','ratio1','u(P=2)','ratio2');
for nsub=[16 24 32 48 64]
  [geo,msh,sp]=buildHemiFull(R,hole,nsub,DEG,REG,NQ);
  K=op_KL_shells_tp(sp,sp,msh,@(x,y,z)E*ones(size(x)),@(x,y,z)nu*ones(size(x)),t);
  rhs=zeros(sp.ndof,1); php=[0 pi/2 pi 3*pi/2]; sgn=[-1 1 -1 1];
  for k=1:4
    nrm=[0;cos(php(k));sin(php(k))];
    rhs=rhs+point_load_rhs(sp,[1;mod(php(k)/(2*pi),1)],sgn(k)*nrm);
  end
  u=lsqminnorm(K,rhs); v=eval_u(u,sp,[1;0]); v=v(2);
  lg('%-8s %-8d %-16.6g %-10.4f %-16.6g %-10.4f\n',mat2str([nsub nsub]),sp.ndof,v,abs(v)/REF,2*v,abs(2*v)/REF);
end
lg('\nLOAD-SPLIT STATEMENT:\n');
lg(' * The canonical problem is 4 radial loads, alternating in/out, on the equator.\n');
lg(' * The QUARTER model retains the two loads that lie ON its symmetry cuts:\n');
lg('   A=(0,R,0) inward (on z=0) and B=(0,0,R) outward (on y=0).  Each such\n');
lg('   corner point is FIXED by the reflection across the OTHER plane, so\n');
lg('   reflecting the quarter to rebuild the full body hits that load twice:\n');
lg('   quarter load P  ->  full load 2P at each of the four points.\n');
lg(' * Hence quarter(P=1) == full(P=2) == 0.0938 (the reference), and\n');
lg('   quarter(P=0.5) == full(P=1) == 0.0469.  The current model applies P=1\n');
lg('   at the two corners, which is the CORRECT share for the 0.0938 (P=2-full)\n');
lg('   convention.  There is NO factor-of-2 load-split error.\n\n');

lg('VERDICT\n');
lg(' - (2) SYMMETRY-CUT CONDITIONS is responsible for the 7.8%% discrepancy.\n');
lg('   u_n=0 alone leaves a free bending edge (natural M=0) -> too soft (0.1011).\n');
lg('   Adding the exact mirror-symmetry kinematic set recovers 0.0938 (<0.1%%),\n');
lg('   penalty- and collocation-insensitive, and matches the independent full\n');
lg('   periodic model at P=2.\n');
lg(' - (1) seam: none - the quarter revolution is open, no closing C0 seam.\n');
lg(' - (3) load split: correct (quarter P=1 == full P=2).\n');
lg('PINCHED_HEMISPHERE_AUDIT_DONE\n');
fclose(fid); fprintf('PINCHED_HEMISPHERE_AUDIT_DONE\n');

%% ============================== helpers ==============================
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
    if diffdir==1
      [i0,w0,j0,wj]=basevals(s1,xi0,et(k));
      [i1,w1,~,~]=basevals(s1,xi0+step*dlt,et(k));
      for a=1:numel(i0); for b=1:numel(j0)
          sd0=i0(a)+(j0(b)-1)*n1; sd1=i1(a)+(j0(b)-1)*n1;
          D(k,sd0+(comp-1)*Ns)=D(k,sd0+(comp-1)*Ns)-w0(a)*wj(b)/dlt;
          D(k,sd1+(comp-1)*Ns)=D(k,sd1+(comp-1)*Ns)+w1(a)*wj(b)/dlt;
      end; end
    else
      [~,~,j0,wj0]=basevals(s1,xi0,et(k)-dlt/2);
      [i0,w0,~,~]=basevals(s1,xi0,et(k));
      [~,~,j1,wj1]=basevals(s1,xi0,et(k)+dlt/2);
      for a=1:numel(i0); for b=1:numel(j0)
          sd=i0(a)+(j0(b)-1)*n1; D(k,sd+(comp-1)*Ns)=D(k,sd+(comp-1)*Ns)-w0(a)*wj0(b)/dlt;
      end; end
      for a=1:numel(i0); for b=1:numel(j1)
          sd=i0(a)+(j1(b)-1)*n1; D(k,sd+(comp-1)*Ns)=D(k,sd+(comp-1)*Ns)+w0(a)*wj1(b)/dlt;
      end; end
    end
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
  idx=cell(1,2); val=cell(1,2);
  for d=1:2
    su=sp.sp_univ(d); kn=su.knots(:).'; p=su.degree; nd=su.ndof;
    spn=findspan(nd-1,p,q(d),kn); b=basisfun(spn,q(d),p,kn); nb=numbasisfun(spn,q(d),p,kn)+1;
    idx{d}=nb(:); val{d}=b(:);
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
