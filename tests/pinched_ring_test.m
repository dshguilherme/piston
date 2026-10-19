% PINCHED_RING_TEST -- cheapest decisive isolation of the BENDING term.
%
% A closed circular band (curved beam strip) of radius R, width b, thickness t,
% pinched by two diametral RADIAL point loads P.  Classical Timoshenko curved-
% beam result for the diametral change:
%       Delta_d = 0.1488 * P R^3 / (E I),   I = b t^3 / 12
% (per load point: half of that).  With R=100, b=2, t=1, E=1e6, nu=0.3, P=1:
%       I = 1/6,  Delta_d = 0.8928,  per-point reference = 0.4464.
%
% The previous attempt was numerically SINGULAR (C0 rational-circle seam plus
% unconstrained rigid-body modes -> RCOND 1.8e-16).  Here the seam is removed by
% GLUING the two parametric ends of the revolution direction (identifying the
% duplicate dofs at theta=0 and theta=2pi), and the 6 rigid-body modes are
% removed by a MINIMAL, non-redundant constraint set at 3 well-separated nodes
% (reported below).  Both fixes are required; each is isolated.
%
% Saves results/pinched_ring.txt.
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
resdir = 'C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(resdir,'dir'), mkdir(resdir); end
fid=fopen(fullfile(resdir,'pinched_ring.txt'),'w');
lg=@(varargin) both(fid,varargin{:});

R=100; b=2; tr=1; E=1e6; nur=0.3; P=1;
I=b*tr^3/12; dD=0.1488*P*R^3/(E*I); dref=dD/2;
lg('==== PINCHED RING (curved beam strip) ====\n');
lg('Timoshenko diametral change = 0.1488 P R^3/(E I) = %.6g  -> per-point ref = %.6g\n', dD, dref);
lg('R=%g b=%g t=%g E=%g nu=%g P=%g  I=%g\n\n', R,b,tr,E,nur,P,I);
lg('%-10s %-8s %-9s %-8s %-16s %-10s\n','variant','nsub','ndof','#glue','u_radial(load)','ratio');

for nsub=[24 48 96]
  nr=nrbcylind(b, R, [0 0 0], 0, 2*pi);         % dir1=theta (C0 seam), dir2=z
  [g,m,sp]=build(nr,[3 3],[2 2],[nsub nsub],[5 5],E,nur,tr);
  rhs=zeros(sp.ndof,1);
  rhs=rhs+point_load_rhs(sp,newton_inv(g,[ R;0;b/2]),[-P;0;0]);
  rhs=rhs+point_load_rhs(sp,newton_inv(g,[-R;0;b/2]),[ P;0;0]);

  % ---- seam glue map: dof(i1=n1,i2,c)  <-  dof(i1=1,i2,c) ----
  s1=sp.scalar_spaces{1}; n1=s1.ndof_dir(1); n2=s1.ndof_dir(2); Ns=n1*n2; n=sp.ndof;
  master=(1:n).';
  for i2=1:n2
    for c=1:3
      master(n1+(i2-1)*n1+(c-1)*Ns) = 1+(i2-1)*n1+(c-1)*Ns;
    end
  end
  nglue=sum(master.'~=(1:n));

  % ---- minimal rigid-body constraint set (6 dofs at 3 non-collinear nodes) ----
  %  P1=(0,-R,b/2) theta=270 : u_x,u_y,u_z   (3 translations)
  %  P2=(R,0,b/2)  theta=0   : u_y,u_z       (in-plane rot about z; one tilt)
  %  P3=(0,R,b/2)  theta=90  : u_z           (the other tilt)
  cset=[dof_at(sp,newton_inv(g,[0;-R;b/2]),1); dof_at(sp,newton_inv(g,[0;-R;b/2]),2); ...
        dof_at(sp,newton_inv(g,[0;-R;b/2]),3); ...
        dof_at(sp,newton_inv(g,[ R; 0;b/2]),2); dof_at(sp,newton_inv(g,[ R; 0;b/2]),3); ...
        dof_at(sp,newton_inv(g,[0; R;b/2]),3)];

  % (a) GLUED + constrained
  u=solve_glue(K_of(sp,m,E,nur,tr),rhs,sp,master,cset);
  q=newton_inv(g,[R;0;b/2]); v=sp_eval(u,sp,g,{q(1),q(2)});
  lg('%-10s %-8s %-9d %-8d %-16.6g %-10.5f\n','glued',mat2str([nsub nsub]),sp.ndof,nglue,v(1),v(1)/dref);

  % (b) UNGLUED + constrained (to exhibit the seam singularity) -- RCOND only
  if nsub==24
    Kt=K_of(sp,m,E,nur,tr); freed=setdiff(1:n,cset); Kf=Kt(freed,freed);
    rc=1/condest(Kf);
    lg('%-10s %-8s %-9d %-8d %-16.3e %-10s\n','unglued',mat2str([nsub nsub]),sp.ndof,0,rc,'(rcond=singular)');
  end
end

lg('\nNOTE: glued ratio should be ~1.0 if the bending term (hence the shell\n');
lg('bending operator and the metric) is correct: the closed-form value isolates\n');
lg('the t^3 bending stiffness with no membrane/inextensional coupling.\n');
lg('PINCHED_RING_DONE\n');
fclose(fid);
fprintf('PINCHED_RING_DONE\n');

%% ============================== helpers ==================================
function K=K_of(sp,m,E,nu,t)
  K=op_KL_shells_tp(sp,sp,m,@(x,y,z)E*ones(size(x)),@(x,y,z)nu*ones(size(x)),t);
end
function u=solve_glue(K,rhs,sp,master,extra_fix)
  n=sp.ndof; slave=find(master(:).'~= (1:n));
  freedofs=setdiff(1:n,[slave(:).' extra_fix(:).']);
  keep=false(n,1); keep(freedofs)=true;
  T=sparse(master(keep),1:sum(keep),1,n,sum(keep));
  Kf=T.'*K*T; ff=T.'*rhs; uf=Kf\ff; u=T*uf;
end
function [geometry,msh,space]=build(nr,degree,regularity,nsub,nquad,E,nu,t) %#ok<INUSD>
  geometry=geo_load(nr);
  degelev=max(degree-(geometry.nurbs.order-1),0);
  nurbs=nrbdegelev(geometry.nurbs,degelev);
  [rk,zeta,nk]=kntrefine(nurbs.knots,nsub-1,nurbs.order-1,regularity);
  nurbs=nrbkntins(nurbs,nk); geometry=geo_load(nurbs);
  rule=msh_gauss_nodes(nquad); [qn,qw]=msh_set_quad_nodes(geometry.nurbs.knots,rule);
  msh=msh_cartesian(geometry.nurbs.knots,qn,qw,geometry);
  sp_scalar=sp_nurbs(geometry.nurbs,msh); space=sp_vector(repmat({sp_scalar},1,3),msh);
end
function gd=dof_at(space,q,c)
  sp=space.scalar_spaces{1}; Ns=prod(sp.ndof_dir); idx=cell(1,2); val=cell(1,2);
  for d=1:2
    su=sp.sp_univ(d); kn=su.knots(:).'; p=su.degree; nd=su.ndof;
    spn=findspan(nd-1,p,q(d),kn); bb=basisfun(spn,q(d),p,kn); nb=numbasisfun(spn,q(d),p,kn)+1;
    idx{d}=nb(:); val{d}=bb(:);
  end
  best=0; sd=1;
  for a=1:numel(idx{1})
    for bb=1:numel(idx{2})
      w=val{1}(a)*val{2}(bb); if w>best, best=w; sd=idx{1}(a)+(idx{2}(bb)-1)*sp.ndof_dir(1); end
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
    spn=findspan(nd-1,p,q(d),kn); bb=basisfun(spn,q(d),p,kn); nb=numbasisfun(spn,q(d),p,kn)+1;
    idx{d}=nb(:); val{d}=bb(:);
  end
  Ns=prod(sp.ndof_dir);
  for a=1:numel(idx{1})
    for bb=1:numel(idx{2})
      sd=idx{1}(a)+(idx{2}(bb)-1)*sp.ndof_dir(1); w=val{1}(a)*val{2}(bb);
      for c=1:3, rhs(sd+(c-1)*Ns)=rhs(sd+(c-1)*Ns)+w*F(c); end
    end
  end
end
function both(fid,varargin)
  fprintf(fid,varargin{:}); fprintf(1,varargin{:});
end
