function ctx = obstacle_build_ctx(kind, spec)
% OBSTACLE_BUILD_CTX  Generic context builder for the Article B OBSTACLE COURSE.
%
% Builds a Kirchhoff-Love curved-shell harmonic context (the SAME field
% contract as tests/articleB_build_ctx.m: Bke, Ske, Me, lm, ne, ndof, Rows,
% Cols, free, dr, F, w, zeta, alpha, beta, E, RHO, t0, nsub, space, msh,
% geometry) but for geometry supplied ANALYTICALLY as a cell-array map (so the
% representation is POLYNOMIAL / non-rational, and the fast-formation machinery
% stays valid).  Closed directions are made PERIODIC with kntunclamp + sp_bspline
% (tests/pinched_c1_closed.m), not with a clamped rational circle.
%
%   kind = 'spherecap'  spherical cap, periodic in theta
%   kind = 'cylshell'   closed cylindrical shell, periodic in theta
%   kind = 'cylpanel'   open cylindrical panel (finite opening)
%   kind = 'scordelis'  Scordelis-Lo roof (delegates to articleB_build_ctx)
%
% spec fields (kind-dependent):
%   common : E NU RHO t0 deg reg nsub nquad free
%   spherecap : R alpha   (alpha = cap half-angle, rad)
%   cylshell  : R L
%   cylpanel  : R L th0   (th0 = half opening angle, rad)
%   load      : loadtype 'pole' (spherecap) or 'point'; for 'point', loaduv
%               (param coords, 1x2) and loaddir (3x1 unit force); loadmag.
%   bc        : cell array of side indices whose ALL components are clamped; or
%               for 'freefree' pass bc={} (no Dirichlet).
%
% C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests/

addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');

if strcmpi(kind,'scordelis')
  ctx = articleB_build_ctx(spec.nsub, spec.E, spec.NU, spec.RHO, spec.t0, spec.freq, spec.zeta);
  ctx.kind = 'scordelis';
  return;
end

deg=spec.deg; reg=spec.reg; nquad=spec.nquad;
ns=spec.nsub; if isscalar(ns), ns=repmat(ns,1,2); end
nsub=ns;
knots0={[0 0 0 0 1 1 1 1],[0 0 0 0 1 1 1 1]};
[knots,zeta]=kntrefine(knots0,ns-1,deg,reg);

switch lower(kind)
  case 'spherecap'
    R=spec.R; al=spec.alpha;
    geo=geo_load({@(P) sph_map(P,R,al), @(P) sph_der(P,R,al), @(P) sph_der2(P,R,al)});
    per_dir=2;
  case 'cylshell'
    R=spec.R; L=spec.L;
    geo=geo_load({@(P) cylc_map(P,R,L), @(P) cylc_der(P,R,L), @(P) cylc_der2(P,R,L)});
    per_dir=1; % theta=2*pi*U; the axial V direction is open
  case 'cylpanel'
    R=spec.R; L=spec.L; th0=spec.th0;
    geo=geo_load({@(P) cylp_map(P,R,L,th0), @(P) cylp_der(P,R,L,th0), @(P) cylp_der2(P,R,L,th0)});
    per_dir=[];   %   direction 2 = axial, clamped
  otherwise
    error('obstacle_build_ctx: unknown kind %s', kind);
end
geo.nurbs=struct('knots',{knots0},'order',deg+1);
rule=msh_gauss_nodes(nquad); [qn,qw]=msh_set_quad_nodes(zeta,rule);
msh=msh_cartesian(zeta,qn,qw,geo);
if isempty(per_dir)
  knotsU=knots;
else
  knotsU=kntunclamp(knots,deg,reg,per_dir);
end
s1=sp_bspline(knotsU,deg,msh,[],per_dir);
space=sp_vector(repmat({s1},1,3),msh);

msh1 = msh_precompute(msh);
sp1  = sp_precompute_param(space,msh1,'value',true,'gradient',true,'hessian',true);
for idim=1:msh1.rdim, x{idim}=reshape(msh1.geo_map(idim,:,:),msh1.nqn,msh1.nel); end
[Bke,Ske,~]=op_KL_shells_elements(sp1,sp1,msh1,spec.E*ones(size(x{1})),spec.NU*ones(size(x{1})));
Me = op_u_v_elements(sp1,sp1,msh1,spec.RHO);
lm = sp1.connectivity.';
ne = msh1.nel; n = size(lm,2); ndof = space.ndof;

Rows = zeros(ne, n*n); Cols = zeros(ne, n*n);
for a=1:n
  for b=1:n
    Rows(:,(a-1)*n+b) = lm(:,a);
    Cols(:,(a-1)*n+b) = lm(:,b);
  end
end

% ---- Dirichlet ----
dr=[];
if isfield(spec,'bc') && ~isempty(spec.bc)
  for s=spec.bc
    for c=1:3
      dr=union(dr, space.boundary(s).dofs(space.boundary(s).comp_dofs{c}));
    end
  end
end
free = setdiff(1:ndof,dr); free=free(:).';

% ---- excitation ----
F = build_load(space, spec);

ctx = struct();
ctx.kind=kind; ctx.Bke=Bke; ctx.Ske=Ske; ctx.Me=Me; ctx.lm=lm; ctx.ne=ne; ctx.ndof=ndof;
ctx.Rows=Rows; ctx.Cols=Cols; ctx.free=free; ctx.dr=dr; ctx.F=F;
ctx.E=spec.E; ctx.RHO=spec.RHO; ctx.t0=spec.t0; ctx.nsub=nsub;
ctx.space=space; ctx.msh=msh1; ctx.geometry=geo;
if isfield(spec,'freq') && ~isempty(spec.freq) && ~isnan(spec.freq)
  ctx.w=2*pi*spec.freq; ctx.zeta=spec.zeta;
  ctx.alpha=2*spec.zeta*ctx.w; ctx.beta=0;
else
  ctx.w=NaN; ctx.zeta=0; ctx.alpha=0; ctx.beta=0;
end
ctx.R0=[];
ctx.spec=spec;
end

% ============================ load builder ================================
function F = build_load(space,spec)
  sp=space.scalar_spaces{1}; Ns=prod(sp.ndof_dir); F=zeros(space.ndof,1);
  mag = 1; if isfield(spec,'loadmag'), mag=spec.loadmag; end
  switch lower(spec.loadtype)
    case 'point'
      q=spec.loaduv(:).'; fd=spec.loaddir(:)/norm(spec.loaddir)*mag;
      [i1,w1]=basfuns(sp.sp_univ(1),q(1),ismember(1,sp.periodic_dir));
      [i2,w2]=basfuns(sp.sp_univ(2),q(2),ismember(2,sp.periodic_dir));
      for a=1:numel(i1); for b=1:numel(i2)
          sd=i1(a)+(i2(b)-1)*sp.ndof_dir(1); w=w1(a)*w2(b);
          for c=1:3, F(sd+(c-1)*Ns)=F(sd+(c-1)*Ns)+w*fd(c); end
      end; end
    case 'pole'
      % axisymmetric point force at the apex (u=0): force shared equally over
      % the theta basis functions of the u=0 row, so the load is symmetric.
      fd=spec.loaddir(:)/norm(spec.loaddir)*mag;
      nv=sp.ndof_dir(2); nu=sp.ndof_dir(1);
      if ismember(2,sp.periodic_dir), nvth=sp.ndof_dir(2); else, nvth=nv; end
      i1=basfuns(sp.sp_univ(1),0,false); i1=i1(:).';  % u=0 basis (single for clamped)
      for b=1:nvth
        sd=i1+(b-1)*nu; w=1/nvth;
        for c=1:3, F(sd+(c-1)*Ns)=F(sd+(c-1)*Ns)+w*fd(c); end
      end
    otherwise
      error('build_load: unknown loadtype');
  end
end

% ============================ helpers =====================================
function [ii,ww]=basfuns(su,q,per)
  kn=su.knots(:).'; p=su.degree; nfull=numel(kn)-p-1;
  spn=findspan(nfull-1,p,q,kn); ww=basisfun(spn,q,p,kn); ii=numbasisfun(spn,q,p,kn)+1;
  if per, ii=mod(ii-1,su.ndof)+1; end
  ii=ii(:).'; ww=ww(:).';
end
function [u,v]=uvp(P)
  if iscell(P), u=P{1}(:).'; v=P{2}(:).'; else, u=P(1,:); v=P(2,:); end
end
function [U,V]=grid(P)
  c=iscell(P);[u,v]=uvp(P);if c,[U,V]=ndgrid(u,v);U=U(:).';V=V(:).';else,U=u;V=v;end
end
% ---- sphere cap (theta periodic) ----
function X=sph_map(P,R,al)
  [U,V]=grid(P); phi=al*U; th=2*pi*V;
  X=[R*sin(phi).*cos(th); R*sin(phi).*sin(th); R*cos(phi)];
end
function J=sph_der(P,R,al)
  [U,V]=grid(P); phi=al*U; th=2*pi*V; n=numel(U); J=zeros(3,2,n);
  J(1,1,:)=R*al*cos(phi).*cos(th); J(2,1,:)=R*al*cos(phi).*sin(th); J(3,1,:)=-R*al*sin(phi);
  J(1,2,:)=-R*2*pi*sin(phi).*sin(th); J(2,2,:)=R*2*pi*sin(phi).*cos(th); J(3,2,:)=0;
end
function H=sph_der2(P,R,al)
  [U,V]=grid(P); phi=al*U; th=2*pi*V; n=numel(U); H=zeros(3,2,2,n);
  H(1,1,1,:)=-R*al^2*sin(phi).*cos(th); H(2,1,1,:)=-R*al^2*sin(phi).*sin(th); H(3,1,1,:)=-R*al^2*cos(phi);
  H(1,2,2,:)=-R*(2*pi)^2*sin(phi).*cos(th); H(2,2,2,:)=-R*(2*pi)^2*sin(phi).*sin(th);
  H(1,1,2,:)=-R*al*2*pi*cos(phi).*sin(th); H(2,1,2,:)=R*al*2*pi*cos(phi).*cos(th);
  H(1,2,1,:)=H(1,1,2,:); H(2,2,1,:)=H(2,1,2,:);
end
% ---- closed cylinder (theta periodic) ----
function X=cylc_map(P,R,L)
  [U,V]=grid(P); th=2*pi*U; X=[R*cos(th); R*sin(th); L*V];
end
function J=cylc_der(P,R,L)
  [U,V]=grid(P); th=2*pi*U; n=numel(U); J=zeros(3,2,n);
  J(1,1,:)=-2*pi*R*sin(th); J(2,1,:)=2*pi*R*cos(th); J(3,2,:)=L*ones(1,n);
end
function H=cylc_der2(P,R,L)
  [U,V]=grid(P); th=2*pi*U; n=numel(U); H=zeros(3,2,2,n);
  H(1,1,1,:)=-(2*pi)^2*R*cos(th); H(2,1,1,:)=-(2*pi)^2*R*sin(th);
end
% ---- open cylindrical panel ----
function X=cylp_map(P,R,L,th0)
  [U,V]=grid(P); ang=th0*(2*U-1); X=[R*sin(ang); R*cos(ang); L*V];
end
function J=cylp_der(P,R,L,th0)
  [U,V]=grid(P); ang=th0*(2*U-1); n=numel(U); J=zeros(3,2,n);
  J(1,1,:)=R*2*th0*cos(ang); J(2,1,:)=-R*2*th0*sin(ang); J(3,2,:)=L*ones(1,n);
end
function H=cylp_der2(P,R,L,th0)
  [U,V]=grid(P); ang=th0*(2*U-1); n=numel(U); H=zeros(3,2,2,n);
  H(1,1,1,:)=-R*(2*th0)^2*sin(ang); H(2,1,1,:)=-R*(2*th0)^2*cos(ang);
end
