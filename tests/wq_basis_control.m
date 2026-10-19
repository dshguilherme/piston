% WQ_BASIS_CONTROL -- fifth directive, item (A).
%
% DECISIVE CONTROL for the owner's hypothesis that the WQ operator's PLAIN
% B-SPLINE space (no rational weighting) is structurally wrong on rational NURBS
% geometry.  Run the SAME manufactured-solution WQ energy test on
%   (i)  a genuinely CURVED but POLYNOMIAL (all weights = 1) B-spline patch, and
%   (ii) a RATIONAL twin: the SAME control net with a NON-CONSTANT weight field.
% If the WQ operator converges on (i) but fails on (ii), the rational weight is
% the cause; if it fails on both, the cause is the varying-coefficient/curvature
% effect already documented, and the rational weight is a separate issue.
%
% Report: convergence histories side by side (E_WQ, E_Gauss, E_ref, rel.err).
%
% Saves results/wq_basis_control.txt.
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
resdir='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(resdir,'dir'), mkdir(resdir); end
fid=fopen(fullfile(resdir,'wq_basis_control.txt'),'w');
lg=@(varargin) both(fid,varargin{:});
E=1.0; NU=0.3; RHO=1.0; D=E/(12*(1-NU^2));

% ---- base (NON-RATIONAL, weights = 1) curved patch ----
n0=nrb4surf([0 0],[1 0],[0 1],[1 1]); n0=nrbdegelev(n0,[2 2]);
n0.coefs(1,2,2)=n0.coefs(1,2,2)+0.15; n0.coefs(2,2,2)=n0.coefs(2,2,2)-0.15;
n0.coefs(1,3,3)=n0.coefs(1,3,3)-0.12; n0.coefs(2,3,3)=n0.coefs(2,3,3)+0.12;
n0.coefs(1,2,3)=n0.coefs(1,2,3)+0.06; n0.coefs(2,3,2)=n0.coefs(2,3,2)-0.06;

% ---- rational twin: same net, NON-CONSTANT weights ----
n1=n0; ni=size(n0.coefs,2); nj=size(n0.coefs,3);
[Ii,Jj]=ndgrid(1:ni,1:nj);
n1.coefs(3,:,:) = reshape(1 + 0.35*((Ii-1)/(ni-1)).*((Jj-1)/(nj-1)), 1, ni, nj);
lg('==== WQ basis control: non-rational vs rational curved patch ====\n');
lg('weight field: min=%.3f max=%.3f  (rational twin)\n', min(n1.coefs(3,:)), max(n1.coefs(3,:)));
lg('%-8s %-12s %6s %13s %13s %13s %9s %9s\n','patch','nsub','N','E_WQ','E_Gauss','E_ref','e_WQ','e_Gauss');

for which={'nonrational','rational'}
  if strcmp(which{1},'nonrational'), nrb=n0; else, nrb=n1; end
  for is=[2 4 8 16]
    nsub=[is is];
    pd=struct('geo_name',nrb); md=struct('degree',[3 3],'regularity',[2 2],'nsub',nsub,'nquad',[5 5]);
    [geometry,msh,space]=buildSpaces(pd,md); te=ones(nsub);
    Hw=harmonic_kd_setup_curved(msh,space,geometry,E,NU,RHO,'cpu','double','plate');
    Kb_w=harmonic_kd_form_curved(Hw,te);
    Kb_g=harmonic_kd_gauss_curved(msh,space,geometry,E,NU,'plate',degree_of(space)+3);
    sp1=space.scalar_spaces{1}.sp_univ(1); sp2=space.scalar_spaces{1}.sp_univ(2);
    B1g=scol(sp1.ndof,sp1.degree,msh.qn{1}(:).',sp1.knots);
    B2g=scol(sp2.ndof,sp2.degree,msh.qn{2}(:).',sp2.knots);
    X1=repmat(msh.qn{1}(:),1,numel(msh.qn{2})); X2=repmat(msh.qn{2}(:).',numel(msh.qn{1}),1);
    W=wfun(X1,X2); c=(B1g\W/B2g.'); c=c(:);
    Ew=real(c.'*Kb_w*c); Eg=real(c.'*Kb_g*c);
    Er=dense_ref(msh,geometry,20,D,NU,'plate');
    lg('%-12s %-12s %6d %13.5e %13.5e %13.5e %9.2e %9.2e\n', which{1}, mat2str(nsub), numel(c), Ew, Eg, Er, ...
       abs(Ew-Er)/abs(Er), abs(Eg-Er)/abs(Er));
  end
end
lg('\nINTERPRETATION: compare the e_WQ column across the two patches.\n');
lg('WQ_BASIS_CONTROL_DONE\n');
fclose(fid);
fprintf('WQ_BASIS_CONTROL_DONE\n');

%% helpers
function d=degree_of(space), d=space.scalar_spaces{1}.sp_univ(1).degree; end
function B=scol(ndof,p,pts,kn)
  iv=findspan(ndof-1,p,pts,kn); b=basisfun(iv,pts,p,kn); nb=numbasisfun(iv,pts,p,kn)+1;
  r=repmat((1:numel(pts)).',1,p+1); B=sparse(r(:),nb(:),b(:),numel(pts),ndof);
end
function w=wfun(xi,eta)
  w=0.5*xi.^3-xi.^2.*eta+0.7*xi.*eta.^2+0.3*eta.^3+xi.*eta+0.2*xi+0.1*eta;
end
function [geometry,msh,space]=buildSpaces(pd,md)
  geometry=geo_load(pd.geo_name);
  degelev=max(md.degree-(geometry.nurbs.order-1),0);
  nurbs=nrbdegelev(geometry.nurbs,degelev);
  [~,zeta,nk]=kntrefine(nurbs.knots,md.nsub-1,nurbs.order-1,md.regularity);
  nurbs=nrbkntins(nurbs,nk); geometry=geo_load(nurbs);
  rule=msh_gauss_nodes(md.nquad); [qn,qw]=msh_set_quad_nodes(zeta,rule);
  msh=msh_cartesian(zeta,qn,qw,geometry);
  space=sp_vector(repmat({sp_nurbs(nurbs,msh)},1,msh.rdim),msh);
end
function E=dense_ref(msh,geometry,ng,D,NU,model)
  rule=msh_gauss_nodes(ng);
  [qd1,qw1]=msh_set_quad_nodes(unique(msh.breaks{1}),rule);
  [qd2,qw2]=msh_set_quad_nodes(unique(msh.breaks{2}),rule);
  E=0;
  for el2=1:size(qd2,2), for el1=1:size(qd1,2), for g2=1:ng, for g1=1:ng
    xi=qd1(g1,el1); eta=qd2(g2,el2); wq=qw1(g1,el1)*qw2(g2,el2);
    [~,J]=geometry.map_der([xi;eta]); [~,~,H]=geometry.map_der2([xi;eta]);
    J=squeeze(J); H=squeeze(H);
    dJ=J(1,1)*J(2,2)-J(1,2)*J(2,1); Ji=inv(J); S=inv_map_der2(Ji,H);
    wxi=1.5*xi^2-2*xi*eta+0.7*eta^2+eta+0.2; weta=-xi^2+1.4*xi*eta+0.9*eta^2+xi+0.1;
    wxixi=3*xi-2*eta; wxieta=-2*xi+1.4*eta+1; wetaeta=1.4*xi+1.8*eta;
    bxx=[Ji(1,1)^2;2*Ji(1,1)*Ji(2,1);Ji(2,1)^2;S(1,1,1);S(2,1,1)];
    byy=[Ji(1,2)^2;2*Ji(1,2)*Ji(2,2);Ji(2,2)^2;S(1,2,2);S(2,2,2)];
    bxy=[Ji(1,1)*Ji(1,2);Ji(1,1)*Ji(2,2)+Ji(2,1)*Ji(1,2);Ji(2,1)*Ji(2,2);S(1,1,2);S(2,1,2)];
    Phi=[wxixi;wxieta;wetaeta;wxi;weta]; wxx=bxx.'*Phi; wyy=byy.'*Phi; wxy=bxy.'*Phi;
    if strcmpi(model,'plate'), dens=wxx^2+wyy^2+2*NU*wxx*wyy+2*(1-NU)*wxy^2;
    else, dens=wxx^2+wyy^2+2*wxy^2; end
    E=E+wq*D*dJ*dens;
  end, end, end, end
end
function S=inv_map_der2(Ji,H)
  S=zeros(2,2,2);
  for al=1:2, for i=1:2, for j=1:2
    s=0; for m=1:2, for b=1:2, for c=1:2, s=s-Ji(al,m)*H(m,b,c)*Ji(b,i)*Ji(c,j); end, end, end
    S(al,i,j)=s;
  end, end, end
end
function both(fid,varargin)
  fprintf(fid,varargin{:}); fprintf(1,varargin{:});
end
