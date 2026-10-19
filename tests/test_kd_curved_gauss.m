% VERIFICATION + DEMONSTRATION of the CONSISTENT REMEDIES for the reduced-WQ
% failure on curvilinear geometry.
%   METHOD 1 (baseline) : harmonic_kd_setup_curved + harmonic_kd_form_curved
%                         (reduced weighted quadrature -- the object that fails)
%   METHOD 2            : harmonic_kd_gauss_curved  (per-element Gauss, metric
%                         and coefficient evaluated at the Gauss points)
%   METHOD 3            : harmonic_kd_setup_cmoment + harmonic_kd_form_curved
%                         (C-weighted / moment-enriched rule, exact for affine C)
%
% PART 0: 1D exactness of both rules for an AFFINE coefficient.
% PART A: affine-map regression: Methods 2,3 vs the existing reduced-WQ baseline.
% PART B: curved NURBS chart (manufactured cubic): convergence + order.
% PART C: cost (setup+form wall time, nominal ops/element) vs achieved accuracy.
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
E=1.0; NU=0.3; RHO=1.0; D=E/(12*(1-NU^2));

%% PART 0 : 1D rule exactness for affine C = 1 + 0.3 x
fprintf('==== PART 0: 1D rule exactness for AFFINE C ====\n');
p=3; nsub=8; kn=[zeros(1,p) linspace(0,1,nsub+1) ones(1,p)];
qn=msh_set_quad_nodes(unique(kn),msh_gauss_nodes(p+1)); sp=sp_bspline_1d_param(kn,p,qn);
Qr=quadrule_bend_fast(sp); Qc=quadrule_bend_cmoment(sp); ndof=sp.ndof;
fprintf('  points/DOF max:  reduced=%d   cmoment=%d\n', max(Qr.nquad_points), max(Qc.nquad_points));
ng=24; dk=unique(kn); [xf,wf]=msh_set_quad_nodes(dk,msh_gauss_nodes(ng)); xf=xf(:).'; wf=wf(:).';
Bf0=scol(ndof,p,xf,kn); [Bf1,Bf2]=dcol(ndof,p,xf,kn,Qr.D1,Qr.D2); Btf={Bf0,Bf1,Bf2};
types=[0 0;1 0;1 1;2 1;2 2];
for q=1:2
  if q==1, Q=Qr; nm='reduced'; else, Q=Qc; nm='cmoment'; end
  er=0;
  for it=1:size(types,1)
    k1=types(it,1); k2=types(it,2);
    for a=1:ndof
      w=Q.(sprintf('quad_weights_%d%d',k1,k2)){a}; if isempty(w), continue; end
      ip=Q.ind_points{a}; xq=Q.all_points(ip); nbrs=unique(sp.connectivity(:,sp.supp{a}));
      for b=nbrs(:).'
        lhs=w*((1+0.3*xq).'.*full(Qcol(Q,k2,ip,b)));
        rhs=wf*(full(Btf{k1+1}(:,a)).*(1+0.3*xf).*full(Btf{k2+1}(:,b)));
        if abs(rhs)>1e-10, er=max(er,abs(lhs-rhs)/abs(rhs)); end
      end
    end
  end
  fprintf('  %-8s rule: max affine-C error = %.3e\n', nm, er);
end

%% PART A : affine regression vs baseline
fprintf('\n==== PART A: affine regression (Gauss / cmoment vs reduced-WQ baseline) ====\n');
geoms={'rect 1.0x0.5 [5 4]',1.0,0.5,[5 4];'square 1x1 [6 6]',1.0,1.0,[6 6]};
for ig=1:size(geoms,1)
 for model={'plate','T1'}
  pd=cantilever_beam(geoms{ig,2},geoms{ig,3});
  md=struct('degree',[3 3],'regularity',[2 2],'nsub',geoms{ig,4},'nquad',[4 4]);
  [geometry,msh,space]=buildSpaces(pd,md); te=ones(md.nsub);
  He=harmonic_kd_setup(msh,space,geometry,E,NU,RHO,'cpu','double',model{1});
  Vv=harmonic_kd_form(He,te,1,0,0,'values');
  Kb_aff=sparse(He.rows,He.cols,Vv.Kb,He.N,He.N); Kb_aff=0.5*(Kb_aff+Kb_aff.');
  Kb_g=harmonic_kd_gauss_curved(msh,space,geometry,E,NU,model{1});
  Hc=harmonic_kd_setup_cmoment(msh,space,geometry,E,NU,RHO,'cpu','double',model{1});
  Kb_c=harmonic_kd_form_curved(Hc,te);
  fprintf('  %-18s %-6s : Gauss %.2e | cmoment %.2e\n', geoms{ig,1}, model{1}, ...
     norm(Kb_g-Kb_aff,'fro')/norm(Kb_aff,'fro'), norm(Kb_c-Kb_aff,'fro')/norm(Kb_aff,'fro'));
 end
end

%% PART B : curved chart convergence
fprintf('\n==== PART B: curved NURBS chart, convergence (analytic ref, ng=20) ====\n');
n0=nrb4surf([0 0],[1 0],[0 1],[1 1]); n0=nrbdegelev(n0,[2 2]);
n0.coefs(1,2,2)=n0.coefs(1,2,2)+0.15; n0.coefs(2,2,2)=n0.coefs(2,2,2)-0.15;
n0.coefs(1,3,3)=n0.coefs(1,3,3)-0.12; n0.coefs(2,3,3)=n0.coefs(2,3,3)+0.12;
n0.coefs(1,2,3)=n0.coefs(1,2,3)+0.06; n0.coefs(2,3,2)=n0.coefs(2,3,2)-0.06;
nsubs=[2 2;4 4;8 8;16 16];
fprintf('  %-8s %6s %13s %13s %13s %9s %9s %7s\n','nsub','N','E_WQ','E_Gauss','E_cmom','e_WQ','e_Gauss','o_Gauss');
fprintf('  %s\n',repmat('-',1,90));
pG=NaN; pC=NaN; ph=NaN;
for is=1:size(nsubs,1)
  pd=struct('geo_name',n0); md=struct('degree',[3 3],'regularity',[2 2],'nsub',nsubs(is,:),'nquad',[5 5]);
  [geometry,msh,space]=buildSpaces(pd,md); te=ones(nsubs(is,:));
  Hw=harmonic_kd_setup_curved(msh,space,geometry,E,NU,RHO,'cpu','double','plate');
  Kb_w=harmonic_kd_form_curved(Hw,te);
  Kb_g=harmonic_kd_gauss_curved(msh,space,geometry,E,NU,'plate',deg_nq(space));
  Hc=harmonic_kd_setup_cmoment(msh,space,geometry,E,NU,RHO,'cpu','double','plate');
  Kb_c=harmonic_kd_form_curved(Hc,te);
  sp1=space.scalar_spaces{1}.sp_univ(1); sp2=space.scalar_spaces{1}.sp_univ(2);
  B1g=scol(sp1.ndof,sp1.degree,msh.qn{1}(:).',sp1.knots);
  B2g=scol(sp2.ndof,sp2.degree,msh.qn{2}(:).',sp2.knots);
  X1=repmat(msh.qn{1}(:),1,numel(msh.qn{2})); X2=repmat(msh.qn{2}(:).',numel(msh.qn{1}),1);
  W=wfun(X1,X2); c=(B1g\W/B2g.'); c=c(:);
  Ew=c.'*Kb_w*c; Eg=c.'*Kb_g*c; Ec=c.'*Kb_c*c;
  Er=dense_ref(msh,geometry,20,D,NU,'plate','analytic',0);
  ew=abs(Ew-Er)/abs(Er); eg=abs(Eg-Er)/abs(Er); ec=abs(Ec-Er)/abs(Er);
  h=1/nsubs(is,1);
  if isnan(pG), og=NaN; oc=NaN; else, og=log(pG/eg)/log(ph/h); oc=log(pC/ec)/log(ph/h); end
  fprintf('  %-8s %6d %13.5e %13.5e %13.5e %9.2e %9.2e %7.2f (cmom o=%.2f)\n', ...
     mat2str(nsubs(is,:)),numel(c),Ew,Eg,Ec,ew,eg,og,oc);
  pG=eg; pC=ec; ph=h;
end

%% PART C : cost + accuracy
fprintf('\n==== PART C: cost (setup+form wall time) vs accuracy ====\n');
fprintf('  %-8s %-10s %9s %12s %10s\n','nsub','method','time_s','rel.err','ops/elem');
for is=2:4
  pd=struct('geo_name',n0); md=struct('degree',[3 3],'regularity',[2 2],'nsub',nsubs(is,:),'nquad',[5 5]);
  [geometry,msh,space]=buildSpaces(pd,md); te=ones(nsubs(is,:));
  sp1=space.scalar_spaces{1}.sp_univ(1); sp2=space.scalar_spaces{1}.sp_univ(2);
  B1g=scol(sp1.ndof,sp1.degree,msh.qn{1}(:).',sp1.knots); B2g=scol(sp2.ndof,sp2.degree,msh.qn{2}(:).',sp2.knots);
  X1=repmat(msh.qn{1}(:),1,numel(msh.qn{2})); X2=repmat(msh.qn{2}(:).',numel(msh.qn{1}),1);
  W=wfun(X1,X2); c=(B1g\W/B2g.'); c=c(:);
  Er=dense_ref(msh,geometry,20,D,NU,'plate','analytic',0);
  t=tic; Hw=harmonic_kd_setup_curved(msh,space,geometry,E,NU,RHO,'cpu','double','plate');
  Kb=harmonic_kd_form_curved(Hw,te); t1=toc(t); e1=abs(c.'*Kb*c-Er)/abs(Er);
  fprintf('  %-8s %-10s %9.3f %12.3e %10d\n', mat2str(nsubs(is,:)),'reducedWQ',t1,e1, 25*numel(Hw.rows));
  nq=deg_nq(space); nloc=(degree_of(space)+1)^2;
  t=tic; Kb=harmonic_kd_gauss_curved(msh,space,geometry,E,NU,'plate',nq); t2=toc(t);
  e2=abs(c.'*Kb*c-Er)/abs(Er);
  fprintf('  %-8s %-10s %9.3f %12.3e %10d\n', mat2str(nsubs(is,:)),'Gauss',t2,e2, nq^2*nloc^2*25);
  t=tic; Hc=harmonic_kd_setup_cmoment(msh,space,geometry,E,NU,RHO,'cpu','double','plate');
  Kb=harmonic_kd_form_curved(Hc,te); t3=toc(t); e3=abs(c.'*Kb*c-Er)/abs(Er);
  fprintf('  %-8s %-10s %9.3f %12.3e %10d\n', mat2str(nsubs(is,:)),'C-weighted',t3,e3, 25*numel(Hc.rows));
end
disp('REMEDY_TEST_DONE');

% ---------------- local functions ----------------
function nq=deg_nq(space), nq=degree_of(space)+3; end
function d=degree_of(space), d=space.scalar_spaces{1}.sp_univ(1).degree; end
function v=Qcol(Q,k2,ip,b)
  switch k2
    case 0, M=Q.B_global; case 1, M=Q.B1_global; otherwise, M=Q.B2_global;
  end
  v=M(ip,b);
end
function B=scol(ndof,p,pts,kn)
  iv=findspan(ndof-1,p,pts,kn); b=basisfun(iv,pts,p,kn); nb=numbasisfun(iv,pts,p,kn)+1;
  r=repmat((1:numel(pts)).',1,p+1); B=sparse(r(:),nb(:),b(:),numel(pts),ndof);
end
function [B1,B2]=dcol(ndof,p,pts,kn,D1,D2)
  kd1=kn(2:end-1); p1=p-1; nd1=ndof-1;
  iv=findspan(nd1-1,p1,pts,kd1); b=basisfun(iv,pts,p1,kd1); nb=numbasisfun(iv,pts,p1,kd1)+1;
  r=repmat((1:numel(pts)).',1,p1+1); Bv1=sparse(r(:),nb(:),b(:),numel(pts),nd1); B1=Bv1*D1;
  kd2=kn(3:end-2); p2=p-2; nd2=ndof-2;
  iv=findspan(nd2-1,p2,pts,kd2); b=basisfun(iv,pts,p2,kd2); nb=numbasisfun(iv,pts,p2,kd2)+1;
  r=repmat((1:numel(pts)).',1,p2+1); Bv2=sparse(r(:),nb(:),b(:),numel(pts),nd2); B2=Bv2*(D2*D1);
end
function w=wfun(xi,eta)
  w=0.5*xi.^3-xi.^2.*eta+0.7*xi.*eta.^2+0.3*eta.^3+xi.*eta+0.2*xi+0.1*eta;
end
function E=dense_ref(msh,geometry,ng,D,NU,model,mode,fdstep)
  rule=msh_gauss_nodes(ng);
  [qd1,qw1]=msh_set_quad_nodes(unique(msh.breaks{1}),rule);
  [qd2,qw2]=msh_set_quad_nodes(unique(msh.breaks{2}),rule);
  nel1=size(qd1,2); nel2=size(qd2,2); E=0;
  for el2=1:nel2, for el1=1:nel1, for g2=1:ng, for g1=1:ng
    xi=qd1(g1,el1); eta=qd2(g2,el2); wq=qw1(g1,el1)*qw2(g2,el2);
    if strcmpi(mode,'analytic')
      [~,J]=geometry.map_der([xi;eta]); [~,~,H]=geometry.map_der2([xi;eta]);
      J=squeeze(J); H=squeeze(H);
    else
      [J,H]=fd_map_derivatives(geometry,xi,eta,fdstep);
    end
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
function [J,H]=fd_map_derivatives(geometry,xi,eta,h0)
  h=min(h0,0.45*min([xi,1-xi,eta,1-eta])); h=max(h,5e-6);
  x0=geometry.map([xi;eta]);
  xp=geometry.map([xi+h;eta]); xm=geometry.map([xi-h;eta]);
  yp=geometry.map([xi;eta+h]); ym=geometry.map([xi;eta-h]);
  xpp=geometry.map([xi+h;eta+h]); xpm=geometry.map([xi+h;eta-h]);
  xmp=geometry.map([xi-h;eta+h]); xmm=geometry.map([xi-h;eta-h]);
  J=[(xp-xm)/(2*h),(yp-ym)/(2*h)]; H=zeros(2,2,2);
  H(:,1,1)=(xp-2*x0+xm)/h^2; H(:,2,2)=(yp-2*x0+ym)/h^2;
  H(:,1,2)=(xpp-xpm-xmp+xmm)/(4*h^2); H(:,2,1)=H(:,1,2);
end
function S=inv_map_der2(Ji,H)
  S=zeros(2,2,2);
  for alpha=1:2, for i=1:2, for j=1:2
    s=0; for m=1:2, for b=1:2, for c=1:2, s=s-Ji(alpha,m)*H(m,b,c)*Ji(b,i)*Ji(c,j); end, end, end
    S(alpha,i,j)=s;
  end, end, end
end
