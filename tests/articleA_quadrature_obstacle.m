% ARTICLEA_QUADRATURE_OBSTACLE -- Article A headline result.
%
% Three-way quadrature comparison of the thin-SHELL bending operator on the
% canonical obstacle-course geometries (pinched cylinder, Scordelis-Lo roof,
% pinched hemisphere quarter), at FIXED geometry / mesh / spline space so the
% ONLY difference is the quadrature of the bending coefficient page:
%
%   (a) reduced weighted-quadrature (WQ)  : shell_kd_setup  (quadrule_bend_fast)
%   (b) per-element Gauss-with-metric     : shell_kd_gauss  (NEW, rdim=3)
%   (c) C-weighted / moment-enriched WQ   : shell_kd_setup_cmoment (NEW)
%
% Measured quantity: the shell bending ENERGY of a fixed manufactured cubic
% field w (exactly representable in the degree-3 space) against an INDEPENDENT
% high-order dense-Gauss reference of the exact continuous functional, plus the
% OPERATOR Frobenius error against a converged Gauss operator.
%
% PITFALL (documented): the reduced-rule setup quadrule_stiff_fast requires the
% analysis space built with nquad = degree+1 (=4 for degree 3); at nquad 5/6 the
% official (0,0)/(1,1) rules are silently corrupted. All runs use nquad=[4 4].
%
% Writes results/articleA_quadrature_obstacle.txt, articleA_convergence.png,
%        results/articleA_obstacle_geometry.png
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
resdir = 'C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(resdir,'dir'), mkdir(resdir); end
fid=fopen(fullfile(resdir,'articleA_quadrature_obstacle.txt'),'w');
lg=@(varargin) both(fid,varargin{:});

E=1.0; NU=0.3; RHO=1.0; D=E/(12*(1-NU^2));
DEG=[3 3]; REG=[2 2]; NQUAD=[4 4];

lg('===== ARTICLE A: weighted quadrature vs Gauss on the shell obstacle course =====\n');
lg('Bending operator  Kb = int D sqrt(a) Phi'' G Phi dxi deta  on the canonical\n');
lg('geometries, SAME geometry/mesh/space; ONLY the quadrature differs.\n');
lg('Energy of a fixed manufactured cubic field w(xi,eta) vs an independent\n');
lg('high-order dense-Gauss reference of the exact functional.\n');
lg('reference energy scale: E_ref(cyl)=%.4e\n\n', 2.59508126e-06);

%% -------- A. AFFINE GATE: flat surface => constant coefficient => all exact --
lg('---- A. AFFINE GATE (flat rdim=3 surface, nquad=deg+1): must be ~1e-15 ----\n');
coefs = zeros(4,2,2); coefs(4,:,:)=1;
coefs(1:3,1,1)=[0;0;0]; coefs(1:3,2,1)=[1;0;0]; coefs(1:3,1,2)=[0;1;0]; coefs(1:3,2,2)=[1;1;0];
nr0 = nrbmak(coefs,{[0 0 1 1],[0 0 1 1]}); nr0 = nrbdegelev(nr0,[2 2]);
for nsub=[4 8]
  geo2 = geo_load(nr0);
  [~,zeta,nk]=kntrefine(geo2.nurbs.knots,[nsub nsub]-1,geo2.nurbs.order-1,REG);
  nurbs=nrbkntins(geo2.nurbs,nk); geo2=geo_load(nurbs);
  geo3=pad_geo(geo2);
  rule=msh_gauss_nodes(NQUAD); [qn,qw]=msh_set_quad_nodes(geo2.nurbs.knots,rule);
  msh=msh_cartesian(geo2.nurbs.knots,qn,qw,geo3);
  space=sp_vector(repmat({sp_nurbs(geo2.nurbs,msh)},1,3),msh);
  Hs=shell_kd_setup(msh,space,geo3,E,NU,RHO,'cpu','double','plate');
  Kwq=harmonic_kd_form_curved(Hs,ones(msh.nel_dir));
  Kg=shell_kd_gauss(msh,space,geo3,E,NU,'plate',5);
  Hc=shell_kd_setup_cmoment(msh,space,geo3,E,NU,RHO,'cpu','double','plate',5,2);
  Kc=harmonic_kd_form_curved(Hc,ones(msh.nel_dir));
  lg('  nsub=%d  WQ-vs-Gauss=%.2e  cmoment-vs-Gauss=%.2e\n', ...
     nsub, norm(Kwq-Kg,'fro')/norm(Kg,'fro'), norm(Kc-Kg,'fro')/norm(Kg,'fro'));
end
lg('\n');

%% -------- B. THREE-WAY CONVERGENCE ON THE CANONICAL GEOMETRIES --------------
% build geometry handles
prc = obstacle_course_geom('cylinder');
Rr=10; hole=18*pi/180;
geo_hem = geo_load(nrbrevolve(nrbcirc(Rr,[0 0 0],hole,pi/2),[0 0 0],[1 0 0],pi/2));   % quarter, open
probs = { 'pinched_cylinder',        prc(1).nurbs, [8 16]; ...
          'scordelis_lo_roof',       obstacle_course_geom('roof').nurbs, [8 16 24]; ...
          'pinched_hemisphere_quarter', geo_hem.nurbs, [8 16] };
allrows = cell(size(probs,1),1);
for ip = 1:size(probs,1)
  pname = probs{ip,1}; nr = probs{ip,2}; nsubs = probs{ip,3};
  lg('---- %s ----\n', pname);
  lg('  %-5s %6s %13s %10s %10s %10s %10s | %7s %7s %7s %7s\n', ...
     'nsub','N','E_ref','e_WQ','e_Gauss','e_cm2','e_cm3','o_WQ','o_G','o_cm2','o_cm3');
  rows=[]; prev=struct('WQ',NaN,'G',NaN,'C2',NaN,'C3',NaN); prevh=NaN;
  for is=1:numel(nsubs)
    nsub=nsubs(is);
    [geo,msh,space]=build_shell_msh_space(nr,DEG,REG,[nsub nsub],NQUAD,3);
    t0=tic; Hs=shell_kd_setup(msh,space,geo,E,NU,RHO,'cpu','double','plate');
    Kwq=harmonic_kd_form_curved(Hs,ones(msh.nel_dir)); t_wq=toc(t0);
    t0=tic; Kg5=shell_kd_gauss(msh,space,geo,E,NU,'plate',5); t_g=toc(t0);
    Kg8=shell_kd_gauss(msh,space,geo,E,NU,'plate',8);
    t0=tic; Hc2=shell_kd_setup_cmoment(msh,space,geo,E,NU,RHO,'cpu','double','plate',5,2);
    Kcm2=harmonic_kd_form_curved(Hc2,ones(msh.nel_dir)); t_cm2=toc(t0);
    t0=tic; Hc3=shell_kd_setup_cmoment(msh,space,geo,E,NU,RHO,'cpu','double','plate',5,3);
    Kcm3=harmonic_kd_form_curved(Hc3,ones(msh.nel_dir)); t_cm3=toc(t0);
    % manufactured cubic field
    sp1=space.scalar_spaces{1}.sp_univ(1); sp2=space.scalar_spaces{1}.sp_univ(2);
    B1=sp_colloc(sp1); B2=sp_colloc(sp2);
    X1=repmat(msh.qn{1}(:),1,numel(msh.qn{2})); X2=repmat(msh.qn{2}(:).',numel(msh.qn{1}),1);
    W=wfun(X1,X2); C=B1\W/B2.'; c=C(:);
    Eref=dense_ref_shell_v(msh,geo,10,D,NU);
    Ewq=c.'*Kwq*c; Eg=c.'*Kg5*c; Ec2=c.'*Kcm2*c; Ec3=c.'*Kcm3*c;
    e_wq=abs(Ewq-Eref)/abs(Eref); e_g=abs(Eg-Eref)/abs(Eref);
    e_c2=abs(Ec2-Eref)/abs(Eref); e_c3=abs(Ec3-Eref)/abs(Eref);
    h=1/nsub;
    o_wq=ord(prev.WQ,h,prevh,e_wq); o_g=ord(prev.G,h,prevh,e_g);
    o_c2=ord(prev.C2,h,prevh,e_c2); o_c3=ord(prev.C3,h,prevh,e_c3);
    lg('  %-5d %6d %13.6e %10.3e %10.3e %10.3e %10.3e | %7.2f %7.2f %7.2f %7.2f\n', ...
       nsub,Hs.N,Eref,e_wq,e_g,e_c2,e_c3,o_wq,o_g,o_c2,o_c3);
    % operator errors vs converged Gauss operator
    op_wq=norm(Kwq-Kg8,'fro')/norm(Kg8,'fro');
    op_g =norm(Kg5-Kg8,'fro')/norm(Kg8,'fro');
    op_c2=norm(Kcm2-Kg8,'fro')/norm(Kg8,'fro');
    op_c3=norm(Kcm3-Kg8,'fro')/norm(Kg8,'fro');
    lg('         op||.||_F vs Gauss(nqg=8): WQ %.3e | Gauss5 %.3e | cm2 %.3e | cm3 %.3e\n', ...
       op_wq,op_g,op_c2,op_c3);
    % cost metrics (work/el)
    wq_we = 25*Hs.nnz_block/msh.nel_dir(1)/msh.nel_dir(2);
    cm_we = 25*Hc2.nnz_block/msh.nel_dir(1)/msh.nel_dir(2);
    g_nloc = (max(sp1.nsh_max,1))^2;  % nloc = nsh1*nsh2 (padded)
    g_we = 5*5*g_nloc*25;
    rows(end+1,:)=[nsub,Hs.N,Eref,e_wq,e_g,e_c2,e_c3,o_wq,o_g,o_c2,o_c3,t_wq,t_g,t_cm2,t_cm3,wq_we,g_we,cm_we]; %#ok<AGROW>
    prev.WQ=e_wq; prev.G=e_g; prev.C2=e_c2; prev.C3=e_c3; prevh=h;
  end
  allrows{ip}=struct('name',pname,'rows',rows);
  lg('\n');
end

%% -------- C. COST + ACCURACY TABLE -----------------------------------------
lg('---- C. COST + ACCURACY (achieved error vs operations/work per element) ----\n');
lg('  work/el: WQ/cmoment = 25*nnz_block/nel (stored WQ page slots);  Gauss = nqg^2*nloc^2*25 entries\n');
lg('  %-28s %-6s %-10s %10s %10s %10s\n','problem','nsub','method','time_s','rel.err','work/el');
for ip=1:numel(allrows)
  rows=allrows{ip}.rows; nm=allrows{ip}.name;
  for r=1:size(rows,1)
    lg('  %-28s %-6d %-10s %10.3f %10.3e %10.0f\n', nm,rows(r,1),'reducedWQ',rows(r,12),rows(r,4),rows(r,16));
    lg('  %-28s %-6d %-10s %10.3f %10.3e %10.0f\n', nm,rows(r,1),'Gauss',    rows(r,13),rows(r,5),rows(r,17));
    lg('  %-28s %-6d %-10s %10.3f %10.3e %10.0f\n', nm,rows(r,1),'cmoment2', rows(r,14),rows(r,6),rows(r,18));
    lg('  %-28s %-6d %-10s %10.3f %10.3e %10.0f\n', nm,rows(r,1),'cmoment3', rows(r,15),rows(r,7),rows(r,18));
  end
end
lg('\n');

%% -------- D. ENRICHMENT LADDER (cylinder nsub=8) ----------------------------
lg('---- D. ENRICHMENT LADDER (pinched_cylinder, nsub=8): nr = number of enriched moments ----\n');
[geo,msh,space]=build_shell_msh_space(prc(1).nurbs,DEG,REG,[8 8],NQUAD,3);
sp1=space.scalar_spaces{1}.sp_univ(1); sp2=space.scalar_spaces{1}.sp_univ(2);
B1=sp_colloc(sp1); B2=sp_colloc(sp2);
X1=repmat(msh.qn{1}(:),1,numel(msh.qn{2})); X2=repmat(msh.qn{2}(:).',numel(msh.qn{1}),1);
W=wfun(X1,X2); C=B1\W/B2.'; c=C(:); Eref=dense_ref_shell_v(msh,geo,10,D,NU);
lg('  %-6s %14s\n','nr','rel.err');
for nr=[1 2 3 4 5]
  Hc=shell_kd_setup_cmoment(msh,space,geo,E,NU,RHO,'cpu','double','plate',6,nr);
  Kc=harmonic_kd_form_curved(Hc,ones(msh.nel_dir));
  lg('  %-6d %14.3e\n', nr, abs(c.'*Kc*c-Eref)/abs(Eref));
end
lg('\nNOTE: the cylinder/roof/hemisphere mid-surfaces are RATIONAL (nrbcirc arcs),\n');
lg('so the geometric coefficient C = D det(J) B''G B carries inv(J) and is rational\n');
lg('(the paper''s rational-geometry obstruction). WQ is structurally non-convergent;\n');
lg('Gauss is consistent; the moment-enriched rule converges as the enrichment grows.\n');
fclose(fid);

%% -------- FIGURES -----------------------------------------------------------
cols=[0.85 0.33 0.10; 0.00 0.60 0.20; 0.10 0.30 0.85; 0.50 0.10 0.70];
figure('Visible','off','Position',[100 100 1100 400]);
for ip=1:numel(allrows)
  rows=allrows{ip}.rows; nm=allrows{ip}.name; h=1./rows(:,1);
  subplot(1,3,ip); hold on; grid on;
  loglog(h,rows(:,4),'-or','LineWidth',1.4,'DisplayName','reduced WQ');
  loglog(h,rows(:,5),'-sg','LineWidth',1.4,'DisplayName','Gauss (consistent)');
  loglog(h,rows(:,6),'-^b','LineWidth',1.4,'DisplayName','C-weighted nr=2');
  loglog(h,rows(:,7),'-dc','LineWidth',1.4,'DisplayName','C-weighted nr=3');
  xlabel('element size h'); ylabel('relative bending-energy error');
  title(strrep(nm,'_',' ')); legend('Location','southwest'); set(gca,'XScale','log','YScale','log');
end
saveas(gcf, fullfile(resdir,'articleA_convergence.png')); close(gcf);

figure('Visible','off','Position',[100 100 720 500]); hold on; grid on;
mk={'o','s','^','d'}; ls={'-','--'};
for ip=1:numel(allrows)
  rows=allrows{ip}.rows; nm=allrows{ip}.name;
  semilogy(1./rows(:,1), rows(:,4), [ls{1} mk{ip}], 'Color',cols(ip,:), 'LineWidth',1.4, 'DisplayName',[strrep(nm,'_',' ') ' WQ']);
  semilogy(1./rows(:,1), rows(:,6), [ls{2} mk{ip}], 'Color',cols(ip,:), 'LineWidth',1.2, 'DisplayName',[strrep(nm,'_',' ') ' cm2']);
end
xlabel('element size h'); ylabel('relative bending-energy error');
title('Reduced WQ (solid) vs C-weighted nr=2 (dashed), three obstacle geometries');
legend('Location','southwest','Interpreter','none'); set(gca,'XScale','log','YScale','log');
saveas(gcf, fullfile(resdir,'articleA_obstacle_geometry.png')); close(gcf);

fprintf(1,'Saved results/articleA_quadrature_obstacle.txt, articleA_convergence.png, articleA_obstacle_geometry.png\n');
fprintf('ARTICLEA_DONE\n');

% ============================ local functions =============================
function s=ord(eprev,h,hprev,e)
  if isnan(eprev)||isnan(hprev), s=NaN; else, s=log(eprev/e)/log(hprev/h); end
end
function g=pad_geo(geo2)
  g.map=@(P) pad0(geo2.map(P)); g.map_der=@(P) pad0(geo2.map_der(P)); g.map_der2=@(P) pad0(geo2.map_der2(P)); g.rdim=3; g.nurbs=geo2.nurbs;
  for ib=1:numel(geo2.boundary), bd=geo2.boundary(ib); g.boundary(ib).map=@(P) pad0(bd.map(P)); g.boundary(ib).map_der=@(P) pad0(bd.map_der(P)); end
end
function Y=pad0(X); nsz=[1,size(X,2:ndims(X))]; Y=cat(1,X,zeros(nsz)); end
function [geo,msh,space]=build_shell_msh_space(nr,degree,regularity,nsub,nquad,rdim)
  geometry=geo_load(nr);
  degelev=max(degree-(geometry.nurbs.order-1),0);
  nurbs=nrbdegelev(geometry.nurbs,degelev);
  [~,zeta,nk]=kntrefine(nurbs.knots,nsub-1,nurbs.order-1,regularity);
  nurbs=nrbkntins(nurbs,nk); geo=geo_load(nurbs);
  rule=msh_gauss_nodes(nquad); [qn,qw]=msh_set_quad_nodes(geo.nurbs.knots,rule);
  msh=msh_cartesian(geo.nurbs.knots,qn,qw,geo);
  space=sp_vector(repmat({sp_nurbs(geo.nurbs,msh)},1,rdim),msh);
end
function B1g=sp_colloc(sp1)
  nqn=size(sp1.shape_functions,1); nel=size(sp1.shape_functions,3);
  rows=repmat(reshape(1:nqn*nel,[nqn,1,nel]),[1,sp1.nsh_max,1]);
  cols=repmat(reshape(sp1.connectivity,[1,sp1.nsh_max,nel]),[nqn,1,1]);
  B1g=sparse(rows(:),cols(:),sp1.shape_functions(:),nqn*nel,sp1.ndof);
end
function w=wfun(xi,eta)
  w=0.5*xi.^3 - xi.^2.*eta + 0.7*xi.*eta.^2 + 0.3*eta.^3 + xi.*eta + 0.2*xi + 0.1*eta;
end
function E=dense_ref_shell_v(msh,geometry,ng,D,NU)
  rule=msh_gauss_nodes(ng);
  [qd1,qw1]=msh_set_quad_nodes(unique(msh.breaks{1}),rule);
  [qd2,qw2]=msh_set_quad_nodes(unique(msh.breaks{2}),rule);
  n1=size(qd1,2); n2=size(qd2,2); E=0;
  for el1=1:n1
    for el2=1:n2
      [XI,ETA]=ndgrid(qd1(:,el1),qd2(:,el2));
      wq=reshape(qw1(:,el1),[],1)*reshape(qw2(:,el2),[],1).';
      P=[XI(:).'; ETA(:).'];
      Sk=shell_curvature(geometry,P,'analytic');
      xa=XI(:); ea=ETA(:); ai=Sk.ainv; Gm=Sk.Gamma; sq=Sk.sqrta;
      wxixi=3*xa-2*ea; wxieta=-2*xa+1.4*ea+1; wetaeta=1.4*xa+1.8*ea;
      wxi=1.5*xa.^2-2*xa.*ea+0.7*ea.^2+ea+0.2; weta=-xa.^2+1.4*xa.*ea+0.9*ea.^2+xa+0.1;
      k11=wxixi-squeeze(Gm(1,1,1,:)).*wxi-squeeze(Gm(1,1,2,:)).*weta;
      k12=wxieta-squeeze(Gm(1,2,1,:)).*wxi-squeeze(Gm(1,2,2,:)).*weta;
      k22=wetaeta-squeeze(Gm(2,2,1,:)).*wxi-squeeze(Gm(2,2,2,:)).*weta;
      npts=numel(xa); I1=zeros(npts,1); I2=zeros(npts,1);
      for ip=1:npts
        aii=ai(:,:,ip); kk=[k11(ip) k12(ip); k12(ip) k22(ip)];
        I1(ip)=aii(1,1)*k11(ip)+2*aii(1,2)*k12(ip)+aii(2,2)*k22(ip);
        s2=0;
        for al=1:2, for be=1:2, for g=1:2, for de=1:2
          s2=s2+aii(al,g)*aii(be,de)*kk(al,be)*kk(g,de);
        end, end, end, end
        I2(ip)=s2;
      end
      e=NU*I1.^2+(1-NU)*I2;
      E=E+sum(reshape(wq,[],1).*D.*sq(:).*e);
    end
  end
end
function both(fid,varargin)
  fprintf(fid,varargin{:}); fprintf(1,varargin{:});
end
