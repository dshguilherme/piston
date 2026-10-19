% TEST_SHELL_VERIFICATION -- verification of the rdim=3 thin-shell kinematics
% and the metric-carrying weighted-quadrature (WQ) bending operator.
%
%   A. REGRESSION GATE: the shell operator on a flat surface embedded in R^3
%      reproduces the existing rdim=2 curved-chart operator
%      (harmonic_kd_setup_curved + harmonic_kd_form_curved) to machine precision.
%   B. POINTWISE KINEMATICS on a cylinder (R=1.5) and a spherical cap:
%      surface metric a_ab, area element sqrt(det a), unit normal,
%      second fundamental form b_ab (vs the Weingarten route -n_,a . x_,b),
%      Christoffel symbols Gamma^g_ab (vs a finite-difference of the metric),
%      the parameterization-invariant curvatures K=det(b)/det(a) and
%      Hm=1/2 a^ab b_ab against the analytic values of the surface, and the
%      covariant bending strain kappa_ab=w_;ab against an independent route
%      (finite differences of w and of the metric, no analytic map_der2).
%   C. CANONICAL OBSTACLE GEOMETRIES: the confirmed blocker (harmonic_kd_setup_curved
%      aborts on a 3x2 Jacobian) is reproduced, and shell_kd_setup RUNS on all
%      three; the pinched-cylinder surface quantities are checked against the
%      analytic cylinder (R=300, K=0, Hm=1/(2R)).
%   D. ENERGY: the WQ shell bending energy vs a dense Gauss reference on curved
%      surfaces (cylinder, spherical cap), reporting the convergence (or NOT)
%      -- this is the reduced-rule characterisation.
%
% Writes results/shell_verification.txt and results/shell_verification.mat and
% results/shell_energy_convergence.png.
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
resdir = 'C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(resdir,'dir'), mkdir(resdir); end
fid = fopen(fullfile(resdir,'shell_verification.txt'),'w');
lg = @(varargin) both(fid, varargin{:});
R_out = struct();

E = 1.0; NU = 0.3; RHO = 1.0;

%% ================= A. REGRESSION GATE (flat embed -> rdim=2) ==============
lg('==== A. REGRESSION GATE: shell operator on a flat surface vs rdim=2 curved operator ====\n');
n0 = nrb4surf([0 0],[1 0],[0 1],[1 1]); n0 = nrbdegelev(n0,[2 2]);
n0.coefs(1,2,2)=n0.coefs(1,2,2)+0.15; n0.coefs(2,2,2)=n0.coefs(2,2,2)-0.15;
n0.coefs(1,3,3)=n0.coefs(1,3,3)-0.12; n0.coefs(2,3,3)=n0.coefs(2,3,3)+0.12;
Aok = true; reg_err = [];
for model = {'plate','T1'}
  for nsubv = [4 8 16]
    [geo2,msh2,sp2] = build_shell_msh_space(n0,[3 3],[2 2],[nsubv nsubv],[4 4],2);
    He = harmonic_kd_setup_curved(msh2,sp2,geo2,E,NU,RHO,'cpu','double',model{1});
    Kb2 = harmonic_kd_form_curved(He,ones(msh2.nel_dir));
    geo3.map=@(P) pad0(geo2.map(P)); geo3.map_der=@(P) pad0(geo2.map_der(P)); geo3.map_der2=@(P) pad0(geo2.map_der2(P));
    for ib=1:numel(geo2.boundary)
      bd=geo2.boundary(ib);
      geo3.boundary(ib).map=@(P) pad0(bd.map(P)); geo3.boundary(ib).map_der=@(P) pad0(bd.map_der(P));
    end
    [qn,qw]=msh_set_quad_nodes(geo2.nurbs.knots, msh_gauss_nodes([4 4]));
    msh3=msh_cartesian(geo2.nurbs.knots,qn,qw,geo3);
    sp3=sp_vector(repmat({sp_nurbs(geo2.nurbs,msh3)},1,3),msh3);
    Hs=shell_kd_setup(msh3,sp3,geo3,E,NU,RHO,'cpu','double',model{1});
    Kb3=harmonic_kd_form_curved(Hs,ones(msh3.nel_dir));
    er=norm(Kb3-Kb2,'fro')/norm(Kb2,'fro');
    reg_err(end+1)=er; %#ok<SAGROW>
    lg('  %-6s nsub=%2d rdim %d->%d  op rel err = %.3e\n',model{1},nsubv,msh2.rdim,msh3.rdim,er);
    Aok = Aok && (er < 1e-13);
  end
end
lg('  A. REGRESSION GATE: %s (max err %.1e)\n\n', ternary(Aok,'OK','FAIL'), max(reg_err));
R_out.reg_max_err = max(reg_err);

%% ================= B. POINTWISE KINEMATICS (cyl + sphere) =================
lg('==== B. POINTWISE surface kinematics: cylinder & spherical cap ====\n');
Rcyl = 1.5; Hh = 2.0;
geo_cyl = geo_load(nrbcylind(Hh, Rcyl, [0 0 0], 0, pi*0.9));
Rsp = 1.0; thcap = pi/4;
arc = nrbcirc(Rsp,[0 0 0],0,thcap); arc = nrbtform(arc, vecrotx(pi/2));
geo_sph = geo_load(nrbrevolve(arc,[0 0 0],[0 0 1],2*pi));
rng(7);
surfaces = {'cylinder', geo_cyl, 1/(2*Rcyl), 0.0; 'sphere', geo_sph, 1/Rsp, 1/Rsp^2};
Btab = zeros(2,6);
for is = 1:2
  nm=surfaces{is,1}; g=surfaces{is,2}; Hm_ref=surfaces{is,3}; K_ref=surfaces{is,4};
  eK=0; eH=0; eb=0; eG=0; ek=0;
  for k=1:60
    xi=0.1+0.8*rand; eta=0.1+0.8*rand; P=[xi;eta];
    Sk = shell_curvature(g,P,'analytic');
    a=Sk.a(:,:,1); ai=Sk.ainv(:,:,1); n=Sk.normal(:,1); b=Sk.b(:,:,1); Gm=Sk.Gamma(:,:,:,1);
    h=1e-4;
    J0 = Jfd(g,xi,eta,h); a0=J0.'*J0; ai0=inv(a0);
    pA=cell(1,2);
    pA{1}=(Afd(g,xi+h,eta,h)-Afd(g,xi-h,eta,h))/(2*h);
    pA{2}=(Afd(g,xi,eta+h,h)-Afd(g,xi,eta-h,h))/(2*h);
    G1=zeros(2,2,2);
    for al=1:2, for be=1:2, for d=1:2
      G1(al,be,d)=0.5*( pA{be}(d,al)+pA{al}(d,be)-pA{d}(al,be) );
    end, end, end
    Gfd=zeros(2,2,2);
    for al=1:2, for be=1:2, for gm=1:2
      s=0; for d=1:2, s=s+ai0(gm,d)*G1(al,be,d); end
      Gfd(al,be,gm)=s;
    end, end, end
    n1=(nfd(g,xi+h,eta,h)-nfd(g,xi-h,eta,h))/(2*h);
    n2=(nfd(g,xi,eta+h,h)-nfd(g,xi,eta-h,h))/(2*h);
    bfd=zeros(2,2); bfd(:,1)=-(n1.'*J0).'; bfd(:,2)=-(n2.'*J0).';
    eK=max(eK, abs(det(b)/det(a)-K_ref));
    eH=max(eH, abs(abs(0.5*trace(ai*b))-Hm_ref));   % |Hm|: normal sign is orientation-dependent
    eb=max(eb, max(abs(b(:)-bfd(:)))/max(abs(b(:))));
    eG=max(eG, max(abs(Gm(:)-Gfd(:)))/max(abs(Gm(:))));
    % covariant bending strain: analytic Gamma+b vs independent (metric-FD Gamma + parametric-FD w)
    w0=wfun(xi,eta);
    wxi=(wfun(xi+h,eta)-wfun(xi-h,eta))/(2*h);
    weta=(wfun(xi,eta+h)-wfun(xi,eta-h))/(2*h);
    wxixi=(wfun(xi+h,eta)-2*w0+wfun(xi-h,eta))/h^2;
    wetaeta=(wfun(xi,eta+h)-2*w0+wfun(xi,eta-h))/h^2;
    wxieta=(wfun(xi+h,eta+h)-wfun(xi+h,eta-h)-wfun(xi-h,eta+h)+wfun(xi-h,eta-h))/(4*h^2);
    k11m=wxixi-Gm(1,1,1)*wxi-Gm(1,1,2)*weta;
    k12m=wxieta-Gm(1,2,1)*wxi-Gm(1,2,2)*weta;
    k22m=wetaeta-Gm(2,2,1)*wxi-Gm(2,2,2)*weta;
    k11f=wxixi-Gfd(1,1,1)*wxi-Gfd(1,1,2)*weta;
    k12f=wxieta-Gfd(1,2,1)*wxi-Gfd(1,2,2)*weta;
    k22f=wetaeta-Gfd(2,2,1)*wxi-Gfd(2,2,2)*weta;
    ek=max(ek, max(abs([k11m-k11f k12m-k12f k22m-k22f]))/max(abs([k11m k12m k22m])) );
  end
  lg('  %-9s  |K-K_ref|=%.2e  |Hm-Hm_ref|=%.2e  |b-weg|=%.2e  |G-Gfd|=%.2e  |k-k_fd|=%.2e\n', ...
     nm, eK, eH, eb, eG, ek);
  Btab(is,:)=[eK eH eb eG ek 0];
end
Bok = all(Btab(:,1) < 1e-8) && all(Btab(:,2) < 1e-6) && all(Btab(:,3) < 1e-5) && all(Btab(:,4) < 1e-5) && all(Btab(:,5) < 1e-3);
lg('  B. POINTWISE: %s\n\n', ternary(Bok,'OK','FAIL'));
R_out.pointwise = Btab;

%% ================= C. CANONICAL OBSTACLE GEOMETRIES =======================
lg('==== C. CANONICAL obstacle-course geometries (rdim=3, ndim=2) ====\n');
probs = obstacle_course_geom();
Ctab = cell(numel(probs),4);
for ip=1:numel(probs)
  pr=probs(ip);
  [geo,msh,space]=build_shell_msh_space(pr.nurbs,[3 3],[2 2],[8 8],[4 4],pr.rdim);
  % (i) confirmed blocker: the rdim=2 operator must ABORT on this 3x2 Jacobian
  aborted=false;
  try
    harmonic_kd_setup_curved(msh,space,geo,pr.E,pr.nu,1,'cpu','double','plate');
  catch
    aborted=true;
  end
  % (ii) the rdim=3 shell operator must RUN and return a finite operator
  ok=false; nrm=NaN;
  try
    Hs=shell_kd_setup(msh,space,geo,pr.E,pr.nu,1,'cpu','double','plate');
    Kb=harmonic_kd_form_curved(Hs,ones(msh.nel_dir));
    nrm=full(norm(Kb,'fro')); ok=isfinite(nrm) && nrm>0;
  catch ME2
    lg('    shell_kd_setup FAILED on %s: %s\n', pr.name, ME2.message);
  end
  lg('  %-24s rdim=%d  rdim2-op:%s  shell-op:%s  ||Kb||_F=%.3e\n', ...
     pr.name, msh.rdim, ternary(aborted,'ABORTS','runs'), ternary(ok,'runs','FAILS'), nrm);
  Ctab{ip,1}=pr.name; Ctab{ip,2}=aborted; Ctab{ip,3}=ok; Ctab{ip,4}=nrm;
end
Cok = all(cellfun(@(x) x, Ctab(:,3)));
lg('  C. CANONICAL: shell operator runs on all three: %s\n\n', ternary(Cok,'OK','FAIL'));
R_out.canonical = Ctab;

%% ================= D. ENERGY: WQ vs dense Gauss (characterisation) ========
lg('==== D. WQ shell bending energy vs dense Gauss reference (curved) ====\n');
lg('   (reduced WQ rules are exact only for a CONSTANT coefficient; on a curved\n');
lg('    surface the shell coefficient varies -> this measures the quadrature error)\n');
D = E/(12*(1-NU^2));
curvs = {'cylinder (R=1.5,H=2, 0..0.9pi)', nrbcylind(Hh,Rcyl,[0 0 0],0,pi*0.9); ...
         'spherical cap (R=1, 45deg)',         sphere_nurbs(Rsp,thcap)};
energy_tab = cell(2,1);
for ic = 1:2
  nm=curvs{ic,1}; nr=curvs{ic,2};
  lg('  -- %s --\n', nm);
  lg('  %-10s %6s %16s %16s %10s %8s\n','nsub','N','E_wq','E_ref','rel.err','order');
  prev=NaN; prevh=NaN; rows=[];
  for nsub=[3 6 12 24]
    [geo,msh,space]=build_shell_msh_space(nr,[3 3],[2 2],[nsub nsub],[5 5],3);
    Hs=shell_kd_setup(msh,space,geo,E,NU,RHO,'cpu','double','plate');
    Kb=harmonic_kd_form_curved(Hs,ones(msh.nel_dir));
    sp1=space.scalar_spaces{1}.sp_univ(1); sp2=space.scalar_spaces{1}.sp_univ(2);
    B1=sp_colloc(sp1); B2=sp_colloc(sp2);
    X1=repmat(msh.qn{1}(:),1,numel(msh.qn{2})); X2=repmat(msh.qn{2}(:).',numel(msh.qn{1}),1);
    W=wfun(X1,X2); C=B1\W/B2.'; c=C(:);
    Ewq=c.'*Kb*c;
    Eref=dense_ref_shell(msh,geo,8,D,NU,'plate');
    rel=abs(Ewq-Eref)/abs(Eref); h=1/nsub;
    if isnan(prev), ord=NaN; else, ord=log(prev/rel)/log(prevh/h); end
    lg('  %-10d %6d %16.8e %16.8e %10.3e %8.2f\n',nsub,Hs.N,Ewq,Eref,rel,ord);
    rows(end+1,:)=[nsub Hs.N Ewq Eref rel ord]; %#ok<SAGROW>
    prev=rel; prevh=h;
  end
  energy_tab{ic}=rows;
end
R_out.energy = energy_tab;
lg('\n  D. NOTE: the reduced WQ rule does NOT converge for the varying curved\n');
lg('     coefficient (order negative at fine meshes). The energy rows above are\n');
lg('     the raw evidence; no convergence claim is made.\n\n');

%% plot
figure('Visible','off','Position',[100 100 700 500]);
hold on; mk=0;
for ic=1:2
  rows=energy_tab{ic};
  mk=mk+1;
  loglog(1./rows(:,1), rows(:,5), '-o', 'DisplayName', curvs{ic,1});
end
xlabel('element size h'); ylabel('relative energy error |E_{WQ}-E_{ref}|/|E_{ref}|');
title('WQ shell bending energy error vs dense-Gauss reference'); legend('Location','northwest');
grid on; set(gca,'XScale','log','YScale','log');
saveas(gcf, fullfile(resdir,'shell_energy_convergence.png'));
close(gcf);

save(fullfile(resdir,'shell_verification.mat'),'R_out','reg_err');
lg('Saved results/shell_verification.txt, shell_verification.mat, shell_energy_convergence.png\n');
fclose(fid);
fprintf('SHELL_VERIFICATION_DONE\n');
disp('SHELL_VERIFICATION_DONE');

% ============================ local functions =============================
function [geo,msh,space]=build_shell_msh_space(nr,degree,regularity,nsub,nquad,rdim)
  geometry=geo_load(nr);
  degelev=max(degree-(geometry.nurbs.order-1),0);
  nurbs=nrbdegelev(geometry.nurbs,degelev);
  [rk,zeta,nk]=kntrefine(nurbs.knots,nsub-1,nurbs.order-1,regularity);
  nurbs=nrbkntins(nurbs,nk); geo=geo_load(nurbs);
  rule=msh_gauss_nodes(nquad); [qn,qw]=msh_set_quad_nodes(geo.nurbs.knots,rule);
  msh=msh_cartesian(geo.nurbs.knots,qn,qw,geo);
  sp_scalar=sp_nurbs(geo.nurbs,msh); space=sp_vector(repmat({sp_scalar},1,rdim),msh);
end
function Y=pad0(X); nsz=[1,size(X,2:ndims(X))]; Y=cat(1,X,zeros(nsz)); end
function nr=sphere_nurbs(R,th)
  arc=nrbcirc(R,[0 0 0],0,th); arc=nrbtform(arc,vecrotx(pi/2)); nr=nrbrevolve(arc,[0 0 0],[0 0 1],2*pi);
end
function J=Jfd(g,xi,eta,h)
  J=[(g.map([xi+h;eta])-g.map([xi-h;eta]))/(2*h), (g.map([xi;eta+h])-g.map([xi;eta-h]))/(2*h)];
end
function A=Afd(g,xi,eta,h); J=Jfd(g,xi,eta,h); A=J.'*J; end
function n=nfd(g,xi,eta,h); J=Jfd(g,xi,eta,h); c=cross(J(:,1),J(:,2)); n=c/norm(c); end
function w=wfun(xi,eta)
  w=0.5*xi.^3 - xi.^2.*eta + 0.7*xi.*eta.^2 + 0.3*eta.^3 + xi.*eta + 0.2*xi + 0.1*eta;
end
function B1g=sp_colloc(sp1)
  nqn=size(sp1.shape_functions,1); nel=size(sp1.shape_functions,3);
  rows=repmat(reshape(1:nqn*nel,[nqn,1,nel]),[1,sp1.nsh_max,1]);
  cols=repmat(reshape(sp1.connectivity,[1,sp1.nsh_max,nel]),[nqn,1,1]);
  B1g=sparse(rows(:),cols(:),sp1.shape_functions(:),nqn*nel,sp1.ndof);
end
function E=dense_ref_shell(msh,geometry,ng,D,NU,model)
  rule=msh_gauss_nodes(ng);
  [qd1,qw1]=msh_set_quad_nodes(unique(msh.breaks{1}),rule);
  [qd2,qw2]=msh_set_quad_nodes(unique(msh.breaks{2}),rule);
  E=0;
  for el2=1:size(qd2,2)
    for el1=1:size(qd1,2)
      for g2=1:ng
        for g1=1:ng
          xi=qd1(g1,el1); eta=qd2(g2,el2); wq=qw1(g1,el1)*qw2(g2,el2);
          Sk=shell_curvature(geometry,[xi;eta],'analytic');
          ai=Sk.ainv(:,:,1); Gm=Sk.Gamma(:,:,:,1); sq=Sk.sqrta(1);
          wxixi=3*xi-2*eta; wxieta=-2*xi+1.4*eta+1; wetaeta=1.4*xi+1.8*eta;
          wxi=1.5*xi^2-2*xi*eta+0.7*eta^2+eta+0.2; weta=-xi^2+1.4*xi*eta+0.9*eta^2+xi+0.1;
          k=[wxixi-Gm(1,1,1)*wxi-Gm(1,1,2)*weta, wxieta-Gm(1,2,1)*wxi-Gm(1,2,2)*weta;
             wxieta-Gm(1,2,1)*wxi-Gm(1,2,2)*weta, wetaeta-Gm(2,2,1)*wxi-Gm(2,2,2)*weta];
          I1=ai(1,1)*k(1,1)+2*ai(1,2)*k(1,2)+ai(2,2)*k(2,2);
          I2=0;
          for al=1:2, for be=1:2, for g=1:2, for de=1:2
            I2=I2+ai(al,g)*ai(be,de)*k(al,be)*k(g,de);
          end, end, end, end
          if strcmpi(model,'plate'), e=NU*I1^2+(1-NU)*I2; else, e=I2; end
          E=E+wq*D*sq*e;
        end
      end
    end
  end
end
function both(fid, varargin)
  fprintf(fid, varargin{:}); fprintf(1, varargin{:});
end
function s=ternary(c,a,b); if c, s=a; else, s=b; end; end
