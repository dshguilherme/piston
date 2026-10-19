% ARTICLEB_R0_WAVECAP_RUN  Phase-cancellation R0 ladder on the MODALLY DENSE
% clamped spherical cap, at an excitation frequency inside the dense band where
% BOTH conditions hold: (A) modal density -> incoherent surface velocity field,
% (B) k*Lmax >~ 1 -> structured sinc kernel.
%
% Reuses tests/obstacle_build_ctx.m ('spherecap': R=10 m, 30 deg half-angle,
% periodic in theta, clamped rim, crown point load) and the src R0 machinery
% (articleB_r0_build/eval).  Everything else (mesh, volume constraint, damping,
% GCMMA settings) matches the existing R0 runs.
%
% Runs the three-operator ladder to convergence and reports, per design:
%   W_rad (R0=3.66), W under R0=I, W under the point-source (ones) operator,
%   pointwise velocity norm, rho_coh, and the DIAGONAL vs OFF-DIAGONAL split of
%   U^H R0 U (the interference term).
%
%   articleB_r0_wavecap_run([128 64], 100, 60)     % mesh, f Hz, nmax
%
% Writes results/articleB_r0_wavecap_*  (new files only).
function articleB_r0_wavecap_run(nsub, freq, nmax, dofigs)
if nargin<1 || isempty(nsub),  nsub=[128 64]; end
if nargin<2 || isempty(freq),  freq=100; end
if nargin<3 || isempty(nmax),  nmax=60; end
if nargin<4 || isempty(dofigs),dofigs=true; end
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(RES,'dir'), mkdir(RES); end
tag=sprintf('%dx%d_f%g',nsub(1),nsub(2),freq);

rho0=1.21; c0=343; zeta=0.05;
spec=struct('R',10,'alpha',deg2rad(30),'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,...
  'deg',[2 2],'reg',[1 1],'nsub',nsub,'nquad',[4 4],'bc',{2},...
  'loadtype','pole','loaddir',[0;0;1],'loadmag',1,'freq',freq,'zeta',zeta);

fprintf('\n==== ARTICLE B R0 WAVE-REGIME ladder, clamped spherical cap, mesh %s, f=%.6g Hz ====\n',mat2str(nsub),freq);
ts=tic;
ctx=obstacle_build_ctx('spherecap',spec);
fprintf('setup %.1f s  N=%d nel=%d  f_exc=%.6g Hz  zeta=%.2f  t0=%.3g\n',toc(ts),ctx.ndof,ctx.ne,freq,zeta,ctx.t0);

% modal content inside the design band
band=[80 130];
Ks=sparse(ctx.Rows(:),ctx.Cols(:),reshape(ctx.Bke.*ctx.t0^3+ctx.Ske.*ctx.t0,[],1),ctx.ndof,ctx.ndof);
M =sparse(ctx.Rows(:),ctx.Cols(:),reshape(ctx.Me.*ctx.t0,[],1),ctx.ndof,ctx.ndof);
[~,D]=eigs(Ks(ctx.free,ctx.free),M(ctx.free,ctx.free),60,(2*pi*105)^2);
freqs=sort(sqrt(abs(real(diag(D)))))/(2*pi);
nband=sum(freqs>=band(1)&freqs<=band(2));
nbelow=sum(freqs<=freq);
Lmax=blockmax(ctx);
fprintf('modal band [%g %g] Hz: %d modes (%.4f modes/Hz); %d modes <= f_exc\n',band(1),band(2),nband,nband/(band(2)-band(1)),nbelow);

R0op=articleB_r0_build(ctx,rho0,c0,128,freq);
kLmax=R0op.k*Lmax;
sv=svd(R0op.KSS); relK1=sv(2)/sv(1);
Kn=R0op.KS.'*pinv((R0op.KSS+R0op.KSS.')/2)*R0op.KS;
relKfro=norm(Kn-ones(ctx.ne),'fro')/norm(Kn,'fro');
mu=measure2(ctx,ctx.t0*ones(ctx.ne,1),R0op);
fprintf('R0 (3.66) at f=%.6g Hz: k=%.4e  Lmax=%.3f m  k*Lmax=%.4f  wavelength=%.3f m\n',freq,R0op.k,Lmax,kLmax,2*pi/R0op.k);
fprintf('  kernel: s2/s1(KSS)=%.4g  ||K-ones||/||K||=%.4g\n',relK1,relKfro);
fprintf('  GATE (uniform design): rho_coh=%.5f  coh_gamma=%.5f  od/tot=%.5f\n',mu.rho_coh,mu.coh2,mu.odfrac);

% ---- quick FD self-check of the gradient on this geometry/frequency -------
fdcheck(ctx,R0op);

opt=struct('tmin',0.25,'tmax',1.0,'volfrac',0.5,'nmax',nmax,'stop',true);
jobs = { 'dyncomp','dyncomp',''; 'inputpower','inputpower',''; ...
         'quadvel-rad','quadvel','rad'; 'quadvel-I','quadvel','I'; 'quadvel-ones','quadvel','ones' };
nJ=size(jobs,1);
L=struct('label',{},'metric',{},'R0mode',{},'tbest',{},'hist',{},'meta',{},'meas',{},'meas0',{});
for j=1:nJ
  lab=jobs{j,1}; met=jobs{j,2}; md=jobs{j,3};
  fprintf('\n---- %d/%d : %s ----\n',j,nJ,lab);
  [tb,h,m]=articleB_r0_opt_one(ctx,met,R0op,md,opt);
  L(j).label=lab; L(j).metric=met; L(j).R0mode=md; L(j).tbest=tb; L(j).hist=h; L(j).meta=m;
  L(j).meas=measure2(ctx,tb,R0op); L(j).meas0=measure2(ctx,opt.volfrac*opt.tmax*ones(ctx.ne,1),R0op);
  fprintf('   %s DONE: Jbest=%.6e iters=%d wall=%.1fs vol=%.4f\n',lab,m.Jbest,m.niter,m.walltime,m.final_vol);
  save(fullfile(RES,sprintf('articleB_r0_wavecap_partial_%s.mat',tag)),'L','j','-v7.3');
end
save(fullfile(RES,sprintf('articleB_r0_wavecap_metrics_%s.mat',tag)),'L','jobs','nsub','nmax','freq','Lmax','kLmax','relK1','relKfro','band','nband','freqs','-v7.3');

%% ================= comparison table =====================================
hdr={'design','Jbest(obj)','obj(dB)','W_rad(3.66)','W_rad/unif','W_I(R0=I)','W_I/unif',...
     'W_ones','W_ones/unif','pwCentre','pwCen/unif','rho_coh','coh_gamma',...
     'diag_rad','offdiag_rad','od/tot'};
data=cell(nJ,numel(hdr));
for j=1:nJ
  Mx=L(j).meas;
  data(j,:)={L(j).label, sprintf('%.6e',L(j).meta.Jbest), sprintf('%.3f',10*log10(max(L(j).meta.Jbest,eps))),...
    sprintf('%.6e',Mx.Wrad), sprintf('%.4f',Mx.Wrad/mu.Wrad),...
    sprintf('%.6e',Mx.WI),   sprintf('%.4f',Mx.WI/mu.WI),...
    sprintf('%.6e',Mx.Wones),sprintf('%.4f',Mx.Wones/mu.Wones),...
    sprintf('%.6e',Mx.pwCen),sprintf('%.4f',Mx.pwCen/mu.pwCen),...
    sprintf('%.4f',Mx.rho_coh), sprintf('%.4f',Mx.coh2),...
    sprintf('%.4e',Mx.diag), sprintf('%.4e',Mx.offdiag), sprintf('%.4f',Mx.odfrac)};
end
fmt='%-13s %13s %9s %13s %11s %13s %9s %13s %11s %11s %11s %8s %8s %12s %12s %8s\n';
fprintf('\n==== ARTICLE B R0 WAVE REGIME (cap) COMPARISON (mesh %s, f=%.6g Hz, kLmax=%.4f) ====\n',mat2str(nsub),freq,kLmax);
fprintf(fmt,hdr{:}); for j=1:nJ, fprintf(fmt,data{j,:}); end
fid=fopen(fullfile(RES,sprintf('articleB_r0_wavecap_table_%s.txt',tag)),'w');
fprintf(fid,'# Article B phase-cancellation ladder, WAVE regime, clamped spherical cap (R=10, 30deg half-angle)\n');
fprintf(fid,'# mesh %s, N=%d, designvars=%d, f_exc=%.6g Hz, zeta=%.2f, volfrac=0.5, t0=%.3g\n',mat2str(nsub),ctx.ndof,ctx.ne,freq,zeta,ctx.t0);
fprintf(fid,'# modal band [%g %g] Hz: %d modes (%.4f modes/Hz); %d modes <= f_exc\n',band(1),band(2),nband,nband/(band(2)-band(1)),nbelow);
fprintf(fid,'# fluid rho0=%.3g c0=%.3g, k=%.4e, Lmax=%.4f m, k*Lmax=%.4f, wavelength=%.3f m\n',rho0,c0,R0op.k,Lmax,kLmax,2*pi/R0op.k);
fprintf(fid,'# kernel rank-1: s2/s1(KSS)=%.4g, ||K-ones||/||K||=%.4g\n',relK1,relKfro);
fprintf(fid,'# W_rad=V^H R0 V (eq 3.66); W_I=w^2||U||^2 (R0=I, phase-BLIND); W_ones=pref*w^2|sum(CU)|^2 (net volume velocity)\n');
fprintf(fid,'# rho_coh=|sum V|/sum|V|; coh_gamma=|sumV|^2/(ne*sum|V|^2); V=(C u)_e centre normal velocity\n');
fprintf(fid,'# diag_rad=sum|z_e|^2 (z=area.*V); offdiag_rad=z^H K z - diag; od/tot=interference fraction\n');
fprintf(fid,fmt,hdr{:}); for j=1:nJ, fprintf(fid,fmt,data{j,:}); end
fclose(fid);
fprintf('table : %s\n',fullfile(RES,sprintf('articleB_r0_wavecap_table_%s.txt',tag)));

%% ============ W_rad vs frequency (uniform) + rank-1 ======================
fscan=logspace(log10(max(1,freq/8)), log10(freq*4), 15);
fid=fopen(fullfile(RES,sprintf('articleB_r0_wavecap_freqscan_%s.txt',tag)),'w');
fprintf(fid,'# W_rad (uniform design) and kernel structure vs frequency, clamped spherical cap, mesh %s, ne=%d\n',mat2str(nsub),ctx.ne);
fprintf(fid,'# f_Hz  k*Lmax  s2/s1(KSS)  rho_coh  W_rad_unif  W_I_unif  od/tot\n');
tu=ctx.t0*ones(ctx.ne,1);
for i=1:numel(fscan)
  ctx.w=2*pi*fscan(i); ctx.alpha=2*zeta*ctx.w; ctx.beta=0;
  Ro=articleB_r0_build(ctx,rho0,c0,128,fscan(i));
  svi=svd(Ro.KSS);
  Mi=measure2(ctx,tu,Ro);
  fprintf(fid,'%12.5g %9.4f %12.5g %10.5f %13.5e %13.5e %9.4f\n',...
    fscan(i),Ro.k*Lmax,svi(2)/svi(1),Mi.rho_coh,Mi.Wrad,Mi.WI,Mi.odfrac);
end
fclose(fid);
fprintf('freqscan: %s\n',fullfile(RES,sprintf('articleB_r0_wavecap_freqscan_%s.txt',tag)));

if dofigs
  articleB_r0_wavecap_figs(tag);
end
fprintf('\nARTICLEB_R0_WAVECAP_RUN_DONE\n');
end

% ============================ local functions =============================
function M=measure2(ctx,t,R0op)
[~,~,aux]=articleB_r0_eval(ctx,t(:),'quadvel',R0op,'I');
u=aux.u; Cu=R0op.C*u; z=R0op.area.*Cu;
Kz=R0op.KS.'*(R0op.KSSinv*(R0op.KS*z));
tot=real(z'*Kz); dg=sum(abs(z).^2); od=tot-dg;
M.WI    = ctx.w^2*real(u'*u);
M.Wones = ctx.w^2*R0op.pref*abs(sum(Cu))^2;
M.Wrad  = ctx.w^2*R0op.pref*tot;
M.pwDOF = ctx.w^2*real(u'*u);
M.pwCen = ctx.w^2*real(Cu'*Cu);
M.rho_coh = abs(sum(Cu))/max(sum(abs(Cu)),eps);
M.coh2  = min(1, abs(sum(Cu))^2/(R0op.ne*max(real(Cu'*Cu),eps)));
M.diag=dg; M.offdiag=od; M.odfrac=od/max(abs(tot),eps);
M.u=u; M.Cu=Cu;
end
function fdcheck(ctx,R0op)
% central-difference verification of the analytic gradient on THIS geometry/freq.
% NOTE the sampled Nystrom kernel KSS is severely ill-conditioned at large
% k*Lmax (cond ~ 1e19); finite-DIFFERENCING the Nystrom operator is therefore
% roundoff-dominated at small h, so we use h=1e-3 (operator-accurate regime).
% The rigorous dense-kernel FD is in tests/articleB_r0_wavecap_fdcheck.m /
% _fddiag.m; here we only need an in-run sanity flag.
ne=ctx.ne; rng(3); pe=randperm(ne,5); pert=zeros(ne,1); pert(pe)=1; h=1e-3;
modes={'I','ones','rad'};
fprintf('  FD self-check (cap, h=1e-3; cond(KSS)=%.1e): ',cond(R0op.KSS));
for im=1:numel(modes)
  [J,dJ]=articleB_r0_eval(ctx,ctx.t0*ones(ne,1),'quadvel',R0op,modes{im});
  Jf=@(tt) Wref(ctx,tt,R0op,modes{im});
  tu=ctx.t0*ones(ne,1);
  dfd=(Jf(tu+h*pert)-Jf(tu-h*pert))/(2*h);
  eu=abs(sum(dJ(pe))-dfd)/max(abs(dfd),eps);
  fprintf('%s:%.1e ',modes{im},eu);
end
fprintf('\n');
end
function W=Wref(ctx,t,R0op,mode)
t=t(:); ne=ctx.ne; lm=ctx.lm; n=size(lm,2); ndof=ctx.ndof;
cK=1+1i*ctx.w*ctx.beta; cM=1i*ctx.w*ctx.alpha-ctx.w^2;
EK=ctx.Bke.*(t.^3)+ctx.Ske.*t; EK=EK(:); EM=ctx.Me.*t; EM=EM(:);
Ks=sparse(ctx.Rows(:),ctx.Cols(:),EK,ndof,ndof); M=sparse(ctx.Rows(:),ctx.Cols(:),EM,ndof,ndof);
KD=cK*Ks+cM*M; u=zeros(ndof,1); u(ctx.free)=KD(ctx.free,ctx.free)\ctx.F(ctx.free);
Cu=R0op.C*u;
switch lower(mode)
  case 'i',    W=ctx.w^2*real(u'*u);
  case 'ones', W=ctx.w^2*R0op.pref*abs(sum(Cu))^2;
  case 'rad',  z=R0op.area.*Cu; Kz=R0op.KS.'*(R0op.KSSinv*(R0op.KS*z));
               W=ctx.w^2*R0op.pref*real(z'*Kz);
end
end
function d=blockmax(ctx)
m=ctx.msh; Jw=m.jacdet.*m.quad_weights; area=sum(Jw,1);
Jw3=reshape(Jw,1,[],ctx.ne); gm=reshape(m.geo_map,3,[],ctx.ne);
X=reshape(squeeze(sum(gm.*Jw3,2))./area,3,ctx.ne);
n=size(X,2); dmax=0; blk=512;
for i=1:blk:n
  ii=i:min(i+blk-1,n);
  DX=X(1,ii).'-X(1,:); DY=X(2,ii).'-X(2,:); DZ=X(3,ii).'-X(3,:);
  DD=DX.^2+DY.^2+DZ.^2; dmax=max(dmax,max(DD(:)));
end
d=sqrt(dmax);
end
