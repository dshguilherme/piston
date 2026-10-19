% ARTICLEB_R0_WAVE_RUN  Phase-cancellation ladder in the acoustic WAVE regime.
%
% Operates the THREE-operator ladder at an excitation frequency chosen so that
% k*Lmax is O(2-3) (acoustic wave regime, the sinc kernel R0 is genuinely
% structured, no longer rank-1) AND the surface-velocity field is incoherent
% (rho_coh clearly < 1) -- the regime where phase cancellation can physically
% exist.  Everything else (geometry, mesh, volume constraint, damping, GCMMA
% settings) is exactly as in articleB_r0_run.
%
% Reports for every design: W_rad (R0=3.66), W under R0=I, W under the
% point-source (ones) operator, pointwise velocity norms, rho_coh, and the
% DIAGONAL vs OFF-DIAGONAL split of U^H R0 U (the interference term).
%
%   articleB_r0_wave_run([128 64], 2.0, 60)     % mesh, freq Hz, nmax
%
% Writes results/articleB_r0_wave_*  (new files only).
function articleB_r0_wave_run(nsub, freq, nmax, dofigs)
if nargin<1 || isempty(nsub),  nsub=[128 64]; end
if nargin<2 || isempty(freq),  freq=2.0; end
if nargin<3 || isempty(nmax),  nmax=60; end
if nargin<4 || isempty(dofigs),dofigs=true; end
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(RES,'dir'), mkdir(RES); end
tag=sprintf('%dx%d_f%g',nsub(1),nsub(2),freq);

E=4.32e8; NU=0.3; RHO=2500; t0=0.25; rho0=1.21; c0=343; zeta=0.05;

fprintf('\n==== ARTICLE B R0 WAVE-REGIME ladder, mesh %s, f=%.6g Hz ====\n',mat2str(nsub),freq);
ts=tic;
ctx=articleB_build_ctx(nsub,E,NU,RHO,t0,freq,zeta);
fprintf('setup %.1f s  N=%d nel=%d  f_exc=%.6g Hz  zeta=%.2f\n',toc(ts),ctx.ndof,ctx.ne,freq,zeta);

R0op=articleB_r0_build(ctx,rho0,c0,128,freq);
Lmax=blockmax(ctx);
kLmax=R0op.k*Lmax;
sv=svd(R0op.KSS); relK1=sv(2)/sv(1);
if ~isempty(R0op.Kfull)
  relKfro=norm(R0op.Kfull-ones(ctx.ne),'fro')/norm(R0op.Kfull,'fro');
else
  Kn=R0op.KS.'*pinv((R0op.KSS+R0op.KSS.')/2)*R0op.KS;
  relKfro=norm(Kn-ones(ctx.ne),'fro')/norm(Kn,'fro');
end
mu=measure2(ctx,t0*ones(ctx.ne,1),R0op);
fprintf('R0 (3.66) at f=%.6g Hz: k=%.4e  Lmax=%.3f m  k*Lmax=%.4f  lambda/k=%g\n',freq,R0op.k,Lmax,kLmax,2*pi/R0op.k);
fprintf('  kernel: s2/s1(KSS)=%.4g  ||K-ones||/||K||=%.4g  (rank-1 => 0)\n',relK1,relKfro);
fprintf('  GATE (uniform design): rho_coh=%.5f  coh_gamma=%.5f  od/tot=%.5f\n',mu.rho_coh,mu.coh2,mu.odfrac);

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
  save(fullfile(RES,sprintf('articleB_r0_wave_partial_%s.mat',tag)),'L','j','-v7.3');
end
save(fullfile(RES,sprintf('articleB_r0_wave_metrics_%s.mat',tag)),'L','jobs','nsub','nmax','freq','Lmax','kLmax','relK1','relKfro','-v7.3');

%% ================= comparison table =====================================
hdr={'design','Jbest(obj)','obj(dB)','W_rad(3.66)','W_rad/unif','W_I(R0=I)','W_I/unif',...
     'W_ones','W_ones/unif','pwCentre','pwCen/unif','rho_coh','coh_gamma',...
     'diag_rad','offdiag_rad','od/tot'};
data=cell(nJ,numel(hdr));
for j=1:nJ
  M=L(j).meas; mu0=mu;
  data(j,:)={L(j).label, sprintf('%.6e',L(j).meta.Jbest), sprintf('%.3f',10*log10(max(L(j).meta.Jbest,eps))),...
    sprintf('%.6e',M.Wrad), sprintf('%.4f',M.Wrad/mu0.Wrad),...
    sprintf('%.6e',M.WI),   sprintf('%.4f',M.WI/mu0.WI),...
    sprintf('%.6e',M.Wones),sprintf('%.4f',M.Wones/mu0.Wones),...
    sprintf('%.6e',M.pwCen),sprintf('%.4f',M.pwCen/mu0.pwCen),...
    sprintf('%.4f',M.rho_coh), sprintf('%.4f',M.coh2),...
    sprintf('%.4e',M.diag), sprintf('%.4e',M.offdiag), sprintf('%.4f',M.odfrac)};
end
fmt='%-13s %13s %9s %13s %11s %13s %9s %13s %11s %11s %11s %8s %8s %12s %12s %8s\n';
fprintf('\n==== ARTICLE B R0 WAVE REGIME COMPARISON (mesh %s, f=%.6g Hz, kLmax=%.4f) ====\n',mat2str(nsub),freq,kLmax);
fprintf(fmt,hdr{:}); for j=1:nJ, fprintf(fmt,data{j,:}); end
fid=fopen(fullfile(RES,sprintf('articleB_r0_wave_table_%s.txt',tag)),'w');
fprintf(fid,'# Article B phase-cancellation ladder, WAVE regime, curved Scordelis-Lo roof\n');
fprintf(fid,'# mesh %s, N=%d, designvars=%d, f_exc=%.6g Hz, zeta=%.2f, volfrac=0.5\n',mat2str(nsub),ctx.ndof,ctx.ne,freq,zeta);
fprintf(fid,'# fluid rho0=%.3g c0=%.3g, k=%.4e, Lmax=%.4f m, k*Lmax=%.4f, wavelength=%.2f m\n',rho0,c0,R0op.k,Lmax,kLmax,2*pi/R0op.k);
fprintf(fid,'# kernel rank-1 indicators: s2/s1(KSS)=%.4g, ||K-ones||/||K||=%.4g\n',relK1,relKfro);
fprintf(fid,'# W_rad=V^H R0 V (eq 3.66); W_I=w^2||U||^2 (R0=I, phase-BLIND); W_ones=pref*w^2|sum(CU)|^2 (net volume velocity)\n');
fprintf(fid,'# rho_coh=|sum V|/sum|V| (1=coherent); coh_gamma=|sumV|^2/(ne*sum|V|^2); V=(C u)_e centre normal velocity\n');
fprintf(fid,'# diag_rad=sum|z_e|^2 (z=area.*V); offdiag_rad=z^H K z - diag; od/tot=interference fraction\n');
fprintf(fid,fmt,hdr{:}); for j=1:nJ, fprintf(fid,fmt,data{j,:}); end
fclose(fid);
fprintf('table : %s\n',fullfile(RES,sprintf('articleB_r0_wave_table_%s.txt',tag)));

%% ============ W_rad vs frequency, uniform design + rank-1 ================
fscan=logspace(log10(freq), log10(20), 13);
fid=fopen(fullfile(RES,sprintf('articleB_r0_wave_freqscan_%s.txt',tag)),'w');
fprintf(fid,'# W_rad (uniform design) and kernel structure vs frequency, mesh %s, ne=%d\n',mat2str(nsub),ctx.ne);
fprintf(fid,'# f_Hz  k*Lmax  s2/s1(KSS)  rho_coh  W_rad_unif  W_I_unif  od/tot\n');
tu=t0*ones(ctx.ne,1);
for i=1:numel(fscan)
  ctx.w=2*pi*fscan(i); ctx.alpha=2*zeta*ctx.w; ctx.beta=0;
  Ro=articleB_r0_build(ctx,rho0,c0,128,fscan(i));
  svi=svd(Ro.KSS);
  Mi=measure2(ctx,tu,Ro);
  fprintf(fid,'%12.5g %9.4f %12.5g %10.5f %13.5e %13.5e %9.4f\n',...
    fscan(i),Ro.k*Lmax,svi(2)/svi(1),Mi.rho_coh,Mi.Wrad,Mi.WI,Mi.odfrac);
end
fclose(fid);
fprintf('freqscan: %s\n',fullfile(RES,sprintf('articleB_r0_wave_freqscan_%s.txt',tag)));

if dofigs
  articleB_r0_wave_figs(tag);
end
fprintf('\nARTICLEB_R0_WAVE_RUN_DONE\n');
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
