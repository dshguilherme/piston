% ARTICLEB_R0_RUN  Phase-cancellation ladder for the quadratic-velocity /
% radiated-power objective on the curved Scordelis-Lo roof shell.
%
% Builds the thesis radiation-resistance R0 (eqs. 3.63-3.66) from knot-span
% centres, FD-verified in tests/articleB_r0_fdcheck.m, and runs the harmonic
% topology optimization for five designs:
%   dyncomp, inputpower           (do NOT use R0; excitation-point metrics)
%   quadvel full R0 (3.66)        ("surface-averaged quadratic velocity" = W_rad)
%   quadvel R0 = I                (phase-BLIND control, thesis R0=I tables)
%   quadvel R0 = rank-1 ones      (net volume velocity; r->0 point-source limit)
% then reports every design under all three operators + the pointwise norms.
%
%   articleB_r0_run([64 32], 60)
%
% Writes results/articleB_r0_*  (figures PNG+PDF, table, .mat). New files only.

function articleB_r0_run(nsub, nmax, dofigs)
if nargin<1 || isempty(nsub), nsub=[256 128]; end
if nargin<2 || isempty(nmax), nmax=60; end
if nargin<3 || isempty(dofigs), dofigs=true; end
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(RES,'dir'), mkdir(RES); end
tag=sprintf('%dx%d',nsub(1),nsub(2));
S=pub_style();

E=4.32e8; NU=0.3; RHO=2500; t0=0.25; rho0=1.21; c0=343;
fprintf('\n==== ARTICLE B R0 phase-cancellation ladder, mesh %s ====\n',tag);
tS=tic;
% f1 from a COARSE mesh (eigenvalue is mesh-robust) to keep the fine setup cheap
ctxc=articleB_build_ctx([32 16],E,NU,RHO,t0,NaN,0);
Kc=sparse(ctxc.Rows(:),ctxc.Cols(:),reshape(ctxc.Bke.*t0^3+ctxc.Ske.*t0,[],1),ctxc.ndof,ctxc.ndof);
Mc=sparse(ctxc.Rows(:),ctxc.Cols(:),reshape(ctxc.Me(:)*t0,[],1),ctxc.ndof,ctxc.ndof);
[~,D]=eigs(Kc(ctxc.free,ctxc.free),Mc(ctxc.free,ctxc.free),6,'sm');
freqs=sort(sqrt(abs(diag(D))))/(2*pi); f1=freqs(1); freq=0.97*f1; zeta=0.05;
ctx=articleB_build_ctx(nsub,E,NU,RHO,t0,freq,zeta);
clear ctxc Kc Mc D
fprintf('setup %.1f s  N=%d nel=%d designvars=%d  f1=%.6g Hz  excitation=%.6g Hz  zeta=%.2f\n',...
        toc(tS),ctx.ndof,ctx.ne,ctx.ne,f1,freq,zeta);

R0op=articleB_r0_build(ctx,rho0,c0,128);
kLmax=R0op.k*maxdist(R0op.centres);
fprintf('R0 (3.66) at f=%.6g Hz: rho0=%.3g c0=%.3g k=%.4e k*Lmax=%.4e  sinc(k*Lmax)=%.6f\n',...
        freq,rho0,c0,R0op.k,kLmax, my_sinc(kLmax));
sv=svd(R0op.KSS);
if ~isempty(R0op.Kfull)
  relRank1=norm(R0op.Kfull-ones(ctx.ne),'fro')/norm(R0op.Kfull,'fro');
  Knys=R0op.KS.'*(R0op.KSSinv*R0op.KS);
  fprintf('  ||K-ones||/||K||=%.4f  Nystrom(m=%d) rel err vs dense=%.2e\n',...
          relRank1,R0op.mrank,norm(Knys-R0op.Kfull,'fro')/norm(R0op.Kfull,'fro'));
else
  fprintf('  rank-1 indicator s2/s1(KSS)=%.3e (>0 => kernel departs from rank-1)\n', sv(2)/sv(1));
end

opt=struct('tmin',0.25,'tmax',1.0,'volfrac',0.5,'nmax',nmax,'stop',true);
jobs = { 'dyncomp','dyncomp',''; 'inputpower','inputpower',''; ...
         'quadvel-rad','quadvel','rad'; 'quadvel-I','quadvel','I'; 'quadvel-ones','quadvel','ones' };
nJ=size(jobs,1);
L=struct('label',{},'metric',{},'R0mode',{},'tbest',{},'hist',{},'meta',{},'meas',{});
for j=1:nJ
  lab=jobs{j,1}; met=jobs{j,2}; md=jobs{j,3};
  fprintf('\n---- %d/%d : %s ----\n',j,nJ,lab);
  [tb,h,m]=articleB_r0_opt_one(ctx,met,R0op,md,opt);
  L(j).label=lab; L(j).metric=met; L(j).R0mode=md; L(j).tbest=tb; L(j).hist=h; L(j).meta=m;
  L(j).meas=measure(ctx,tb,R0op); L(j).meas0=measure(ctx,opt.volfrac*opt.tmax*ones(ctx.ne,1),R0op);
  fprintf('   %s DONE: Jbest=%.6e iters=%d wall=%.1fs vol=%.4f stir=%.1f\n',...
     lab,m.Jbest,m.niter,m.walltime,m.final_vol,m.tailspread10);
  save(fullfile(RES,sprintf('articleB_r0_partial_%s.mat',tag)),'L','j','-v7.3');  % checkpoint
end
uniform=L(1).meas0;
save(fullfile(RES,sprintf('articleB_r0_metrics_%s.mat',tag)),'L','jobs','nsub','nmax','freq','f1','freqs','-v7.3');

%% ================= comparison table =====================================
hdr={'design','Jbest(obj)','obj(dB)','W_rad(3.66)','W_rad/unif','W_I(R0=I)','W_I/unif',...
     'W_ones','W_ones/unif','pwDOF w^2||U||^2','pwDOF/unif','pwCentre','coh_gamma','cancel'};
data=cell(nJ,numel(hdr));
for j=1:nJ
  M=L(j).meas; M0=L(j).meas0; mu=L(1).meas0;
  data(j,:)={L(j).label, sprintf('%.6e',L(j).meta.Jbest), sprintf('%.3f',10*log10(max(L(j).meta.Jbest,eps))),...
    sprintf('%.6e',M.Wrad), sprintf('%.4f',M.Wrad/mu.Wrad),...
    sprintf('%.6e',M.WI),   sprintf('%.4f',M.WI/mu.WI),...
    sprintf('%.6e',M.Wones),sprintf('%.4f',M.Wones/mu.Wones),...
    sprintf('%.6e',M.pwDOF), sprintf('%.4f',M.pwDOF/mu.pwDOF),...
    sprintf('%.6e',M.pwCen), sprintf('%.4f',M.coh2), sprintf('%.4f',1-M.coh2)};
end
fmt='%-13s %13s %8s %13s %9s %13s %9s %13s %11s %13s %11s %13s %8s %8s\n';
fprintf('\n==== ARTICLE B R0 DESIGN COMPARISON (mesh %s, f=%.6g Hz) ====\n',tag,freq);
fprintf(fmt,hdr{:});
for j=1:nJ, fprintf(fmt,data{j,:}); end
fid=fopen(fullfile(RES,sprintf('articleB_r0_table_%s.txt',tag)),'w');
fprintf(fid,'# Article B phase-cancellation ladder, curved Scordelis-Lo roof, mesh %s, N=%d, designvars=%d\n',tag,ctx.ndof,ctx.ne);
fprintf(fid,'# f_exc=%.6g Hz, zeta=%.2f, volfrac=0.5, fluid rho0=%.3g c0=%.3g, k=%.4e, k*Lmax=%.4e\n',freq,zeta,rho0,c0,R0op.k,kLmax);
fprintf(fid,'# W_rad=V^H R0 V (eq 3.66); W_I=w^2||U||^2 (R0=I, phase-blind); W_ones=pref*w^2|sum(CU)|^2 (net volume velocity)\n');
fprintf(fid,'# coh_gamma=|sum(CU)|^2/(ne*sum|CU|^2) in [0,1] (1=fully coherent,0=full cancellation); cancel=1-coh_gamma\n');
fprintf(fid,fmt,hdr{:});
for j=1:nJ, fprintf(fid,fmt,data{j,:}); end
fclose(fid);
fprintf('table : %s\n',fullfile(RES,sprintf('articleB_r0_table_%s.txt',tag)));

%% ================= frequency dependence of the R0 kernel ================
fscan = logspace(log10(freq), log10(freq*80), 9);
fid=fopen(fullfile(RES,sprintf('articleB_r0_freqscan_%s.txt',tag)),'w');
fprintf(fid,'# R0 (3.66) sinc kernel structure vs assembly frequency (mesh %s, ne=%d)\n',tag,ctx.ne);
fprintf(fid,'# f_Hz  k*Lmax  sinc(kLmax)  ||K-ones||/||K|| (Nystrom, m=128)\n');
ctxs=ctx;
for i=1:numel(fscan)
  ctxs.w=2*pi*fscan(i); Ro=articleB_r0_build(ctxs,rho0,c0,128);
  kl=Ro.k*maxdist(Ro.centres);
  sv=svd(Ro.KSS); rel=sv(2)/sv(1);   % 2nd/1st singular value = rank-1 indicator
  fprintf(fid,'%12.5g %12.5g %12.6f %12.5g\n',fscan(i),kl,my_sinc(kl),rel);
end
fclose(fid);
fprintf('freqscan: %s\n',fullfile(RES,sprintf('articleB_r0_freqscan_%s.txt',tag)));

%% ================= FIGURES (white background, PNG + PDF) ================
if dofigs
lab3={'dyncomp','inputpower','quadvel-rad'};
idx3=[1 2 3]; idx5=1:5;
alllab={L.label};

% ---- F1: five designs, same colour scale ---------------------------------
fig=figure('Units','inches','Position',[1 1 12.5 2.6],'Color','w');
tl=tiledlayout(1,nJ,'TileSpacing','compact','Padding','compact');
for j=idx5
  nexttile; smooth_field(reshape(L(j).tbest,nsub(1),nsub(2))); colormap(gca,parula(256));
  caxis([opt.tmin opt.tmax]); cb=colorbar; cb.Label.String='t (m)';
  xlabel('i (\xi)'); ylabel('j (\eta)');
  title(sprintf('%s (%.2f dB)',L(j).label,10*log10(max(L(j).meta.Jbest,eps))),'Interpreter','none','FontWeight','normal');
end
sgtitle(sprintf('Optimised thickness fields, curved Scordelis-Lo roof shell, mesh %s, N=%d, f=%.4g Hz',tag,ctx.ndof,freq),...
        'Interpreter','none');
save_fig(fig, fullfile(RES,sprintf('articleB_r0_designs_%s',tag))); close(fig);

% ---- F2: three primary designs side by side, same colour scale -----------
fig=figure('Units','inches','Position',[1 1 8.0 2.6],'Color','w');
tl=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for jj=1:3
  j=idx3(jj); nexttile; smooth_field(reshape(L(j).tbest,nsub(1),nsub(2))); colormap(gca,parula(256));
  caxis([opt.tmin opt.tmax]); cb=colorbar; cb.Label.String='t (m)';
  xlabel('i (\xi)'); ylabel('j (\eta)');
  title(sprintf('%s',L(j).label),'Interpreter','none','FontWeight','normal');
end
sgtitle(sprintf('Three objectives, same colour scale (%s)',tag),'Interpreter','none');
save_fig(fig, fullfile(RES,sprintf('articleB_r0_compare_%s',tag))); close(fig);

% ---- F3: bar chart W_rad (full R0) vs pointwise norm vs R0=I -------------
fig=figure('Units','inches','Position',[1 1 8.6 3.4],'Color','w');
Rr=zeros(3,3);
for jj=1:3, j=idx3(jj);
  Rr(jj,:)=[L(j).meas.Wrad/mu.Wrad, L(j).meas.WI/mu.WI, L(j).meas.Wones/mu.Wones]; end
hold on; grid on; box on;
b=bar(1:3, Rr, 0.85,'grouped');
cols=[0.122 0.403 0.678; 0.850 0.470 0.090; 0.176 0.615 0.271];
for c=1:3, b(c).FaceColor=cols(c,:); end
for c=1:3, for jj=1:3
  text(b(c).XEndPoints(jj),Rr(jj,c),sprintf('%.3f',Rr(jj,c)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
end, end
set(gca,'XTick',1:3,'XTickLabel',lab3,'TickLabelInterpreter','none');
ylabel('W / W_{uniform}  (lower = better)');
legend({'W_{rad} (full R0, eq.3.66)','W with R0=I (phase-BLIND)','W with R0=rank1 (net vol. vel.)'},...
       'Location','northoutside','Orientation','horizontal');
title(sprintf('Cancellation ladder on optimised designs (%s)',tag),'Interpreter','none');
save_fig(fig, fullfile(RES,sprintf('articleB_r0_bars_%s',tag))); close(fig);

% ---- F4: phase/coherence field for the qv-rad design ---------------------
jrad=3; u=L(jrad).meas.u; Cu=L(jrad).meas.Cu;
gi = reshape(L(jrad).meas0.Cu, nsub(1), nsub(2));   % uniform-design centre vel
gr = reshape(Cu, nsub(1), nsub(2));
fig=figure('Units','inches','Position',[1 1 10.5 3.2],'Color','w');
tl=tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
% (a) phase map of the optimised design
nexttile; smooth_field(angle(gr)); colormap(gca,hsv(256)); caxis([-pi pi]); cb=colorbar; cb.Label.String='arg(V) [rad]';
xlabel('i (\xi)'); ylabel('j (\eta)'); title('phase of centre velocity, qv-rad','Interpreter','none','FontWeight','normal');
% (b) complex-plane sum: individual phasors vs their resultant (cancellation)
nexttile; hold on; box on; grid on;
y=Cu; y=y/ max(abs(y));
plot(real(y),imag(y),'.','Color',[0.55 0.55 0.55]);
plot([0 real(sum(y)/max(abs(y)))],[0 imag(sum(y)/max(abs(y)))],'-o','Color',[0.85 0.15 0.15],'LineWidth',2.2);
plot(cos(linspace(0,2*pi,200)),sin(linspace(0,2*pi,200)),'k:');
axis equal; xlabel('Re V (normalised)'); ylabel('Im V');
title(sprintf('phasors sum: |\\SigmaV|/\\Sigma|V| = %.3f',abs(sum(y))/sum(abs(y))),'Interpreter','tex','FontWeight','normal');
% (c) per-element contribution to the coherent sum, cumulative magnitude
nexttile; hold on; box on; grid on;
cum=abs(cumsum(sort(abs(y),'descend')));
plot(cum,'LineWidth',1.8,'Color',[0.122 0.403 0.678]);
yline(abs(sum(y))/max(abs(y))*numel(y),'r--',sprintf('|\\SigmaV| scaled'));
xlabel('elements (sorted by |V|)'); ylabel('cumulative |V|');
title('incoherent sum >> coherent sum','Interpreter','none','FontWeight','normal');
sgtitle(sprintf('Phase/coherence field, quadvel-rad design (mesh %s)',tag),'Interpreter','none');
save_fig(fig, fullfile(RES,sprintf('articleB_r0_phase_%s',tag))); close(fig);
end

fprintf('\nALL R0 FIGURES:\n');
d=dir(fullfile(RES,sprintf('articleB_r0_*%s*',tag)));
for i=1:numel(d), fprintf('  %s\n',fullfile(RES,d(i).name)); end
disp('ARTICLEB_R0_RUN_DONE');
end

% ============================ local functions =============================
function M=measure(ctx,t,R0op)
[~,~,aux]=articleB_r0_eval(ctx,t(:),'quadvel',R0op,'I');
u=aux.u; Cu=R0op.C*u; z=R0op.area.*Cu;
Kz=R0op.KS.'*(R0op.KSSinv*(R0op.KS*z));
M.WI    = ctx.w^2*real(u'*u);
M.Wones = ctx.w^2*R0op.pref*abs(sum(Cu))^2;
M.Wrad  = ctx.w^2*R0op.pref*real(z'*Kz);
M.pwDOF = ctx.w^2*real(u'*u);
M.pwCen = ctx.w^2*real(Cu'*Cu);
M.coh2  = min(1, abs(sum(Cu))^2/(R0op.ne*max(real(Cu'*Cu),eps)));
M.u=u; M.Cu=Cu;
end
function smooth_field(tg)
[n1,n2]=size(tg);
Xi=linspace(1,n1,6*n1); Yi=linspace(1,n2,6*n2);
Zi=interp2(1:n1,1:n2,tg.',Xi,Yi.','linear');
imagesc(Xi,Yi,Zi); axis image; set(gca,'YDir','normal');
end
function save_fig(fig, base)
set(fig,'Color','w','InvertHardcopy','off');
try
  print(fig,[base '.png'],'-dpng','-r300');
  print(fig,[base '.pdf'],'-dpdf','-vector');
catch err
  warning('save_fig failed for %s: %s',base,err.message);
end
end
function d=maxdist(X)
% memory-safe blockwise max pairwise distance (never forms ne x ne)
n=size(X,2); dmax=0; blk=512;
for i=1:blk:n
  ii=i:min(i+blk-1,n);
  DX=X(1,ii).'-X(1,:); DY=X(2,ii).'-X(2,:); DZ=X(3,ii).'-X(3,:);
  DD=DX.^2+DY.^2+DZ.^2; dmax=max(dmax,max(DD(:)));
end
d=sqrt(dmax);
end
function s=my_sinc(z)
s=ones(size(z)); idx=z~=0; s(idx)=sin(z(idx))./z(idx);
end
