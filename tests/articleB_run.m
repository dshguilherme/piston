% ARTICLEB_RUN  Article B benchmark: harmonic topology optimization of a
% CURVED shell (Scordelis-Lo roof) for the three objectives (dynamic
% compliance, active input power, quadratic velocity), with an explicit
% NON-DEGENERACY audit and readable figures.
%
%   articleB_run([128 64], 100)
%
% Writes results/articleB_*  (figures, .mat histories, .txt table).
% Does NOT modify src/ (uses the verified src/gcmma_solver.m) and does NOT
% touch geopdes/immersed-iga sources.

function articleB_run(nsub, nmax)
if nargin<1 || isempty(nsub), nsub=[128 64]; end
if nargin<2 || isempty(nmax), nmax=100; end
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');

RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(RES,'dir'), mkdir(RES); end
tag=sprintf('%dx%d',nsub(1),nsub(2));

% ---- benchmark definition (curved Scordelis-Lo roof shell) -----------------
E=4.32e8; NU=0.3; RHO=2500; t0=0.25;
fprintf('\n==== ARTICLE B curved-shell benchmark, mesh %s ====\n',tag);
tS=tic; ctx0=articleB_build_ctx(nsub,E,NU,RHO,t0,NaN,0);
% first elastic eigenvalue => excitation frequency near resonance
Ks0=sparse(ctx0.Rows(:),ctx0.Cols(:),reshape(ctx0.Bke.*t0^3+ctx0.Ske.*t0,[],1),ctx0.ndof,ctx0.ndof);
M0 =sparse(ctx0.Rows(:),ctx0.Cols(:),reshape(ctx0.Me(:)*t0,[],1),ctx0.ndof,ctx0.ndof);
[~,D]=eigs(Ks0(ctx0.free,ctx0.free),M0(ctx0.free,ctx0.free),6,'sm');
freqs=sort(sqrt(abs(diag(D))))/(2*pi);
f1=freqs(1); freq=0.97*f1; zeta=0.05;
ctx=articleB_build_ctx(nsub,E,NU,RHO,t0,freq,zeta);
fprintf('setup %.1f s  N=%d nel=%d  f1=%.4f Hz  excitation=%.4f Hz  zeta=%.2f  alpha=%.5g\n',...
        toc(tS),ctx.ndof,ctx.ne,f1,freq,zeta,ctx.alpha);
fprintf('natural freqs (Hz): %s\n',mat2str(freqs.',6));

opt=struct('tmin',0.25,'tmax',1.0,'volfrac',0.5,'nmax',nmax,'stop',true);
metrics={'dyncomp','inputpower','quadvel'};
nM=numel(metrics);
R=struct('metric',{},'tbest',{},'hist',{},'meta',{},'deg',{});
for k=1:nM
  fprintf('\n---- metric %d/%d : %s ----\n',k,nM,metrics{k});
  [tb,h,m]=articleB_shell_opt_one(ctx,metrics{k},opt);
  R(k).metric=metrics{k}; R(k).tbest=tb; R(k).hist=h; R(k).meta=m;
  R(k).deg=degeneracy_metrics(tb,nsub(1),nsub(2),opt.tmin,opt.tmax);
  fprintf('  %s DONE: Jbest=%.6e dB=%.3f iters=%d evals=%d wall=%.1fs vol=%.4f tail10=%.4f (%s)\n',...
    metrics{k},m.Jbest,10*log10(max(m.Jbest,eps)),m.niter,m.numeval,m.walltime,m.final_vol,m.tailspread10,m.stopReason);
end

save(fullfile(RES,sprintf('articleB_shell_metrics_%s.mat',tag)),'R','metrics','nsub','nmax','freq','f1','freqs','ctx');

% ---- per-metric figures ----------------------------------------------------
paths={};
for k=1:nM
  hh=R(k).hist; dg=R(k).deg;
  fig=figure('Position',[40 40 1750 560],'Color','w');
  subplot(1,3,1);
  tg=reshape(R(k).tbest,nsub(1),nsub(2));
  imagesc(tg.'); axis image; set(gca,'YDir','normal');
  cb=colorbar; cb.Label.String='thickness t (m)'; caxis([opt.tmin opt.tmax]);
  xlabel('element i  (circumferential \xi)'); ylabel('element j  (axial \eta)');
  title(sprintf('%s: best thickness field\nJ_{best}=%.4e (%.2f dB)  CV=%.2f  frac@mid=%.2f',...
        R(k).metric,R(k).meta.Jbest,10*log10(max(R(k).meta.Jbest,eps)),dg.CV,dg.frac_mid),'Interpreter','none');
  subplot(1,3,2);
  ii=1:numel(hh.J);
  semilogy(ii,hh.J,'.','MarkerSize',9,'DisplayName','running J'); hold on;
  semilogy(ii,hh.Jbest,'-','LineWidth',2,'DisplayName','best J'); grid on;
  xlabel('outer iteration'); ylabel('objective J  (raw, log scale)');
  legend('Location','northeast');
  title(sprintf('%s: objective trajectory',R(k).metric),'Interpreter','none');
  subplot(1,3,3);
  plot(ii,hh.vol,'-o','LineWidth',1.2,'MarkerSize',3); hold on;
  yline(opt.volfrac,'k--','volfrac'); grid on; ylim([0.45 0.55]);
  xlabel('outer iteration'); ylabel('volume fraction \Sigma t / \Sigma t_{max}');
  title(sprintf('%s: volume-fraction history',R(k).metric),'Interpreter','none');
  fp=fullfile(RES,sprintf('articleB_design_%s_%s.png',R(k).metric,tag));
  print(fig,'-dpng','-r140',fp); close(fig); paths{end+1}=fp; %#ok<AGROW>
  fprintf('figure: %s\n',fp);
end

% ---- side-by-side comparison ----------------------------------------------
fig=figure('Position',[40 40 1750 560],'Color','w');
for k=1:nM
  subplot(1,nM,k);
  imagesc(reshape(R(k).tbest,nsub(1),nsub(2)).'); axis image; set(gca,'YDir','normal');
  caxis([opt.tmin opt.tmax]); cb=colorbar; cb.Label.String='thickness t (m)';
  xlabel('element i  (\xi)'); ylabel('element j  (\eta)');
  title(sprintf('%s\nJ_{best}=%.4e (%.2f dB)\nCV=%.2f  frac@mid=%.2f',R(k).metric,...
        R(k).meta.Jbest,10*log10(max(R(k).meta.Jbest,eps)),R(k).deg.CV,R(k).deg.frac_mid),'Interpreter','none');
end
sgtitle(sprintf('Best thickness fields, curved Scordelis-Lo roof shell, mesh %s, N=%d, f=%.3f Hz',tag,ctx.ndof,freq));
fp=fullfile(RES,sprintf('articleB_compare_%s.png',tag)); print(fig,'-dpng','-r140',fp); close(fig);
paths{end+1}=fp; fprintf('figure: %s\n',fp);

% ---- degeneracy histograms -------------------------------------------------
fig=figure('Position',[40 40 1300 400],'Color','w');
colors=lines(nM);
for k=1:nM
  subplot(1,nM,k); hold on;
  histogram(R(k).tbest,20,'FaceColor',colors(k,:));
  xline(opt.tmin,'k--'); xline(opt.tmax,'k--'); xline(mean(R(k).tbest),'r-');
  xlabel('thickness t'); ylabel('count'); xlim([opt.tmin-0.02 opt.tmax+0.02]);
  title(sprintf('%s: thickness histogram\nstd=%.3f  CV=%.2f  conc_{1%%}=%.2f',...
        R(k).metric,R(k).deg.std,R(k).deg.CV,R(k).deg.conc1pct),'Interpreter','none');
end
fp=fullfile(RES,sprintf('articleB_degeneracy_%s.png',tag)); print(fig,'-dpng','-r140',fp); close(fig);
paths{end+1}=fp; fprintf('figure: %s\n',fp);

% ---- table -----------------------------------------------------------------
hdr={'metric','best J (raw)','best J (dB)','iters','GCMMA evals','wall time(s)',...
     'final volfrac','tail spread(10)','stop','-','frac@tmin','frac@tmax','frac@mid',...
     'std(t)','CV','span(hot)','conc_1pct'};
data=cell(nM,numel(hdr));
for k=1:nM
  dg=R(k).deg;
  data(k,:)={R(k).metric, sprintf('%.6e',R(k).meta.Jbest), sprintf('%.3f',10*log10(max(R(k).meta.Jbest,eps))),...
    sprintf('%d',R(k).meta.niter), sprintf('%d',R(k).meta.numeval), sprintf('%.1f',R(k).meta.walltime),...
    sprintf('%.4f',R(k).meta.final_vol), sprintf('%.4f',R(k).meta.tailspread10), R(k).meta.stopReason, '|',...
    sprintf('%.3f',dg.frac_tmin), sprintf('%.3f',dg.frac_tmax), sprintf('%.3f',dg.frac_mid),...
    sprintf('%.4f',dg.std), sprintf('%.3f',dg.CV), sprintf('%.3f',dg.hot_span), sprintf('%.3f',dg.conc1pct)};
end
fprintf('\n==== ARTICLE B COMPARISON TABLE (mesh %s) ====\n',tag);
for k=1:nM, fprintf('%-11s %14s %11s %5s %11s %11s %11s %13s %9s %1s %8s %8s %8s %8s %6s %9s %9s\n',data{k,:}); end
fid=fopen(fullfile(RES,sprintf('articleB_metrics_table_%s.txt',tag)),'w');
fprintf(fid,'# Article B: curved Scordelis-Lo roof shell, mesh %s, N=%d, f=%.4f Hz, zeta=%.2f\n',tag,ctx.ndof,freq,zeta);
fprintf(fid,'%-11s %14s %11s %5s %11s %11s %11s %13s %9s %1s %8s %8s %8s %8s %6s %9s %9s\n',hdr{:});
for k=1:nM, fprintf(fid,'%-11s %14s %11s %5s %11s %11s %11s %13s %9s %1s %8s %8s %8s %8s %6s %9s %9s\n',data{k,:}); end
fclose(fid);
fprintf('table : %s\n',fullfile(RES,sprintf('articleB_metrics_table_%s.txt',tag)));
fprintf('\nALL FIGURES:\n'); for i=1:numel(paths), fprintf('  %s\n',paths{i}); end
disp('ARTICLEB_RUN_DONE');
end

% =========================== local functions ===============================
function dg=degeneracy_metrics(t,nel1,nel2,tmin,tmax)
t=t(:); n=numel(t); m=mean(t); s=std(t);
dev=abs(t-m); [~,ord]=sort(dev,'descend');
n1=max(1,round(0.01*n)); conc1pct=sum(dev(ord(1:n1)))/max(sum(dev),eps);
% hot set: t > mean+std ; connected components (4-neighbour) + spatial span
tg=reshape(t,nel1,nel2); hot=tg> (m+s);
lab=zeros(nel1,nel2); cur=0;
for i=1:nel1, for j=1:nel2
  if hot(i,j) && lab(i,j)==0, cur=cur+1; stack=[i j];
    while ~isempty(stack)
      p=stack(end,:); stack(end,:)=[]; lab(p(1),p(2))=cur;
      dirs=[-1 0;1 0;0 -1;0 1];
      for di=1:4
        q=p+dirs(di,:); if q(1)>=1&&q(1)<=nel1&&q(2)>=1&&q(2)<=nel2&&hot(q(1),q(2))&&lab(q(1),q(2))==0
          stack(end+1,:)=q; end %#ok<AGROW>
      end
    end
  end
end, end
[I,J]=find(hot);
if numel(I)>=2
  hot_span=sqrt(var(I/(nel1))+var(J/(nel2)));   % normalised spatial spread
else, hot_span=0; end
dg=struct('mean',m,'std',s,'CV',s/max(m,eps),'range',[min(t) max(t)],...
  'frac_tmin',mean(t<=tmin+1e-6),'frac_tmax',mean(t>=tmax-1e-6),...
  'frac_mid',mean(t>tmin+1e-6 & t<tmax-1e-6),...
  'conc1pct',conc1pct,'nblobs',cur,'nhot',numel(I),'hot_span',hot_span);
end
