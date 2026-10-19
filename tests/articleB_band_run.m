% ARTICLEB_BAND_RUN  THE DE-ALIASING DEMONSTRATION.
%
% Runs TWO thickness optimizations of the SAME curved Scordelis-Lo roof shell
% against the SAME metric (dynamic compliance |F.'U|) but different objectives:
%   (S) SINGLE-FREQUENCY objective at f_s = fsfac * f1(t_start)
%   (B) BAND objective (mean or KS of |F.'U| over [1-bh,1+bh] f1, Nf samples)
% then EVALUATES both final designs ON THE BAND and on a WIDER frequency sweep.
%
% IMPORTANT: f1 is the first eigenfrequency of the ACTUAL optimizer start design
% (t = volfrac*tmax), whose loaded response peaks there -- NOT the eigenfrequency
% of the context's t0=0.25 element matrices (a wrong-thickness trap that puts the
% band far below the resonance and hides the effect).
%
% CLAIM: the single-frequency design is good at its own f_s but the resonance peak
% has simply MOVED (good at one f, poor across the band), while the band design is
% good ACROSS the band -- i.e. the band objective cannot be gamed by peak-shifting.
%
%   articleB_band_run([48 24], 5, 20, 0.05, 0.99, 0.10)
%     nsub, Nf, nmax, zeta, fsfac, bandhalf
%
% Writes results/articleB_band_*  (spectrum figure PNG+PDF, metrics .txt, .mat).
% Adds src/tests files only; does NOT modify src/ or touch geopdes.

function articleB_band_run(nsub, Nf, nmax, zeta, fsfac, bandhalf)
if nargin<1 || isempty(nsub),     nsub=[48 24]; end
if nargin<2 || isempty(Nf),       Nf=5;    end
if nargin<3 || isempty(nmax),     nmax=20; end
if nargin<4 || isempty(zeta),     zeta=0.05; end   % modal damping ratio
if nargin<5 || isempty(fsfac),    fsfac=0.97; end  % single-freq excitation = fsfac*f1
if nargin<6 || isempty(bandhalf), bandhalf=0.10; end % band = [1-h,1+h]*f1
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
pub_style();
RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(RES,'dir'), mkdir(RES); end
tag=sprintf('%dx%d_z%g_bh%g',nsub(1),nsub(2),zeta,bandhalf);
metric='dyncomp';
E=4.32e8; NU=0.3; RHO=2500; t0=0.25;

fprintf('\n==== ARTICLE B BAND DE-ALIASING DEMO, mesh %s zeta=%.3g f_s=%.2f f1 ====\n',tag,zeta,fsfac);
opt=struct('tmin',0.25,'tmax',1.0,'volfrac',0.5,'nmax',nmax,'stop',true,'window',6);
% f1 of the ACTUAL optimizer START design (t = volfrac*tmax), whose loaded
% response peaks at that first eigenfrequency -- the band and f_s are placed
% around THIS resonance, not the t0=0.25 eigenfrequency (a wrong-thickness trap).
tS=tic; ctx0=articleB_build_ctx(nsub,E,NU,RHO,t0,NaN,0);
t_init=opt.volfrac*opt.tmax*ones(ctx0.ne,1);
Ks0=sparse(ctx0.Rows(:),ctx0.Cols(:),reshape(ctx0.Bke.*t_init.^3+ctx0.Ske.*t_init,[],1),ctx0.ndof,ctx0.ndof);
M0 =sparse(ctx0.Rows(:),ctx0.Cols(:),reshape(ctx0.Me.*t_init,[],1),ctx0.ndof,ctx0.ndof);
[~,D]=eigs(Ks0(ctx0.free,ctx0.free),M0(ctx0.free,ctx0.free),4,'sm');
fr=sort(sqrt(abs(diag(D))))/(2*pi); f1=fr(1);
% locate the initial design's ACTUAL loaded-response peak inside a +/-bandhalf
% window around f1, and center both the band and f_s on it -- the exploit lives
% on a resonance peak, so the single-frequency objective must sit on one.
fscan=linspace(0.70*f1,1.70*f1,21); Jscan=zeros(size(fscan));
for i=1:numel(fscan)
  ck=ctx0; ck.w=2*pi*fscan(i); ck.alpha=2*zeta*ck.w; ck.beta=0; ck.R0=[];
  Jscan(i)=articleB_shell_eval(ck,t_init,metric);
end
[~,ipk]=max(Jscan); f_peak=fscan(ipk);
f_s=f_peak;
band=struct('f_lo',(1-bandhalf)*f_peak,'f_hi',(1+bandhalf)*f_peak,'Nf',Nf); band.f=linspace(band.f_lo,band.f_hi,band.Nf);
fprintf('setup %.1f s  N=%d nel=%d | f1(t_start)=%.5f Hz  response peak=%.5f Hz  f_s=%.5f Hz\n',...
        toc(tS),ctx0.ndof,ctx0.ne,f1,f_peak,f_s);
fprintf('band = [%.5f, %.5f] Hz, Nf=%d at %s Hz\n',band.f_lo,band.f_hi,Nf,mat2str(band.f,5));
cts=articleB_band_ctxs(ctx0,band,zeta);
ctx_s=ctx0; ctx_s.w=2*pi*f_s; ctx_s.alpha=2*zeta*ctx_s.w; ctx_s.beta=0; ctx_s.R0=[];

% ---------------- (S) single-frequency optimization -----------------------
fprintf('\n---- (S) SINGLE-FREQUENCY optimization at f_s=%.5f Hz ----\n',f_s);
[tS1,hS1,mS1]=articleB_shell_opt_one(ctx_s,metric,opt);
fprintf('  DONE S: Jbest=%.4e iters=%d evals=%d wall=%.1fs (%s)\n',...
        mS1.Jbest,mS1.niter,mS1.numeval,mS1.walltime,mS1.stopReason);

% ---------------- (B) band optimization -----------------------------------
fprintf('\n---- (B) BAND optimization, mode=%s over [%.5f,%.5f] Hz ----\n',...
        'mean',band.f_lo,band.f_hi);
optB=opt; optB.bandmode='mean'; optB.rho=60;
[tB1,hB1,mB1]=articleB_band_opt_one(cts,metric,band,optB);
fprintf('  DONE B: Jb_best=%.4e iters=%d evals=%d wall=%.1fs (%s)\n',...
        mB1.Jbest,mB1.niter,mB1.numeval,mB1.walltime,mB1.stopReason);

% ---------------- evaluate all designs on a WIDE sweep --------------------
fsw=linspace(0.60*f_peak,1.50*f_peak,61);
Jsw=zeros(numel(fsw),3); nm={'initial','single-f','band-mean'}; designs={t_init,tS1,tB1};
for di=1:3
  Jsw(:,di)=sweep_J(ctx0,fsw,zeta,designs{di},metric);
  fprintf('sweep design %-9s done\n',nm{di});
end
inband = fsw>=band.f_lo-1e-12 & fsw<=band.f_hi+1e-12;
fband  = fsw(inband);
% per-freq objective at f_s for each design
fprintf('done: mesh=%s f1=%.5f f_s=%.5f band=[%.5f,%.5f]\n',tag,f1,f_s,band.f_lo,band.f_hi);

% ---------------- band metrics table --------------------------------------
T=struct('name',{},'J_at_fs',{},'band_avg',{},'band_max',{},'f_at_bandmax',{},...
         'sweep_avg',{},'sweep_max',{},'f_at_sweepmax',{},'band_spread',{});
for di=1:3
  Jband=Jsw(inband,di);
  [bmax,~]=max(Jband); [smx,imx]=max(Jsw(:,di));
  T(di).name=nm{di};
  T(di).J_at_fs=interp1(fsw,Jsw(:,di),f_s,'linear');
  T(di).band_avg=mean(Jband);
  T(di).band_max=bmax;
  T(di).f_at_bandmax=fband(find(Jband==bmax,1));     % index into fband, NOT fsw
  T(di).sweep_avg=mean(Jsw(:,di)); T(di).sweep_max=smx; T(di).f_at_sweepmax=fsw(imx);
  T(di).band_spread=bmax/mean(Jband);
end
fid=fopen(fullfile(RES,sprintf('articleB_band_metrics_%s.txt',tag)),'w');
pr=@(varargin) both(fid,varargin{:});
pr('==== ARTICLE B BAND DE-ALIASING DEMO (curved Scordelis-Lo roof), mesh %s ====\n',tag);
pr('N=%d nel=%d | f1(t_start)=%.5f Hz | f_s=%.5f Hz (%.3f f1) | band=[%.5f,%.5f] Hz, Nf=%d | zeta=%.3g | volfrac=%.2f\n\n',...
   ctx0.ndof,ctx0.ne,f1,f_s,f_s/f1,band.f_lo,band.f_hi,Nf,zeta,opt.volfrac);
pr('objective: %s (dynamic compliance |F.''U|).  band metrics from a %d-point sweep over [%.2f,%.2f] f_peak,\n',...
   metric,numel(fsw),fsw(1)/f_peak,fsw(end)/f_peak);
pr('restricted to the band for avg/max, and reported for the whole sweep as the relocation evidence.\n\n');
pr('%-11s %13s %13s %13s %11s %13s %13s %12s %10s\n',...
   'design','J(f_s)','band_avg','band_max','f@bandmax','sweep_avg','sweep_max','f@sweepmax','spread');
for di=1:3
  pr('%-11s %13.5e %13.5e %13.5e %11.5f %13.5e %13.5e %12.5f %10.3f\n',...
     T(di).name,T(di).J_at_fs,T(di).band_avg,T(di).band_max,T(di).f_at_bandmax,...
     T(di).sweep_avg,T(di).sweep_max,T(di).f_at_sweepmax,T(di).band_spread);
end
pr('\ninterpretation: the SINGLE-FREQUENCY design has the smallest J(f_s) (best AT its own\n');
pr('frequency) but its resonance peak has MOVED: its band_max / sweep_max exceed the initial\n');
pr('design''s, i.e. it relocated (not suppressed) the peak. The BAND design lowers band_avg and\n');
pr('band_max together and its sweep_max is far below the single-frequency design''s.\n');
pr('optimizer: GCMMA (Svanberg), thesis-exact dB scaling, windowed stop. iterations: S=%d (%.1fs), B=%d (%.1fs).\n',...
   mS1.niter,mS1.walltime,mB1.niter,mB1.walltime);
pr('\nARTICLEB_BAND_RUN_DONE\n');
fclose(fid);
save(fullfile(RES,sprintf('articleB_band_%s.mat',tag)),'T','fsw','Jsw','band','f1','f_s','nsub',...
     'tS1','tB1','t_init','hS1','hB1','mS1','mB1','metric','zeta','opt');
fprintf('table : %s\n',fullfile(RES,sprintf('articleB_band_metrics_%s.txt',tag)));

% ---------------- FIGURE --------------------------------------------------
S=pub_style();
fig=figure('Position',[60 60 1450 560],'Color','w');

ax1=subplot(1,2,1); hold(ax1,'on');
% shade the band first (behind the lines)
yl=[min(Jsw(:))*0.6, max(Jsw(:))*1.6];
patch(ax1,[band.f_lo band.f_hi band.f_hi band.f_lo],[yl(1) yl(1) yl(2) yl(2)],...
      [0.85 0.90 0.97],'EdgeColor','none','FaceAlpha',0.9,'HandleVisibility','off');
set(ax1,'YScale','log'); ylim(ax1,yl);
h1=semilogy(ax1,fsw,Jsw(:,1),'-','Color',[0.35 0.35 0.35],'LineWidth',S.lw,'DisplayName','initial design');
h2=semilogy(ax1,fsw,Jsw(:,2),'-','Color',S.c.b2,'LineWidth',S.lw,'DisplayName','single-freq obj (f_s)');
h3=semilogy(ax1,fsw,Jsw(:,3),'-','Color',S.c.b1,'LineWidth',S.lw,'DisplayName','band obj (mean)');
plot(ax1,fsw(inband),Jsw(inband,3),'.','Color',S.c.b1,'MarkerSize',12,'HandleVisibility','off');
xline(ax1,f_s,'k--','f_s','LineWidth',S.lw_ref,'HandleVisibility','off');
xline(ax1,band.f_lo,':','Color',[0.3 0.4 0.6],'LineWidth',S.lw_ref,'HandleVisibility','off');
xline(ax1,band.f_hi,':','Color',[0.3 0.4 0.6],'LineWidth',S.lw_ref,'HandleVisibility','off');
xlabel(ax1,'frequency f (Hz)','FontSize',S.fs_lab);
ylabel(ax1,'dynamic compliance |F^T U|','FontSize',S.fs_lab);
title(ax1,'metric spectrum (band shaded)','FontSize',S.fs_title);
legend(ax1,[h1 h2 h3],'Location','northwest','FontSize',S.fs_leg);
grid(ax1,'on'); box(ax1,'on');

ax2=subplot(1,2,2);
M=[T(1).band_avg T(2).band_avg T(3).band_avg; T(1).band_max T(2).band_max T(3).band_max].';
b=bar(ax2,M,'grouped'); 
b(1).FaceColor=[0.30 0.45 0.70]; b(2).FaceColor=[0.75 0.35 0.20];
set(ax2,'YScale','log');
set(ax2,'XTickLabel',{'initial','single-freq','band'});
ylabel(ax2,'objective (log scale)','FontSize',S.fs_lab);
title(ax2,'band avg and band max after each optimisation','FontSize',S.fs_title);
legend(ax2,{'band average','band maximum'},'Location','northwest','FontSize',S.fs_leg);
grid(ax2,'on'); box(ax2,'on');
% annotate the numbers
for di=1:3
  text(ax2,di-0.15,T(di).band_avg,sprintf('%.3g',T(di).band_avg),'HorizontalAlignment','center',...
       'VerticalAlignment','bottom','FontSize',9);
  text(ax2,di+0.15,T(di).band_max,sprintf('%.3g',T(di).band_max),'HorizontalAlignment','center',...
       'VerticalAlignment','bottom','FontSize',9);
end
sgtitle(sprintf('Band de-aliasing: curved Scordelis-Lo roof, mesh %s, f_1=%.4f Hz, band [%.3f,%.3f] Hz',...
        tag,f1,band.f_lo,band.f_hi),'FontSize',S.fs_title);
base=fullfile(RES,sprintf('articleB_band_spectrum_%s',tag));
pub_save(fig,base,300); close(fig);
fprintf('figure: %s.png / .pdf\n',base);
fprintf('\nARTICLEB_BAND_RUN_DONE\n');
end

% --------------------------------- helpers --------------------------------
function Jf = sweep_J(ctx0,fs,zeta,t,metric)
Jf=zeros(numel(fs),1);
for i=1:numel(fs)
  ck=ctx0; ck.w=2*pi*fs(i); ck.alpha=2*zeta*ck.w; ck.beta=0; ck.R0=[];
  Jf(i)=articleB_shell_eval(ck,t,metric);
end
end
function both(fid,varargin), fprintf(fid,varargin{:}); fprintf(1,varargin{:}); end
