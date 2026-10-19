function articleB_obstacle_run(cases, nmax)
% ARTICLEB_OBSTACLE_RUN  Article B OBSTACLE COURSE for harmonic topology
% optimization of curved shells.
%
%   articleB_obstacle_run([], 30)
%
% Builds a NEW benchmark suite of genuinely DIFFICULT curved-shell harmonic
% topology-optimization problems, MEASURES each difficulty criterion, optimizes
% the element thickness with the verified GCMMA (src/gcmma_solver.m) under a
% BAND-averaged harmonic objective, audits non-degeneracy, and writes readable
% WHITE-background figures (PNG + PDF) + a table to results/articleB_obstacle_*.
%
% Benchmarks (kind in obstacle_build_ctx.m):
%   'spherecap'  clamped spherical cap, periodic theta, crown point load
%                DIFFICULTY: modal density (a CLUSTER of modes inside the band)
%   'cylpanel'   open cylindrical panel, clamped all round
%                DIFFICULTY: mode VEERING (avoided crossings across a param sweep)
%   'scordelis'  Scordelis-Lo roof (classical CONTROL case)
%
% NEW files only; src/, results/ (existing), geopdes and immersed-iga untouched.

if nargin<1 || isempty(cases), cases = {'spherecap','cylpanel','scordelis'}; end
if nargin<2 || isempty(nmax), nmax = 30; end
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(RES,'dir'), mkdir(RES); end
fid=fopen(fullfile(RES,'articleB_obstacle_log.txt'),'w');
LG=@(varargin) both(fid,varargin{:});
LG('==== ARTICLE B OBSTACLE COURSE (harmonic TopOpt of curved shells) ====\n');
LG('started %s\n\n', datestr(now));

zeta=0.05;
metrics={'dyncomp','inputpower','quadvel'};
opt=struct('tmin',0.25,'tmax',1.0,'volfrac',0.5,'nmax',nmax,'stop',true,'window',10);

TBL=struct('case',{},'crit',{},'critval',{},'nband',{},'band',{},'metric',{},...
  'Jbest',{},'iters',{},'wall',{},'ndof',{},'nel',{},'CV',{},'frac_mid',{},'conc1pct',{},'verdict',{});

for ci=1:numel(cases)
  kind=cases{ci};
  LG('\n==================== BENCHMARK %d : %s ====================\n',ci,kind);
  try
  tS=tic;
  switch lower(kind)
    case 'spherecap'
      spec=struct('R',10,'alpha',deg2rad(30),'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,...
        'deg',[2 2],'reg',[1 1],'nsub',[128 64],'nquad',[4 4],'bc',{2},...
        'loadtype','pole','loaddir',[0;0;1],'loadmag',1);
      ctx=obstacle_build_ctx('spherecap',spec);
      band=[80 130];
    case 'cylpanel'
      spec=struct('R',25,'L',50,'th0',deg2rad(40),'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,...
        'deg',[2 2],'reg',[1 1],'nsub',[80 48],'nquad',[4 4],'bc',{[1 2 3 4]},...
        'loadtype','point','loaduv',[0.5 0.5],'loaddir',[0;1;0],'loadmag',1);
      ctx=obstacle_build_ctx('cylpanel',spec);
      band=[];
    case 'scordelis'
      ctx=articleB_build_ctx([48 24],4.32e8,0.3,2500,0.25,NaN,0);
      ctx.kind='scordelis'; spec=struct('nsub',[48 24]);
      band=[];
  end
  LG('built %s: ndof=%d nel=%d (%.1fs)\n',kind,ctx.ndof,ctx.ne,toc(tS));

  % ---- spectrum + difficulty measurement -----------------------------------
  tS=tic;
  if strcmp(kind,'cylpanel')
    f0=eigband(ctx,4,[]); f1=f0(1); band=[f1 2.2*f1];
    [freqs_all,Vall]=eigband(ctx,60,1.5*f1);
  elseif strcmp(kind,'scordelis')
    [freqs_all,Vall]=eigband(ctx,60,0.12); f1=freqs_all(1); band=[0.9*f1 1.6*f1];
  else
    [freqs_all,Vall]=eigband(ctx,60,105); f1=freqs_all(1); band=[80 130];
  end
  LG('eigband %.1fs\n',toc(tS));
  nband=sum(freqs_all>=band(1)&freqs_all<=band(2));
  LG('f1=%.4f Hz ; band [%.3f %.3f] Hz ; %d modes in band ; modal density %.4f modes/Hz\n',...
     f1,band(1),band(2),nband,nband/(band(2)-band(1)));
  LG('freqs<=band-top: %s\n',mat2str(freqs_all(freqs_all<=band(2)+1e-9).',5));

  % ---- veering measurement (cylpanel only) ---------------------------------
  veer=struct('Ls',[],'Fs',[],'pair',[],'gap',[],'Lstar',[]);
  if strcmp(kind,'cylpanel')
    veer=measure_veering(spec,LG);
  end

  % ---- band objective frequencies ------------------------------------------
  nf=3; fr=linspace(band(1)+0.05*(band(2)-band(1)), band(2)-0.05*(band(2)-band(1)), nf);
  LG('band objective averaged over %d frequencies: %s Hz\n',nf,mat2str(fr,6));

  % ---- optimize ------------------------------------------------------------
  for k=1:numel(metrics)
    LG('\n---- %s / %s ----\n',kind,metrics{k});
    [tb,h,m]=obstacle_opt_band(ctx,fr,metrics{k},zeta,opt,LG);
    dg=degeneracy_metrics(tb,ctx.nsub(1),ctx.nsub(2),opt.tmin,opt.tmax);
    verdict=nondeg_verdict(dg);
    LG('  Jbest=%.6e (%.3f dB) iters=%d wall=%.1fs vol=%.4f CV=%.2f frac@mid=%.2f conc1pct=%.3f -> %s\n',...
       m.Jbest,10*log10(max(m.Jbest,eps)),m.niter,m.walltime,m.final_vol,dg.CV,dg.frac_mid,dg.conc1pct,verdict);
    R(k)=struct('metric',metrics{k},'tbest',tb,'hist',h,'meta',m,'deg',dg,'verdict',verdict); %#ok<AGROW>
    TBL(end+1)=struct('case',kind,'crit',critname(kind),'critval',critval(kind,band,nband,f1,veer),...
      'nband',nband,'band',band,'metric',metrics{k},'Jbest',m.Jbest,'iters',m.niter,...
      'wall',m.walltime,'ndof',ctx.ndof,'nel',ctx.ne,'CV',dg.CV,'frac_mid',dg.frac_mid,...
      'conc1pct',dg.conc1pct,'verdict',verdict); %#ok<AGROW>
  end
  save(fullfile(RES,sprintf('articleB_obstacle_%s.mat',kind)),'R','band','freqs_all','nband','veer','ctx','opt','fr');

  % ---- figures (each non-fatal: a failed export must not lose the data) ----
  tryfig(@() fig_designs(ctx,R,opt,kind),'designs',LG);
  tryfig(@() fig_spectrum(ctx,band,freqs_all,kind),'spectrum',LG);
  if strcmp(kind,'cylpanel'), tryfig(@() fig_veering(veer,kind),'veering',LG); end
  tryfig(@() fig_geometry(ctx,band,kind),'geometry',LG);
  tryfig(@() fig_traj(R,kind),'traj',LG);

  LG('\nbenchmark %s done (%.1fs total)\n',kind,toc(tS));
  clear R
  catch ME
    LG('\n*** ERROR in %s: %s\n',kind,ME.message);
    LG('%s\n',getReport(ME,'extended','hyperlinks','off'));
  end
end

% ---- table ----------------------------------------------------------------
write_table(fullfile(RES,'articleB_obstacle_table.txt'),TBL);
LG('\nTABLE written. OBSTACLE_DONE\n');
fclose(fid);
fprintf('OBSTACLE_DONE\n');
end

% =====================================================================
function [f,V]=eigband(ctx,k,fc)
  Ks=sparse(ctx.Rows(:),ctx.Cols(:),reshape(ctx.Bke.*ctx.t0^3+ctx.Ske.*ctx.t0,[],1),ctx.ndof,ctx.ndof);
  M =sparse(ctx.Rows(:),ctx.Cols(:),reshape(ctx.Me.*ctx.t0,[],1),ctx.ndof,ctx.ndof);
  if nargin<3||isempty(fc)
    [VV,D]=eigs(Ks(ctx.free,ctx.free),M(ctx.free,ctx.free),k,'smallestabs');
  else
    [VV,D]=eigs(Ks(ctx.free,ctx.free),M(ctx.free,ctx.free),k,(2*pi*fc)^2);
  end
  d=real(diag(D)); [f,ix]=sort(sqrt(abs(d))); f=f/(2*pi); V=VV(:,ix);
end

function [J,dJ]=obstacle_band_eval(ctx0,freqs,t,metric,zeta)
  t=t(:); ne=ctx0.ne; J=0; dJ=zeros(ne,1);
  for f=freqs
    c=ctx0; c.w=2*pi*f; c.alpha=2*zeta*c.w; c.beta=0;
    [Jf,dJf]=articleB_shell_eval(c,t,metric);
    J=J+Jf; dJ=dJ+dJf;
  end
  J=J/numel(freqs); dJ=dJ/numel(freqs);
end

function [tbest,hist,meta]=obstacle_opt_band(ctx,freqs,metric,zeta,opt,LG)
  ne=ctx.ne;
  t=opt.volfrac*opt.tmax*ones(ne,1);
  tmin_e=opt.tmin*ones(ne,1); tmax_e=opt.tmax*ones(ne,1);
  V0=opt.volfrac*sum(tmax_e); Abar=100; dbfac=10;
  st=[]; Jbest=inf; tbest=t; p0=[]; nholdcnt=0;
  hist=struct('J',[],'Jbest',[],'f0',[],'vol',[],'chg',[],'innerit',[],'wall',[]);
  gopt=opt; gopt.kmax=20; gopt.growth=1.30; gopt.raa0eps=1e-3; gopt.raaeps=1e-3;
  stopReason='nmax'; ts=tic;
  for it=1:opt.nmax
    t0=tic;
    [J,dJ]=obstacle_band_eval(ctx,freqs,t,metric,zeta);
    gv=sum(t)/V0-1; dgv=(1/V0)*ones(ne,1);
    Jdb=Abar+dbfac*log10(max(J,eps));
    if isempty(p0), p0=Jdb; end
    f0s=100*Jdb/p0; df0s=100*(dbfac/log(10))*dJ/max(J,eps)/p0;
    hist.f0(end+1)=f0s;
    gopt.eval=@(tq) deal(obj_scaled(tq,ctx,freqs,metric,zeta,V0,p0), sum(tq)/V0-1);
    [tnew,st,info]=gcmma_solver(t,tmin_e,tmax_e,f0s,df0s,gv,dgv,st,it,gopt);
    hist.innerit(it)=info.innerit;
    chg=max(abs(tnew(:)-t));
    if J<Jbest, Jbest=J; tbest=t; end
    hist.J(end+1)=J; hist.Jbest(end+1)=Jbest; hist.vol(end+1)=sum(t)/sum(tmax_e);
    hist.chg(end+1)=chg; hist.wall(end+1)=toc(t0);
    LG('  [%-11s] it %3d: J=%.5e  f0(dB)=%+.4f  Jbest=%.5e  vol=%.4f  chg=%.2e  inner=%d  (%.1fs)\n',...
       metric,it,J,f0s,Jbest,hist.vol(end),chg,info.innerit,toc(t0));
    t=tnew(:);
    if opt.stop && numel(hist.J)>=2*opt.window+1
      W=opt.window; fo=hist.f0;
      fstat=abs(mean(fo(end-W+1:end))-mean(fo(end-2*W+1:end-W)))/max(abs(mean(fo(end-2*W+1:end-W))),eps);
      xstat=chg/max(opt.tmax-opt.tmin,eps);
      if xstat<1e-4 && fstat<1e-4 && gv<=1e-6, nholdcnt=nholdcnt+1; else, nholdcnt=0; end
      if nholdcnt>=3, stopReason='converged'; LG('  [%s] STOP it %d converged\n',metric,it); break; end
    end
  end
  hist.stopReason=stopReason; hist.walltime=toc(ts); hist.numeval=sum(hist.innerit);
  hist.tailspread10=max(hist.J(max(1,end-9):end))/max(min(hist.J(max(1,end-9):end)),eps);
  meta=struct('metric',metric,'niter',numel(hist.J),'numeval',sum(hist.innerit),...
    'walltime',hist.walltime,'final_vol',sum(tbest)/sum(tmax_e),'stopReason',stopReason,'Jbest',Jbest);
end
function f0q=obj_scaled(tq,ctx,freqs,metric,zeta,V0,p0)
  [Jq,~]=obstacle_band_eval(ctx,freqs,tq(:),metric,zeta);
  f0q=100*(100+10*log10(max(Jq,eps)))/p0;
end

% ---------------- veering measurement -------------------------------------
function veer=measure_veering(spec,LG)
  Ls=linspace(30,90,9); Fs=[]; Vs={};
  LG('  veering sweep over L=%s ...\n',mat2str([Ls(1) Ls(end)],3));
  for i=1:numel(Ls)
    s=spec; s.L=Ls(i); s.nsub=[32 32];
    c=obstacle_build_ctx('cylpanel',s);
    [f,V]=eigband(c,8); Fs(i,:)=f(:).'; Vs{i}=V;
  end
  gaps=diff(Fs,1,2); [gmin,idx]=min(gaps(:)); [ii,jj]=ind2sub(size(gaps),idx);
  % MAC correlation of branch jj across the sweep to show the shape swap
  mac=zeros(numel(Ls)-1,1);
  for i=1:numel(Ls)-1
    a=Vs{i}(:,jj); b=Vs{i+1}(:,jj);
    mac(i)=abs(a'*b)^2/((a'*a)*(b'*b));
  end
  veer=struct('Ls',Ls,'Fs',Fs,'pair',[jj jj+1],'gap',gmin,'Lstar',Ls(ii),'mac',mac);
  LG('  VEERING: branches %d & %d, MIN GAP = %.4f Hz at L=%.2f m\n',jj,jj+1,gmin,Ls(ii));
  LG('   branch%d = %s\n',jj,mat2str(Fs(:,jj).',4));
  LG('   branch%d = %s\n',jj+1,mat2str(Fs(:,jj+1).',4));
  LG('   MAC(branch %d) across sweep = %s\n',jj,mat2str(mac.',3));
end

% ---------------- degeneracy + verdict ------------------------------------
function dg=degeneracy_metrics(t,nel1,nel2,tmin,tmax)
  t=t(:); n=numel(t); m=mean(t); s=std(t);
  dev=abs(t-m); [~,ord]=sort(dev,'descend'); n1=max(1,round(0.01*n));
  conc1pct=sum(dev(ord(1:n1)))/max(sum(dev),eps);
  tg=reshape(t,nel1,nel2); hot=tg>(m+s);
  [I,J]=find(hot);
  if numel(I)>=2, hot_span=sqrt(var(I/(nel1))+var(J/(nel2))); else, hot_span=0; end
  dg=struct('mean',m,'std',s,'CV',s/max(m,eps),'range',[min(t) max(t)],...
    'frac_tmin',mean(t<=tmin+1e-6),'frac_tmax',mean(t>=tmax-1e-6),...
    'frac_mid',mean(t>tmin+1e-6 & t<tmax-1e-6),'conc1pct',conc1pct,'hot_span',hot_span);
end
function v=nondeg_verdict(dg)
  % Per Article B's established metric set, the DEGENERATE signatures are:
  %   (i)  uniform collapse  -> thickness CV ~ 0 (design loses its shape);
  %   (ii) localized spot     -> conc1pct -> 1 (deviation in the hottest 1%).
  % NOTE: frac@tmin/frac@tmax being 0 is NORMAL for this graded design family
  % (the Scordelis-Lo baseline design is fully interior, frac@mid=1, yet
  %  non-degenerate with CV=0.68, conc1pct=0.016), so bounds are NOT used.
  if dg.CV>0.15 && dg.conc1pct<0.30
    v='NON-DEGENERATE';
  elseif dg.CV<0.05 || dg.conc1pct>0.60
    v='DEGENERATE';
  else
    v='MARGINAL';
  end
end
function s=critname(kind)
  switch kind, case 'spherecap', s='modal density'; case 'cylpanel', s='mode veering';
    otherwise, s='classical control'; end
end
function s=critval(kind,band,nband,f1,veer)
  switch kind
    case 'spherecap', s=sprintf('%d modes in [%.0f,%.0f] Hz (%.3f/Hz)',nband,band(1),band(2),nband/(band(2)-band(1)));
    case 'cylpanel', s=sprintf('min branch gap %.3f Hz at L=%.1f (f1=%.2f)',veer.gap,veer.Lstar,f1);
    otherwise, s=sprintf('f1=%.4f Hz',f1);
  end
end

% ---------------- figures --------------------------------------------------
function fig_designs(ctx,R,opt,kind)
  RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
  S=pub_style(); n1=ctx.nsub(1); n2=ctx.nsub(2); nM=numel(R);
  fig=figure('Position',[30 30 1750 620],'Color','w');
  for k=1:nM
    tg=reshape(R(k).tbest,n1,n2); [Xf,Yf,Tf]=smooth_field(tg,4);
    subplot(1,nM,k);
    surf(Xf,Yf,Tf,'EdgeColor','none','FaceColor','interp'); view(2); axis tight;
    cb=colorbar; cb.Label.String='thickness t (m)'; caxis([opt.tmin opt.tmax]); colormap(parula);
    xlabel('element i  (\xi)'); ylabel('element j  (\eta)');
    title(sprintf('%s:  J_{best}=%.4e (%.2f dB)\nCV=%.2f  conc_{1%%}=%.3f  %s',...
      R(k).metric,R(k).meta.Jbest,10*log10(max(R(k).meta.Jbest,eps)),R(k).deg.CV,R(k).deg.conc1pct,R(k).verdict),...
      'Interpreter','none');
  end
  sgtitle(sprintf('%s optimal thickness fields  (mesh %dx%d, %d elems, %d dof, band-averaged objectives)',...
    kind,n1,n2,ctx.ne,ctx.ndof),'Interpreter','none');
  base=fullfile(RES,sprintf('articleB_obstacle_design_%s',kind));
  pub_save(fig,base,300); close(fig);
end

function fig_spectrum(ctx,band,freqs_all,kind)
  RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
  S=pub_style();
  fin=freqs_all(freqs_all<=band(2)); n=numel(fin);
  fig=figure('Position',[40 40 900 500],'Color','w');
  stem(1:n,fin,'filled','Color',S.c.b1,'MarkerFaceColor',S.c.b1); hold on; grid on;
  xl=xlim; patch([xl(1) xl(2) xl(2) xl(1)],[band(1) band(1) band(2) band(2)],[0.85 0.9 1],...
    'FaceAlpha',0.35,'EdgeColor','none');
  yline(band(1),'k--','band low'); yline(band(2),'k--','band high');
  xlabel('mode index in band'); ylabel('natural frequency (Hz)');
  title(sprintf('%s: mode spectrum in the design band [%.1f, %.1f] Hz  --  %d modes',kind,band(1),band(2),n),'Interpreter','none');
  base=fullfile(RES,sprintf('articleB_obstacle_spectrum_%s',kind));
  pub_save(fig,base,300); close(fig);
end

function fig_veering(veer,kind)
  RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
  S=pub_style();
  fig=figure('Position',[40 40 900 560],'Color','w');
  hold on; grid on;
  plot(veer.Ls,veer.Fs,'-','Color',[0.6 0.6 0.6],'LineWidth',1.0);
  j=veer.pair(1);
  plot(veer.Ls,veer.Fs(:,j),'-o','Color',S.c.b1,'LineWidth',2,'MarkerSize',4);
  plot(veer.Ls,veer.Fs(:,j+1),'-s','Color',S.c.b2,'LineWidth',2,'MarkerSize',4);
  xline(veer.Lstar,'k--','LineWidth',1.2);
  plot(veer.Lstar,veer.Fs(veer.Ls==veer.Lstar,j),'kp','MarkerSize',14,'MarkerFaceColor','y');
  xlabel('panel length L (m)'); ylabel('natural frequency (Hz)');
  title(sprintf('%s: frequency-branch VEERING -- branches %d & %d, minimum gap %.3f Hz at L=%.1f m',...
    kind,j,j+1,veer.gap,veer.Lstar),'Interpreter','none');
  legend({'all branches',sprintf('branch %d',j),sprintf('branch %d',j+1),'L^*'},'Location','best');
  base=fullfile(RES,sprintf('articleB_obstacle_veering_%s',kind));
  pub_save(fig,base,300); close(fig);
end

function fig_geometry(ctx,band,kind)
  RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
  S=pub_style();
  n1=ctx.nsub(1); n2=ctx.nsub(2);
  [U,V]=meshgrid(linspace(0,1,60),linspace(0,1,60));
  P=[U(:).'; V(:).'];
  X=ctx.geometry.map(P);
  Xg=reshape(X(1,:),60,60); Yg=reshape(X(2,:),60,60); Zg=reshape(X(3,:),60,60);
  fig=figure('Position',[40 40 1000 520],'Color','w');
  surf(Xg,Yg,Zg,zeros(60),'EdgeColor','none','FaceAlpha',0.85); hold on; axis equal;
  colormap(parula); view(-35,25); grid on;
  % excitation arrow
  q0=[0.5 0.5]; if strcmp(kind,'spherecap'), q0=[0 0.5]; end
  xq=ctx.geometry.map(q0.'); 
  switch kind
    case 'spherecap', dvec=[0;0;-1]; txt='crown point load';
    case 'cylpanel',  dvec=[0;1;0]; txt='point load';
    otherwise,        dvec=[0;0;-1]; txt='body load';
  end
  quiver3(xq(1),xq(2),xq(3),dvec(1),dvec(2),dvec(3),0.5*max(max(Xg(:))-min(Xg(:)),max(Zg(:))-min(Zg(:))),...
    'r','LineWidth',2.2,'MaxHeadSize',0.8);
  text(xq(1),xq(2),xq(3),['  ' txt],'Color','r','FontWeight','bold');
  xlabel('x'); ylabel('y'); zlabel('z');
  title(sprintf('%s benchmark geometry; design band [%.1f, %.1f] Hz',kind,band(1),band(2)),'Interpreter','none');
  base=fullfile(RES,sprintf('articleB_obstacle_geometry_%s',kind));
  pub_save(fig,base,300); close(fig);
end

function fig_traj(R,kind)
  RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
  S=pub_style(); nM=numel(R);
  fig=figure('Position',[40 40 1600 460],'Color','w');
  cols={S.c.b1,S.c.b2,S.c.b3};
  subplot(1,2,1); hold on; grid on;
  for k=1:nM, semilogy(1:numel(R(k).hist.J),R(k).hist.J,'-','Color',cols{k},'DisplayName',R(k).metric); end
  xlabel('outer iteration'); ylabel('objective J (raw, log)'); legend('Location','best');
  title(sprintf('%s: band-averaged objective trajectories',kind),'Interpreter','none');
  subplot(1,2,2); hold on; grid on;
  for k=1:nM, plot(1:numel(R(k).hist.vol),R(k).hist.vol,'-','Color',cols{k},'DisplayName',R(k).metric); end
  yline(0.5,'k--'); xlabel('outer iteration'); ylabel('volume fraction'); ylim([0.45 0.55]); legend('Location','best');
  title(sprintf('%s: volume-fraction history',kind),'Interpreter','none');
  base=fullfile(RES,sprintf('articleB_obstacle_traj_%s',kind));
  pub_save(fig,base,300); close(fig);
end

function [Xf,Yf,Tf]=smooth_field(tg,up)
  [n1,n2]=size(tg);
  xi=1:n1; et=1:n2; xf=linspace(1,n1,up*n1); ef=linspace(1,n2,up*n2);
  Tf=interp2(xi,et,tg.',xf,ef.','spline');
  [Xf,Yf]=meshgrid(xf,ef);
end

function write_table(path,T)
  % Merge: keep rows already in the file whose benchmark is NOT re-run here,
  % so re-running a single benchmark does not clobber the others.
  old={};
  if exist(path,'file')
    lines=readlines(path); old=lines(:);
  end
  names=unique({T.case});
  fid=fopen(path,'w');
  fprintf(fid,'# Article B OBSTACLE COURSE: harmonic TopOpt of curved shells\n');
  fprintf(fid,'# difficulty criterion with MEASURED value; band-averaged objective; non-degeneracy verdict\n');
  hdr=sprintf('%-11s %-16s %-42s %6s %10s %10s %11s %6s %9s %7s %8s %8s %9s %-22s',...
    'benchmark','criterion','measured value','#band','band(Hz)','metric','J_best','iters','wall(s)','ndof','nel','CV','frac@mid','verdict');
  fprintf(fid,'%s\n',hdr);
  for i=1:numel(old)
    ln=old{i};
    if isempty(ln) || startsWith(strtrim(ln),'#') || startsWith(ln,hdr(1:min(11,end))), continue; end
    keep=true; for n=1:numel(names), if startsWith(ln,names{n}), keep=false; end, end
    if keep && ~isempty(strtrim(ln)), fprintf(fid,'%s\n',ln); end
  end
  for i=1:numel(T)
    fprintf(fid,'%-11s %-16s %-42s %6d %10s %10s %11.4e %6d %9.1f %7d %8d %8.3f %8.3f %-22s\n',...
      T(i).case,T(i).crit,T(i).critval,T(i).nband,sprintf('[%.3g %.3g]',T(i).band(1),T(i).band(2)),...
      T(i).metric,T(i).Jbest,T(i).iters,T(i).wall,T(i).ndof,T(i).nel,T(i).CV,T(i).frac_mid,T(i).verdict);
  end
  fclose(fid);
  fprintf('TABLE %s\n',path);
end

function tryfig(fn,tag,LG)
  for attempt=1:2
    try, fn(); return; catch ME, LG('  [%s] figure attempt %d failed: %s\n',tag,attempt,ME.message); end
  end
  LG('  [%s] figure FAILED after retries\n',tag);
end
function both(fid,varargin), fprintf(fid,varargin{:}); fprintf(1,varargin{:}); end
