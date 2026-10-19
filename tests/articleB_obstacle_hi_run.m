function articleB_obstacle_hi_run(kind, nmax, res)
% ARTICLEB_OBSTACLE_HI_RUN  HIGH-RESOLUTION re-run of the Article B obstacle
% course, with the dB-scaling guard (src/objective_scale_guard.m) wired in.
%
%   articleB_obstacle_hi_run('spherecap', 30, [256 128])
%   articleB_obstacle_hi_run('cylpanel',  30, [256 128])
%   articleB_obstacle_hi_run('scordelis', 30, [128 64])
%
% Same benchmark/difficulty/degeneracy pipeline as tests/articleB_obstacle_run.m
% but (i) MESH is raised to the requested nsub, (ii) the design is developed with
% nmax >= 25-30 GCMMA outer iterations (design variables == element count, so a
% finer mesh buys a readable design field), (iii) the objective scaling goes
% through objective_scale_guard so the descent direction stays correct even when
% J falls below 1e-10 (the defect fixed for Article B).  Designs are rendered
% with spline interpolation onto a 4x finer display grid and surf(...,'FaceColor'
% 'interp'); every figure is saved through pub_save, which FORCES a white
% background and reads the PNG corners back as white-check evidence.
%
% Writes results/articleB_obstacle_hi_*.  ADD-only file.

if nargin<1 || isempty(kind), kind='spherecap'; end
if nargin<2 || isempty(nmax), nmax=30; end
if nargin<3 || isempty(res),  res=[256 128]; end
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
fid=fopen(fullfile(RES,sprintf('articleB_obstacle_hi_log_%s.txt',kind)),'w');
LG=@(varargin) both(fid,varargin{:});
LG('==== ARTICLE B OBSTACLE COURSE (HI-RES, guarded scaling) : %s ====\n',kind);
LG('started %s ; nmax=%d ; mesh=%dx%d\n\n', datestr(now), nmax, res(1), res(2));

zeta=0.05;
metrics={'dyncomp','inputpower','quadvel'};
opt=struct('tmin',0.25,'tmax',1.0,'volfrac',0.5,'nmax',nmax,'stop',true,'window',10);

TBL=struct('case',{},'crit',{},'critval',{},'nband',{},'band',{},'metric',{},...
  'Jbest',{},'iters',{},'wall',{},'ndof',{},'nel',{},'CV',{},'frac_mid',{},'conc1pct',{},'verdict',{});

try
  tS=tic;
  switch lower(kind)
    case 'spherecap'
      spec=struct('R',10,'alpha',deg2rad(30),'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,...
        'deg',[2 2],'reg',[1 1],'nsub',res,'nquad',[4 4],'bc',{2},...
        'loadtype','pole','loaddir',[0;0;1],'loadmag',1);
      ctx=obstacle_build_ctx('spherecap',spec);
      band=[80 130];
    case 'cylpanel'
      spec=struct('R',25,'L',50,'th0',deg2rad(40),'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,...
        'deg',[2 2],'reg',[1 1],'nsub',res,'nquad',[4 4],'bc',{[1 2 3 4]},...
        'loadtype','point','loaduv',[0.5 0.5],'loaddir',[0;1;0],'loadmag',1);
      ctx=obstacle_build_ctx('cylpanel',spec);
      band=[];
    case 'scordelis'
      ctx=articleB_build_ctx(res,4.32e8,0.3,2500,0.25,NaN,0);
      ctx.kind='scordelis'; spec=struct('nsub',res);
      band=[];
    otherwise
      error('unknown kind %s',kind);
  end
  LG('built %s: ndof=%d nel=%d (%.1fs)\n',kind,ctx.ndof,ctx.ne,toc(tS));

  % ---- spectrum + difficulty measurement -----------------------------------
  tS=tic;
  if strcmp(kind,'cylpanel')
    f0=eigband(ctx,4,[]); f1=f0(1); band=[f1 2.2*f1];
    [freqs_all,Vall]=eigband(ctx,60,1.5*f1); %#ok<NASGU>
  elseif strcmp(kind,'scordelis')
    [freqs_all,Vall]=eigband(ctx,60,0.12); f1=freqs_all(1); band=[0.9*f1 1.6*f1]; %#ok<NASGU>
  else
    [freqs_all,Vall]=eigband(ctx,60,105); f1=freqs_all(1); band=[80 130]; %#ok<NASGU>
  end
  LG('eigband %.1fs\n',toc(tS));
  nband=sum(freqs_all>=band(1)&freqs_all<=band(2));
  LG('f1=%.4f Hz ; band [%.3f %.3f] Hz ; %d modes in band ; modal density %.4f modes/Hz\n',...
     f1,band(1),band(2),nband,nband/(band(2)-band(1)));

  % ---- veering measurement (cylpanel only) ---------------------------------
  veer=struct('Ls',[],'Fs',[],'pair',[],'gap',[],'Lstar',[]);
  if strcmp(kind,'cylpanel')
    veer=measure_veering(spec,res,LG);
  end

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
  % strip the big operator arrays before saving (keep the designs + spectra)
  ctx_s=rmfield(ctx, intersect(fieldnames(ctx),{'Bke','Ske','Me','Rows','Cols','space','msh','geometry'}));
  save(fullfile(RES,sprintf('articleB_obstacle_hi_%s.mat',kind)),'R','band','freqs_all','nband','veer','ctx_s','opt','fr');

  tryfig(@() fig_designs(ctx,R,opt,kind),'designs',LG);
  tryfig(@() fig_spectrum(ctx,band,freqs_all,kind),'spectrum',LG);
  if strcmp(kind,'cylpanel'), tryfig(@() fig_veering(veer,kind),'veering',LG); end
  tryfig(@() fig_traj(R,kind),'traj',LG);

  write_table(fullfile(RES,'articleB_obstacle_hi_table.txt'),TBL);
  LG('\nbenchmark %s done (%.1fs total)\n',kind,toc(tS));
catch ME
  LG('\n*** ERROR in %s: %s\n',kind,ME.message);
  LG('%s\n',getReport(ME,'extended','hyperlinks','off'));
end
LG('\nHI_DONE\n');
fclose(fid);
fprintf('HI_DONE %s\n',kind);
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
  V0=opt.volfrac*sum(tmax_e);
  st=[]; Jbest=inf; tbest=t; p0=[]; nholdcnt=0;
  hist=struct('J',[],'Jbest',[],'f0',[],'vol',[],'chg',[],'innerit',[],'wall',[],'branch',{{}});
  gopt=opt; gopt.kmax=20; gopt.growth=1.30; gopt.raa0eps=1e-3; gopt.raaeps=1e-3;
  stopReason='nmax'; ts=tic;
  for it=1:opt.nmax
    t0=tic;
    [J,dJ]=obstacle_band_eval(ctx,freqs,t,metric,zeta);
    gv=sum(t)/V0-1; dgv=(1/V0)*ones(ne,1);
    % ---- GUARDED SCALING (replaces Jdb=100+10*log10(J)) -------------------
    [f0s,df0s,g]=objective_scale_guard(J,dJ,p0,struct());
    p0=g.p0; hist.branch{end+1}=g.branch;
    hist.f0(end+1)=f0s;
    gopt.eval=@(tq) deal(obj_guard(tq,ctx,freqs,metric,zeta,V0,p0), sum(tq)/V0-1);
    [tnew,st,info]=gcmma_solver(t,tmin_e,tmax_e,f0s,df0s,gv,dgv,st,it,gopt);
    hist.innerit(it)=info.innerit;
    chg=max(abs(tnew(:)-t));
    if J<Jbest, Jbest=J; tbest=t; end
    hist.J(end+1)=J; hist.Jbest(end+1)=Jbest; hist.vol(end+1)=sum(t)/sum(tmax_e);
    hist.chg(end+1)=chg; hist.wall(end+1)=toc(t0);
    LG('  [%-11s] it %3d: J=%.5e  f0(guarded)=%.4f  psi=%.4f  br=%s  Jbest=%.5e  vol=%.4f  chg=%.2e  inner=%d  (%.1fs)\n',...
       metric,it,J,f0s,g.psi,g.branch,Jbest,hist.vol(end),chg,info.innerit,toc(t0));
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
  meta=struct('metric',metric,'niter',numel(hist.J),'numeval',sum(hist.innerit),...
    'walltime',hist.walltime,'final_vol',sum(tbest)/sum(tmax_e),'stopReason',stopReason,'Jbest',Jbest);
end
function f0q=obj_guard(tq,ctx,freqs,metric,zeta,V0,p0)
  [Jq,~]=obstacle_band_eval(ctx,freqs,tq(:),metric,zeta);
  [f0q,~]=objective_scale_guard(Jq,0,p0,struct());
end

function veer=measure_veering(spec,~,LG)
  Ls=linspace(30,90,9); Fs=[]; Vs={};
  LG('  veering sweep over L=%s ...\n',mat2str([Ls(1) Ls(end)],3));
  for i=1:numel(Ls)
    s=spec; s.L=Ls(i); s.nsub=[32 32];
    c=obstacle_build_ctx('cylpanel',s);
    [f,V]=eigband(c,8); Fs(i,:)=f(:).'; Vs{i}=V;
  end
  gaps=diff(Fs,1,2); [gmin,idx]=min(gaps(:)); [ii,jj]=ind2sub(size(gaps),idx);
  mac=zeros(numel(Ls)-1,1);
  for i=1:numel(Ls)-1
    a=Vs{i}(:,jj); b=Vs{i+1}(:,jj);
    mac(i)=abs(a'*b)^2/((a'*a)*(b'*b));
  end
  veer=struct('Ls',Ls,'Fs',Fs,'pair',[jj jj+1],'gap',gmin,'Lstar',Ls(ii),'mac',mac);
  LG('  VEERING: branches %d & %d, MIN GAP = %.4f Hz at L=%.2f m\n',jj,jj+1,gmin,Ls(ii));
  LG('   MAC(branch %d) across sweep = %s\n',jj,mat2str(mac.',3));
end

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
  n1=ctx.nsub(1); n2=ctx.nsub(2); nM=numel(R);
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
  sgtitle(sprintf('%s optimal thickness fields  (mesh %dx%d, %d elems, %d dof, guarded scaling)',...
    kind,n1,n2,ctx.ne,ctx.ndof),'Interpreter','none');
  base=fullfile(RES,sprintf('articleB_obstacle_hi_design_%s',kind));
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
  title(sprintf('%s (hi-res): mode spectrum in band [%.1f, %.1f] Hz  --  %d modes',kind,band(1),band(2),n),'Interpreter','none');
  base=fullfile(RES,sprintf('articleB_obstacle_hi_spectrum_%s',kind));
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
  xlabel('panel length L (m)'); ylabel('natural frequency (Hz)');
  title(sprintf('%s (hi-res): frequency-branch VEERING -- branches %d & %d, min gap %.3f Hz at L=%.1f m',...
    kind,j,j+1,veer.gap,veer.Lstar),'Interpreter','none');
  legend({'all branches',sprintf('branch %d',j),sprintf('branch %d',j+1),'L^*'},'Location','best');
  base=fullfile(RES,sprintf('articleB_obstacle_hi_veering_%s',kind));
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
  title(sprintf('%s (hi-res): band-averaged objective trajectories',kind),'Interpreter','none');
  subplot(1,2,2); hold on; grid on;
  for k=1:nM, plot(1:numel(R(k).hist.vol),R(k).hist.vol,'-','Color',cols{k},'DisplayName',R(k).metric); end
  yline(0.5,'k--'); xlabel('outer iteration'); ylabel('volume fraction'); ylim([0.45 0.55]); legend('Location','best');
  title(sprintf('%s (hi-res): volume-fraction history',kind),'Interpreter','none');
  base=fullfile(RES,sprintf('articleB_obstacle_hi_traj_%s',kind));
  pub_save(fig,base,300); close(fig);
end

function [Xf,Yf,Tf]=smooth_field(tg,up)
  [n1,n2]=size(tg);
  xi=1:n1; et=1:n2; xf=linspace(1,n1,up*n1); ef=linspace(1,n2,up*n2);
  Tf=interp2(xi,et,tg.',xf,ef.','spline');
  [Xf,Yf]=meshgrid(xf,ef);
end

function write_table(path,T)
  old={};
  if exist(path,'file'), lines=readlines(path); old=lines(:); end
  names=unique({T.case});
  fid=fopen(path,'w');
  fprintf(fid,'# Article B OBSTACLE COURSE (HI-RES, guarded scaling): harmonic TopOpt of curved shells\n');
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
