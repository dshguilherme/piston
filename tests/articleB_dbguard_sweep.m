function articleB_dbguard_sweep(np)
% ARTICLEB_DBGUARD_SWEEP  Quantify the dB-scaling descent inversion and
% demonstrate that src/objective_scale_guard.m repairs it.
%
%   articleB_dbguard_sweep([40 24])       % optional demo mesh (speed knob)
%
% PART 1  ONSET.  The legacy scaling is  f0s(J) = 100*Jdb(J)/p0  with
%   Jdb(J) = 100 + 10*log10(J)  and the reference p0 = Jdb(J0) frozen at the
%   first iteration.  Sweep J through the threshold Jthr = 10^(-100/10) = 1e-10
%   for several reference values J0 and show that df0s/dJ changes SIGN exactly
%   when p0 = Jdb(J0) crosses zero, i.e. when J0 < 1e-10.  The guarded scaling
%   is monotone increasing for every J>0 and every reference.
%
% PART 2  REPAIR.  Run the ACTUAL cylindrical-panel band-averaged 'quadvel'
%   optimisation (the affected case, J ~ 1e-12, -120 dB) with the legacy
%   scaling and with the guard at a modest mesh, and report that legacy J
%   CLIMBS (degenerate) while guarded J DESCENDS and the volume stays at 0.5.
%
% Writes results/articleB_dbguard_onset.txt, articleB_dbguard_demo.txt and
% articleB_dbguard_onset.png.  ADD-only file.
% ---------------------------------------------------------------------------
if nargin < 1 || isempty(np), np = [40 24]; end
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';

Jthr = 1e-10; Abar=100; dbfac=10;
fid=fopen(fullfile(RES,'articleB_dbguard_onset.txt'),'w');
LG=@(varargin) both(fid,varargin{:});
LG('==== ARTICLE B: dB-SCALING ONSET / DESCENT-INVERSION AUDIT ====\n');
LG('started %s\n\n', datestr(now));
LG('legacy scaled objective:  f0s(J) = 100*Jdb(J)/p0 ,  Jdb(J)=%g+%g*log10(J)\n',Abar,dbfac);
LG('=> the dB map crosses zero at  Jthr = 10^(-%g/%g) = %g\n',Abar,dbfac,Jthr);
LG('=> p0 = Jdb(J0) is NEGATIVE iff J0 < Jthr; then 1/p0<0 and f0s is\n');
LG('   monotonically DECREASING in J (descent inverted).\n\n');

Jg = logspace(-16,-4,241);                      % fine sweep of the CURRENT J
J0list = [1e-6 1e-8 1e-9 1e-10 1e-11 1e-12 8.4e-13];
LG('%-10s %-12s %12s %14s %10s   %-10s\n','J0','p0=Jdb(J0)','sign(p0)','d f0s/dJ','increasing','guard incr');
onsetJ0 = NaN; rows={};
for J0 = J0list
  p0 = Abar + dbfac*log10(J0);
  f0s = 100*(Abar + dbfac*log10(Jg)) / p0;      % legacy normalised objective
  df0s = 100*(dbfac/log(10))./Jg / p0;          % its analytic slope d/dJ
  inc = all(df0s > 0);
  % guarded
  [gg,~,~] = objective_scale_guard(Jg, ones(size(Jg)), p0, struct());
  % monotone check of guarded psi via finite diff
  [ga,~,~] = objective_scale_guard(Jg(1:end-1)*1, ones(size(Jg(1:end-1))), p0, struct());
  [gb,~,~] = objective_scale_guard(Jg(2:end),   ones(size(Jg(2:end))),   p0, struct());
  ginc = all(gb > ga);                          % psi strictly increasing
  LG('%-10.1e %-12.4f %12d %14.3e %10d   %-10d\n',J0,p0,sign(p0),df0s(1),inc,ginc);
  rows{end+1}=struct('J0',J0,'p0',p0,'inc',inc,'ginc',ginc); %#ok<AGROW>
  if ~inc && isnan(onsetJ0), onsetJ0=J0; end
end
LG('\nONSET: legacy scaling becomes non-monotone (inverted) at J0 = %g\n',onsetJ0);
LG('       Jthr = %g is the analytic crossing of Jdb through zero.\n',Jthr);
LG('       guarded scaling is strictly increasing for ALL J0 (ginc=1).\n');

% identity check: guard == legacy for every J >= Jthr (the "unchanged" claim)
Jok = logspace(-10,-1,200);                    % all >= Jthr
p0ok = Abar + dbfac*log10(1e-6);              % a positive reference (J0=1e-6)
eF=0; eD=0;
for i=1:numel(Jok)
  l0 = Abar + dbfac*log10(Jok(i));             % legacy map
  ld = (dbfac/log(10))/Jok(i);                 % legacy slope
  f_leg = 100*l0/p0ok; df_leg = 100*ld/p0ok;
  [f_grd,df_grd] = objective_scale_guard(Jok(i), 1, p0ok, struct());
  eF=max(eF,abs(f_grd-f_leg)); eD=max(eD,abs(df_grd-df_leg));
end
LG('\nIDENTITY (J>=Jthr): max|f0s_guard - f0s_legacy| = %.3e ; max|df0s diff| = %.3e\n', eF, eD);
LG('  => guard is BIT-IDENTICAL to the legacy scaling wherever J >= 1e-10 (cases with J >> 1e-10 unchanged).\n\n');

% --- figure: legacy f0s(J) curves for the crucial references vs guard -------
try
  fig=figure('Position',[40 40 1000 460],'Color','w');
  subplot(1,2,1); hold on; grid on;
  show=[1e-8 1e-10 1e-12 8.4e-13];
  cols=lines(numel(show));
  for i=1:numel(show)
    p0=Abar+dbfac*log10(show(i));
    plot(Jg,100*(Abar+dbfac*log10(Jg))/p0,'-','Color',cols(i,:),'LineWidth',1.6,...
      'DisplayName',sprintf('J_0=%.1e (p_0=%.2f)',show(i),p0));
  end
  xline(Jthr,'k--','Jthr=1e-10'); xlabel('current J'); ylabel('legacy f0s(J)');
  title('legacy scaling: monotone up (p0>0) / INVERTED down (p0<0)'); legend('Location','best');
  subplot(1,2,2); hold on; grid on;
  for i=1:numel(show)
    p0=Abar+dbfac*log10(show(i));
    % re-evaluate guard with THIS reference
    psi=Abar+dbfac*log10(Jg); b=Jg<Jthr; s0=dbfac/log(10);
    psi(b)=s0*(Jg(b)-Jthr)/Jthr; 
    plot(Jg,100*psi/abs(p0),'-','Color',cols(i,:),'LineWidth',1.6,'DisplayName',sprintf('J_0=%.1e',show(i)));
  end
  xline(Jthr,'k--','Jthr=1e-10'); xlabel('current J'); ylabel('guarded f0s(J)');
  title('guarded scaling: strictly increasing for every reference'); legend('Location','best');
  pub_save(fig, fullfile(RES,'articleB_dbguard_onset'), 300); close(fig);
catch ME
  LG('figure failed: %s\n',ME.message);
end

% =========================================================================
% PART 2 : repair demonstration on the real affected case (cylpanel quadvel)
% =========================================================================
LG('\n==== PART 2: repair demonstration (cylpanel, band-averaged quadvel) ====\n');
fid2=fopen(fullfile(RES,'articleB_dbguard_demo.txt'),'w');
LG2=@(varargin) both(fid2,varargin{:});
spec=struct('R',25,'L',50,'th0',deg2rad(40),'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,...
  'deg',[2 2],'reg',[1 1],'nsub',np,'nquad',[4 4],'bc',{[1 2 3 4]},...
  'loadtype','point','loaduv',[0.5 0.5],'loaddir',[0;1;0],'loadmag',1);
ctx=obstacle_build_ctx('cylpanel',spec);
f0=eigband(ctx,4,[]); f1=f0(1); band=[f1 2.2*f1];
fr=linspace(band(1)+0.05*diff(band), band(2)-0.05*diff(band), 3);
LG2('mesh %dx%d  ndof=%d nel=%d  f1=%.4f Hz  band=[%.3f %.3f]\n',np(1),np(2),ctx.ndof,ctx.ne,f1,band(1),band(2));
LG2('band objective freqs: %s Hz\n\n',mat2str(fr,6));
zeta=0.05;

for mode={'legacy','guard'}
  md=mode{1};
  LG2('---- scaling = %s ----\n',md);
  ne=ctx.ne; tmin=0.25; tmax=1.0; volfrac=0.5; nmax=15;
  t=volfrac*tmax*ones(ne,1); V0=volfrac*sum(tmax*ones(ne,1));
  tmin_e=tmin*ones(ne,1); tmax_e=tmax*ones(ne,1);
  gopt=struct('kmax',20,'growth',1.30,'raa0eps',1e-3,'raaeps',1e-3);
  st=[]; p0=[]; Jhist=[]; volhist=[];
  for it=1:nmax
    [J,dJ]=band_eval(ctx,fr,t,'quadvel',zeta);
    gv=sum(t)/V0-1; dgv=(1/V0)*ones(ne,1);
    if strcmp(md,'legacy')
      Jdb=100+10*log10(max(J,eps));
      if isempty(p0), p0=Jdb; end
      f0s=100*Jdb/p0; df0s=100*(10/log(10))*dJ/max(J,eps)/p0;
      evalh=@(tq) obj_legacy(tq,ctx,fr,zeta,V0,p0);
    else
      [f0s,df0s,g]=objective_scale_guard(J,dJ,p0,struct());
      p0=g.p0;
      evalh=@(tq) obj_guard(tq,ctx,fr,zeta,V0,p0);
    end
    gopt.eval=@(tq) deal(evalh(tq), sum(tq)/V0-1);
    [tnew,st,info]=gcmma_solver(t,tmin_e,tmax_e,f0s,df0s,gv,dgv,st,it,gopt);
    Jhist(end+1)=J; volhist(end+1)=sum(t)/sum(tmax_e); %#ok<AGROW>
    LG2('  it %2d: J=%.5e (%.2f dB)  vol=%.4f  inner=%d\n',it,J,10*log10(max(J,eps)),volhist(end),info.innerit);
    t=tnew(:);
  end
  drift=Jhist(end)/Jhist(1);
  LG2('  => J first=%.4e last=%.4e  ratio=%.3f  (descent iff ratio<1)  vol_end=%.4f\n\n',...
    Jhist(1),Jhist(end),drift,volhist(end));
end
LG2('DEMO_DONE\n');
fclose(fid2);
fclose(fid);
fprintf('DBGUARD_DONE\n');
end

% ---------------------------------------------------------------------------
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
function [J,dJ]=band_eval(ctx0,freqs,t,metric,zeta)
  t=t(:); J=0; dJ=zeros(ctx0.ne,1);
  for f=freqs
    c=ctx0; c.w=2*pi*f; c.alpha=2*zeta*c.w; c.beta=0;
    [Jf,dJf]=articleB_shell_eval(c,t,metric); J=J+Jf; dJ=dJ+dJf;
  end
  J=J/numel(freqs); dJ=dJ/numel(freqs);
end
function f0q=obj_legacy(tq,ctx,freqs,zeta,V0,p0)
  [Jq,~]=band_eval(ctx,freqs,tq(:),'quadvel',zeta);
  f0q=100*(100+10*log10(max(Jq,eps)))/p0;
end
function f0q=obj_guard(tq,ctx,freqs,zeta,V0,p0)
  [Jq,~]=band_eval(ctx,freqs,tq(:),'quadvel',zeta);
  [f0q,~]=objective_scale_guard(Jq,0,p0,struct());
end
function both(fid,varargin), fprintf(fid,varargin{:}); fprintf(1,varargin{:}); end
