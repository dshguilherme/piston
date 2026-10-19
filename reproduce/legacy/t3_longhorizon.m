function t3_longhorizon(deadline, nmax, ckpt_every, outbase)
% T3-LH: continue the V2 (quadvel-rad, warm-from-AIP) panel optimization from
% its 600-iteration endpoint (tvb in results/t2_panel_quadvelrad_warm_aip.mat)
% until opt.deadline (posixtime of an unzoned local datetime, same convention
% as posixtime(datetime('now')) in articleB_r0_opt_one). Checkpointed.
if nargin<1 || isempty(deadline),   deadline = inf; end
if nargin<2 || isempty(nmax),       nmax = 100000; end
if nargin<3 || isempty(ckpt_every), ckpt_every = 250; end
if nargin<4 || isempty(outbase),    outbase = 't3_longhorizon'; end
coderoot=fileparts(fileparts(mfilename('fullpath')));
RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fullfile(coderoot,'src')); addpath(fullfile(coderoot,'tests'));
fprintf('driver: %s\n',which('articleB_r0_opt_one'));
rho_air=1.21; c_air=343; zeta_ref=0.05;
name='cylpanel'; f=6.53547; ns=[40,24];
spec=struct('R',25,'L',50,'th0',deg2rad(40),'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,...
  'deg',[2 2],'reg',[1 1],'nsub',[40 24],'nquad',[4 4],'bc',{[1 2 3 4]},...
  'loadtype','point','loaduv',[0.5 0.5],'loaddir',[0;1;0],'loadmag',1,'freq',f,'zeta',zeta_ref);
ctx=obstacle_build_ctx(name,spec); ctx.alpha=0; ctx.beta=0.1/ctx.w; ctx.zeta=zeta_ref;
Rop=articleB_r0_build(ctx,rho_air,c_air,128,f);
fprintf('\n=== T3-LH PANEL f=%.6g Hz mesh=%dx%d alpha=%.6g beta=%.8g ===\n',f,ns(1),ns(2),ctx.alpha,ctx.beta);
opt=struct('tmin',0.025,'tmax',0.10,'volfrac',0.5,'nmax',600,'stop',true,...
  'ftol',1e-4,'xchgtol',1e-4,'nhold',3,'window',10,'feastol',1e-6,...
  'Abar',100,'dbfac',10);
S=load(fullfile(RES,'t2_panel_quadvelrad_warm_aip.mat'),'tvb');
opt.tinit=S.tvb;
opt.nmax=nmax; opt.stop=true; opt.jfloor=0; opt.jfloor_factor=0;
opt.deadline=deadline; opt.ckpt_every=ckpt_every;
opt.ckpt_file=fullfile(RES,[outbase '_ckpt.mat']);
fprintf('nmax=%d ckpt_every=%d deadline=%.0f (now=%.0f, %.1f h left) ckpt=%s\n',...
  nmax,ckpt_every,deadline,posixtime(datetime('now')),(deadline-posixtime(datetime('now')))/3600,opt.ckpt_file);
[tbest,hist,meta]=articleB_r0_opt_one(ctx,'quadvel',Rop,'rad',opt);
M=measure_local(ctx,tbest,Rop);
meta.deadline=deadline; meta.source='t2_panel_quadvelrad_warm_aip.mat:tvb'; meta.f=f; meta.ns=ns;
n=numel(hist.J);
if n>0, fprintf('J(1)=%.10e\n',hist.J(1)); end
save(fullfile(RES,[outbase '.mat']),'tbest','hist','meta','M','-v7.3');
fid=fopen(fullfile(RES,[outbase '.txt']),'w');
fprintf(fid,'T3-LH quadvel-rad warm continuation from t2 V2 endpoint\n');
fprintf(fid,'stopReason: %s\n',meta.stopReason);
fprintf(fid,'iterations: %d\n',n);
fprintf(fid,'walltime_s: %.1f\n',hist.walltime);
fprintf(fid,'Jbest: %.10e\n',meta.Jbest);
if n>0, fprintf(fid,'J(1): %.10e\n',hist.J(1)); end
fprintf(fid,'Wrad: %.10e\n',M.Wrad);
fprintf(fid,'o/d: %.6f\n',M.o_over_d);
fprintf(fid,'Jbest every 500 its:\n');
for k=unique([500:500:n n])
  if k>=1, fprintf(fid,'  it %6d  Jbest %.10e\n',k,hist.Jbest(k)); end
end
if n>=2
  idx=max(1,n-1999):n;
  p=polyfit(idx,log10(hist.Jbest(idx)),1);
  fprintf(fid,'slope log10(Jbest) per 1000 its over last %d its: %.6e\n',numel(idx),1000*p(1));
else
  fprintf(fid,'slope log10(Jbest) per 1000 its: n/a (fewer than 2 its)\n');
end
fclose(fid);
fprintf('RESULT T3-LH: it=%d stop=%s Jbest=%.10e Wrad=%.8e o/d=%.4f\n',n,meta.stopReason,meta.Jbest,M.Wrad,M.o_over_d);
fprintf('T3_LONGHORIZON_DONE\n');
end

function M=measure_local(ctx,t,Rop)
[~,~,aux]=articleB_r0_eval(ctx,t,'quadvel',Rop,'rad'); u=aux.u; CU=Rop.C*u; z=Rop.area.*CU;
if isempty(Rop.Kfull),error('Requires exact Kfull'); end
Kz=Rop.Kfull*z; q=real(z'*Kz); d=sum(abs(z).^2); od=q-d;
M.Wrad=ctx.w^2*Rop.pref*q; M.WI=ctx.w^2*real(u'*u); M.Wones=ctx.w^2*Rop.pref*abs(sum(CU))^2;
M.rho_coh=abs(sum(CU))/max(sum(abs(CU)),realmin); M.diag=d; M.offdiag=od;
M.odfrac=od/max(abs(q),realmin); M.o_over_d=od/max(abs(d),realmin);
end
