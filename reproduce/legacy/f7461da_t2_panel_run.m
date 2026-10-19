function t2_panel_run()
% T2 PANEL: characterise dyncomp non-convergence under thesis damping.
% Runs panel dyncomp (cold-start) with enhanced checkpointing + diagnostics,
% then warm-start probe: dyncomp from AIP seed, and quadvel-rad from AIP seed.
%
% MATLAB SEQUENTIAL constraint: all runs one at a time, checkpointed between.
root='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course';
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
RES=fullfile(root,'results');
rho_air=1.21; c_air=343; zeta_ref=0.05;
name='cylpanel'; f=6.53547; ns=[40,24];
spec=struct('R',25,'L',50,'th0',deg2rad(40),'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,...
  'deg',[2 2],'reg',[1 1],'nsub',[40 24],'nquad',[4 4],'bc',{[1 2 3 4]},...
  'loadtype','point','loaduv',[0.5 0.5],'loaddir',[0;1;0],'loadmag',1,'freq',f,'zeta',zeta_ref);
ctx=obstacle_build_ctx(name,spec); ctx.alpha=0; ctx.beta=0.1/ctx.w; ctx.zeta=zeta_ref;
Rop=articleB_r0_build(ctx,rho_air,c_air,128,f);
fprintf('\n=== T2 PANEL f=%.6g Hz mesh=%dx%d alpha=%.6g beta=%.8g ===\n',f,ns(1),ns(2),ctx.alpha,ctx.beta);

% --- RUN 1: dyncomp cold-start (no tinit) ---
opt=struct('tmin',0.025,'tmax',0.10,'volfrac',0.5,'nmax',600,'stop',true,...
  'ftol',1e-4,'xchgtol',1e-4,'nhold',3,'window',10,'feastol',1e-6,...
  'Abar',100,'dbfac',10);
if isfield(opt,'tinit'), opt=rmfield(opt,'tinit'); end
fprintf('\n--- RUN 1: dyncomp (cold-start, thesis damping) ---\n');
[tb,h,m]=articleB_r0_opt_one(ctx,'dyncomp',Rop,[],opt);
M=measure_local(ctx,tb,Rop);
fprintf('RESULT dyncomp-cold: it=%d stop=%s eval=%d J=%.8e Wrad=%.8e o/tot=%.5f o/d=%.4f time=%.1fs\n',...
  m.niter,m.stopReason,m.numeval,m.Jbest,M.Wrad,M.odfrac,M.o_over_d,toc(tic));

% Save this run's history + design
matfile_full=fullfile(RES,'t2_panel_dyncomp_cold.mat');
save(matfile_full,'tb','h','m','M','ns','f','ctx','-v7.3');
fprintf('Saved: %s\n',matfile_full);

% --- Extract AIP converged design for warm-start ----
% AIP = inputpower. Re-run from cold if needed (should converge by ~it 306).
fprintf('\n--- RUN 2: inputpower (AIP, cold-start, for warm-start seed) ---\n');
[ap,ha,ma]=articleB_r0_opt_one(ctx,'inputpower',Rop,[],opt);
Ma=measure_local(ctx,ap,Rop);
fprintf('RESULT inputpower: it=%d stop=%s eval=%d J=%.8e Wrad=%.8e time=%.1fs\n',...
  ma.niter,ma.stopReason,ma.numeval,ma.Jbest,Ma.Wrad,toc(tic));
matfile_aip=fullfile(RES,'t2_panel_inputpower_seed.mat');
save(matfile_aip,'ap','ha','ma','Ma','ns','f','-v7.3');
fprintf('Saved: %s\n',matfile_aip);

% --- RUN 3: dyncomp warm-started from AIP seed ---
fprintf('\n--- RUN 3: dyncomp (warm-start from AIP seed) ---\n');
opt_ws=opt; opt_ws.tinit=ap;
[twb,hwb,mwb]=articleB_r0_opt_one(ctx,'dyncomp',Rop,[],opt_ws);
Mwb=measure_local(ctx,twb,Rop);
fprintf('RESULT dyncomp-warm(AIP): it=%d stop=%s eval=%d J=%.8e Wrad=%.8e o/tot=%.5f o/d=%.4f time=%.1fs\n',...
  mwb.niter,mwb.stopReason,mwb.numeval,mwb.Jbest,Mwb.Wrad,Mwb.odfrac,Mwb.o_over_d,toc(tic));
matfile_ws=fullfile(RES,'t2_panel_dyncomp_warm_aip.mat');
save(matfile_ws,'twb','hwb','mwb','Mwb','ns','f','-v7.3');
fprintf('Saved: %s\n',matfile_ws);

% --- RUN 4: quadvel-rad warm-started from AIP seed (V2 phase-aware) ---
fprintf('\n--- RUN 4: quadvel-rad (warm-start from AIP seed, V2) ---\n');
opt_v2=opt; opt_v2.tinit=ap;
[tvb,hvb,mvb]=articleB_r0_opt_one(ctx,'quadvel',Rop,'rad',opt_v2);
Mvb=measure_local(ctx,tvb,Rop);
fprintf('RESULT quadvel-rad-warm(AIP): it=%d stop=%s eval=%d J=%.8e Wrad=%.8e o/tot=%.5f o/d=%.4f time=%.1fs\n',...
  mvb.niter,mvb.stopReason,mvb.numeval,mvb.Jbest,Mvb.Wrad,Mvb.odfrac,Mvb.o_over_d,toc(tic));
matfile_v2=fullfile(RES,'t2_panel_quadvelrad_warm_aip.mat');
save(matfile_v2,'tvb','hvb','mvb','Mvb','ns','f','-v7.3');
fprintf('Saved: %s\n',matfile_v2);

fprintf('\nT2_PANEL_DONE\n');
end

function M=measure_local(ctx,t,Rop)
[~,~,aux]=articleB_r0_eval(ctx,t,'quadvel',Rop,'rad'); u=aux.u; CU=Rop.C*u; z=Rop.area.*CU;
if isempty(Rop.Kfull),error('Requires exact Kfull'); end
Kz=Rop.Kfull*z; q=real(z'*Kz); d=sum(abs(z).^2); od=q-d;
M.Wrad=ctx.w^2*Rop.pref*q; M.WI=ctx.w^2*real(u'*u); M.Wones=ctx.w^2*Rop.pref*abs(sum(CU))^2;
M.rho_coh=abs(sum(CU))/max(sum(abs(CU)),realmin); M.diag=d; M.offdiag=od;
M.odfrac=od/max(abs(q),realmin); M.o_over_d=od/max(abs(d),realmin);
end
