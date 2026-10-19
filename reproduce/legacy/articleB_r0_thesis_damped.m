function articleB_r0_thesis_damped()
% R0 ladder with the thesis Eq. 6.118 stiffness-proportional damping.
% For each fixed-tone run: alpha=0, beta=0.1/omega (equivalent 5% modal
% damping at the excitation frequency); exact Kfull is used on these meshes.
root='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course';
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master')); addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation'); addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
RES=fullfile(root,'results'); rho_air=1.21; c_air=343; zeta_ref=0.05;
jobs={'dyncomp','dyncomp','';'inputpower','inputpower','';'quadvel-rad','quadvel','rad';'quadvel-I','quadvel','I';'quadvel-ones','quadvel','ones'};
cases={
 struct('name','spherecap','freq',100,'spec',struct('R',10,'alpha',deg2rad(30),'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,'deg',[2 2],'reg',[1 1],'nsub',[40 24],'nquad',[4 4],'bc',{2},'loadtype','pole','loaddir',[0;0;1],'loadmag',1,'freq',100,'zeta',zeta_ref),'tmin',0.025,'tmax',0.10),...
 struct('name','cylpanel','freq',6.53547,'spec',struct('R',25,'L',50,'th0',deg2rad(40),'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,'deg',[2 2],'reg',[1 1],'nsub',[40 24],'nquad',[4 4],'bc',{[1 2 3 4]},'loadtype','point','loaduv',[0.5 0.5],'loaddir',[0;1;0],'loadmag',1,'freq',6.53547,'zeta',zeta_ref),'tmin',0.025,'tmax',0.10)};
volfrac=0.5; nmax=600;
for ic=1:numel(cases)
 C=cases{ic}; name=C.name; f=C.freq; ns=C.spec.nsub;
 ctx=obstacle_build_ctx(name,C.spec); ctx.alpha=0; ctx.beta=0.1/ctx.w; ctx.zeta=zeta_ref; alpha=ctx.alpha; beta=ctx.beta;
 Rop=articleB_r0_build(ctx,rho_air,c_air,128,f);
 fprintf('\n=== %s f=%.6g Hz mesh=%dx%d alpha=%.6g beta=%.8g (thesis Eq.6.118) ===\n',name,f,ns(1),ns(2),alpha,beta);
 opt=struct('tmin',C.tmin,'tmax',C.tmax,'volfrac',volfrac,'nmax',nmax,'stop',true,'ftol',1e-4,'xchgtol',1e-4,'nhold',3,'window',10,'feastol',1e-6,'objscale','log','Abar',100,'dbfac',10);
 L=struct('label',{},'metric',{},'R0mode',{},'tbest',{},'hist',{},'meta',{},'meas',{});
 ftag=strrep(sprintf('%.8g',f),'.','p'); matfile=fullfile(RES,sprintf('articleB_r0_thesisdamp_%s_%dx%d_f%s_partial.mat',name,ns(1),ns(2),ftag));
 if exist(matfile,'file')
  S=load(matfile,'L'); L=S.L;
  for k=1:numel(L), if ~isfield(L(k).meta,'total_iter'), L(k).meta.total_iter=L(k).meta.niter; end, end
  fprintf('Resume: loaded %d saved designs from %s\n',numel(L),matfile);
 end
 for j=1:size(jobs,1)
  if j<=numel(L) && strcmp(L(j).meta.stopReason,'converged')
   fprintf('SKIP %s already converged\n',L(j).label);
   % Re-measure with current measure() to ensure all metric columns (incl. o_over_d) are present
   L(j).meas=measure(ctx,L(j).tbest,Rop);
   save(matfile,'L','ns','f','alpha','beta','volfrac','opt','-v7.3');
   continue;
  end
  if j<=numel(L)
   opt.tinit=L(j).tbest;
   if isfield(L(j).meta,'total_iter'), prev_total=L(j).meta.total_iter; else, prev_total=L(j).meta.niter; end
  else
   if isfield(opt,'tinit'), opt=rmfield(opt,'tinit'); end
   prev_total=0;
  end
  fprintf('-- %d/%d %s --\n',j,size(jobs,1),jobs{j,1}); ticjob=tic;
  [tb,h,m]=articleB_r0_opt_one(ctx,jobs{j,2},Rop,jobs{j,3},opt); M=measure(ctx,tb,Rop);
  m.total_iter=prev_total+m.niter;
  L(j)=struct('label',jobs{j,1},'metric',jobs{j,2},'R0mode',jobs{j,3},'tbest',tb,'hist',h,'meta',m,'meas',M);
  save(matfile,'L','ns','f','alpha','beta','volfrac','opt','-v7.3');
  fprintf('RESULT %s it=%d stop=%s eval=%d J=%.5e Wrad=%.5e WI=%.5e rho=%.5g od/tot=%.5g od/diag=%.5g time=%.1fs\n',...
    L(j).label,m.niter,m.stopReason,m.numeval,m.Jbest,M.Wrad,M.WI,M.rho_coh,M.odfrac,M.o_over_d,toc(ticjob));
 end
 save(matfile,'L','ns','f','alpha','beta','volfrac','opt','-v7.3');
 txt=fullfile(RES,sprintf('articleB_r0_thesisdamp_%s_%dx%d_f%s.txt',name,ns(1),ns(2),ftag)); fid=fopen(txt,'w');
 fprintf(fid,'R0 ladder using thesis Eq. 6.118 stiffness-proportional damping\ncase=%s mesh=%dx%d frequency=%.8g Hz alpha=0 beta=0.1/omega=%.12g s (zeta=0.05 at f) t=[%.4f,%.4f] volfrac<=%.2f nmax=%d\n',name,ns(1),ns(2),f,beta,C.tmin,C.tmax,volfrac,nmax);
 fprintf(fid,'design outer_iter stopReason evals Jbest Wrad WI Wones rho_coh CV vol offdiag_over_total o_over_d wall_s\n');
 for j=1:numel(L),M=L(j).meas;fprintf(fid,'%s %d %s %d %.8e %.8e %.8e %.8e %.6g %.5g %.5g %.7g %.7g %.2f\n',L(j).label,L(j).meta.total_iter,L(j).meta.stopReason,L(j).meta.numeval,L(j).meta.Jbest,M.Wrad,M.WI,M.Wones,M.rho_coh,std(L(j).tbest)/mean(L(j).tbest),L(j).meta.final_vol,M.odfrac,M.o_over_d,L(j).meta.walltime);end
 fclose(fid); fprintf('saved %s\n',txt);
end
fprintf('THESIS_DAMPED_R0_DONE\n');
end
function M=measure(ctx,t,Rop)
[~,~,aux]=articleB_r0_eval(ctx,t,'quadvel',Rop,'rad'); u=aux.u; CU=Rop.C*u; z=Rop.area.*CU;
if isempty(Rop.Kfull),error('This diagnostic ladder requires exact Kfull on the 40x24 mesh.');end
Kz=Rop.Kfull*z; q=real(z'*Kz); d=sum(abs(z).^2); od=q-d;
M.Wrad=ctx.w^2*Rop.pref*q; M.WI=ctx.w^2*real(u'*u); M.Wones=ctx.w^2*Rop.pref*abs(sum(CU))^2;
M.rho_coh=abs(sum(CU))/max(sum(abs(CU)),realmin); M.diag=d; M.offdiag=od; M.odfrac=od/max(abs(q),realmin); M.o_over_d=od/max(abs(d),realmin);
end
