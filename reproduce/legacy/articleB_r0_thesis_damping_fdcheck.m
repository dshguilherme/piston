function articleB_r0_thesis_damping_fdcheck()
root='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course';addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');addpath(fullfile(root,'src'));addpath(fullfile(root,'tests'));
s=struct('R',25,'L',50,'th0',deg2rad(40),'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,'deg',[2 2],'reg',[1 1],'nsub',[12 8],'nquad',[4 4],'bc',{[1 2 3 4]},'loadtype','point','loaduv',[0.5 0.5],'loaddir',[0;1;0],'loadmag',1,'freq',6.53547,'zeta',0.05);
ctx=obstacle_build_ctx('cylpanel',s);ctx.alpha=0;ctx.beta=0.1/ctx.w;Rop=articleB_r0_build(ctx,1.21,343,128,s.freq);t=0.05*ones(ctx.ne,1);e=round(ctx.ne/2);h=2e-6;
metrics={'dyncomp','inputpower','quadvel'};modes={'','';'','';'rad','I'};
for k=1:3
 for j=1:size(modes,2)
  if strcmp(metrics{k},'quadvel'), mode=modes{k,j};else,mode='';end
  [J,dJ]=articleB_r0_eval(ctx,t,metrics{k},Rop,mode);tp=t;tm=t;tp(e)=tp(e)+h;tm(e)=tm(e)-h;
  Jp=articleB_r0_eval(ctx,tp,metrics{k},Rop,mode);Jm=articleB_r0_eval(ctx,tm,metrics{k},Rop,mode);fd=(Jp-Jm)/(2*h); rel=abs(fd-dJ(e))/max([abs(fd),abs(dJ(e)),realmin]);
  fprintf('%s|%s J=%.6e dJ=%.6e fd=%.6e rel=%.3e alpha=%.3g beta=%.7g\n',metrics{k},mode,J,dJ(e),fd,rel,ctx.alpha,ctx.beta);
  if rel>5e-4,error('Thesis damping FD check failed for %s/%s.',metrics{k},mode);end
 end
end
fprintf('THESIS_DAMPING_FD_OK\n');
end
