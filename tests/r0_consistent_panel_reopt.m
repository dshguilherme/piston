function r0_consistent_panel_reopt(nmax,use_old_seed)
% Re-optimize the legacy panel geometry under the spline-consistent sinc model.
% All outputs are for the existing displacement-fixed boundary condition; this
% does not claim a full Kirchhoff-Love rotational clamp or a BEM exterior model.
if nargin<1,nmax=250;end
if nargin<2,use_old_seed=false;end
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
out=fullfile(root,'results','r0_revision'); if ~exist(out,'dir'),mkdir(out);end
spec=struct('R',25,'L',50,'th0',deg2rad(40),'E',2e11,'NU',0.3,'RHO',7850, ...
 't0',0.1,'deg',[2 2],'reg',[1 1],'nsub',[40 24],'nquad',[4 4], ...
 'bc',[1 2 3 4],'loadtype','point','loaduv',[0.5 0.5], ...
 'loaddir',[0;1;0],'loadmag',1,'freq',6.53547,'zeta',0.05);
ctx=obstacle_build_ctx('cylpanel',spec); ctx.alpha=0; ctx.beta=0.1/ctx.w;
Rop=r0_build_consistent(ctx,1.21,343,[4 4],32,64);
fprintf('resolved evaluator: %s\n',which('articleB_r0_eval'));
fprintf('resolved builder: %s\n',which('obstacle_build_ctx'));
fprintf('design variables %d, FE dofs %d, acoustic angular factor %d x %d\n',ctx.ne,ctx.ndof,size(Rop.factor,1),size(Rop.factor,2));
old=load(fullfile(root,'results','t3_longhorizon.mat'),'tbest'); told=old.tbest(:);
opt=struct('tmin',0.025,'tmax',0.1,'volfrac',0.5,'nmax',nmax,'stop',true, ...
 'vweights',Rop.element_area,'window',12,'nhold',4,'ftol',2e-5,'xchgtol',2e-5);
if use_old_seed, opt.tinit=told; end
[tbest,hist,meta]=articleB_r0_opt_one(ctx,'quadvel',Rop,'consistent',opt);
function [W,U]=evaldesign(t)
    KE=ctx.Bke.*(t.^3)+ctx.Ske.*t; ME=ctx.Me.*t;
    KD=(1+1i*ctx.w*ctx.beta)*sparse(ctx.Rows(:),ctx.Cols(:),KE(:),ctx.ndof,ctx.ndof)- ...
       ctx.w^2*sparse(ctx.Rows(:),ctx.Cols(:),ME(:),ctx.ndof,ctx.ndof);
    U=zeros(ctx.ndof,1); U(ctx.free)=KD(ctx.free,ctx.free)\ctx.F(ctx.free);
    W=ctx.w^2*norm(Rop.factor*U)^2;
end
[Wnew,Unew]=evaldesign(tbest); [Wold,Uold]=evaldesign(told);
fprintf('Consistent-model reevaluation: old design %.12e W, new design %.12e W\n',Wold,Wnew);
fprintf('weighted volume old %.8f, new %.8f; stop=%s, iter=%d\n', ...
 sum(Rop.element_area.*told)/sum(Rop.element_area)/opt.tmax, ...
 sum(Rop.element_area.*tbest)/sum(Rop.element_area)/opt.tmax,meta.stopReason,meta.niter);
run_name='consistent_panel_reoptimization'; if use_old_seed, run_name='consistent_panel_reoptimization_warmstart'; end
save(fullfile(out,[run_name '.mat']),'tbest','told','hist','meta','Wnew','Wold','Unew','Uold','spec','opt','Rop','-v7.3');
fid=fopen(fullfile(out,[run_name '.txt']),'w');
fprintf(fid,'Consistent-spline prescribed-source sinc reoptimization, cylindrical panel.\n');
fprintf(fid,'Boundary condition: existing displacement-fixed edges; not a verified KL rotation clamp.\n');
fprintf(fid,'Acoustic model: baffled-Rayleigh-type sinc source operator on midsurface; not an unbaffled BEM.\n');
fprintf(fid,'W old design = %.12e W; W consistent reoptimized = %.12e W\n',Wold,Wnew);
fprintf(fid,'Weighted volume old=%.9g new=%.9g stop=%s iterations=%d\n', ...
 sum(Rop.element_area.*told)/sum(Rop.element_area)/opt.tmax, ...
 sum(Rop.element_area.*tbest)/sum(Rop.element_area)/opt.tmax,meta.stopReason,meta.niter);
fclose(fid); fprintf('R0_CONSISTENT_PANEL_REOPT_DONE\n');
end
