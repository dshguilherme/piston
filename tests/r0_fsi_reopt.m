function r0_fsi_reopt(cell, zfile, startcode, nmax, tag)
% Experiment 7, optional step: re-optimise one pilot cell with the exact two-sided radiation
% operator AND fluid loading in the state solve (FSI in the loop).
%   cell      : pilot cell name, e.g. 'c30_k3'
%   zfile     : .mat with iwZ = i w s Z_plus on the free dofs (tests/r0_fsi_export_z.py)
%   startcode : 'U' (uniform t = 0.05) or a pilot run code ('XU','XS','SU',...) whose tbest is the start
%   nmax      : GCMMA iteration cap;  tag : output tag
% Optimiser, bounds, volume constraint, dB scaling, stopping rule, normalisation J(start) = Jref
% and the exact factor Fx are those of the pilot (tests/r0_pilot_runs.m); only the state
% solve differs (r0_fsi_eval).  Writes results/r0_revision/fsi/reopt_<cell>_<tag>.mat
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fileparts(root)); addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
pil=fullfile(root,'results','r0_revision','pilot'); out=fullfile(root,'results','r0_revision','fsi');
maxNumCompThreads(4);
th=sscanf(cell,'c%d_'); Jref=5.9464e-7;
ctx=pilot_ctx(th);
F=load(fullfile(pil,['factor_' cell '.mat']),'Fx'); Fx=F.Fx; clear F
Z=load(zfile,'iwZ','s'); s=Z.s;
R0=struct('factor',Fx,'iwZ',Z.iwZ); clear Z
R0vac=struct('factor',Fx);
opt=struct('tmin',0.025,'tmax',0.10,'volfrac',0.5,'nmax',nmax,'stop',true,'ftol',1e-4, ...
    'xchgtol',1e-4,'nhold',3,'window',10,'feastol',1e-6,'Abar',100,'dbfac',10);
if strcmp(startcode,'U'), tstart=0.05*ones(ctx.ne,1);
else, P=load(fullfile(pil,['run_' cell '_' startcode '.mat']),'tbest'); tstart=P.tbest; end
opt.tinit=tstart;
% physical (unnormalised) powers of the start, fluid-loaded and in vacuo (power in the medium of s: x max(s,1) is
% applied in the report, not here; J below is w^2 ||Fx U||^2 with the rho0 c0 factor)
J0=r0_fsi_eval(ctx,tstart,R0); J0v=r0_pilot_eval(ctx,tstart,R0vac);
Ju=r0_fsi_eval(ctx,0.05*ones(ctx.ne,1),R0); Juv=r0_pilot_eval(ctx,0.05*ones(ctx.ne,1),R0vac);
fprintf('== FSI re-opt %s start %s s=%.5f: J_fsi(start) %.5e  J_vac(start) %.5e  J_fsi(uniform) %.5e  J_vac(uniform) %.5e ==\n', ...
    cell,startcode,s,J0,J0v,Ju,Juv);
% FD check of the FSI adjoint at the start
rng(517); d=randn(ctx.ne,1); d=d/norm(d); h=1e-6;
[Jc,g]=r0_fsi_eval(ctx,tstart,R0); fd=(r0_fsi_eval(ctx,tstart+h*d,R0)-r0_fsi_eval(ctx,tstart-h*d,R0))/(2*h);
fprintf('  FD check FSI adjoint: J %.5e fd %+.9e adj %+.9e rel %.2e\n',Jc,fd,g'*d,abs(fd-g'*d)/abs(g'*d));
cstart=Jref/J0; R0.factor=R0.factor*sqrt(cstart);
tic; [tbest,hist,meta]=r0_fsi_track(ctx,'quadvel',R0,'consistent',opt); wt=toc;
R0.factor=Fx;
Jb=r0_fsi_eval(ctx,tbest,R0); Jbv=r0_pilot_eval(ctx,tbest,R0vac);
fprintf('RESULT %s %s: it=%d stop=%s J_fsi(best) %.5e (%.2f dB vs fsi uniform, %.2f dB vs start)  J_vac(best) %.5e  vol %.6f  time %.0f s\n', ...
    cell,tag,meta.niter,meta.stopReason,Jb,10*log10(Jb/Ju),10*log10(Jb/J0),Jbv,meta.final_vol,wt);
save(fullfile(out,['reopt_' cell '_' tag '.mat']),'tbest','hist','meta','opt','cstart','s','J0','J0v','Ju','Juv','Jb','Jbv','startcode','-v7');
fprintf('R0_FSI_REOPT_DONE\n');
end

function ctx=pilot_ctx(th)
spec=struct('R',25,'L',50,'th0',deg2rad(40),'alpha',deg2rad(30), ...
    'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,'deg',[2 2],'reg',[1 1], ...
    'nsub',[40 24],'nquad',[4 4],'bc',[1 2 3 4],'loadtype','point', ...
    'loaduv',[0.5 0.5],'loaddir',[0;1;0],'loadmag',1,'freq',6.53547,'zeta',0.05);
spec.R=10; spec.bc=2; spec.loadtype='pole'; spec.loaddir=[0;0;1]; spec.freq=100;
spec.alpha=deg2rad(th);
ctx=obstacle_build_ctx('spherecap',spec); ctx.alpha=0; ctx.beta=0.1/ctx.w;
end
