% REPRODUCE-PACKAGE COPY. Do not edit numerics here; edit the original and re-copy.
% Original: paper_jsv/reproduce_plan/rescued_scratch/prof_eval.m  (sha256 d971ec8021672597)
% Changed: original lines 2-4 (root=..., addpath(genpath(geopdes)), addpath(fileparts(root)), src, tests) replaced by
%          addpath(<package folder>); c = config(true); root = c.ROOT;  (config.m puts GeoPDEs, ROOT, src and tests on the path).
%          Everything else unchanged. Nothing is written (timings are printed to stdout only).
% Profiling helper only: not used for any paper number. Loads results/r0_revision/pilot/factor_c30_k3.mat and
%          results/r0_revision/pilot/ops_c30_k3.mat (release assets, not in git).
maxNumCompThreads(2);
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..')); c = config(true); root = c.ROOT;
out=fullfile(root,'results','r0_revision','pilot');
spec=struct('R',25,'L',50,'th0',deg2rad(40),'alpha',deg2rad(30), ...
    'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,'deg',[2 2],'reg',[1 1], ...
    'nsub',[40 24],'nquad',[4 4],'bc',[1 2 3 4],'loadtype','point', ...
    'loaduv',[0.5 0.5],'loaddir',[0;1;0],'loadmag',1,'freq',6.53547,'zeta',0.05);
spec.R=10; spec.bc=2; spec.loadtype='pole'; spec.loaddir=[0;0;1]; spec.freq=100; spec.alpha=deg2rad(30);
tic; ctx=obstacle_build_ctx('spherecap',spec); ctx.alpha=0; ctx.beta=0.1/ctx.w; fprintf('ctx %.1fs\n',toc);
F=load(fullfile(out,'factor_c30_k3.mat'),'Fx'); R0=struct('factor',F.Fx);
P=load(fullfile(out,'ops_c30_k3.mat'),'Fs_c'); fprintf('Fx %dx%d  Fs_c %dx%d\n',size(F.Fx),size(P.Fs_c));
whos
t=0.05*ones(ctx.ne,1)+0.01*randn(ctx.ne,1); t=min(max(t,0.025),0.1);
r0_pilot_eval(ctx,t,R0);
tic; for i=1:10, J=r0_pilot_eval(ctx,t,R0); end; fprintf('J only: %.4f s\n',toc/10);
tic; for i=1:10, [J,g]=r0_pilot_eval(ctx,t,R0); end; fprintf('J+grad: %.4f s\n',toc/10);
cK = 1 + 1i*ctx.w*ctx.beta; cM = 1i*ctx.w*ctx.alpha - ctx.w^2;
tic; for i=1:10, EmatK = ctx.Bke.*(t.^3) + ctx.Ske.*t; EmatM = ctx.Me.*t; KD = sparse(ctx.Rows(:), ctx.Cols(:), cK*EmatK(:) + cM*EmatM(:), ctx.ndof, ctx.ndof); end; fprintf('assembly: %.4f s (nnz entries %d)\n',toc/10,numel(ctx.Rows));
A=KD(ctx.free,ctx.free);
tic; for i=1:10, [L,U,PP,Q]=lu(A); end; fprintf('lu: %.4f s, nnz(L+U) %d\n',toc/10,nnz(L)+nnz(U));
u=zeros(ctx.ndof,1); u(ctx.free)=Q*(U\(L\(PP*ctx.F(ctx.free))));
tic; for i=1:10, R0U=R0.factor'*(R0.factor*u); end; fprintf('factor matvec: %.4f s\n',toc/10);
fprintf('ne %d ndof %d nfree %d lm %dx%d\n',ctx.ne,ctx.ndof,numel(ctx.free),size(ctx.lm));
