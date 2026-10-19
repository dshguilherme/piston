% REPRODUCE-PACKAGE COPY. Do not edit numerics here; edit the original and re-copy.
% Original: paper_jsv/reproduce_plan/rescued_scratch/prof_lu.m  (sha256 f4d7aa50b63b99d5)
% Changed: original lines 2-4 (root=..., addpath(genpath(geopdes)), addpath(fileparts(root)), src, tests) replaced by
%          addpath(<package folder>); c = config(true); root = c.ROOT;  (config.m puts GeoPDEs, ROOT, src and tests on the path).
%          Everything else unchanged. Nothing is read from or written to results (timings are printed to stdout only).
% Profiling helper only: not used for any paper number. (This script builds its own context; the sibling prof_eval.m
%          loads results/r0_revision/pilot/factor_c30_k3.mat and ops_c30_k3.mat, release assets.)
maxNumCompThreads(2);
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..')); c = config(true); root = c.ROOT;
spec=struct('R',25,'L',50,'th0',deg2rad(40),'alpha',deg2rad(30), ...
    'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,'deg',[2 2],'reg',[1 1], ...
    'nsub',[40 24],'nquad',[4 4],'bc',[1 2 3 4],'loadtype','point', ...
    'loaduv',[0.5 0.5],'loaddir',[0;1;0],'loadmag',1,'freq',6.53547,'zeta',0.05);
spec.R=10; spec.bc=2; spec.loadtype='pole'; spec.loaddir=[0;0;1]; spec.freq=100; spec.alpha=deg2rad(30);
ctx=obstacle_build_ctx('spherecap',spec); ctx.alpha=0; ctx.beta=0.1/ctx.w;
t=0.05*ones(ctx.ne,1)+0.01*randn(ctx.ne,1); t=min(max(t,0.025),0.1);
cK = 1 + 1i*ctx.w*ctx.beta; cM = - ctx.w^2;
EmatK = ctx.Bke.*(t.^3) + ctx.Ske.*t; EmatM = ctx.Me.*t;
KD = sparse(ctx.Rows(:), ctx.Cols(:), cK*EmatK(:) + cM*EmatM(:), ctx.ndof, ctx.ndof);
A=KD(ctx.free,ctx.free); b=ctx.F(ctx.free);
[L,U,P,Q]=lu(A); x0=Q*(U\(L\(P*b)));
tic; for i=1:5, [L,U,P,Q]=lu(A); end; fprintf('umf 4out: %.4f\n',toc/5);
tic; for i=1:5, x=Q*(U\(L\(P*b))); end; fprintf('  solve: %.4f\n',toc/5);
tic; for i=1:5, x=A\b; end; fprintf('backslash: %.4f relx %.1e\n',toc/5,norm(x-x0)/norm(x0));
for nm={'amd','symamd','symrcm'}
  p=feval(nm{1},A); Ap=A(p,p);
  tic; for i=1:5, [L,U,P]=lu(Ap); end; tt=toc/5; y=U\(L\(P*b(p))); x=zeros(size(b)); x(p)=y;
  fprintf('%s + lu3: %.4f nnz %d relx %.1e\n',nm{1},tt,nnz(L)+nnz(U),norm(x-x0)/norm(x0));
  tic; for i=1:5, [L,U,P]=lu(Ap,0.01); end; tt=toc/5; y=U\(L\(P*b(p))); x(p)=y;
  fprintf('%s + lu3 thr0.01: %.4f nnz %d relx %.1e\n',nm{1},tt,nnz(L)+nnz(U),norm(x-x0)/norm(x0));
end
tic; for i=1:5, [L,U,P,Q]=lu(A,[0.1 0.001]); end; fprintf('umf thr: %.4f\n',toc/5);
tic; for i=1:5, d=decomposition(A,'lu'); end; fprintf('decomp lu: %.4f\n',toc/5);
tic; for i=1:5, d=decomposition(A,'ldl'); end; x=d\b; fprintf('decomp ldl: %.4f relx %.1e\n',toc/5,norm(x-x0)/norm(x0));
Af=full(A); tic; [L,U,P]=lu(Af); fprintf('dense lu: %.4f\n',toc);
