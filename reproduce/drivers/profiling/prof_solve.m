% REPRODUCE-PACKAGE COPY. Do not edit numerics here; edit the original and re-copy.
% Original: paper_jsv/reproduce_plan/rescued_scratch/prof_solve.m  (sha256 958a6deabacfd4aa)
% Changed: original lines 3-5 (root=..., addpath(genpath(geopdes)), addpath(fileparts(root)), src, tests) replaced, inside the function, by
%          addpath(<package folder>); c = config(true); root = c.ROOT;  (config.m puts GeoPDEs, ROOT, src and tests on the path).
%          Everything else unchanged. Nothing is read from or written to results (timings are printed to stdout only).
% Profiling helper only: not used for any paper number. (Builds its own context; the sibling prof_eval.m loads
%          results/r0_revision/pilot/factor_c30_k3.mat and ops_c30_k3.mat, release assets.)
function prof_solve(nthr)
maxNumCompThreads(nthr);
addpath(fullfile(fileparts(mfilename('fullpath')),'..','..')); c = config(true); root = c.ROOT;
spec=struct('R',10,'L',50,'th0',deg2rad(40),'alpha',deg2rad(30), ...
    'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,'deg',[2 2],'reg',[1 1], ...
    'nsub',[40 24],'nquad',[4 4],'bc',2,'loadtype','pole', ...
    'loaduv',[0.5 0.5],'loaddir',[0;0;1],'loadmag',1,'freq',100,'zeta',0.05);
ctx=obstacle_build_ctx('spherecap',spec); ctx.alpha=0; ctx.beta=0.1/ctx.w;
t=0.05*ones(ctx.ne,1);
EmatK = ctx.Bke.*(t.^3) + ctx.Ske.*t; Ks = sparse(ctx.Rows(:), ctx.Cols(:), EmatK(:), ctx.ndof, ctx.ndof);
EmatM = ctx.Me.*t; M = sparse(ctx.Rows(:), ctx.Cols(:), EmatM(:), ctx.ndof, ctx.ndof);
KD=(1+1i*ctx.w*ctx.beta)*Ks-ctx.w^2*M; A=KD(ctx.free,ctx.free); b=ctx.F(ctx.free);
fprintf('n %d nnz %d\n',size(A,1),nnz(A));
tic; for k=1:3, x=A\b; end; fprintf('backslash %.3f s\n',toc/3);
tic; for k=1:3, dA=decomposition(A,'lu'); x=dA\b; end; fprintf('decomposition lu %.3f s\n',toc/3);
p=symamd(A); tic; for k=1:3, x2=zeros(size(b)); x2(p)=A(p,p)\b(p); end; fprintf('symamd+backslash %.3f s rel %.1e\n',toc/3,norm(x2-x)/norm(x));
tic; for k=1:3, [L,U,P,Q]=lu(A); x3=Q*(U\(L\(P*b))); end; fprintf('lu 4-out %.3f s rel %.1e\n',toc/3,norm(x3-x)/norm(x));
tic; for k=1:3, [L,D,P]=ldl(A); end; fprintf('ldl %.3f s\n',toc/3);
end
