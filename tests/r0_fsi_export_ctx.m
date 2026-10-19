function r0_fsi_export_ctx()
% R0 FSI check (Experiment 7): export the pilot's structural model (unchanged) for the
% Python FSI solver.  For theta0 = 30, 90, 120 deg the pilot context (pilot_ctx, copied
% verbatim from tests/r0_pilot_runs.m) is built and its element arrays are saved:
%   K(t) = sum_e Bke t_e^3 + Ske t_e,  M(t) = sum_e Me t_e  (sparse via Rows/Cols),
%   K_D = (1 + i w beta) K - w^2 M  (alpha = 0),  free dofs, load F, w, beta.
% Also saves the in-vacuo uniform state (t = 0.05) computed here with r0_pilot_eval's
% assembly, for the Python assembly check.  Writes results/r0_revision/fsi/ctx_th<deg>.mat
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fileparts(root)); addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
out=fullfile(root,'results','r0_revision','fsi'); if ~exist(out,'dir'), mkdir(out); end
maxNumCompThreads(2);
for th=[30 90 120]
    ctx=pilot_ctx(th);
    Bke=ctx.Bke; Ske=ctx.Ske; Me=ctx.Me; Rows=ctx.Rows; Cols=ctx.Cols; %#ok<NASGU>
    free=ctx.free(:); F=full(ctx.F(:)); w=ctx.w; beta=ctx.beta; alpha=ctx.alpha; ndof=ctx.ndof; ne=ctx.ne; lm=ctx.lm; %#ok<NASGU>
    t=0.05*ones(ctx.ne,1);
    cK=1+1i*ctx.w*ctx.beta; cM=1i*ctx.w*ctx.alpha-ctx.w^2;
    KD=sparse(ctx.Rows(:),ctx.Cols(:),cK*reshape(ctx.Bke.*(t.^3)+ctx.Ske.*t,[],1)+cM*reshape(ctx.Me.*t,[],1),ndof,ndof);
    Uu=complex(zeros(ndof,1)); Uu(ctx.free)=KD(ctx.free,ctx.free)\ctx.F(ctx.free); %#ok<NASGU>
    save(fullfile(out,sprintf('ctx_th%d.mat',th)),'Bke','Ske','Me','Rows','Cols','free','F','w','beta','alpha', ...
        'ndof','ne','lm','Uu','-v7');
    fprintf('theta0=%d: ndof %d, free %d, ne %d, Bke %dx%d, w %.6f, beta %.6e, |F| %.3f\n',th,ndof,numel(free),ne, ...
        size(Bke,1),size(Bke,2),w,beta,norm(F));
end
fprintf('R0_FSI_EXPORT_CTX_DONE\n');
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
