function r0_nystrom_rank_sweep()
% Test Nyström compression separately from quadrature and acoustic model error.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
out=fullfile(root,'results','r0_revision');
rows={};
for kind={'cylpanel','spherecap'}
    name=kind{1};
    spec=struct('R',25,'L',50,'th0',deg2rad(40),'alpha',deg2rad(30), ...
      'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,'deg',[2 2],'reg',[1 1], ...
      'nsub',[40 24],'nquad',[4 4],'bc',[1 2 3 4],'loadtype','point', ...
      'loaduv',[0.5 0.5],'loaddir',[0;1;0],'loadmag',1,'freq',6.53547,'zeta',0.05);
    if strcmp(name,'spherecap'), spec.R=10; spec.bc=2; spec.loadtype='pole'; spec.loaddir=[0;0;1]; spec.freq=100; end
    ctx=obstacle_build_ctx(name,spec); ctx.alpha=0; ctx.beta=0.1/ctx.w;
    if strcmp(name,'spherecap')
        D=load(fullfile(root,'results','articleB_r0_thesisdamp_spherecap_40x24_f100_partial.mat'),'L');
        picks=[]; tags={};
        for j=1:numel(D.L)
            if any(strcmp(D.L(j).label,{'quadvel-rad','quadvel-I','quadvel-ones'}))
                picks(:,end+1)=D.L(j).tbest(:); tags{end+1}=D.L(j).label; %#ok<AGROW>
            end
        end
    else
        D=load(fullfile(root,'results','t3_longhorizon.mat'),'tbest'); picks=D.tbest(:); tags={'legacy_long_horizon'};
    end
    U=zeros(ctx.ndof,size(picks,2));
    % Assemble each structural state independently of the acoustic compression.
    for j=1:size(picks,2)
      t=picks(:,j); ke=ctx.Bke.*t.^3+ctx.Ske.*t; me=ctx.Me.*t;
      A=sparse(ctx.Rows(:),ctx.Cols(:),((1+1i*ctx.w*ctx.beta)*ke(:)-ctx.w^2*me(:)),ctx.ndof,ctx.ndof);
      U(ctx.free,j)=A(ctx.free,ctx.free)\ctx.F(ctx.free);
    end
    full=articleB_r0_build(ctx,1.21,343,128); Z=full.area.*(full.C*U);
    qref=real(sum(conj(Z).*(full.Kfull*Z),1));
    ranks=[16 32 48 64 96 128 192 256 384 512 768 960];
    for m=ranks
      op=articleB_r0_build(ctx,1.21,343,m);
      for j=1:size(picks,2)
        za=op.area.*(op.C*U(:,j)); ka=op.KS.'*(op.KSSinv*(op.KS*za));
        q=real(za'*ka); err=(q-qref(j))/max(abs(qref(j)),realmin);
        rows(end+1,:)={name,tags{j},m,qref(j),q,err,norm(full.Kfull-(op.KS.'*op.KSSinv*op.KS),'fro')/norm(full.Kfull,'fro')}; %#ok<AGROW>
      end
    end
end
T=cell2table(rows,'VariableNames',{'geometry','field','rank','dense_quadratic','compressed_quadratic', ...
 'relative_quadratic_error','relative_kernel_frobenius_error'});
writetable(T,fullfile(out,'nystrom_rank_sweep.csv')); save(fullfile(out,'nystrom_rank_sweep.mat'),'T');
fprintf('R0_NYSTROM_RANK_SWEEP_DONE\n'); disp(T);
end
