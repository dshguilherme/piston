function r0_frozen_field_audit()
% Separate acoustic quadrature changes from structural-state changes.
% Preserves legacy BCs and designs explicitly; not a validation of those BCs.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fileparts(root)); addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
out=fullfile(root,'results','r0_revision'); if ~exist(out,'dir'), mkdir(out); end
fprintf('Builder: %s\n',which('obstacle_build_ctx'));
fprintf('Legacy radiation builder: %s\n',which('articleB_r0_build'));
rows={};
for geometry={'cylpanel','spherecap'}
    kind=geometry{1};
    spec=struct('R',25,'L',50,'th0',deg2rad(40),'alpha',deg2rad(30), ...
        'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,'deg',[2 2],'reg',[1 1], ...
        'nsub',[40 24],'nquad',[4 4],'bc',[1 2 3 4],'loadtype','point', ...
        'loaduv',[0.5 0.5],'loaddir',[0;1;0],'loadmag',1,'freq',6.53547,'zeta',0.05);
    if strcmp(kind,'spherecap')
        spec.R=10; spec.bc=2; spec.loadtype='pole'; spec.loaddir=[0;0;1]; spec.freq=100;
    end
    ctx=obstacle_build_ctx(kind,spec); ctx.alpha=0; ctx.beta=0.1/ctx.w;
    designs=0.05*ones(ctx.ne,1); labels={'uniform'};
    if strcmp(kind,'cylpanel')
        T=load(fullfile(root,'results','t3_longhorizon.mat'),'tbest');
        designs(:,end+1)=T.tbest(:); labels{end+1}='legacy_long_horizon';
    else
        S=load(fullfile(root,'results','articleB_r0_thesisdamp_spherecap_40x24_f100_partial.mat'));
        if isfield(S,'L')
            for j=1:numel(S.L)
                if isfield(S.L(j),'tbest') && numel(S.L(j).tbest)==ctx.ne
                    designs(:,end+1)=S.L(j).tbest(:); labels{end+1}=S.L(j).label;
                end
            end
        end
    end
    U=complex(zeros(ctx.ndof,size(designs,2)));
    for j=1:size(designs,2)
        t=designs(:,j); KE=ctx.Bke.*(t.^3)+ctx.Ske.*t; ME=ctx.Me.*t;
        A=sparse(ctx.Rows(:),ctx.Cols(:),(1+1i*ctx.w*ctx.beta)*KE(:)-ctx.w^2*ME(:),ctx.ndof,ctx.ndof);
        U(ctx.free,j)=A(ctx.free,ctx.free)\ctx.F(ctx.free);
    end
    legacy=articleB_r0_build(ctx,1.21,343,128);
    Z=legacy.area.*(legacy.C*U);
    legacy_power=ctx.w^2*legacy.pref*real(sum(conj(Z).*(legacy.Kfull*Z),1));
    for nq=[4,6,8]
        samples=r0_surface_samples(ctx,U,[nq nq]);
        Y=samples.normal_displacement;
        if nq==4
            direct=zeros(size(Y)); qn=ctx.msh.nqn;
            spOrig=sp_precompute_param(ctx.space,ctx.msh,'value',true);
            for ee=1:ctx.ne
                for cc=1:3
                    ids=ctx.lm(ee,(cc-1)*size(ctx.lm,2)/3+(1:size(ctx.lm,2)/3));
                    loc=(cc-1)*numel(ids)+(1:numel(ids));
                    vals=reshape(spOrig.shape_functions(cc,:,loc,ee),qn,[]);
                    nn=reshape(ctx.msh.normal(cc,:,ee),qn,1);
                    direct((ee-1)*qn+(1:qn),:)=direct((ee-1)*qn+(1:qn),:)+nn.*(vals*U(ids,:));
                end
            end
            trace_err=norm(direct-Y,'fro')/max(norm(direct,'fro'),eps);
            fprintf('trace pullback check %s relative error %.3e\\n',kind,trace_err);
            assert(trace_err<1e-11,'Refreshed quadrature trace does not reproduce original context.');
        end
        joint=zeros(ctx.ne,size(U,2));
        for j=1:size(U,2)
            joint(:,j)=sum(reshape(samples.weights.*Y(:,j),samples.nq,ctx.ne),1)'./samples.area;
        end
        joint_Z=samples.area.*joint;
        joint_power=ctx.w^2*legacy.pref*real(sum(conj(joint_Z).*(legacy.Kfull*joint_Z),1));
        constant=kron(joint,ones(samples.nq,1));
        for na=[24,40,64]
            op=r0_angular_gram(samples.x,samples.weights,[constant,Y],ctx.w,1.21,343,na,2*na);
            power=ctx.w^2*sum(abs(op.factor).^2,1);
            nf=size(U,2);
            for j=1:nf
                rows(end+1,:)={kind,labels{j},nq,na,legacy_power(j),joint_power(j),power(j),power(nf+j)}; %#ok<AGROW>
                fprintf('%s %s q=%d angular=%d: legacy=%.9e joint=%.9e pair=%.9e consistent=%.9e\n', ...
                    kind,labels{j},nq,na,legacy_power(j),joint_power(j),power(j),power(nf+j));
            end
        end
    end
    save(fullfile(out,[kind '_frozen_states.mat']),'U','designs','labels','spec');
end
T=cell2table(rows,'VariableNames',{'geometry','design','surface_order','angular_order', ...
    'legacy_centroid_W','joint_average_centroid_W','finite_patch_W','consistent_spline_W'});
writetable(T,fullfile(out,'frozen_field_quadrature.csv'));
save(fullfile(out,'frozen_field_quadrature.mat'),'T');
fprintf('R0_FROZEN_FIELD_AUDIT_DONE\n');
end
