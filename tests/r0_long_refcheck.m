function r0_long_refcheck(cell)
% R0_LONG_REFCHECK  Reference sensitivity of the Experiment-8 optima of one cell.
% For every checkpoint ck_*_<cell>_*.mat in results/r0_revision/longruns, evaluate the run's
% own-objective pick (part A: segment best; part B: best exact iterate) and its best exact iterate
% under: the exact factor used in the runs (Richardson L 4000/8000, N = 48), the L = 8000 factor,
% the N = 32 factor, the consistent sinc used in the runs and the spectral (Funk-Hecke) sinc.
% Same quantities and same factor files as the pilot 'eval' stage (exp6, pilot_reference.csv).
% Output: longruns/refcheck_<cell>.csv (dB vs the uniform design of the same factor).
maxNumCompThreads(2); warning('off','backtrace');
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fileparts(root)); addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
pilot=fullfile(root,'results','r0_revision','pilot'); out=fullfile(root,'results','r0_revision','longruns');
th=sscanf(cell,'c%d_');
spec=struct('R',25,'L',50,'th0',deg2rad(40),'alpha',deg2rad(30), ...
    'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,'deg',[2 2],'reg',[1 1], ...
    'nsub',[40 24],'nquad',[4 4],'bc',[1 2 3 4],'loadtype','point', ...
    'loaduv',[0.5 0.5],'loaddir',[0;1;0],'loadmag',1,'freq',6.53547,'zeta',0.05);
spec.R=10; spec.bc=2; spec.loadtype='pole'; spec.loaddir=[0;0;1]; spec.freq=100; spec.alpha=deg2rad(th);
ctx=obstacle_build_ctx('spherecap',spec); ctx.alpha=0; ctx.beta=0.1/ctx.w;
F=load(fullfile(pilot,['factor_' cell '.mat']),'Fx','Fx8000','Fx32','Fs');
P=load(fullfile(pilot,['ops_' cell '.mat']),'Fs_c');
files=dir(fullfile(out,['ck_*_' cell '_*.mat']));
files=files(~contains({files.name},'_tmp'));
T=0.05*ones(ctx.ne,1); lab={'uniform'}; its=0;
for i=1:numel(files)
    L=load(fullfile(out,files(i).name),'S'); S=L.S;
    nm=erase(erase(files(i).name,'ck_'),'.mat');
    if startsWith(nm,'A_') && ~isempty(S.tbest_seg)
        [~,ib]=min(S.H.Jown); T(:,end+1)=S.tbest_seg; lab{end+1}=[nm ':pick']; its(end+1)=S.H.it(ib); %#ok<AGROW>
        if isfield(S,'tfbest_seg') && ~isempty(S.tfbest_seg)
            T(:,end+1)=S.tfbest_seg; lab{end+1}=[nm ':pick_feasible']; its(end+1)=S.itfbest_seg; %#ok<AGROW>
        end
    end
    if ~isempty(S.txbest)
        T(:,end+1)=S.txbest; lab{end+1}=[nm ':bestexact']; its(end+1)=S.itxbest; %#ok<AGROW>
        if mean(S.txbest)/0.10>0.5*(1+1e-5)
            T(:,end+1)=project_volume(S.txbest); lab{end+1}=[nm ':bestexact_scaled_to_vol0.5']; its(end+1)=S.itxbest; %#ok<AGROW>
        end
    end
    if isfield(S,'txfbest') && ~isempty(S.txfbest)
        T(:,end+1)=S.txfbest; lab{end+1}=[nm ':bestexact_feasible']; its(end+1)=S.itxfbest; %#ok<AGROW>
    elseif ~isempty(S.Tit)
        % checkpoints written before feasible tracking existed: best feasible among the stored
        % every-10th designs
        [ok,loc]=ismember(S.Tit,S.H.it);
        cand=find(ok & mean(S.T,1)/0.10<=0.5*(1+1e-5));
        if ~isempty(cand)
            [~,j]=min(S.H.Jx(loc(cand))); j=cand(j);
            T(:,end+1)=S.T(:,j); lab{end+1}=[nm ':bestexact_feasible_stored10']; its(end+1)=S.Tit(j); %#ok<AGROW>
        end
    end
end
w2=ctx.w^2; nd=size(T,2); R=zeros(nd,5);
for j=1:nd
    t=T(:,j); KE=ctx.Bke.*(t.^3)+ctx.Ske.*t; ME=ctx.Me.*t;
    A=sparse(ctx.Rows(:),ctx.Cols(:),(1+1i*ctx.w*ctx.beta)*KE(:)-ctx.w^2*ME(:),ctx.ndof,ctx.ndof);
    u=complex(zeros(ctx.ndof,1)); u(ctx.free)=A(ctx.free,ctx.free)\ctx.F(ctx.free); us=single(u);
    R(j,:)=w2*[sum(abs(F.Fx*u).^2), sum(abs(double(F.Fx8000*us)).^2), sum(abs(double(F.Fx32*us)).^2), ...
        sum(abs(P.Fs_c*u).^2), sum(abs(double(F.Fs*us)).^2)];
end
f=fopen(fullfile(out,['refcheck_' cell '.csv']),'w');
fprintf(f,'cell,design,iter,vol,exact_dB,exact_L8000_dB,exact_N32_dB,sinc_dB,sinc_spec_dB,eL_rel,eN_rel\n');
for j=1:nd
    d=10*log10(R(j,:)./R(1,:));
    fprintf(f,'%s,%s,%d,%.6f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4e,%.4e\n',cell,lab{j},its(j),mean(T(:,j))/0.10,d, ...
        abs(R(j,2)-R(j,1))/R(j,1),abs(R(j,3)-R(j,1))/R(j,1));
end
fclose(f);
fprintf('R0_LONG_refcheck_%s_DONE (%d designs)\n',cell,nd);
end

function t=project_volume(t)
% uniform thinning to volume fraction 0.5 within the bounds [0.025, 0.10] (post-hoc feasibility
% check of designs on which GCMMA's soft volume constraint was violated)
for k=1:50
    v=mean(t)/0.10;
    if v<=0.5, break; end
    fr=t>0.025;
    t(fr)=t(fr)-(v-0.5)*0.10*numel(t)/nnz(fr);
    t=max(t,0.025);
end
end
