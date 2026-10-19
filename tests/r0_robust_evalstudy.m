function r0_robust_evalstudy(cell, zdir, nsamp, sigmas, tagin, fsisig, onlyrobust)
% R0_ROBUST_EVALSTUDY  Experiment 10, task 3: fragility of existing optimised designs under random
% relative thickness errors h_e = t_e (1 + sigma xi_e), xi_e ~ N(0,1) independent per element.
% No optimisation.  For every design and sigma, nsamp fresh samples (seed 7777, the SAME samples for
% every design and sigma, so comparisons between designs are paired) are solved; each state gives
% the exact power (Fx, Richardson), the exact power with the L = 8000 factor (reference check) and
% the sinc power; with zdir given, the exact power with fluid loading in the state solve (device
% s = 1 and, for the scaled cells, physical air density) is also computed, at the sigmas listed in
% fsisig (default: all; NaN in the CSV elsewhere).  onlyrobust = true evaluates only the uniform
% design and the robust-run designs (ck_<cell>_*.mat in results/r0_revision/robust).
% Designs: uniform; anchor c30_k18: Experiment-5 S0, SX, X0, XS (+ SR, XR); other cells: pilot XU,
% XS, SU, SX; every Experiment-8 long-run checkpoint of the cell (best feasible exact design,
% 'L<run>'); for c30_k3 the Experiment-7 FSI-in-the-loop optima ('Fdev', 'Fphys').
% Output: results/r0_revision/robust/evalstudy_<cell>.csv (one row per design, sigma, sample;
% sigma = 0 rows are the nominal designs), evalstudy_<cell>.log
if nargin<3||isempty(nsamp), nsamp=200; end
if nargin<4||isempty(sigmas), sigmas=[0.001 0.003 0.01 0.03]; end
if nargin<5, tagin=''; end
if nargin<6||isempty(fsisig), fsisig=[0 sigmas]; end
maxNumCompThreads(2); warning('off','backtrace');
[ctx,ops,G]=r0_robust_setup(cell);
tag=cell; if ~isempty(tagin), tag=[cell '_' tagin]; end
lf=fopen(fullfile(G.out,['evalstudy_' tag '.log']),'a');
lg=@(varargin) logboth(lf,varargin{:});
lg('== r0_robust_evalstudy %s %s: nsamp %d, sigmas %s ==\n',tag,datestr(now,'yyyy-mm-ddTHH:MM:SS'),nsamp,mat2str(sigmas));
[T,names]=designs(cell,G,ctx.ne);
if nargin>=7 && onlyrobust
    keepd=strcmp(names,'uniform') | startsWith(names,'R');
    T=T(:,keepd); names=names(keepd);
end
F8=load(fullfile(G.pilot,['factor_' cell '.mat']),'Fx8000'); F8=double(F8.Fx8000);
extra={ops.Sraw.factor, F8}; clear F8
variants={struct('name','vac','iwZ',[],'s',0)};
if ~isempty(zdir)
    for v={'dev','phys'}
        fn=fullfile(zdir,['iwz_' cell '_' v{1} '.mat']);
        if exist(fn,'file'), Z=load(fn,'iwZ','s'); variants{end+1}=struct('name',v{1},'iwZ',Z.iwZ,'s',Z.s); end %#ok<AGROW>
    end
    clear Z
end
rng(7777); Xi=randn(ctx.ne,nsamp);
fn=fullfile(G.out,['evalstudy_' tag '.csv']);
f=fopen(fn,'w');
fprintf(f,'cell,design,sigma,sample,vol,P_exact,P_exact_L8000,P_sinc');
for iv=2:numel(variants), fprintf(f,',P_exact_fsi_%s',variants{iv}.name); end
fprintf(f,'\n');
for id=1:numel(names)
    t=T(:,id); tc=tic;
    for sg=[0 sigmas]
        if sg==0, X=zeros(ctx.ne,1); else, X=Xi; end
        R=struct('sigma',sg,'Xi',X,'solver','gmres','extra',{extra});
        if sg==0, R.solver='direct'; end
        [~,~,in]=evalJ(ctx,t,ops.Xraw,R);
        Pv=[in.P in.Pext(:,2) in.Pext(:,1)];
        nfb=in.nfallback;
        for iv=2:numel(variants)
            if ~any(abs(fsisig-sg)<1e-12), Pv=[Pv nan(size(X,2),1)]; continue; end %#ok<AGROW>
            Ff=struct('factor',ops.Xraw.factor,'iwZ',variants{iv}.iwZ);
            Rf=R; Rf.extra={}; Rf.solver='gmres';
            [~,~,inf2]=evalJ(ctx,t,Ff,Rf);
            Pv=[Pv inf2.P]; nfb=nfb+inf2.nfallback; %#ok<AGROW>
        end
        for j=1:size(X,2)
            fprintf(f,'%s,%s,%.4g,%d,%.6f',cell,names{id},sg,j,mean(t.*(1+sg*X(:,j)))/0.10);
            fprintf(f,',%.8e',Pv(j,:)); fprintf(f,'\n');
        end
        lg('  %-8s sigma %.3f: P_exact mean %.4e median %.4e [p5 %.4e, p95 %.4e] (nominal-design ref below), fsi cols %d, fallbacks %d\n', ...
            names{id},sg,mean(Pv(:,1)),median(Pv(:,1)),prctile(Pv(:,1),5),prctile(Pv(:,1),95),numel(variants)-1,nfb);
    end
    lg('  design %s done in %.0f s\n',names{id},toc(tc));
end
fclose(f);
lg('R0_ROBUST_EVALSTUDY_%s_DONE\n',tag);
fclose(lf);
end

function [J,dJ,in]=evalJ(ctx,t,op,R)
clear r0_robust_eval
[J,~,in]=r0_robust_eval(ctx,t,op,R); dJ=[];
end

function [T,names]=designs(cell,G,ne)
T=0.05*ones(ne,1); names={'uniform'};
if strcmp(cell,'c30_k18')
    for c={'S0','SX','X0','XS','SR','XR'}
        P=load(fullfile(G.exp5,['run_' c{1} '.mat']),'tbest'); T(:,end+1)=P.tbest(:); names{end+1}=c{1}; %#ok<AGROW>
    end
else
    for c={'XU','XS','SU','SX'}
        P=load(fullfile(G.pilot,['run_' cell '_' c{1} '.mat']),'tbest'); T(:,end+1)=P.tbest(:); names{end+1}=c{1}; %#ok<AGROW>
    end
end
L=dir(fullfile(G.long,['ck_A_' cell '_*.mat']));
for k=1:numel(L)
    if contains(L(k).name,'_tmp'), continue; end
    S=load(fullfile(L(k).folder,L(k).name),'S'); S=S.S;
    nm=regexprep(L(k).name,['ck_A_' cell '_(\w+)\.mat'],'$1');
    if isfield(S,'txfbest') && ~isempty(S.txfbest), tb=S.txfbest; else, tb=S.tbest_seg; end
    T(:,end+1)=tb(:); names{end+1}=sprintf('L%s_it%d',nm,S.it); %#ok<AGROW>
end
Rk=dir(fullfile(G.out,['ck_' cell '_*.mat']));
for k=1:numel(Rk)
    if contains(Rk(k).name,'_tmp'), continue; end
    S=load(fullfile(Rk(k).folder,Rk(k).name),'S'); S=S.S;
    nm=regexprep(Rk(k).name,['ck_' cell '_(.+)\.mat'],'$1');
    if ~isempty(S.trobbest), tb=S.trobbest; else, tb=S.t; end
    T(:,end+1)=tb(:); names{end+1}=sprintf('R%s_it%d',nm,S.itrobbest); %#ok<AGROW>
end
if strcmp(cell,'c30_k3')
    for c={'dev','phys'}
        fn=fullfile(G.fsi,['reopt_c30_k3_' c{1} '_XU.mat']);
        if exist(fn,'file'), P=load(fn,'tbest'); T(:,end+1)=P.tbest(:); names{end+1}=['F' c{1}]; end %#ok<AGROW>
    end
end
end

function logboth(f,varargin)
s=sprintf(varargin{:}); fprintf('%s',s); fprintf(f,'%s',s);
end
