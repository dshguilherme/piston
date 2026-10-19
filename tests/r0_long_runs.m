function r0_long_runs(job)
% R0_LONG_RUNS  Experiment 8: long optimisation runs and multi-step strategies on the regime-pilot
% cells (exp6). One job per MATLAB session; resumable from checkpoints.
%
%   r0_long_runs('validate')
%       r0_long_eval vs r0_pilot_eval (bitwise), r0_long_core vs the stored pilot run c30_k3 XU
%       (first 30 iterations), and chunked (15 + 15 with a save/load in between) vs straight 30.
%   r0_long_runs('A <cell> <run> <target>')       run in {XU, XS, SU, SX}
%       Part A, long run: warm restart from the stored pilot optimum (run_<cell>_<run>.mat tbest;
%       anchor c30_k18 uses the Experiment-5 run files) and continue with the same rung up to
%       <target> iterations IN TOTAL (pilot iterations included). At the pilot-to-long restart the
%       GCMMA state is RESET (the pilot did not store it) and the start is rescaled to J = Jref
%       (pilot warm-start rule); inside the long run every chunk CARRIES the full optimiser state.
%   r0_long_runs('B <cell> <strategy> <budget>')
%       Part B, multi-step strategy from the uniform design with a total budget of <budget>
%       outer iterations (stages defined in strategy_plan below, fixed before any B result).
%       Every stage change is a warm restart under the pilot rule (state reset, J(start) = Jref).
%
% Records (results/r0_revision/longruns/): ck_<kind>_<cell>_<name>.mat (checkpoint = full state
% and history), hist_<kind>_<cell>_<name>.csv (one row per outer iteration with exact, sinc and
% ERP powers of the target cell in dB vs uniform), chunks_<kind>.csv and chunks.log (one line per
% chunk). A file longruns/STOP makes a running job stop cleanly at the next 25-iteration mark.
maxNumCompThreads(2); warning('off','backtrace');
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fileparts(root)); addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
G.pilot=fullfile(root,'results','r0_revision','pilot');
G.exp5=fullfile(root,'results','r0_revision','exact_opt');
G.out=fullfile(root,'results','r0_revision','longruns');
if ~exist(G.out,'dir'), mkdir(G.out); end
G.Jref=5.9464e-7;
G.opt=struct('tmin',0.025,'tmax',0.10,'volfrac',0.5,'nmax',inf,'stop',true,'ftol',1e-4, ...
    'xchgtol',1e-4,'nhold',3,'window',10,'feastol',1e-6,'Abar',100,'dbfac',10);
tok=strsplit(strtrim(job));
if strcmp(tok{1},'validate'), validate(G); return; end
if strcmp(tok{1},'R'), r0_long_refcheck(tok{2}); return; end
kind=tok{1}; cell=tok{2}; name=tok{3}; target=str2double(tok{4});
[th,kR]=cell_def(cell);
ctx=pilot_ctx(th); ctx.cachekey=sprintf('th%d',th);
G.opt.vweights=ones(ctx.ne,1);
ops=load_ops(G,cell);
tag=sprintf('%s_%s_%s',kind,cell,name);
ckfile=fullfile(G.out,['ck_' tag '.mat']);
stopfile=fullfile(G.out,'STOP');
switch kind
    case 'A'
        plan=struct('rung',name(1),'src',cell,'n',inf,'move',0,'switch',[]);
        chunk=600;
    case 'B'
        plan=strategy_plan(name,cell,target);
        chunk=300;
    otherwise, error('unknown job kind %s',kind);
end
if exist(ckfile,'file')
    L=load(ckfile,'S'); S=L.S; clear L
    fprintf('== resume %s from checkpoint: it %d, stage %d (stage it %d), done %d ==\n',tag,S.it,S.stage,S.stage_it,S.done);
else
    S=new_state(ctx,G);
    if strcmp(kind,'A')
        fn=runfile(cell,name,G.pilot,G.exp5);
        P=load(fn,'tbest','hist'); S.t=P.tbest(:); S.it=numel(P.hist.J); S.uniform_start=false;
        S.source=fn;
    else
        S.t=0.05*ones(ctx.ne,1); S.it=0; S.uniform_start=true; S.source='uniform';
    end
    S.it0=S.it; S.stage=1;
    fprintf('== new %s (theta0 %d, kR %.4f), start %s at it %d, target %d ==\n',tag,th,kR,S.source,S.it,target);
end
srcops=struct();
while ~S.done && S.it<target
    st=plan(S.stage);
    [Fbase,srcops]=stage_factor(G,cell,ops,st,srcops);
    if S.stage_it==0 && ~S.stage_ready
        S.st=[]; S.git=0; S.p0=[]; S.nhold=0; S.f0seg=[]; S.Jbest_seg=inf; S.tbest_seg=S.t; S.swcount=0;
        S.Jfbest_seg=inf; S.tfbest_seg=[]; S.itfbest_seg=0;
        if S.uniform_start && S.it==0
            S.cstart=1;
        else
            S.cstart=G.Jref/r0_long_eval(ctx,S.t,Fbase);
        end
        S.stage_ready=true;
        S.stage_log{S.stage}=sprintf('stage %d: rung %c on %s from it %d, cstart %.5e, move %g',S.stage,st.rung,st.src,S.it,S.cstart,st.move);
        fprintf('  %s\n',S.stage_log{S.stage});
    end
    if S.cstart==1, Fobj=Fbase; else, Fobj=struct('factor',Fbase.factor*sqrt(S.cstart)); end
    nstage=st.n-S.stage_it;
    nchunk=chunk-mod(S.it-S.it0,chunk);
    n=min([nstage, nchunk, target-S.it]);
    cfg=struct('move',st.move,'switch',st.switch,'stopfile',stopfile,'rung',st.rung,'src',st.src);
    tc=tic; [S,reason]=r0_long_core(ctx,S,Fobj,ops,G.opt,n,cfg); tchunk=toc(tc);
    stage_end = S.stage_it>=st.n || any(strcmp(reason,{'converged','switch'}));
    if stage_end
        S.stage_status{S.stage}=reason;
        if strcmp(reason,'switch'), S.stage_log{S.stage}=[S.stage_log{S.stage} ', switched: ' S.swreason]; end
        if S.stage==numel(plan)
            if strcmp(reason,'converged'), S.done=true; S.status='converged'; end
        else
            S.stage=S.stage+1; S.stage_it=0; S.stage_ready=false;
        end
    end
    if S.it>=target && ~S.done, S.status='budget'; end
    if ~S.done && S.it<target, S.status='running'; end
    write_records(G,S,tag,kind,cell,name,ckfile,tchunk,reason,st);
    if strcmp(reason,'stopfile'), fprintf('STOP file seen at it %d\n',S.it); break; end
end
fprintf('R0_LONG_%s_DONE it=%d status=%s walltime=%.0fs warn=%d\n',tag,S.it,S.status,S.walltime,S.nwarn);
end

function S=new_state(ctx,G)
S=struct('t',[],'st',[],'p0',[],'git',0,'nhold',0,'f0seg',[],'cstart',1,'Jbest_seg',inf,'tbest_seg',[], ...
    'it',0,'it0',0,'stage',1,'stage_it',0,'stage_ready',false,'done',false,'status','running', ...
    'walltime',0,'nwarn',0,'swcount',0,'swreason','','Jref',G.Jref,'Jxbest',inf,'txbest',[],'itxbest',0, ...
    'uniform_start',true,'source','');
S.stage_log={}; S.stage_status={};
S.H=struct('it',[],'stage',[],'rung',char([]),'Jown',[],'Jx',[],'Js',[],'Je',[],'vol',[],'chg',[], ...
    'inner',[],'cos',[],'f0',[],'time',[]);
S.T=zeros(ctx.ne,0); S.Tit=[];
end

function plan=strategy_plan(name,cell,budget)
% Multi-step strategies (fixed before any part-B run; nothing tuned afterwards).
% Stage lengths are outer iterations; the last stage runs to the budget.
mk=@(r,src,n,mv,sw) struct('rung',r,'src',src,'n',n,'move',mv,'switch',sw);
sw=struct('cosmin',0.9,'ncons',3,'lag',10,'maxits',round(budget/2));
switch name
    case 'P'    % plain exact from uniform, uninterrupted
        plan=mk('X',cell,inf,0,[]);
    case 'S1'   % staged models ERP -> sinc -> exact
        plan=[mk('E',cell,round(budget/8),0,[]), mk('S',cell,round(budget*7/24),0,[]), mk('X',cell,inf,0,[])];
    case 'S2'   % sinc descent, then exact polish
        plan=[mk('S',cell,round(budget/4),0,[]), mk('X',cell,inf,0,[])];
    case 'S3'   % sinc descent, then exact with a trust region (move limit 0.005 = 1/15 of the box)
        plan=[mk('S',cell,round(budget/4),0,[]), mk('X',cell,inf,0.005,[])];
    case 'S4'   % continuation in kR: sinc at the anchor kR (where sinc works), then exact down in kR
        switch cell
            case 'c30_k3', plan=[mk('S','c30_k18',round(budget/4),0,[]), mk('X',cell,inf,0,[])];
            case 'c30_k1', plan=[mk('S','c30_k18',round(budget/6),0,[]), mk('X','c30_k3',round(budget/4),0,[]), mk('X',cell,inf,0,[])];
            otherwise, error('S4 defined for c30 cells only');
        end
    case 'S5'   % sinc while its projected gradient agrees with exact and exact power keeps falling
        plan=[mk('S',cell,inf,0,sw), mk('X',cell,inf,0,[])];
    otherwise, error('unknown strategy %s',name);
end
end

function [F,srcops]=stage_factor(G,cell,ops,st,srcops)
if strcmp(st.src,cell), F=ops.(st.rung); return; end
key=['c_' st.src];
if ~isfield(srcops,key), srcops.(key)=load_ops(G,st.src); end
F=srcops.(key).(st.rung);
end

function ops=load_ops(G,cell)
% normalised factors, built exactly as in r0_pilot_runs (J(uniform) = Jref for every rung)
O=load(fullfile(G.pilot,['ops_' cell '.mat']),'cn','Fs_c','Fe_c');
F=load(fullfile(G.pilot,['factor_' cell '.mat']),'Fx');
ops.X=struct('factor',F.Fx*sqrt(O.cn(3)));
ops.S=struct('factor',O.Fs_c*sqrt(O.cn(2)));
ops.E=struct('factor',O.Fe_c*sqrt(O.cn(1)));
ops.cn=O.cn;
end

function write_records(G,S,tag,kind,cell,name,ckfile,tchunk,reason,st)
H=S.H; db=@(J) 10*log10(J/G.Jref);
n=numel(H.Jx);
if n==0, return; end
if strcmp(kind,'A')
    [~,ib]=min(H.Jown);
else
    [~,ib]=min(H.Jx);
end
[~,ix]=min(H.Jx);
c0=find(H.it>S.it-600,1); if isempty(c0), c0=1; end
gain=db(min(H.Jown(1:max(c0-1,1))))-db(min(H.Jown));
if c0==1, gain=NaN; end
fn=fullfile(G.out,['chunks_' kind '.csv']);
if ~exist(fn,'file')
    f=fopen(fn,'w');
    fprintf(f,'job,kind,cell,name,stage,rung,src,reason,status,iter_total,iter_new,stage_iter,vol_cur,own_dB_cur,exact_dB_cur,sinc_dB_cur,erp_dB_cur,pick_iter,pick_own_dB,pick_exact_dB,pick_sinc_dB,pick_vol,best_exact_dB_path,best_exact_iter,own_gain_last600_dB,walltime_cum_s,chunk_walltime_s,inner_mean_last600,singular_warn_cum,timestamp\n');
    fclose(f);
end
f=fopen(fn,'a');
fprintf(f,'%s,%s,%s,%s,%d,%c,%s,%s,%s,%d,%d,%d,%.6f,%.4f,%.4f,%.4f,%.4f,%d,%.4f,%.4f,%.4f,%.6f,%.4f,%d,%.4f,%.1f,%.1f,%.3f,%d,%s\n', ...
    tag,kind,cell,name,S.stage,st.rung,st.src,reason,S.status,S.it,S.it-S.it0,S.stage_it,H.vol(n),db(H.Jown(n)),db(H.Jx(n)),db(H.Js(n)),db(H.Je(n)), ...
    H.it(ib),db(H.Jown(ib)),db(H.Jx(ib)),db(H.Js(ib)),H.vol(ib),db(H.Jx(ix)),H.it(ix),gain,S.walltime,tchunk,mean(H.inner(c0:end)),S.nwarn, ...
    datestr(now,'yyyy-mm-ddTHH:MM:SS'));
fclose(f);
f=fopen(fullfile(G.out,'chunks.log'),'a');
fprintf(f,'%s %s it=%d stage=%d(%c@%s) reason=%s status=%s exact=%.3f sinc=%.3f own=%.3f dB (pick: exact %.3f, sinc %.3f) best_exact=%.3f vol=%.5f chunk=%.0fs warn=%d\n', ...
    datestr(now,'yyyy-mm-ddTHH:MM:SS'),tag,S.it,S.stage,st.rung,st.src,reason,S.status,db(H.Jx(n)),db(H.Js(n)),db(H.Jown(n)), ...
    db(H.Jx(ib)),db(H.Js(ib)),db(H.Jx(ix)),H.vol(ib),tchunk,S.nwarn);
fclose(f);
tmp=[ckfile(1:end-4) '_tmp.mat'];
save(tmp,'S','-v7'); movefile(tmp,ckfile,'f');
hf=fullfile(G.out,['hist_' tag '.csv']); f=fopen(hf,'w');
fprintf(f,'it,stage,rung,own_dB,exact_dB,sinc_dB,erp_dB,vol,chg,inner,cos_sinc_exact,time_s\n');
for j=1:n
    fprintf(f,'%d,%d,%c,%.6f,%.6f,%.6f,%.6f,%.7f,%.4e,%d,%.6f,%.4f\n',H.it(j),H.stage(j),H.rung(j),db(H.Jown(j)),db(H.Jx(j)),db(H.Js(j)),db(H.Je(j)),H.vol(j),H.chg(j),H.inner(j),H.cos(j),H.time(j));
end
fclose(f);
end

function validate(G)
cell='c30_k3'; [th,~]=cell_def(cell);
ctx=pilot_ctx(th); ctx.cachekey=sprintf('th%d',th); G.opt.vweights=ones(ctx.ne,1);
ops=load_ops(G,cell);
rng(11); tr=0.025+0.075*rand(ctx.ne,1);
[J1,g1]=r0_pilot_eval(ctx,tr,ops.X); J0=r0_long_eval(ctx,tr,ops.X); [J2,g2]=r0_long_eval(ctx,tr,ops.X);
fprintf('eval: pilot vs long J diff %.3e (rel %.2e), J(cache hit) diff %.3e, grad max diff %.3e\n',J2-J1,abs(J2-J1)/J1,J0-J1,max(abs(g2-g1)));
[Js1,gs1]=r0_pilot_eval(ctx,tr,ops.S); [Js2,gs2]=r0_long_eval(ctx,tr,ops.S);
fprintf('eval sinc (cache hit on a second factor): J diff %.3e, grad max diff %.3e\n',Js2-Js1,max(abs(gs2-gs1)));
P=load(fullfile(G.pilot,'run_c30_k3_XU.mat'),'hist');
cfg=struct('move',0,'switch',[],'stopfile','','rung','X','src',cell);
S=new_state(ctx,G); S.t=0.05*ones(ctx.ne,1); S.cstart=1; S.Jbest_seg=inf;
tc=tic; [Sa,~]=r0_long_core(ctx,S,ops.X,ops,G.opt,30,cfg); ta=toc(tc);
d=abs(Sa.H.Jown-P.hist.J(1:30))./P.hist.J(1:30);
fprintf('core vs pilot run c30_k3 XU, its 1..30: max rel diff %.3e (it 1 %.1e, it 10 %.1e, it 30 %.1e); %.2f s/it here vs %.2f s/it pilot average\n', ...
    max(d),d(1),d(10),d(30),ta/30,650.5/600);
dT=max(max(abs(Sa.T-P.hist.T(:,[10 20 30]))));
fprintf('  designs at its 10/20/30: max abs diff %.3e\n',dT);
[Sb,~]=r0_long_core(ctx,S,ops.X,ops,G.opt,15,cfg);
tmp=fullfile(tempdir,'r0_long_validate.mat'); save(tmp,'Sb','-v7'); L=load(tmp); Sb=L.Sb; delete(tmp);
[Sb,~]=r0_long_core(ctx,Sb,ops.X,ops,G.opt,15,cfg);
fprintf('chunked 15+15 (save/load) vs straight 30: max |dJ| %.3e, max |dt| %.3e, isequal(st) %d\n', ...
    max(abs(Sb.H.Jown-Sa.H.Jown)),max(abs(Sb.t-Sa.t)),isequaln(Sb.st,Sa.st));
warning('off','MATLAB:nearlySingularMatrix'); lastwarn('',''); x=[1 1;1 1+1e-15]\[1;1]; [m,id]=lastwarn; %#ok<ASGLU>
fprintf('lastwarn with the warning display off: id "%s"\n',id); warning('on','MATLAB:nearlySingularMatrix');
fprintf('R0_LONG_validate_DONE\n');
end

function fn=runfile(cell,code,out,exp5)
map=struct('XU','X0','SU','S0','SX','SX','XS','XS');
if strcmp(cell,'c30_k18') && isfield(map,code)
    fn=fullfile(exp5,['run_' map.(code) '.mat']);
else
    fn=fullfile(out,['run_' cell '_' code '.mat']);
end
end

function [th,kR]=cell_def(cell)
switch cell
    case 'c30_k1',    th=30;  kR=1.0;
    case 'c30_k3',    th=30;  kR=3.0;
    case 'c90_k1p5',  th=90;  kR=1.5;
    case 'c90_k7',    th=90;  kR=7.0;
    case 'c120_k1p5', th=120; kR=1.5;
    case 'c120_k7',   th=120; kR=7.0;
    case 'c30_k18',   th=30;  kR=2*pi*100/343*10;
    otherwise, error('unknown cell %s',cell);
end
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
