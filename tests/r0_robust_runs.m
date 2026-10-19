function r0_robust_runs(job)
% R0_ROBUST_RUNS  Experiment 10: robust shell optimisation under random relative thickness errors.
% One job per MATLAB session; resumable from its checkpoint (same pattern as r0_long_runs).
%
%   r0_robust_runs('<cell> <rung> <start> <N> <sigma> <mode> <target> [<fsi>]')
%     cell   : regime-pilot cell (c30_k18 = the paper case, c30_k3, ...)
%     rung   : X (exact two-sided, Fx) | S (spline-consistent sinc)
%     start  : U (uniform t = 0.05) | X (nominal exact optimum: anchor Experiment-5 X0, else pilot XU)
%              | S (nominal sinc optimum: anchor S0, else pilot SU) | SX | XS (other stored runs)
%     N      : number of fixed samples (common random numbers: rng(4242), randn(ne,32), first N columns)
%     sigma  : relative thickness standard deviation (0.01 = 1 %)
%     mode   : mean | meanstd1 (mean + 1 std) | meanstd2
%     target : total outer iterations
%     fsi    : optional dev | phys: fluid loading in the state solve (iwZ from
%              tests/r0_robust_zexport.py in the directory given by env R0_ROBUST_ZDIR)
% Optimiser: the pilot's GCMMA, bounds [0.025, 0.10], volume fraction <= 0.5, thesis dB scaling,
% windowed stopping rule, through r0_robust_core.  Normalised factors (J(uniform) = Jref per rung)
% and the start rescale J_rob(start) = Jref (pilot warm-start rule) for every start.
% Records in results/r0_revision/robust/: ck_<tag>.mat (state + history), hist_<tag>.csv,
% chunks.csv and chunks.log (one line per 50-iteration chunk).  A file robust/STOP stops a running
% job cleanly at the next 25-iteration mark.
maxNumCompThreads(2); warning('off','backtrace');
tok=strsplit(strtrim(job));
cell=tok{1}; rung=tok{2}; start=tok{3}; N=str2double(tok{4}); sigma=str2double(tok{5});
mode=tok{6}; target=str2double(tok{7}); fsi=''; if numel(tok)>=8, fsi=tok{8}; end
[ctx,ops,G]=r0_robust_setup(cell);
tag=sprintf('%s_%s%s_N%d_s%g_%s',cell,rung,start,N,100*sigma,mode);
if ~isempty(fsi), tag=[tag '_fsi' fsi]; end
ckfile=fullfile(G.out,['ck_' tag '.mat']);
stopfile=fullfile(G.out,'STOP');
rng(4242); Xall=randn(ctx.ne,32); Xi=Xall(:,1:N);
R=struct('sigma',sigma,'Xi',Xi,'solver','gmres');
switch mode
    case 'mean', R.mode='mean'; R.k=0;
    case 'meanstd1', R.mode='meanstd'; R.k=1;
    case 'meanstd2', R.mode='meanstd'; R.k=2;
    otherwise, error('unknown mode %s',mode);
end
Fbase=ops.(rung);
if ~isempty(fsi)
    Z=load(fullfile(getenv('R0_ROBUST_ZDIR'),['iwz_' cell '_' fsi '.mat']),'iwZ');
    Fbase.iwZ=Z.iwZ; clear Z
end
if exist(ckfile,'file')
    L=load(ckfile,'S'); S=L.S; clear L
    assert(isequal(S.Xi,Xi),'checkpoint samples differ');
    fprintf('== resume %s from checkpoint: it %d ==\n',tag,S.it);
else
    S=new_state(ctx,G);
    if strcmp(start,'U'), S.t=0.05*ones(ctx.ne,1); S.source='uniform';
    else
        fn=runfile(cell,start,G); P=load(fn,'tbest'); S.t=P.tbest(:); S.source=fn;
    end
    S.t0=S.t; S.Xi=Xi; S.job=job;
    Rb=R; S.cstart=G.Jref/r0_robust_eval(ctx,S.t,Fbase,Rb);
    fprintf('== new %s (theta0 %d, kR %.4f): start %s, N %d, sigma %g, mode %s, cstart %.6e, target %d ==\n', ...
        tag,G.th,G.kR,S.source,N,sigma,mode,S.cstart,target);
end
Fobj=Fbase; Fobj.factor=Fbase.factor*sqrt(S.cstart);
chunk=50;
while ~S.done && S.it<target
    n=min(chunk-mod(S.it,chunk),target-S.it);
    cfg=struct('move',0,'stopfile',stopfile,'rung',rung,'src',cell,'R',R);
    tc=tic; [S,reason]=r0_robust_core(ctx,S,Fobj,ops,G.opt,n,cfg); tchunk=toc(tc);
    if strcmp(reason,'converged'), S.done=true; S.status='converged';
    elseif S.it>=target, S.status='budget';
    else, S.status='running'; end
    write_records(G,S,tag,ckfile,tchunk,reason);
    if strcmp(reason,'stopfile'), fprintf('STOP file seen at it %d\n',S.it); break; end
end
fprintf('R0_ROBUST_%s_DONE it=%d status=%s walltime=%.0fs warn=%d\n',tag,S.it,S.status,S.walltime,S.nwarn);
end

function S=new_state(ctx,G)
S=struct('t',[],'st',[],'p0',[],'git',0,'nhold',0,'f0seg',[],'cstart',1,'it',0,'done',false, ...
    'status','running','walltime',0,'nwarn',0,'Jref',G.Jref,'Jrobbest',inf,'trobbest',[],'itrobbest',0,'source','');
S.H=struct('it',[],'Jrob',[],'Jstd',[],'Jx',[],'Js',[],'vol',[],'chg',[],'inner',[],'f0',[],'time',[],'gmres',[],'nfb',[]);
S.T=zeros(ctx.ne,0); S.Tit=[];
end

function fn=runfile(cell,code,G)
map=struct('X','X0','S','S0','SX','SX','XS','XS');
if strcmp(cell,'c30_k18')
    fn=fullfile(G.exp5,['run_' map.(code) '.mat']);
else
    map=struct('X','XU','S','SU','SX','SX','XS','XS');
    fn=fullfile(G.pilot,['run_' cell '_' map.(code) '.mat']);
end
end

function write_records(G,S,tag,ckfile,tchunk,reason)
H=S.H; db=@(J) 10*log10(J/G.Jref); n=numel(H.Jrob);
if n==0, return; end
c0=max(1,n-49);
fn=fullfile(G.out,'chunks.csv');
if ~exist(fn,'file')
    f=fopen(fn,'w');
    fprintf(f,'job,reason,status,iter,rob_dB_cur,rob_std_rel_cur,nom_exact_dB_cur,nom_sinc_dB_cur,vol_cur,best_rob_dB,best_rob_iter,walltime_cum_s,chunk_walltime_s,s_per_it_chunk,inner_mean_chunk,gmres_it_mean_chunk,fallbacks_chunk,singular_warn_cum,timestamp\n');
    fclose(f);
end
f=fopen(fn,'a');
fprintf(f,'%s,%s,%s,%d,%.4f,%.4e,%.4f,%.4f,%.6f,%.4f,%d,%.1f,%.1f,%.3f,%.3f,%.2f,%d,%d,%s\n',tag,reason,S.status,S.it, ...
    db(H.Jrob(n)),H.Jstd(n)/H.Jrob(n),db(H.Jx(n)),db(H.Js(n)),H.vol(n),db(S.Jrobbest),S.itrobbest,S.walltime,tchunk, ...
    mean(H.time(c0:n)),mean(H.inner(c0:n)),mean(H.gmres(c0:n)),sum(H.nfb(c0:n)),S.nwarn,datestr(now,'yyyy-mm-ddTHH:MM:SS'));
fclose(f);
f=fopen(fullfile(G.out,'chunks.log'),'a');
fprintf(f,'%s %s it=%d reason=%s status=%s robust=%.3f dB (std/mean %.3f) nominal exact=%.3f sinc=%.3f vol=%.5f best_robust=%.3f@%d chunk=%.0fs (%.2f s/it) warn=%d\n', ...
    datestr(now,'yyyy-mm-ddTHH:MM:SS'),tag,S.it,reason,S.status,db(H.Jrob(n)),H.Jstd(n)/H.Jrob(n),db(H.Jx(n)),db(H.Js(n)),H.vol(n), ...
    db(S.Jrobbest),S.itrobbest,tchunk,mean(H.time(c0:n)),S.nwarn);
fclose(f);
tmp=[ckfile(1:end-4) '_tmp.mat'];
save(tmp,'S','-v7'); movefile(tmp,ckfile,'f');
hf=fullfile(G.out,['hist_' tag '.csv']); f=fopen(hf,'w');
fprintf(f,'it,rob_dB,rob_std_over_mean,nom_exact_dB,nom_sinc_dB,vol,chg,inner,gmres_it,time_s\n');
for j=1:n
    fprintf(f,'%d,%.6f,%.5e,%.6f,%.6f,%.7f,%.4e,%d,%.2f,%.4f\n',H.it(j),db(H.Jrob(j)),H.Jstd(j)/H.Jrob(j),db(H.Jx(j)),db(H.Js(j)),H.vol(j),H.chg(j),H.inner(j),H.gmres(j),H.time(j));
end
fclose(f);
end
