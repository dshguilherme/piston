function r0_pilot_runs(cell, stage, nmax)
% R0 regime pilot: re-optimise the cap thickness under each radiation rung, then evaluate
% every stored design under every rung.  The structural model, adjoint (articleB_r0_eval),
% GCMMA, bounds, volume constraint, dB scaling and stopping rule are those of the paper's
% cap ladder, used through r0_pilot_track (= r0_opt_one_track with the computationally cheaper but
% numerically equivalent r0_pilot_eval; equivalence checked in stage 'ops' and 'replicate').
%
% Rungs (all through the 'consistent' interface J = w^2 ||factor U||^2):
%   E  R0  ERP, P = (rho c / 2) int |v_n|^2 dS on the acoustic q8 samples (phase-blind)
%   S  R1  the paper's spline-consistent sinc (r0_build_consistent, q4, 32 x 64 directions)
%   X  R4  exact spectral two-sided cap operator (tests/r0_pilot_factor.py, Fx)
%   R2 = 2 x R1 is not run: after the per-rung normalisation below it is the same objective.
% Normalisation: factor_r <- factor_r * sqrt(Jref / J_r(uniform)), Jref = 5.9464e-7 W (exact
% uniform power of the paper cap at 100 Hz, so the paper's dB scaling, Abar = 100, stays in
% the range it was used in).  Argmin and predicted dB reductions are unaffected; every rung
% starts from the same scaled objective value (f0 = 100 at uniform).
% Acoustic medium: rho' c' = 1.21 * 343 and k = kR / R (c' = w R / kR), as in r0_pilot_factor.py.
%
%   stage 'ops'    : build/cache the rung factors and normalisation, operator and adjoint checks
%   stage 'XU','SU','EU' : rung X/S/E from the uniform design t = 0.05
%   stage 'SX','EX' : sinc / ERP warm-started at the exact optimum (XU tbest)
%   stage 'XS'     : exact warm-started at the sinc optimum (SU tbest)
%   stage 'SUn'    : (anchor only) normalised sinc from uniform; the Experiment-5 sinc runs are unnormalised
%   stage 'replicate' : (anchor) exact rung from uniform for 40 its with r0_pilot_track vs Experiment-5 X0
%   stage 'eval'   : states + powers of uniform, 40 random smooth designs, every 5th path design and
%                    every run optimum under every rung -> pilot/eval_<cell>.csv, runs_<cell>.csv
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fileparts(root)); addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
out=fullfile(root,'results','r0_revision','pilot');
exp5=fullfile(root,'results','r0_revision','exact_opt');
[th,kR]=cell_def(cell);
ctx=pilot_ctx(th); R=10; rhoc=1.21*343; Jref=5.9464e-7;
cp=ctx.w*R/kR;                               % effective sound speed, rho' = rhoc / cp
opsfile=fullfile(out,['ops_' cell '.mat']);
if nargin<3||isempty(nmax), nmax=600; end
opt=struct('tmin',0.025,'tmax',0.10,'volfrac',0.5,'nmax',nmax,'stop',true,'ftol',1e-4, ...
    'xchgtol',1e-4,'nhold',3,'window',10,'feastol',1e-6,'Abar',100,'dbfac',10);
if strcmp(stage,'ops')
    S=load(fullfile(out,sprintf('samples_th%d.mat',th)),'w','B','Uu');
    F=load(fullfile(out,['factor_' cell '.mat']),'Fx');
    Rs=r0_build_consistent(ctx,rhoc/cp,cp,[4 4],32,64);
    Fe=spdiags(sqrt(0.5*rhoc*S.w(:)),0,numel(S.w),numel(S.w))*S.B;
    t=0.05*ones(ctx.ne,1);
    Ju=[articleB_r0_eval(ctx,t,'quadvel',struct('factor',Fe),'consistent'), ...
        articleB_r0_eval(ctx,t,'quadvel',Rs,'consistent'), ...
        articleB_r0_eval(ctx,t,'quadvel',struct('factor',F.Fx),'consistent')];
    cn=Jref./Ju;                                 % [E S X]
    Fs_c=Rs.factor; Fe_c=Fe; Fx_c=F.Fx;           %#ok<NASGU>
    fprintf('%s: theta0 %d kR %.4f c'' %.3f m/s; uniform J: ERP %.5e sinc %.5e exact %.5e; sinc/exact %.4f, sigma_exact %.4e\n',...
        cell,th,kR,cp,Ju(1),Ju(2),Ju(3),Ju(2)/Ju(3),Ju(3)/Ju(1));
    fprintf('  normalisation c = Jref/J(uniform): E %.5e  S %.5e  X %.5e\n',cn);
    % uniform state check against the stored samples state
    [~,~,aux]=articleB_r0_eval(ctx,t,'quadvel',struct('factor',Fe),'consistent');
    fprintf('  uniform state vs samples file: rel diff %.2e\n',norm(aux.u-S.Uu)/norm(S.Uu));
    % adjoint check of each normalised rung
    rng(517); d=randn(ctx.ne,1); d=d/norm(d); h=1e-5;
    facs={Fe*sqrt(cn(1)),Rs.factor*sqrt(cn(2)),F.Fx*sqrt(cn(3))}; nm='ESX';
    for r=1:3
        R0=struct('factor',facs{r});
        [J,g]=articleB_r0_eval(ctx,t,'quadvel',R0,'consistent');
        Jp=articleB_r0_eval(ctx,t+h*d,'quadvel',R0,'consistent'); Jm=articleB_r0_eval(ctx,t-h*d,'quadvel',R0,'consistent');
        fd=(Jp-Jm)/(2*h); an=g'*d;
        fprintf('  FD check rung %s: J=%.5e fd=%+.9e adj=%+.9e rel=%.2e\n',nm(r),J,fd,an,abs(fd-an)/abs(an));
    end
    % r0_pilot_eval vs articleB_r0_eval on a random design, every rung
    tr=0.025+0.075*rand(ctx.ne,1);
    for r=1:3
        R0=struct('factor',facs{r});
        [J1,g1]=articleB_r0_eval(ctx,tr,'quadvel',R0,'consistent'); [J2,g2]=r0_pilot_eval(ctx,tr,R0); J3=r0_pilot_eval(ctx,tr,R0);
        fprintf('  r0_pilot_eval vs articleB_r0_eval rung %s: J rel %.2e, J(no grad) rel %.2e, grad rel %.2e\n',nm(r),abs(J2-J1)/J1,abs(J3-J1)/J1,norm(g2-g1)/norm(g1));
    end
    save(opsfile,'Fs_c','Fe_c','cn','Ju','Jref','cp','th','kR','-v7');
    fprintf('R0_PILOT_RUNS_%s_ops_DONE\n',cell); return
end
O=load(opsfile,'cn');
switch stage
  case {'XU','SU','EU','SX','EX','XS','SUn'}
    switch stage(1)
        case 'X', F=load(fullfile(out,['factor_' cell '.mat']),'Fx'); R0=struct('factor',F.Fx*sqrt(O.cn(3))); clear F
        case 'S', P=load(opsfile,'Fs_c'); R0=struct('factor',P.Fs_c*sqrt(O.cn(2))); clear P
        case 'E', P=load(opsfile,'Fe_c'); R0=struct('factor',P.Fe_c*sqrt(O.cn(1))); clear P
    end
    o=opt;
    switch stage(2)
        case 'X', P=load(runfile(cell,'XU',out,exp5),'tbest'); o.tinit=P.tbest;
        case 'S', P=load(runfile(cell,'SU',out,exp5),'tbest'); o.tinit=P.tbest;
    end
    % Warm starts: rescale so that J(start) = Jref as well. The paper's scaling
    % f0 = 100 (Abar + 10 log10 J) / p0 needs p0 = Abar + 10 log10 J(start) > 0, i.e. J(start) > 1e-10;
    % warm starts at deep optima violate this (p0 < 0 flips the sign of the objective).
    cstart=1;
    if isfield(o,'tinit')
        J0=r0_pilot_eval(ctx,o.tinit,R0); cstart=Jref/J0; R0.factor=R0.factor*sqrt(cstart);
    end
    fprintf('== pilot run %s %s (theta0 %d, kR %.4f), normalisation %.5e, start rescale %.5e ==\n',cell,stage,th,kR,O.cn('ESX'==stage(1)),cstart);
    [tbest,hist,meta]=r0_pilot_track(ctx,'quadvel',R0,'consistent',o);
    fprintf('RESULT %s %s: it=%d stop=%s evals=%d Jbest=%.6e vol=%.6f time=%.0fs\n',cell,stage,meta.niter, ...
        meta.stopReason,meta.numeval,meta.Jbest,meta.final_vol,meta.walltime);
    save(fullfile(out,['run_' cell '_' stage '.mat']),'tbest','hist','meta','o','cstart','-v7');
  case 'grad'
    % gradient test (course correction 3): at uniform and at the run end points, compare the
    % feasible-projected gradients of log P (ERP, sinc vs exact), and take a short feasible
    % descent step generated by each rung; report the change of each rung's own and the exact power.
    F=load(fullfile(out,['factor_' cell '.mat']),'Fx'); P=load(opsfile,'Fs_c','Fe_c');
    ops={struct('factor',P.Fe_c),struct('factor',P.Fs_c),struct('factor',F.Fx)}; nm='ESX'; clear F P
    pts={'uniform','XU','SU','EU'}; fg=fopen(fullfile(out,['grad_' cell '.csv']),'w');
    fprintf(fg,'cell,point,rung,cos_angle_vs_exact,angle_deg,relerr_vs_exact,step_maxdt,own_dB_change,exact_dB_change,own_firstorder_dB,exact_firstorder_dB\n');
    for ip=1:numel(pts)
        if ip==1, t=0.05*ones(ctx.ne,1);
        else, fn=runfile(cell,pts{ip},out,exp5); if ~exist(fn,'file'), continue; end; Q=load(fn,'tbest'); t=Q.tbest; end
        fr=(t>opt.tmin+1e-6)&(t<opt.tmax-1e-6);
        J=zeros(1,3); G=zeros(ctx.ne,3);
        for r=1:3
            [J(r),g]=r0_pilot_eval(ctx,t,ops{r}); g=g/J(r); g(~fr)=0; g(fr)=g(fr)-mean(g(fr)); G(:,r)=g;
        end
        for r=1:3
            c=G(:,r)'*G(:,3)/(norm(G(:,r))*norm(G(:,3))); re=norm(G(:,r)-G(:,3))/norm(G(:,3));
            d=-G(:,r)/max(abs(G(:,r)));
            for h=[5e-4 2e-3]
                tn=t+h*d; tn=min(max(tn,opt.tmin),opt.tmax); tn(fr)=tn(fr)-mean(tn(fr)-t(fr));
                tn=min(max(tn,opt.tmin),opt.tmax);
                Jn=[r0_pilot_eval(ctx,tn,ops{r}), r0_pilot_eval(ctx,tn,ops{3})];
                dt=tn-t; fo=10/log(10)*[G(:,r)'*dt, G(:,3)'*dt];
                fprintf(fg,'%s,%s,%s,%.6f,%.3f,%.4e,%.1e,%.6e,%.6e,%.6e,%.6e\n',cell,pts{ip},nm(r),c,acosd(min(1,max(-1,c))),re,h, ...
                    10*log10(Jn(1)/J(r)),10*log10(Jn(2)/J(3)),fo(1),fo(2));
            end
        end
        fprintf('%s grad at %s: cos(S,X)=%.4f cos(E,X)=%.4f\n',cell,pts{ip},G(:,2)'*G(:,3)/norm(G(:,2))/norm(G(:,3)),G(:,1)'*G(:,3)/norm(G(:,1))/norm(G(:,3)));
    end
    fclose(fg);
  case 'replicate'
    assert(strcmp(cell,'c30_k18'));
    F=load(fullfile(out,['factor_' cell '.mat']),'Fx'); R0=struct('factor',F.Fx*sqrt(O.cn(3))); clear F
    o=opt; o.nmax=40; o.stop=false;
    [~,hist,meta]=r0_pilot_track(ctx,'quadvel',R0,'consistent',o);
    X=load(fullfile(exp5,'run_X0.mat'),'hist');
    c=O.cn(3); n=numel(hist.J);
    for it=[1 2 5 10 20 30 40], if it<=n, fprintf('  it %2d: pilot J/c %.9e  exp5 X0 J %.9e  rel %.2e\n',it,hist.J(it)/c,X.hist.J(it),abs(hist.J(it)/c-X.hist.J(it))/X.hist.J(it)); end, end
    fprintf('  replicate: %d its in %.0f s (%.2f s/it, %d inner evals); normalisation c_X = %.8f\n',n,meta.walltime,meta.walltime/n,meta.numeval,c);
  case 'eval'
    S=load(fullfile(out,sprintf('samples_th%d.mat',th)),'w','B','Trand');
    F=load(fullfile(out,['factor_' cell '.mat']),'Fx','Fx8000','Fx32','Fs');
    P=load(opsfile,'Fs_c','Fe_c');
    runs={'XU','SU','EU','SX','EX','XS','SUn'};
    T=0.05*ones(ctx.ne,1); dset={'uniform'}; run={''}; iter=0;
    T=[T S.Trand]; dset=[dset repmat({'random'},1,size(S.Trand,2))]; run=[run repmat({''},1,size(S.Trand,2))]; iter=[iter 1:size(S.Trand,2)];
    fr=fopen(fullfile(out,['runs_' cell '.csv']),'w');
    fprintf(fr,'cell,run,rung,start,niter,stop,numeval,walltime_s,final_vol,Jbest_own_normalised,tail100_dB,tail50_dB,source\n');
    for r=1:numel(runs)
        fn=runfile(cell,runs{r},out,exp5); if ~exist(fn,'file'), continue; end
        Q=load(fn,'tbest','hist','meta'); n=size(Q.hist.T,2);
        sel=5:5:n; T=[T Q.hist.T(:,sel) Q.tbest]; %#ok<AGROW>
        dset=[dset repmat({'path'},1,numel(sel)) {'optimum'}]; run=[run repmat(runs(r),1,numel(sel)+1)]; iter=[iter sel 0]; %#ok<AGROW>
        jb=cummin(Q.hist.J(:)); t100=NaN; t50=NaN;
        if n>101, t100=10*log10(jb(end-100)/jb(end)); end
        if n>51, t50=10*log10(jb(end-50)/jb(end)); end
        src='pilot'; if contains(fn,'exact_opt'), src='exp5'; end
        fprintf(fr,'%s,%s,%s,%s,%d,%s,%d,%.1f,%.6f,%.6e,%.5f,%.5f,%s\n',cell,runs{r},runs{r}(1),runs{r}(2:end),n, ...
            Q.meta.stopReason,Q.meta.numeval,Q.meta.walltime,Q.meta.final_vol,Q.meta.Jbest,t100,t50,src);
    end
    fclose(fr);
    nd=size(T,2); fprintf('%s eval: %d designs\n',cell,nd);
    Pe=zeros(nd,1); Pi=Pe; Ps=Pe; Pss=Pe; Px=Pe; Px8=Pe; Px32=Pe; vol=Pe;
    w2=ctx.w^2;
    for j=1:nd
        t=T(:,j); KE=ctx.Bke.*(t.^3)+ctx.Ske.*t; ME=ctx.Me.*t;
        A=sparse(ctx.Rows(:),ctx.Cols(:),(1+1i*ctx.w*ctx.beta)*KE(:)-ctx.w^2*ME(:),ctx.ndof,ctx.ndof);
        u=complex(zeros(ctx.ndof,1)); u(ctx.free)=A(ctx.free,ctx.free)\ctx.F(ctx.free);
        Pe(j)=w2*sum(abs(P.Fe_c*u).^2); Pi(j)=w2*sum(abs(u).^2); Ps(j)=w2*sum(abs(P.Fs_c*u).^2);
        us=single(u); Pss(j)=w2*sum(abs(double(F.Fs*us)).^2); Px(j)=w2*sum(abs(F.Fx*u).^2);
        Px8(j)=w2*sum(abs(double(F.Fx8000*us)).^2); Px32(j)=w2*sum(abs(double(F.Fx32*us)).^2);
        vol(j)=mean(t)/0.10;
    end
    fe=fopen(fullfile(out,['eval_' cell '.csv']),'w');
    fprintf(fe,'cell,theta0_deg,kR,set,run,iter,vol,P_erp,P_I,P_sinc,P_sinc_spec,P_exact,P_exact_L8000,P_exact_N32\n');
    for j=1:nd
        fprintf(fe,'%s,%d,%.6f,%s,%s,%d,%.6f,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n',cell,th,kR,dset{j},run{j},iter(j),vol(j), ...
            Pe(j),Pi(j),Ps(j),Pss(j),Px(j),Px8(j),Px32(j));
    end
    fclose(fe);
    save(fullfile(out,['eval_' cell '_designs.mat']),'T','dset','run','iter','-v7');
end
fprintf('R0_PILOT_RUNS_%s_%s_DONE\n',cell,stage);
end

function fn=runfile(cell,code,out,exp5)
% anchor cell: reuse the Experiment-5 runs (same structure, same exact factor, unnormalised sinc)
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
