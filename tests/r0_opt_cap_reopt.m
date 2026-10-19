function r0_opt_cap_reopt(stage)
% Re-optimise the 30-deg spherical-cap thickness design at 100 Hz with the radiated
% power evaluated by (a) the paper's spline-consistent sinc operator and (b) the EXACT
% spectral two-sided operator reduced to the structural coefficients
% (tests/r0_opt_exact_factor.py -> results/r0_revision/exact_opt/cap_exact_factor_q8.mat).
% Both operators enter articleB_r0_eval through the same 'consistent' interface
% (J = w^2 ||factor*U||^2) so the structural model, adjoint, GCMMA, bounds, volume
% constraint, dB scaling and stopping rule are identical to the paper's cap ladder
% (tests/articleB_r0_thesis_damped.m: tmin 0.025, tmax 0.10, volfrac 0.5, nmax 600).
%   stage 'check'  : operator sanity + adjoint FD check for the exact operator
%   stage 'S0','X0': sinc / exact from the uniform design t = 0.05
%   stage 'XS','SX': exact warm-started at the sinc optimum, sinc warm-started at the exact optimum
%   stage 'SR','XR': sinc / exact warm-started at the paper's quadvel-rad (legacy Nystrom) design
%   stage 'states' : structural states U of every design on every stored path (for r0_opt_cap_eval.py)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fileparts(root)); addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
out=fullfile(root,'results','r0_revision','exact_opt');
ctx=cap_ctx();
F=load(fullfile(out,'cap_exact_factor_q8.mat'),'Fx');
Rx=struct('factor',F.Fx);
Rs=r0_build_consistent(ctx,1.21,343,[4 4],32,64);
fz=load(fullfile(root,'results','r0_revision','spherecap_frozen_states.mat'),'designs','labels');
opt=struct('tmin',0.025,'tmax',0.10,'volfrac',0.5,'nmax',600,'stop',true,'ftol',1e-4, ...
    'xchgtol',1e-4,'nhold',3,'window',10,'feastol',1e-6,'Abar',100,'dbfac',10);
switch stage
  case 'check'
    for j=1:size(fz.designs,2)
        Js=articleB_r0_eval(ctx,fz.designs(:,j),'quadvel',Rs,'consistent');
        Jx=articleB_r0_eval(ctx,fz.designs(:,j),'quadvel',Rx,'consistent');
        fprintf('%-13s J_sinc(consistent q4/32x64)=%.6e  J_exact=%.6e  ratio %.4f\n',fz.labels{j},Js,Jx,Js/Jx);
    end
    t=0.05*ones(ctx.ne,1); rng(517); d=randn(ctx.ne,1); d=d/norm(d);
    for R0={Rx,Rs}
        [J,g]=articleB_r0_eval(ctx,t,'quadvel',R0{1},'consistent');
        for h=[1e-4 1e-5 1e-6]
            Jp=articleB_r0_eval(ctx,t+h*d,'quadvel',R0{1},'consistent');
            Jm=articleB_r0_eval(ctx,t-h*d,'quadvel',R0{1},'consistent');
            fd=(Jp-Jm)/(2*h); an=g'*d;
            fprintf('FD check J=%.4e h=%g fd=%+.9e adj=%+.9e rel=%.2e\n',J,h,fd,an,abs(fd-an)/abs(an));
        end
    end
  case {'S0','X0','XS','SX','SR','XR'}
    switch stage(1), case 'S', R0=Rs; tag='sinc'; case 'X', R0=Rx; tag='exact'; end
    switch stage(2)
        case '0', o=opt;
        case 'S', P=load(fullfile(out,'run_S0.mat'),'tbest'); o=opt; o.tinit=P.tbest;
        case 'X', P=load(fullfile(out,'run_X0.mat'),'tbest'); o=opt; o.tinit=P.tbest;
        case 'R', o=opt; o.tinit=fz.designs(:,strcmp(fz.labels,'quadvel-rad'));
    end
    fprintf('== run %s: operator %s, start %s ==\n',stage,tag,stage(2));
    [tbest,hist,meta]=r0_opt_one_track(ctx,'quadvel',R0,'consistent',o);
    fprintf('RESULT %s: it=%d stop=%s evals=%d Jbest=%.6e vol=%.6f time=%.0fs\n',stage,meta.niter, ...
        meta.stopReason,meta.numeval,meta.Jbest,meta.final_vol,meta.walltime);
    save(fullfile(out,['run_' stage '.mat']),'tbest','hist','meta','o','tag','-v7');
  case 'states'
    runs={'S0','X0','XS','SX','SR','XR'}; T=[]; run=[]; iter=[];
    for r=1:numel(runs)
        fn=fullfile(out,['run_' runs{r} '.mat']); if ~exist(fn,'file'), continue; end
        P=load(fn,'hist','tbest'); n=size(P.hist.T,2);
        T=[T P.hist.T P.tbest]; run=[run r*ones(1,n+1)]; iter=[iter 1:n 0]; %#ok<AGROW>  iter 0 = tbest
    end
    T=[T fz.designs]; run=[run zeros(1,size(fz.designs,2))]; iter=[iter -(1:size(fz.designs,2))];
    U=complex(zeros(ctx.ndof,size(T,2)));
    for j=1:size(T,2)
        t=T(:,j); KE=ctx.Bke.*(t.^3)+ctx.Ske.*t; ME=ctx.Me.*t;
        A=sparse(ctx.Rows(:),ctx.Cols(:),(1+1i*ctx.w*ctx.beta)*KE(:)-ctx.w^2*ME(:),ctx.ndof,ctx.ndof);
        U(ctx.free,j)=A(ctx.free,ctx.free)\ctx.F(ctx.free);
    end
    Psinc_q4=ctx.w^2*sum(abs(Rs.factor*U).^2,1); w0=ctx.w; frozen_labels=fz.labels; %#ok<NASGU>
    save(fullfile(out,'cap_path_states.mat'),'T','U','run','iter','runs','Psinc_q4','w0','frozen_labels','-v7');
    fprintf('states: %d designs saved\n',size(T,2));
end
fprintf('R0_OPT_CAP_REOPT_%s_DONE\n',stage);
end

function ctx=cap_ctx()
spec=struct('R',25,'L',50,'th0',deg2rad(40),'alpha',deg2rad(30), ...
    'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,'deg',[2 2],'reg',[1 1], ...
    'nsub',[40 24],'nquad',[4 4],'bc',[1 2 3 4],'loadtype','point', ...
    'loaduv',[0.5 0.5],'loaddir',[0;1;0],'loadmag',1,'freq',6.53547,'zeta',0.05);
spec.R=10; spec.bc=2; spec.loadtype='pole'; spec.loaddir=[0;0;1]; spec.freq=100;
ctx=obstacle_build_ctx('spherecap',spec); ctx.alpha=0; ctx.beta=0.1/ctx.w;
end
