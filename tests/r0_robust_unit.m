function r0_robust_unit(cell, zdir)
% R0_ROBUST_UNIT  Checks of the robust objective r0_robust_eval (Experiment 10) on one cell.
%  (1) sigma = 0 reproduces the nominal objective and gradient of r0_pilot_eval (exact and sinc)
%  (2) sigma -> 0 (1e-6, N = 16): J_rob and its gradient converge to the nominal ones
%  (3) central finite-difference check of dJ_rob along a random unit direction, sigma = 1%, N = 8,
%      modes 'mean' and 'meanstd' (k = 1), exact and sinc, at the uniform design and the exact optimum
%  (4) GMRES (nominal-LU preconditioner) vs one direct LU per sample, sigma = 1% and 3%
%  (5) fluid loading: sigma = 0 vs tests/r0_fsi_eval (dense LU), GMRES vs dense direct, FD check
%  (6) cost: seconds per evaluation (J only / J + gradient) vs N, GMRES iterations, memory
% Writes results/r0_revision/robust/unit_<cell>.log and unit_<cell>.csv (timing rows).
maxNumCompThreads(2); warning('off','backtrace');
[ctx,ops,G]=r0_robust_setup(cell);
if nargin<2, zdir=''; end
lf=fopen(fullfile(G.out,['unit_' cell '.log']),'a');
lg=@(varargin) logboth(lf,varargin{:});
lg('== r0_robust_unit %s %s (theta0 %d, kR %.4f), threads %d ==\n',cell,datestr(now,'yyyy-mm-ddTHH:MM:SS'),G.th,G.kR,maxNumCompThreads);
ne=ctx.ne; tu=0.05*ones(ne,1);
switch cell
    case 'c30_k18', P=load(fullfile(G.exp5,'run_X0.mat'),'tbest');
    otherwise,      P=load(fullfile(G.pilot,['run_' cell '_XU.mat']),'tbest');
end
tx=P.tbest(:);
pts={tu,tx}; pnm={'uniform','exactopt'};
opn={'X','S'};
% (1) sigma = 0
for ip=1:2
    for r=1:2
        F=ops.(opn{r});
        [J1,g1]=r0_pilot_eval(ctx,pts{ip},F);
        R=struct('sigma',0,'Xi',[],'solver','direct');
        [J2,g2]=r0_robust_eval(ctx,pts{ip},F,R);
        lg('(1) sigma=0 %-8s %s: J %.6e, |dJ|/J %.3e, grad rel %.3e\n',pnm{ip},opn{r},J1,abs(J2-J1)/J1,norm(g2-g1)/norm(g1));
    end
end
% (2) sigma -> 0
rng(1001); Xi16=randn(ne,16);
for ip=1:2
    F=ops.X; [J1,g1]=r0_pilot_eval(ctx,pts{ip},F);
    for sg=[1e-2 1e-3 1e-4 1e-6]
        R=struct('sigma',sg,'Xi',Xi16,'solver','gmres');
        [J2,g2,inf2]=r0_robust_eval(ctx,pts{ip},F,R);
        lg('(2) %-8s X sigma=%.0e N=16: J_rob/J_nom-1 %+.3e, grad rel %.3e, cos %.8f, gmres its %d-%d\n',pnm{ip},sg,J2/J1-1, ...
            norm(g2-g1)/norm(g1),g1'*g2/norm(g1)/norm(g2),min(inf2.gmres_it),max(inf2.gmres_it));
    end
end
% (3) FD check
rng(2002); d=randn(ne,1); d=d/norm(d); Xi8=randn(ne,8);
for ip=1:2
    R=struct('sigma',0,'Xi',[],'solver','direct'); [J,g]=r0_robust_eval(ctx,pts{ip},ops.X,R);
    for h=[1e-3 1e-4 1e-5]
        fd=(r0_robust_eval(ctx,pts{ip}+h*d,ops.X,R)-r0_robust_eval(ctx,pts{ip}-h*d,ops.X,R))/(2*h);
        lg('(3) FD %-8s X nominal (sigma=0) h=%.0e: J %.6e fd %+.9e adj %+.9e rel %.2e\n',pnm{ip},h,J,fd,g'*d,abs(fd-g'*d)/abs(g'*d));
    end
    for r=1:2
        for md={'mean','meanstd'}
          for solver={'direct','gmres'}
            R=struct('sigma',0.01,'Xi',Xi8,'mode',md{1},'k',1,'solver',solver{1});
            F=ops.(opn{r});
            [J,g]=r0_robust_eval(ctx,pts{ip},F,R);
            hs=[1e-3 1e-4 1e-5]; if strcmp(solver{1},'gmres'), hs=[1e-4 1e-5]; end
            for h=hs
                Jp=r0_robust_eval(ctx,pts{ip}+h*d,F,R); Jm=r0_robust_eval(ctx,pts{ip}-h*d,F,R);
                fd=(Jp-Jm)/(2*h); an=g'*d;
                lg('(3) FD %-8s %s %-7s %-6s sigma=1%% N=8 h=%.0e: J %.6e fd %+.9e adj %+.9e rel %.2e\n',pnm{ip},opn{r},md{1},solver{1},h,J,fd,an,abs(fd-an)/abs(an));
            end
          end
        end
    end
end
% (4) GMRES vs direct
for ip=1:2
    for sg=[0.01 0.03]
        Rg=struct('sigma',sg,'Xi',Xi8,'solver','gmres'); Rd=Rg; Rd.solver='direct';
        [Jg,gg,ig]=r0_robust_eval(ctx,pts{ip},ops.X,Rg); [Jd,gd,id]=r0_robust_eval(ctx,pts{ip},ops.X,Rd);
        lg('(4) gmres vs direct %-8s X sigma=%.0e N=8: J rel %.2e, per-sample P max rel %.2e, grad rel %.2e, gmres its %d-%d, fallbacks %d\n', ...
            pnm{ip},sg,abs(Jg-Jd)/Jd,max(abs(ig.P-id.P)./id.P),norm(gg-gd)/norm(gd),min(ig.gmres_it),max(ig.gmres_it),ig.nfallback);
    end
end
% (5) fluid loading
if ~isempty(zdir)
    Z=load(fullfile(zdir,['iwz_' cell '_dev.mat']),'iwZ','s');
    for ip=1:2
        Ff=struct('factor',ops.X.factor,'iwZ',Z.iwZ);
        [J1,g1]=r0_fsi_eval(ctx,pts{ip},Ff);
        [J2,g2]=r0_robust_eval(ctx,pts{ip},Ff,struct('sigma',0,'Xi',[],'solver','gmres'));
        J0=r0_pilot_eval(ctx,pts{ip},ops.X);
        lg('(5) FSI dev %-8s sigma=0: r0_fsi_eval J %.6e (P_fsi/P_vac %.4f), gmres J rel %.2e, grad rel %.2e\n',pnm{ip},J1,J1/J0,abs(J2-J1)/J1,norm(g2-g1)/norm(g1));
        Rg=struct('sigma',0.01,'Xi',Xi8(:,1:4),'solver','gmres'); Rd=Rg; Rd.solver='direct';
        [Jg,gg,ig]=r0_robust_eval(ctx,pts{ip},Ff,Rg); [Jd,gd]=r0_robust_eval(ctx,pts{ip},Ff,Rd);
        lg('(5) FSI dev %-8s sigma=1%% N=4: gmres vs dense direct J rel %.2e grad rel %.2e, gmres its %d-%d, fallbacks %d\n',pnm{ip},abs(Jg-Jd)/Jd,norm(gg-gd)/norm(gd),min(ig.gmres_it),max(ig.gmres_it),ig.nfallback);
        h=1e-6; Jp=r0_robust_eval(ctx,pts{ip}+h*d,Ff,Rg); Jm=r0_robust_eval(ctx,pts{ip}-h*d,Ff,Rg);
        fd=(Jp-Jm)/(2*h); an=gg'*d;
        lg('(5) FSI dev %-8s sigma=1%% N=4 FD h=1e-6: fd %+.9e adj %+.9e rel %.2e\n',pnm{ip},fd,an,abs(fd-an)/abs(an));
    end
    clear Z Ff
end
% (6) cost
fc=fopen(fullfile(G.out,['unit_' cell '.csv']),'a');
fprintf(fc,'cell,point,solver,N,s_per_eval_J,s_per_eval_Jgrad,gmres_it_mean,mem_matlab_MB\n');
rng(3003); XiT=randn(ne,32);
for ip=1:2
    for solver={'gmres','direct'}
        for N=[1 8 16 32]
            if strcmp(solver{1},'direct') && N>8, continue; end
            R=struct('sigma',0.01,'Xi',XiT(:,1:N),'solver',solver{1});
            t1=pts{ip}*(1+1e-9); clear r0_robust_eval
            tc=tic; r0_robust_eval(ctx,t1,ops.X,R); tJ=toc(tc);
            tc=tic; [~,~,in1]=r0_robust_eval(ctx,t1,ops.X,R); tG=toc(tc);
            t2=pts{ip}*(1+2e-9); clear r0_robust_eval
            tc=tic; [~,~]=r0_robust_eval(ctx,t2,ops.X,R); tJG=toc(tc);
            mm=NaN; try, M=memory; mm=M.MemUsedMATLAB/2^20; catch, end
            lg('(6) cost %-8s %-6s N=%2d: J only %.3f s, gradient after cached J %.3f s, J+grad cold %.3f s, gmres its mean %.1f, MATLAB memory %.0f MB\n', ...
                pnm{ip},solver{1},N,tJ,tG,tJG,mean(in1.gmres_it),mm);
            fprintf(fc,'%s,%s,%s,%d,%.4f,%.4f,%.2f,%.0f\n',cell,pnm{ip},solver{1},N,tJ,tJG,mean(in1.gmres_it),mm);
        end
    end
end
fclose(fc);
lg('R0_ROBUST_UNIT_%s_DONE\n',cell);
fclose(lf);
end

function logboth(f,varargin)
s=sprintf(varargin{:}); fprintf('%s',s); fprintf(f,'%s',s);
end
