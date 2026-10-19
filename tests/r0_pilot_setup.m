function r0_pilot_setup(stage)
% R0 regime pilot, structural set-up (no acoustic solve, no BEM).
%   stage 'modes'   : natural frequencies of the uniform caps (theta0 = 30, 90, 120 deg;
%                     t = 0.025, 0.05, 0.10) and a 1-200 Hz forced-response sweep of the
%                     uniform t = 0.05 design with the paper's damping -> pilot/modes_*.csv
%   stage 'samples' : acoustic surface samples (q8 Gauss points per element) of each cap's
%                     spline space, B: U -> normal displacement, plus element centroids
%                     and 40 seeded random smooth thickness fields   -> pilot/samples_th<deg>.mat
%   stage 'timing'  : cost of one articleB_r0_eval call with a dense 1730 x ndof factor
% Structure = the paper's cap model (tests/r0_opt_export_samples.m cap_ctx): steel, R = 10 m,
% nsub [40 24], quadratic C1 splines, clamped edge (bc = 2), axisymmetric pole load,
% 100 Hz, stiffness-proportional damping beta = 0.1/w (loss factor 0.1 at 100 Hz).
% Only the half-angle alpha changes between geometries.
if nargin<1, stage='modes'; end
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fileparts(root)); addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
out=fullfile(root,'results','r0_revision','pilot'); if ~exist(out,'dir'), mkdir(out); end
TH=[30 90 120];
switch stage
  case 'modes'
    fid=fopen(fullfile(out,'modes_natural.csv'),'w');
    fprintf(fid,'theta0_deg,t_uniform,mode,f_Hz,Omega,pole_participation\n');
    fs=fopen(fullfile(out,'modes_sweep.csv'),'w');
    fprintf(fs,'theta0_deg,f_Hz,loss_factor,abs_drive_point,kinetic_energy_J\n');
    for th=TH
        ctx=pilot_ctx(th);
        cL=sqrt(ctx.E/(ctx.RHO*(1-0.3^2)));
        for t0=[0.025 0.05 0.10]
            t=t0*ones(ctx.ne,1);
            [K,M]=km(ctx,t); fr=ctx.free;
            Kf=K(fr,fr); Mf=M(fr,fr); Kf=(Kf+Kf.')/2; Mf=(Mf+Mf.')/2;
            [V,D]=eigs(Kf,Mf,60,'smallestabs');
            [lam,ix]=sort(real(diag(D))); V=V(:,ix);
            f=sqrt(max(lam,0))/2/pi;
            mm=real(sum(conj(V).*(Mf*V),1)).';
            pf=abs(V.'*ctx.F(fr))./sqrt(mm); pf=pf/max(pf);
            for n=1:numel(f)
                fprintf(fid,'%d,%.3f,%d,%.6f,%.6f,%.4e\n',th,t0,n,f(n),2*pi*f(n)*10/cL,pf(n));
            end
            fprintf('theta0=%d t=%.3f: f1=%.3f Hz, f60=%.3f Hz, modes below 100 Hz: %d (pole-participating >1e-3: %d)\n',...
                th,t0,f(1),f(end),sum(f<100),sum(f<100 & pf>1e-3));
        end
        t=0.05*ones(ctx.ne,1); [K,M]=km(ctx,t); fr=ctx.free;
        for fq=1:1:200
            w=2*pi*fq; A=(1+1i*w*ctx.beta)*K-w^2*M;
            u=zeros(ctx.ndof,1); u(fr)=A(fr,fr)\ctx.F(fr);
            ke=0.25*w^2*real(u'*M*u);
            fprintf(fs,'%d,%d,%.4f,%.6e,%.6e\n',th,fq,w*ctx.beta,abs(ctx.F.'*u),ke);
        end
    end
    fclose(fid); fclose(fs);
  case 'samples'
    for th=TH
        ctx=pilot_ctx(th);
        s=r0_surface_samples(ctx,zeros(ctx.ndof,0),[8 8]);
        x=s.x; w=s.weights; B=s.basis; w0=ctx.w; ndof=ctx.ndof; ne=ctx.ne; theta0=deg2rad(th); %#ok<NASGU>
        r=sqrt(sum(x.^2,1));
        fprintf('theta0=%d: %d points, |x| in [%.10f, %.10f], area %.8f (exact %.8f), max polar %.6f deg\n',...
            th,numel(w),min(r),max(r),sum(w),2*pi*100*(1-cosd(th)),max(acosd(x(3,:)./r)));
        % element centroids (from the structural quadrature)
        g=ctx.msh.geo_map; c=squeeze(mean(g,2));               % 3 x ne
        rc=sqrt(sum(c.^2,1)); ecth=acos(c(3,:)./rc); ecph=atan2(c(2,:),c(1,:));
        % uniform state at 100 Hz
        t=0.05*ones(ctx.ne,1); [K,M]=km(ctx,t); A=(1+1i*ctx.w*ctx.beta)*K-ctx.w^2*M;
        Uu=complex(zeros(ctx.ndof,1)); Uu(ctx.free)=A(ctx.free,ctx.free)\ctx.F(ctx.free);
        % 40 seeded random smooth thickness fields, mean 0.05 (volume fraction 0.5), within [0.025, 0.10]
        rng(1000+th); nr=40; Trand=zeros(ctx.ne,nr); u=ecth(:)/theta0; v=ecph(:)/(2*pi);
        for j=1:nr
            g=zeros(ctx.ne,1);
            for a=0:3
                for b=0:3
                    g=g+randn/(1+a+b)*cos(a*pi*u).*cos(2*pi*b*v);
                    if b>0, g=g+randn/(1+a+b)*cos(a*pi*u).*sin(2*pi*b*v); end
                end
            end
            g=g-mean(g); g=g/max(abs(g)); A=0.005+0.020*rand;
            tj=0.05+A*g;
            for k=1:50, tj=min(max(tj,0.025),0.10); tj=tj-mean(tj)+0.05; end
            tj=min(max(tj,0.025),0.10);
            Trand(:,j)=tj;
        end
        save(fullfile(out,sprintf('samples_th%d.mat',th)),'x','w','B','w0','ndof','ne','theta0', ...
            'ecth','ecph','Uu','Trand','-v7');
    end
  case 'timing'
    ctx=pilot_ctx(30); rng(1);
    Fd=complex(randn(1730,ctx.ndof),randn(1730,ctx.ndof))*1e-3; R0=struct('factor',Fd);
    t=0.05*ones(ctx.ne,1); tic; for k=1:10, [J,g]=articleB_r0_eval(ctx,t,'quadvel',R0,'consistent'); end %#ok<ASGLU>
    fprintf('articleB_r0_eval, dense 1730 x %d factor: %.3f s per call (threads %d)\n',ctx.ndof,toc/10,maxNumCompThreads);
end
fprintf('R0_PILOT_SETUP_%s_DONE\n',stage);
end

function [K,M]=km(ctx,t)
KE=ctx.Bke.*(t.^3)+ctx.Ske.*t; ME=ctx.Me.*t;
K=sparse(ctx.Rows(:),ctx.Cols(:),KE(:),ctx.ndof,ctx.ndof);
M=sparse(ctx.Rows(:),ctx.Cols(:),ME(:),ctx.ndof,ctx.ndof);
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
