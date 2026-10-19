function r0_viz_sinc_freq_sweep()
% Prescribed-source sinc radiated power of FROZEN spherical-cap fields vs frequency.
% ALGEBRA (frozen velocity). U = displacement at w0=2*pi*100 (ctx.w = w0). Normal velocity
% v = i*w0*u_n is held fixed. Rayleigh/sinc power
%   P(f) = (rho*omega_f^2/(4*pi*c)) * v^H [W K(k_f) W] v ,  k_f = omega_f/c,  |v|^2 = w0^2|u_n|^2.
% Legacy builder articleB_r0_build(ctx,rho,c,m,f): pref = omega_f^2*rho/(4*pi*c) (acoustic
% Rayleigh factor, uses queried omega) and Kfull uses k_f; the displacement->velocity factor is
% the separate leading ctx.w^2 = w0^2.  So P = w0^2 * pref(f) * real(Z'*Kfull(f)*Z): the code
% formula of the audit is already the frozen-velocity form once fquery is passed (ctx.w stays w0).
% Consistent: r0_angular_gram(...,omega_f,...) puts omega_f^2 in pref and k_f in the phase;
% P = w0^2 * sum|factor|^2.  No re-solve; no existing file is written.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fileparts(root)); addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
out=fullfile(root,'results','r0_revision','viz'); if ~exist(out,'dir'), mkdir(out); end
logf=fopen(fullfile(out,'cap_sinc_freq_sweep.log'),'w');
lg=@(varargin) [fprintf(varargin{:}); fprintf(logf,varargin{:})];
t0=tic;
spec=struct('R',25,'L',50,'th0',deg2rad(40),'alpha',deg2rad(30), ...
    'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,'deg',[2 2],'reg',[1 1], ...
    'nsub',[40 24],'nquad',[4 4],'bc',[1 2 3 4],'loadtype','point', ...
    'loaduv',[0.5 0.5],'loaddir',[0;1;0],'loadmag',1,'freq',6.53547,'zeta',0.05);
spec.R=10; spec.bc=2; spec.loadtype='pole'; spec.loaddir=[0;0;1]; spec.freq=100;
ctx=obstacle_build_ctx('spherecap',spec); ctx.alpha=0; ctx.beta=0.1/ctx.w;
S=load(fullfile(root,'results','r0_revision','spherecap_frozen_states.mat'),'U','labels');
U=S.U; labels=S.labels; nf=size(U,2);
lg('labels: %s\n',strjoin(labels,', ')); lg('ctx.w=%.9g (2*pi*100=%.9g)\n',ctx.w,2*pi*100);
samples=r0_surface_samples(ctx,U,[8 8]); Y=samples.normal_displacement;
rho=1.21; c=343; w0=ctx.w; R=10;
kR=[0.1 0.3 1 2 4 6 8 10 12 14 16 18.3183 21];
kRd=logspace(log10(0.1),log10(21),60);
rows={};
% listed points: na=64 (all), na=96 for the highest two
for i=1:numel(kR)
    f=kR(i)*c/(2*pi*R);
    [Pc,Pn]=evalf(f,64);
    for j=1:nf, rows(end+1,:)={kR(i),f,labels{j},Pc(j),Pn(j),64,0}; end %#ok<AGROW>
    if i>=numel(kR)-1
        [~,Pn96]=evalf(f,96);
        for j=1:nf, rows(end+1,:)={kR(i),f,labels{j},Pc(j),Pn96(j),96,0}; end %#ok<AGROW>
        lg('kR=%.4f f=%.4f angular convergence (consistent) na=64 vs 96, max rel diff %.3e\n',kR(i),f,max(abs(Pn96-Pn)./abs(Pn)));
        for j=1:nf, lg('   %-14s P64=%.9e P96=%.9e rel=%.2e\n',labels{j},Pn(j),Pn96(j),abs(Pn96(j)-Pn(j))/abs(Pn(j))); end
    end
    if abs(kR(i)-18.3183)<1e-9
        % validation gate at 100 Hz (use exact f=100 rather than kR rounding)
        [Pc100,Pn100]=evalf(100,64);
        T=readtable(fullfile(root,'results','r0_revision','frozen_field_quadrature.csv'));
        worst=0;
        for j=1:nf
            m=strcmp(T.geometry,'spherecap')&strcmp(T.design,labels{j})&T.surface_order==8&T.angular_order==64;
            e1=abs(Pc100(j)-T.legacy_centroid_W(m))/abs(T.legacy_centroid_W(m));
            e2=abs(Pn100(j)-T.consistent_spline_W(m))/abs(T.consistent_spline_W(m));
            worst=max([worst e1 e2]);
            lg('GATE %-14s centroid %.6e vs %.6e (rel %.2e) | consistent %.6e vs %.6e (rel %.2e)\n',labels{j},Pc100(j),T.legacy_centroid_W(m),e1,Pn100(j),T.consistent_spline_W(m),e2);
        end
        lg('VALIDATION GATE worst relative error %.3e -> %s\n',worst,string(worst<1e-8));
        assert(worst<1e-8,'Validation gate failed');
    end
end
for i=1:numel(kRd)
    f=kRd(i)*c/(2*pi*R);
    [Pc,Pn]=evalf(f,64);
    for j=1:nf, rows(end+1,:)={kRd(i),f,labels{j},Pc(j),Pn(j),64,1}; end %#ok<AGROW>
end
T=cell2table(rows,'VariableNames',{'kR','freq_Hz','design','P_centroid','P_consistent','angular_order','dense_flag'});
writetable(T,fullfile(out,'cap_sinc_freq_sweep.csv'));
lg('runtime %.1f s\nDONE\n',toc(t0)); fclose(logf);

    function [Pc,Pn]=evalf(f,na)
        legacy=articleB_r0_build(ctx,rho,c,128,f);
        Z=legacy.area.*(legacy.C*U);
        Pc=(w0^2*legacy.pref*real(sum(conj(Z).*(legacy.Kfull*Z),1)));
        op=r0_angular_gram(samples.x,samples.weights,Y,2*pi*f,rho,c,na,2*na);
        Pn=w0^2*sum(abs(op.factor).^2,1);
    end
end
