function r0_consistent_gradient_test()
% Gradient test for consistent angular-Gram objective, independent of optimizer.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
spec=struct('R',25,'L',50,'th0',deg2rad(40),'E',2e11,'NU',0.3,'RHO',7850, ...
 't0',0.1,'deg',[2 2],'reg',[1 1],'nsub',[8 6],'nquad',[4 4], ...
 'bc',[1 2 3 4],'loadtype','point','loaduv',[0.5 0.5], ...
 'loaddir',[0;1;0],'loadmag',1,'freq',6.53547,'zeta',0.05);
ctx=obstacle_build_ctx('cylpanel',spec); ctx.alpha=0; ctx.beta=0.1/ctx.w;
Rop=r0_build_consistent(ctx,1.21,343,[4 4],20,40);
t=0.05*ones(ctx.ne,1); rng(413); d=randn(ctx.ne,1); d=d/norm(d);
[J,g]=articleB_r0_eval(ctx,t,'quadvel',Rop,'consistent');
for h=[1e-3 3e-4 1e-4 3e-5]
    Jp=articleB_r0_eval(ctx,t+h*d,'quadvel',Rop,'consistent');
    Jm=articleB_r0_eval(ctx,t-h*d,'quadvel',Rop,'consistent');
    fd=(Jp-Jm)/(2*h); an=g'*d; rel=abs(fd-an)/max([abs(fd),abs(an),eps]);
    fprintf('h=%g fd=%+.9e adj=%+.9e rel=%g\n',h,fd,an,rel);
end
fprintf('R0_CONSISTENT_GRADIENT_TEST_DONE\n');
end
