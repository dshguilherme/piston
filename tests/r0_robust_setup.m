function [ctx, ops, G] = r0_robust_setup(cell)
% R0_ROBUST_SETUP  Structural context and radiation factors of a regime-pilot cell for Experiment 10.
% ctx  : the pilot's cap model (pilot_ctx of tests/r0_pilot_runs.m, copied verbatim)
% ops.X, ops.S : exact and sinc factors normalised so that J(uniform) = Jref (as r0_long_runs load_ops)
% ops.Xraw, ops.Sraw : unnormalised factors (physical power in W, pilot medium)
% G    : paths, Jref, the pilot optimiser settings, and the cell definition (theta0, kR)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fileparts(root)); addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
G.root=root;
G.pilot=fullfile(root,'results','r0_revision','pilot');
G.exp5=fullfile(root,'results','r0_revision','exact_opt');
G.fsi=fullfile(root,'results','r0_revision','fsi');
G.long=fullfile(root,'results','r0_revision','longruns');
G.out=fullfile(root,'results','r0_revision','robust');
if ~exist(G.out,'dir'), mkdir(G.out); end
G.Jref=5.9464e-7;
G.opt=struct('tmin',0.025,'tmax',0.10,'volfrac',0.5,'nmax',inf,'stop',true,'ftol',1e-4, ...
    'xchgtol',1e-4,'nhold',3,'window',10,'feastol',1e-6,'Abar',100,'dbfac',10);
[G.th,G.kR]=cell_def(cell); G.cell=cell;
ctx=pilot_ctx(G.th); ctx.cachekey=sprintf('th%d',G.th);
G.opt.vweights=ones(ctx.ne,1);
O=load(fullfile(G.pilot,['ops_' cell '.mat']),'cn','Fs_c');
F=load(fullfile(G.pilot,['factor_' cell '.mat']),'Fx');
ops.Xraw=struct('factor',F.Fx); ops.Sraw=struct('factor',O.Fs_c);
ops.X=struct('factor',F.Fx*sqrt(O.cn(3)));
ops.S=struct('factor',O.Fs_c*sqrt(O.cn(2)));
ops.cn=O.cn;
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
