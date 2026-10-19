function r0_opt_export_samples(nqlist)
% Export acoustic surface samples of the spherical-cap spline space so that the
% exact spectral two-sided operator (tests/r0_nov_capspec.py) can be reduced to a
% Hermitian PSD form on the STRUCTURAL coefficient space (r0_opt_exact_factor.py).
% Writes results/r0_revision/exact_opt/cap_samples_q<nq>.mat with
%   x (3 x np), w (np x 1), B (np x ndof sparse: U -> normal displacement),
%   U (frozen states of the six paper designs), designs, labels, ctx frequency.
% No structural or acoustic solve beyond the six frozen states.
if nargin<1, nqlist=[8 12]; end
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fileparts(root)); addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
out=fullfile(root,'results','r0_revision','exact_opt'); if ~exist(out,'dir'), mkdir(out); end
ctx=cap_ctx();
S=load(fullfile(root,'results','r0_revision','spherecap_frozen_states.mat'),'U','designs','labels');
% recompute the frozen states from the stored designs (consistency gate)
U=complex(zeros(ctx.ndof,size(S.designs,2)));
for j=1:size(S.designs,2)
    t=S.designs(:,j); KE=ctx.Bke.*(t.^3)+ctx.Ske.*t; ME=ctx.Me.*t;
    A=sparse(ctx.Rows(:),ctx.Cols(:),(1+1i*ctx.w*ctx.beta)*KE(:)-ctx.w^2*ME(:),ctx.ndof,ctx.ndof);
    U(ctx.free,j)=A(ctx.free,ctx.free)\ctx.F(ctx.free);
end
gate=norm(U-S.U,'fro')/norm(S.U,'fro');
fprintf('frozen-state recompute rel diff %.3e\n',gate); assert(gate<1e-8);   % sparse-direct round-off (thread count) gives ~3e-10
designs=S.designs; labels=S.labels; w0=ctx.w; ndof=ctx.ndof; ne=ctx.ne; %#ok<NASGU>
for nq=nqlist
    s=r0_surface_samples(ctx,zeros(ndof,0),[nq nq]);
    x=s.x; w=s.weights; B=s.basis; area=s.area; %#ok<NASGU>
    r=sqrt(sum(x.^2,1));
    fprintf('nq=%d: %d points, |x| in [%.12f, %.12f], total area %.10f (sphere cap 2piR^2(1-cos30)=%.10f)\n',...
        nq,numel(w),min(r),max(r),sum(w),2*pi*100*(1-cos(pi/6)));
    save(fullfile(out,sprintf('cap_samples_q%d.mat',nq)),'x','w','B','area','U','designs','labels','w0','ndof','ne','-v7');
end
fprintf('R0_OPT_EXPORT_DONE\n');
end

function ctx=cap_ctx()
spec=struct('R',25,'L',50,'th0',deg2rad(40),'alpha',deg2rad(30), ...
    'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,'deg',[2 2],'reg',[1 1], ...
    'nsub',[40 24],'nquad',[4 4],'bc',[1 2 3 4],'loadtype','point', ...
    'loaduv',[0.5 0.5],'loaddir',[0;1;0],'loadmag',1,'freq',6.53547,'zeta',0.05);
spec.R=10; spec.bc=2; spec.loadtype='pole'; spec.loaddir=[0;0;1]; spec.freq=100;
ctx=obstacle_build_ctx('spherecap',spec); ctx.alpha=0; ctx.beta=0.1/ctx.w;
end
