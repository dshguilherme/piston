function r0_seam_trace_tests()
% Direct asymmetric spline trace checks for periodic sphere-cap and cylinder seams.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
base=struct('R',10,'L',20,'alpha',deg2rad(30),'E',2e11,'NU',0.3,'RHO',7850, ...
 't0',0.1,'deg',[2 2],'reg',[1 1],'nsub',[12 8],'nquad',[4 4], ...
 'bc',[],'loadtype','point','loaduv',[0.37 0.42],'loaddir',[0;1;0], ...
 'loadmag',1,'freq',100,'zeta',0.05);
for kind={'spherecap','cylshell'}
    k=kind{1}; spec=base;
    if strcmp(k,'cylshell'), spec.R=2; spec.L=20; expected=1; seam_dim=1;
    else, expected=2; seam_dim=2; end
    ctx=obstacle_build_ctx(k,spec);
    periodic=ctx.space.scalar_spaces{1}.periodic_dir;
    assert(ismember(expected,periodic),'Wrong periodic direction in %s.',k);
    rng(100+seam_dim); u=randn(ctx.ndof,1);
    q={linspace(0.13,0.87,7),linspace(0.19,0.81,5)};
    qa=q; qb=q; qa{seam_dim}=0; qb{seam_dim}=1;
    [ea,Fa]=sp_eval(u,ctx.space,ctx.geometry,qa,{'value','gradient'});
    [eb,Fb]=sp_eval(u,ctx.space,ctx.geometry,qb,{'value','gradient'});
    n=numel(q{3-seam_dim}); va=ea{1}; ga=ea{2}; vb=eb{1}; gb=eb{2};
    % Since seam coordinate is singleton, outputs align with the other coordinate.
    va=reshape(va,3,[]); vb=reshape(vb,3,[]); Fa=reshape(Fa,3,[]); Fb=reshape(Fb,3,[]);
    va=va(:,1:n); vb=vb(:,1:n); Fa=Fa(:,1:n); Fb=Fb(:,1:n);
    ga=reshape(ga,3,3,[]); gb=reshape(gb,3,3,[]);
    valerr=norm(va-vb,'fro')/max(norm(va,'fro'),1);
    xerr=norm(Fa-Fb,'fro')/max(norm(Fa,'fro'),1);
    graderr=norm(ga-gb,'fro')/max(norm(ga,'fro'),1);
    fprintf('%s seam value %.3e geometry %.3e physical-gradient %.3e\n',k,valerr,xerr,graderr);
    assert(valerr<1e-10 && xerr<1e-10 && graderr<1e-8,'Seam trace failure in %s.',k);
end
fprintf('R0_SEAM_TRACE_TESTS_PASS\n');
end
