function op = r0_build_consistent(ctx, rho0, c0, nquad, nmu, nphi)
%R0_BUILD_CONSISTENT Build a consistent spline-pullback sinc power factor.
% This is a refinement of the prescribed-source sinc model, not a general
% curved-shell exterior acoustic solution. Uses peak displacement coefficients:
% P = omega^2 ||L U||^2. Angular sampling controls positive Gram compression;
% surface Gauss sampling and angular sampling must be refined independently.
if nargin<2||isempty(rho0),rho0=1.21;end
if nargin<3||isempty(c0),c0=343;end
if nargin<4||isempty(nquad),nquad=ctx.spec.nquad;end
if nargin<5||isempty(nmu),nmu=32;end
if nargin<6||isempty(nphi),nphi=64;end
s=r0_surface_samples(ctx,zeros(ctx.ndof,0),nquad);
g=r0_angular_gram(s.x,s.weights,s.basis,ctx.w,rho0,c0,nmu,nphi,'rayleigh');
op=g; op.surface_order=nquad; op.surface_points=numel(s.weights);
op.element_area=s.area; op.total_area=sum(s.area);
op.dof=ctx.ndof; op.geometry=ctx.kind; op.power_variable='peak displacement';
end
