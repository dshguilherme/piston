function out = render_shell_thickness(panels, opt)
%RENDER_SHELL_THICKNESS  Render optimised shells as PHYSICAL 3D geometry whose
%thickness genuinely varies over the surface, with colour only as support.
%
%   out = render_shell_thickness(panels, opt)
%
% THE SMOOTHING STEP ("blending").  The design field produced by the optimiser
% is ELEMENT-WISE CONSTANT (one thickness per knot span), which is why direct
% pictures of it look blocky.  This renderer SMOOTHS it by projecting the
% element-centre values onto a genuinely COARSER tensor-product B-spline space
% (clamped uniform, degree opt.smooth_deg, with opt.smooth_frac of the elements
% in each direction) by least squares.  This IS the blending that removes the
% blockiness: a plain L2 projection onto the ANALYSIS knot vector reproduces the
% element-to-element jumps (Gibbs-like spikes and overshoot beyond the design
% bounds), whereas the coarse space yields a smooth t_hat(xi,eta).  The
% mid-surface is then offset along its unit normal by t_hat/2 to build the real
% solid sandwich (outer + inner surfaces + the four rim walls).  Colour is the
% SECONDARY channel; the 3D form carries the geometry.
%
% PANELS (struct array).  Fields read here:
%   geometry   geopdes geo struct; must carry .nurbs (NURBS toolbox surface)
%   nsub       [n1 n2] number of analysis elements per direction
%   t_elem     ne x 1 element-wise field, physical units (i fastest, then j),
%              ne = n1*n2,  OR
%   t_mat      n1 x n2 element grid
%   title      panel title (char)
%   cgrid      OPTIONAL grid to colour the surface by (instead of t_hat), e.g.
%              a RESPONSE field (mode 2).  Resampled if its size differs.
%   cmap/clim/clab   per-panel colour overrides (else opt defaults)
%   (extra fields such as msh/sp_scalar are tolerated and ignored)
%
% OPT (all optional):
%   smooth_frac  coarse-space resolution as a fraction of the element count in
%                each direction (default 0.20; smaller = smoother).
%   smooth_deg   spline degree of the smoothing space (default 2).
%   texag        thickness display EXAGGERATION (default 8): offset = texag*t/2.
%   ng           fine render-grid resolution (default 200).
%   cmap         default colormap (parula -- perceptually uniform).
%   clim         default [lo hi] shared colour limits (auto from data if empty).
%   clab         default colorbar label ('thickness t (m)').
%   view [az el], figsize, outbase (saves png/pdf/eps), sgtitle.
%
% OUT: .fig, .tmin/.tmax (physical thickness over panels), .clim, .texag,
%   .files, .checks.

addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');

if nargin < 2, opt = struct(); end
def = struct('texag',8, 'ng',200, 'cmap',parula(256), 'clim',[], ...
             'clab','thickness t (m)', 'view',[-37 22], 'figsize',[], ...
             'outbase','', 'sgtitle','', 'smooth_frac',0.20, 'smooth_deg',2, ...
             'blur',3);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opt, fn{i}) || isempty(opt.(fn{i})), opt.(fn{i}) = def.(fn{i}); end
end

nM = numel(panels);
tP = tic;

% ---- 1. smooth-project every field, build fine geometry -------------------
P = struct('Tf',{},'X',{},'Y',{},'Z',{},'Nx',{},'Ny',{},'Nz',{});
rawLo = zeros(1,nM); rawHi = zeros(1,nM);
for k = 1:nM
    pk = panels(k);
    if isfield(pk,'t_mat') && ~isempty(pk.t_mat)
        tg = double(pk.t_mat);
    else
        tg = reshape(double(pk.t_elem(:)), pk.nsub(1), pk.nsub(2));
    end
    rawLo(k) = min(tg(:)); rawHi(k) = max(tg(:));
    n1 = pk.nsub(1); n2 = pk.nsub(2);
    ngx = opt.ng; ngy = opt.ng;
    ar = n2 / n1;                       % keep fine quads roughly square
    if ar > 1, ngy = round(opt.ng*sqrt(ar)); else, ngx = round(opt.ng/sqrt(ar)); end

    % --- blending: de-spike (Gaussian blur of the element values, which stays
    %     inside the design bounds) then least-squares fit onto a COARSER
    %     tensor-product B-spline space -> smooth t_hat without ringing ------
    tg = gauss_blur2(tg, opt.blur);
    [spu1, spu2] = coarse_spaces(n1, n2, opt.smooth_frac, opt.smooth_deg);
    xi  = ((1:n1)-0.5)/n1;  eta = ((1:n2)-0.5)/n2;      % element centres
    Nu  = bspline_basis_row(spu1, xi);                  % n1 x m1
    Nv  = bspline_basis_row(spu2, eta);                 % n2 x m2
    Cmat = (Nu.'*Nu) \ (Nu.' * tg * Nv) / (Nv.'*Nv);    % m1 x m2 control net

    % --- fine grid: geometry + unit normal + smoothed t_hat ---------------
    U = linspace(0,1,ngx); V = linspace(0,1,ngy);
    [X,Y,Z,Nx,Ny,Nz] = nurbs_surface_grid(pk.geometry.nurbs, U, V);
    Tf = eval_field_grid(Cmat, spu1, spu2, U, V);       % ngy x ngx

    P(k) = struct('Tf',Tf,'X',X,'Y',Y,'Z',Z,'Nx',Nx,'Ny',Ny,'Nz',Nz);
end

tmin = min(cellfun(@(A) min(A(:)), {P.Tf}));
tmax = max(cellfun(@(A) max(A(:)), {P.Tf}));
fprintf(['  [render] blended+rendered %d panels in %.1fs | raw t %.3f..%.3f | ' ...
         'blended t %.3f..%.3f\n'], nM, toc(tP), min(rawLo), max(rawHi), tmin, tmax);
if isempty(opt.clim), climv = [min(rawLo) max(rawHi)]; else, climv = opt.clim; end

% ---- 2. draw ---------------------------------------------------------------
if isempty(opt.figsize), opt.figsize = [3.6*nM + 0.9, 4.7]; end
fig = figure('Units','inches','Position',[0.6 0.6 opt.figsize(1) opt.figsize(2)], ...
             'Color','w');
tl = tiledlayout(1, nM, 'TileSpacing','compact','Padding','loose');

files = {}; checks = struct();
for k = 1:nM
    pk = panels(k); Q = P(k);
    ax = nexttile; hold(ax,'on'); set(ax,'Color','w');

    e = opt.texag;
    Xo = Q.X + (e/2)*Q.Tf.*Q.Nx;  Yo = Q.Y + (e/2)*Q.Tf.*Q.Ny;  Zo = Q.Z + (e/2)*Q.Tf.*Q.Nz;
    Xi = Q.X - (e/2)*Q.Tf.*Q.Nx;  Yi = Q.Y - (e/2)*Q.Tf.*Q.Ny;  Zi = Q.Z - (e/2)*Q.Tf.*Q.Nz;

    % per-panel colour: grid + map + limits + label
    if isfield(pk,'cgrid') && ~isempty(pk.cgrid), Cg = pk.cgrid; else, Cg = Q.Tf; end
    if ~isequal(size(Cg), size(Q.X))       % resample a supplied grid to ours
        xq = linspace(0,1,size(Q.X,2)); yq = linspace(0,1,size(Q.X,1));
        x0 = linspace(0,1,size(Cg,2));  y0 = linspace(0,1,size(Cg,1));
        [XQ,YQ] = meshgrid(xq,yq); [X0,Y0] = meshgrid(x0,y0);
        Cg = interp2(X0, Y0, Cg, XQ, YQ, 'linear');
    end
    if isfield(pk,'cmap') && ~isempty(pk.cmap), cm = pk.cmap; else, cm = opt.cmap; end
    if isfield(pk,'clim') && ~isempty(pk.clim), cl = pk.clim; else, cl = climv; end
    if isfield(pk,'clab') && ~isempty(pk.clab), lab = pk.clab; else, lab = opt.clab; end

    % inner surface: pale grey (the underside of the solid)
    surf(ax, Xi, Yi, Zi, 'EdgeColor','none','FaceColor',[0.83 0.845 0.865], ...
         'FaceLighting','gouraud','FaceAlpha',1);
    draw_rim(ax, Xo,Yo,Zo, Xi,Yi,Zi);   % lofted rim closes the solid
    % outer surface: colour = value (secondary channel)
    surf(ax, Xo, Yo, Zo, Cg, 'EdgeColor','none','FaceColor','interp','FaceAlpha',1);
    shading(ax,'interp');

    colormap(ax, cm); caxis(ax, cl);
    material(ax, [0.45 0.62 0.16 9 0.35]);
    camlight(ax,'headlight'); camlight(ax,'left'); lighting(ax,'gouraud');
    axis(ax,'equal'); axis(ax,'tight'); view(ax, opt.view);
    grid(ax,'off'); set(ax,'Box','on');
    set(ax,'LineWidth',0.7,'XColor',[0 0 0],'YColor',[0 0 0],'ZColor',[0 0 0]);
    xlabel(ax,'x','Color',[0 0 0]); ylabel(ax,'y','Color',[0 0 0]); zlabel(ax,'z','Color',[0 0 0]);
    cb = colorbar(ax); cb.Color = [0 0 0];
    cb.Label.String = lab; cb.Label.Color = [0 0 0];
    if isfield(pk,'title') && ~isempty(pk.title)
        title(ax, pk.title, 'Interpreter','none','Color',[0 0 0],'FontSize',10);
    end
    text(ax, 0.02, -0.15, sprintf('t = %.3f - %.3f m', min(Q.Tf(:)), max(Q.Tf(:))), ...
        'Units','normalized','Color',[0.15 0.15 0.15],'FontSize',8,'Interpreter','none');
end

if ~isempty(opt.sgtitle)
    annotation(fig, 'textbox', [0.02 0.005 0.96 0.055], 'String', opt.sgtitle, ...
        'EdgeColor','none','HorizontalAlignment','center','VerticalAlignment','middle', ...
        'Color',[0 0 0],'FontSize',9,'Interpreter','none','FitBoxToText','off');
end

% ---- 3. save (PNG/PDF/EPS, white background) ------------------------------
out = struct('fig',fig,'tmin',tmin,'tmax',tmax,'clim',climv,'texag',opt.texag, ...
             'files',{files},'checks',checks);
if ~isempty(opt.outbase)
    drawnow;
    out.files = pub_save_all(fig, opt.outbase, 300);
end
end

% ===========================================================================
function A = gauss_blur2(A, sigma)
% Separable Gaussian blur of a 2D grid with replicate padding, used to de-spike
% the element-constant design BEFORE the spline fit.  A weighted average cannot
% leave the [min,max] range, so the design bounds are preserved.  (Hand-rolled
% replicate padding: padarray needs the Image Processing Toolbox.)
if sigma <= 0, return; end
r = max(1, ceil(3*sigma));
k = exp(-((-r:r).^2)/(2*sigma^2)); k = k(:)/sum(k);
A = conv2([repmat(A(1,:),r,1); A; repmat(A(end,:),r,1)], k, 'valid');
A = conv2([repmat(A(:,1),1,r), A, repmat(A(:,end),1,r)], k.', 'valid');
end

function [spu1, spu2] = coarse_spaces(n1, n2, frac, deg)
% Clamped uniform tensor-product spline spaces with ~frac*n elements per
% direction (at least 4).  This is the blending resolution: coarser -> smoother.
nc1 = max(4, round(n1*frac)); nc2 = max(4, round(n2*frac));
spu1 = struct('knots', unif_clamped(nc1,deg), 'degree',deg, 'ndof', nc1+deg);
spu2 = struct('knots', unif_clamped(nc2,deg), 'degree',deg, 'ndof', nc2+deg);
end

function kn = unif_clamped(nc, deg)
kn = [zeros(1,deg), linspace(0,1,nc+1), ones(1,deg)];
end

function Tf = eval_field_grid(Cmat, spu1, spu2, U, V)
% t_hat(xi,eta) = sum Nv(eta) Cmat Nu(xi);  ndgrid order Tf(ky,kx).
Nu = bspline_basis_row(spu1, U);       % nU x m1
Nv = bspline_basis_row(spu2, V);       % nV x m2
Tf = Nv * Cmat.' * Nu.';               % nV x nU
end

function N = bspline_basis_row(spu, q)
% nq x ndof matrix of univariate B-spline basis values at the points q.
kn = spu.knots(:).'; p = spu.degree; nbf = spu.ndof;
q = q(:); nq = numel(q); N = zeros(nq, nbf);
for i = 1:nq
    spn = findspan(nbf-1, p, q(i), kn);
    idx = numbasisfun(spn, q(i), p, kn) + 1;
    N(i, idx) = basisfun(spn, q(i), p, kn);
end
end

% ===========================================================================
function draw_rim(ax, Xo,Yo,Zo, Xi,Yi,Zi)
% Close the solid along the four parametric edges (the lofted rim) so local
% thickness reads as a real wall.  Each edge is a 2-row strip [outer; inner].
g = [0.55 0.58 0.62];
ox = {Xo(1,:), Xo(end,:), Xo(:,1).', Xo(:,end).'};
oy = {Yo(1,:), Yo(end,:), Yo(:,1).', Yo(:,end).'};
oz = {Zo(1,:), Zo(end,:), Zo(:,1).', Zo(:,end).'};
ix = {Xi(1,:), Xi(end,:), Xi(:,1).', Xi(:,end).'};
iy = {Yi(1,:), Yi(end,:), Yi(:,1).', Yi(:,end).'};
iz = {Zi(1,:), Zi(end,:), Zi(:,1).', Zi(:,end).'};
for i = 1:4
    surf(ax, [ox{i}; ix{i}], [oy{i}; iy{i}], [oz{i}; iz{i}], ...
         'EdgeColor','none','FaceColor',g,'FaceLighting','gouraud','FaceAlpha',1);
end
end

% ===========================================================================
function files = pub_save_all(fig, base, res)
% Save PNG + PDF + EPS (white background).  Heavily shaded 3D surfaces: use
% pub_save's 'raster' mode -- a high-resolution raster inside the PDF/EPS
% wrapper is the honest format (painters cannot keep a shaded 3D surface vector
% and would try, very slowly, to vectorise thousands of polygons).
files = {};
if exist('pub_save','file')
    pub_save(fig, base, res, 'raster');
    files = {[base '.png']};
    if exist([base '.pdf'],'file'), files{end+1} = [base '.pdf']; end
    if exist([base '.eps'],'file'), files{end+1} = [base '.eps']; end
else
    set(fig,'Color','w');
    exportgraphics(fig, [base '.png'], 'Resolution', res, 'BackgroundColor','white');
    exportgraphics(fig, [base '.pdf'], 'ContentType','image','Resolution',res,'BackgroundColor','white');
    exportgraphics(fig, [base '.eps'], 'ContentType','image','Resolution',res,'BackgroundColor','white');
    files = {[base '.png'],[base '.pdf'],[base '.eps']};
end
end
