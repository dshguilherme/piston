function vis3d_make_figs
% VIS3D_MAKE_FIGS  Apply the physical-geometry shell renderer to the existing
% Article B design fields.  Produces (all under results/, prefix vis3d_):
%   vis3d_SL_designs       3 Scordelis-Lo designs (dyncomp/inputpower/quadvel),
%                          128x64, physical thickness, shared colour scale
%   vis3d_SL_beforeafter   uniform (equal-volume) shell vs optimised design
%   vis3d_obstacle_SL_designs  3 obstacle-course Scordelis-Lo designs (48x24)
%   vis3d_SL_response      RESPONSE mode: |normal velocity| and signed (phase)
%                          normal velocity riding on the SAME physical geometry
%
% Every figure is written as .png (300 dpi screen), .pdf (pdflatex build) and
% .eps (journal submission vector sibling) with a forced WHITE background.
%
% NOTE ON THE SAVED CONTEXT.  The .mat files were written with geopdes class
% objects (msh_cartesian / sp_scalar) whose class definitions do not match on
% load, so those fields come back as plain structures.  Rather than depend on
% them, this driver REBUILDS the scalar B-spline space and the mesh from the
% geometry's own NURBS (ctx.geometry.nurbs), which reproduces the analysis knot
% vector exactly -> the L2 projection is onto the analysis space.
%
% This file ADDS only; it does not modify src/, the existing result files,
% geopdes or immersed-iga.

repo = 'C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course';
RES  = fullfile(repo,'results');
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath(fullfile(repo,'src')); addpath(fullfile(repo,'tests'));
pub_style();

fprintf('==== VIS3D: physical-geometry shell rendering ====\n');
made = {};

%% ---------------------------------------------------------------- (a) 128x64
L = load(fullfile(RES,'articleB_shell_metrics_128x64.mat'),'R','nsub','ctx','freq');
nsub = L.nsub; ctx = L.ctx;
[sps128, msh128] = build_space_mesh(ctx);
lab = {'dynamic compliance','input power','quadratic velocity'};
P = panel_struct();
for k = 1:3
    P(k) = mk_panel(ctx, sps128, msh128, L.R(k).tbest, nsub, ...
        sprintf('%s  (J=%.2f dB)', lab{k}, 10*log10(L.R(k).meta.Jbest)));
end
o = struct('texag',10, 'ng',200, 'clab','surface thickness t (m)', ...
    'sgtitle', sprintf(['Scordelis-Lo roof, %dx%d: optimised shells as physical 3D ' ...
    'geometry (spline-blended thickness, display x10)'], nsub(1), nsub(2)), ...
    'outbase', fullfile(RES,'vis3d_SL_designs'));
r = render_shell_thickness(P, o); close(r.fig);
made = note(made, o.outbase);

%% ------------------------------------------------------------- (c) before/after
tm = mean(L.R(1).tbest);
P2 = panel_struct();
P2(1) = mk_panel(ctx, sps128, msh128, tm*ones(numel(L.R(1).tbest),1), nsub, ...
    sprintf('uniform shell, t = %.3f m', tm));
P2(2) = mk_panel(ctx, sps128, msh128, L.R(1).tbest, nsub, ...
    'optimised (dynamic compliance)');
o = struct('texag',10, 'ng',200, 'clab','surface thickness t (m)', ...
    'sgtitle', ['Equal-volume comparison: the same curved shell as physical 3D ' ...
    'geometry - thickness varies over the surface (display x10)'], ...
    'outbase', fullfile(RES,'vis3d_SL_beforeafter'));
r = render_shell_thickness(P2, o); close(r.fig);
made = note(made, o.outbase);

%% --------------------------------------------------- (b) obstacle-course 48x24
O = load(fullfile(RES,'articleB_obstacle_scordelis.mat'),'R','ctx','band','fr');
obsub = O.ctx.nsub; octx = O.ctx;
[osps, omsh] = build_space_mesh(octx);
P3 = panel_struct();
for k = 1:3
    P3(k) = mk_panel(octx, osps, omsh, O.R(k).tbest, obsub, ...
        sprintf('%s  (J=%.2f dB)', lab{k}, 10*log10(O.R(k).meta.Jbest)));
end
o = struct('texag',10, 'ng',200, 'clab','surface thickness t (m)', ...
    'sgtitle', sprintf(['Obstacle-course Scordelis-Lo design, %dx%d: optimised shells ' ...
    'as physical 3D geometry (spline-blended thickness, display x10)'], ...
    obsub(1), obsub(2)), ...
    'outbase', fullfile(RES,'vis3d_obstacle_SL_designs'));
r = render_shell_thickness(P3, o); close(r.fig);
made = note(made, o.outbase);

%% ------------------------------------------------ (d) RESPONSE mode (phase)
% Take the obstacle-course dyncomp design; evaluate the harmonic response at the
% MIDDLE band frequency fr(2) (zeta = 0.05, matching the run) and ride the
% normal-velocity response on the SAME physical geometry.
fm = O.fr(2); zeta = 0.05;
[cdn, cabs] = resp_normal_velocity(octx, osps, O.R(1).tbest, fm, zeta, 200);
P4 = panel_struct();
p = mk_panel(octx, osps, omsh, O.R(1).tbest, obsub, 'thickness t');
P4(1) = p;
P4(2) = p; P4(2).title = sprintf('|normal velocity|  (f = %.4f Hz)', fm);
P4(2).cgrid = cabs; P4(2).cmap = parula(256); P4(2).clab = '|v_n| (normalised)';
P4(3) = p; P4(3).title = sprintf('signed normal velocity  (f = %.4f Hz)', fm);
P4(3).cgrid = cdn;  P4(3).cmap = div_cmap(256);
P4(3).clim = [-1 1]*max(abs(cdn(:))); P4(3).clab = 'Re(v_n), normalised';
o = struct('texag',10, 'ng',200, 'clab','surface thickness t (m)', ...
    'sgtitle', ['Response mode: the optimised shell as physical geometry; colour is ' ...
    'the support channel (diverging map shows antiphase regions)'], ...
    'outbase', fullfile(RES,'vis3d_SL_response'));
% keep the thickness panel on a shared thickness colour bar
r = render_shell_thickness(P4, o); close(r.fig);
made = note(made, o.outbase);

%% --------------------------------------------------------------- report
fprintf('\n---- vis3d figure files ----\n');
for i = 1:numel(made)
    for ext = {'.png','.pdf','.eps'}
        f = [made{i} ext{1}]; d = dir(f);
        if ~isempty(d)
            fprintf('  %s  %10d bytes\n', f, d.bytes);
        else
            fprintf('  %s  (MISSING)\n', f);
        end
    end
end
fprintf('VIS3D_DONE\n');
end

% ===========================================================================
function [sps, msh1] = build_space_mesh(ctx)
% Rebuild the scalar B-spline space and the precomputed mesh from the geometry's
% own NURBS (the refined knots match the analysis knot vector exactly).
nrb   = ctx.geometry.nurbs;
knots = nrb.knots;
rule  = msh_gauss_nodes([3 3]);
[qn,qw] = msh_set_quad_nodes(knots, rule);
msh0  = msh_cartesian(knots, qn, qw, ctx.geometry);
msh1  = msh_precompute(msh0);
sps   = sp_nurbs(nrb, msh0);
end

function P = panel_struct()
P = struct('geometry',{},'msh',{},'sp_scalar',{},'nsub',{},'t_elem',{}, ...
           't_mat',{},'title',{},'cgrid',{},'cmap',{},'clim',{},'clab',{});
end

function p = mk_panel(ctx, sps, msh, t, nsub, ttl)
p = struct('geometry',ctx.geometry,'msh',msh,'sp_scalar',sps,'nsub',nsub, ...
           't_elem',t(:),'t_mat',[],'title',ttl,'cgrid',[],'cmap',[], ...
           'clim',[],'clab',[]);
end

function made = note(made, base)
made{end+1} = base;
end

% ===========================================================================
function [cdn, cabs] = resp_normal_velocity(ctx, sps, t, freq, zeta, ng)
% Harmonic response of the optimised design at one frequency; return the normal
% velocity on a fine parameter grid: signed (Re -> antiphase) and magnitude.
w = 2*pi*freq;
cK = 1; cM = 1i*w*(2*zeta*w) - w^2;
t  = t(:);
Ks = sparse(ctx.Rows(:), ctx.Cols(:), ctx.Bke.*(t.^3) + ctx.Ske.*t, ctx.ndof, ctx.ndof);
M  = sparse(ctx.Rows(:), ctx.Cols(:), ctx.Me.*t,                   ctx.ndof, ctx.ndof);
KD = cK*Ks + cM*M;
u  = zeros(ctx.ndof,1);
u(ctx.free) = KD(ctx.free,ctx.free) \ ctx.F(ctx.free);

Ns = sps.ndof; nu = sps.ndof_dir(1); nv = sps.ndof_dir(2);
% unit normal on the fine grid (vectorised NURBS evaluation, see
% nurbs_surface_grid.m -- geometry.map_der is far too slow on a fine grid)
U = linspace(0,1,ng); V = linspace(0,1,ng);
[~,~,~,Nx,Ny,Nz] = nurbs_surface_grid(ctx.geometry.nurbs, U, V);
% displacement components (sp_vector stores component dofs in contiguous blocks)
Dx = eval_grid(reshape(u(0*Ns+(1:Ns)),nu,nv), sps, U, V);
Dy = eval_grid(reshape(u(1*Ns+(1:Ns)),nu,nv), sps, U, V);
Dz = eval_grid(reshape(u(2*Ns+(1:Ns)),nu,nv), sps, U, V);
un = Dx.*Nx + Dy.*Ny + Dz.*Nz;          % complex normal displacement
vn = 1i*w*un;                            % normal velocity (complex)
cdn  = real(vn);                         % signed -> antiphase map
cabs = abs(vn);                          % magnitude
mx = max(abs(cdn(:))); if mx>0, cdn = cdn/mx; end
mx = max(cabs(:));     if mx>0, cabs = cabs/mx; end
end

function Tf = eval_grid(Cmat, sps, U, V)
Nu = bbas(sps.sp_univ(1), U); Nv = bbas(sps.sp_univ(2), V);
Tf = Nv * Cmat.' * Nu.';               % ndgrid order (numel(V) x numel(U))
end

function N = bbas(spu, q)
kn = spu.knots(:).'; p = spu.degree; nbf = spu.ndof;
q = q(:); nq = numel(q); N = zeros(nq, nbf);
for i=1:nq
    spn = findspan(nbf-1, p, q(i), kn);
    idx = numbasisfun(spn, q(i), p, kn) + 1;
    N(i, idx) = basisfun(spn, q(i), p, kn);
end
end

function cm = div_cmap(n)
b = [0.23 0.30 0.75]; w = [1 1 1]; r = [0.70 0.11 0.15];
cm = zeros(n,3); tt = linspace(0,1,n).';
for i=1:n
    if tt(i) < 0.5, a = tt(i)/0.5; cm(i,:) = b + (w-b)*a;
    else,           a = (tt(i)-0.5)/0.5; cm(i,:) = w + (r-w)*a; end
end
end
