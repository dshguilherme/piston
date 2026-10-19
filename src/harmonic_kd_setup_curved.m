function S = harmonic_kd_setup_curved(msh, space, geometry, YOUNG, POISSON, RHO, device, precision, bend_model)
% HARMONIC_KD_SETUP_CURVED One-time setup of the metric-carrying bending
% operator for a GENERAL (curved / non-affine) 2D map. This is the curved
% counterpart of harmonic_kd_setup.m: it builds the identical weighted-
% quadrature PAGES (WBp / WBe, the same rule types (k1,k2) in {0,1,2}), but the
% bending coefficient page is the 5-component curved page of
% quadrule_bend_curved.m (25 slots) instead of the affine 3-component page
% (9 slots), and the parametric component table is
%     xi xi = (2,0), xi eta = (1,1), eta eta = (0,2), xi = (1,0), eta = (0,1).
%
% The bending operator it feeds (harmonic_kd_form_curved.m) is
%     Kb = int D det(J) Phi' G Phi dx ,  Phi = [w_xixi; w_xieta; w_etaeta; w_xi; w_eta]
% which is the EXACT second-order chain rule of the physical Hessian
% (Christoffel / inverse-map second-derivative terms included). For an affine
% map the two extra components carry no energy and Kb reproduces
% harmonic_kd_setup exactly.
%
% Only the BENDING page is built here (that is the object the task concerns);
% the mass and membrane pages of harmonic_kd_setup.m are affine-safe and
% unchanged. Everything here is a NEW file: no existing src/ file is modified.
%
% INPUTS are the same as harmonic_kd_setup.m (RHO accepted for signature parity
% but unused, since the mass page is not built).
%   bend_model: 'T1' (default) or 'plate'.

if nargin < 7 || isempty(device),     device     = 'cpu'; end
if nargin < 8 || isempty(precision),  precision  = 'double'; end
if nargin < 9 || isempty(bend_model), bend_model = 'T1'; end
assert(msh.ndim == 2, 'harmonic_kd_setup_curved: only 2D charts implemented.');

t_setup = tic;
nsd = 2;
sp_s = space.scalar_spaces{1};
n = sp_s.ndof_dir;
N = prod(n);

%% 1. WQ rules (identical to harmonic_kd_setup.m)
for d = 1:nsd
    sp1d = sp_s.sp_univ(d);
    nb{d} = cellfun(@(x) unique(sp1d.connectivity(:, x)).', sp1d.supp, 'UniformOutput', false);
    Q2{d} = quadrule_bend_fast(sp1d);
    QR(d) = Q2{d};
    brk{d} = [sp_s.knots{d}(1), sp_s.knots{d}(end)];
    qn{d} = QR(d).all_points';
end
new_msh = msh_cartesian(brk, qn, [], geometry);
space_wq = space.constructor(new_msh);
nqtot = cellfun(@numel, qn);

%% 2. Pages per direction (identical construction to harmonic_kd_setup.m)
for d = 1:nsd
    sp1d = space_wq.scalar_spaces{1}.sp_univ(d);
    nd = n(d);
    nnmax(d) = max(cellfun(@numel, nb{d}));
    nqmax(d) = max(QR(d).nquad_points);
    P{d}   = ones(nqmax(d), nd);
    J{d}   = ones(nnmax(d), nd);
    len{d} = zeros(1, nd);
    lself{d} = zeros(1, nd);
    for k1 = 0:2
        for k2 = 0:2
            WBp{d}{k1+1, k2+1} = zeros(nnmax(d), nqmax(d), nd); %#ok<AGROW>
        end
    end
    for ii = 1:nd
        pts = QR(d).ind_points{ii};
        ng  = nb{d}{ii};
        nq = numel(pts); nn = numel(ng);
        P{d}(1:nq, ii) = pts;
        J{d}(1:nn, ii) = ng;
        len{d}(ii) = nn;
        lself{d}(ii) = find(ng == ii, 1);
        trial = {sp1d.shape_functions(pts, ng).', ...           % k2 = 0
                 sp1d.shape_function_gradients(pts, ng).', ...  % k2 = 1
                 full(Q2{d}.B2_global(pts, ng)).'};             % k2 = 2
        for k1 = 0:2
            for k2 = 0:2
                wgt = QR(d).(sprintf('quad_weights_%d%d', k1, k2)){ii};
                WBp{d}{k1+1, k2+1}(1:nn, 1:nq, ii) = wgt .* trial{k2+1};
            end
        end
    end
end

[I1, I2] = ind2sub(n, 1:N);

%% 3. Curved bending coefficient page (25 slots)
[bendC, cinfo] = quadrule_bend_curved(geometry, qn, YOUNG, POISSON, bend_model, 'analytic');
bend = struct('model', lower(bend_model), 'D', cinfo.D, 'comps', cinfo.comps, ...
              'detJ', cinfo.detJ, 'mapmode', cinfo.mapmode);

%% 4. Gather index and pre-gathered coefficient pages
LIN = reshape(P{1}(:, I1), nqmax(1), 1, N) + ...
      (reshape(P{2}(:, I2), 1, nqmax(2), N) - 1) * nqtot(1);

benchk = max(abs(bendC(:)));
isconst = all(max(abs(bendC - bendC(1,:)), [], 1) <= 1e-12 * max(benchk, eps));
bendKp = cell(1, 25);
if isconst
    for k = 1:25, bendKp{k} = reshape(bendC(1,k), 1, 1); end
else
    for k = 1:25, bendKp{k} = reshape(bendC(LIN, k), nqmax(1), nqmax(2), N); end
end

%% 5. Sparsity (identical to harmonic_kd_setup.m)
mask = (reshape((1:nnmax(1))', nnmax(1), 1, 1) <= reshape(len{1}(I1), 1, 1, N)) & ...
       (reshape(1:nnmax(2), 1, nnmax(2), 1) <= reshape(len{2}(I2), 1, 1, N));
ROWp = repmat(reshape(1:N, 1, 1, N), nnmax(1), nnmax(2), 1);
COLp = reshape(J{1}(:, I1), nnmax(1), 1, N) + (reshape(J{2}(:, I2), 1, nnmax(2), N) - 1) * n(1);
rows = ROWp(mask);
cols = COLp(mask);

% Expanded per-row pages: WBe{d}{k1+1,k2+1} is [nnmax x nqmax x N]
WBe = cell(1, nsd);
for d = 1:nsd
    Id = I1; if d == 2, Id = I2; end
    for k1 = 0:2
        for k2 = 0:2
            WBe{d}{k1+1, k2+1} = WBp{d}{k1+1, k2+1}(:, :, Id);
        end
    end
end

%% 6. Density-to-WQ-point element maps (for optional thickness pages)
knots_u1 = unique(sp_s.knots{1});
knots_u2 = unique(sp_s.knots{2});
e1 = discretize(qn{1}, knots_u1);
e2 = discretize(qn{2}, knots_u2);

%% 7. Assemble the setup struct (CPU/double; the GPU path is out of scope here)
S = struct();
S.device = lower(device); S.precision = precision;
S.N = N; S.n = n; S.nqtot = nqtot; S.nnmax = nnmax; S.nqmax = nqmax;
S.WBe = WBe; S.bendKp = bendKp; S.bend = bend; S.bendK_const = isconst;
S.comps = cinfo.comps;
S.lself = {lself{1}(I1), lself{2}(I2)};
S.LIN = uint32(LIN); S.mask = mask; S.rows = rows; S.cols = cols;
S.e1 = e1; S.e2 = e2;
S.YOUNG = YOUNG; S.POISSON = POISSON;
S.nnz_block = numel(rows);
S.nel_dir = [size(sp_s.sp_univ(1).connectivity,2), size(sp_s.sp_univ(2).connectivity,2)];
S.QR = QR;
S.t_setup = toc(t_setup);
end
