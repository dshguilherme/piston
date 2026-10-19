function S = harmonic_kd_setup(msh, space, geometry, YOUNG, POISSON, RHO, device, precision, bend_model)
% HARMONIC_KD_SETUP One-time setup for the matrix-free harmonic K_D operator
% (flat Kirchhoff-Love plate in PHYSICAL coordinates, decoupled transverse
% field, N-DOF scalar operator used by the T1 validation suite).
%
%   K_D = (t^3 Kb + t Km) + j*w*(alpha*M(t) + beta*(t^3 Kb + t Km)) - w^2 M(t)
%
% BENDING IN PHYSICAL COORDINATES (metric-carrying)
% The physical Hessian is expanded in the parametric Hessian components
% h = [w_xi xi ; w_xi eta ; w_eta eta] through the inverse map Jacobian,
%     w_xx = a'h,  w_yy = b'h,  w_xy = c'h,
%     a = [Ji11^2; 2 Ji11 Ji21; Ji21^2]
%     b = [Ji12^2; 2 Ji12 Ji22; Ji22^2]
%     c = [Ji11 Ji12; Ji11 Ji22 + Ji12 Ji21; Ji21 Ji22],   Ji = inv(J)
% so that the energy density is D*det(J) * h' M h with the symmetric 3x3
%     M = a a' + b b' + nu (a b' + b a') + 2(1-nu) c c'      ('plate')
%     M = a a' + b b' + 2 c c'                                ('T1')
% The nine entries M_{alpha,beta} are stored as a [nq1 x nq2 x 9] page, so the
% bending contraction runs over the nine (alpha,beta) terms, each pairing
%   direction 1 : rule (i,k) ,  direction 2 : rule (j,l)
% for components alpha = (i,j), beta = (k,l) with the component map
%   xi xi = (2,0) ,  xi eta = (1,1) ,  eta eta = (0,2).
% This is exact for AFFINE maps (the second-derivative-of-map term vanishes);
% for curved shells the Christoffel terms belong to the shell kinematics module.
%
% bend_model: 'T1' (default) or 'plate' (full KL plate with Poisson coupling).

if nargin < 7 || isempty(device), device = 'cpu'; end
if nargin < 8 || isempty(precision), precision = 'double'; end
if nargin < 9 || isempty(bend_model), bend_model = 'T1'; end
assert(msh.ndim == 2, 'harmonic_kd_setup: only 2D (flat plate) implemented.');
assert(nargin >= 6 && ~isempty(RHO), 'harmonic_kd_setup: RHO (mass density) required.');

t_setup = tic;
nsd = 2;
sp_s = space.scalar_spaces{1};
n = sp_s.ndof_dir;
N = prod(n);

%% 1. WQ rules (stiffness family + all bending rule types)
for d = 1:nsd
    sp1d = sp_s.sp_univ(d);
    nb{d} = cellfun(@(x) unique(sp1d.connectivity(:, x)).', sp1d.supp, 'UniformOutput', false);
    Q2{d} = quadrule_bend_fast(sp1d); %#ok<AGROW>
    QR(d) = Q2{d};                    %#ok<AGROW>
    brk{d} = [sp_s.knots{d}(1), sp_s.knots{d}(end)];
    qn{d} = QR(d).all_points';
end
new_msh = msh_cartesian(brk, qn, [], geometry);
space_wq = space.constructor(new_msh);
nqtot = cellfun(@numel, qn);

%% 2. Pages per direction: type (k1,k2) -> weight .* trial basis of order k2
% trial derivative order k2 = 0 -> Bv, 1 -> Bd, 2 -> Bxx
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
        lself{d}(ii) = find(ng == ii, 1);   % the dof is always its own neighbour
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

%% 3. Constitutive / metric pages
jac = msh.map_der(qn);

% 3a. Membrane (plane stress D, 3x3) -- identical to wq_setup.m
E = zeros(3, 2, 2);
E(:, :, 1) = [1 0; 0 0; 0 1];
E(:, :, 2) = [0 0; 0 1; 1 0];
lambda = YOUNG * POISSON / ((1 + POISSON) * (1 - 2 * POISSON));
mu = YOUNG / (2 * (1 + POISSON));
D0 = blkdiag(lambda * ones(nsd) + 2 * mu * eye(nsd), mu);
Cm = C_ijkl(E, jac, D0, nsd, nqtot);

% 3b. Metric-carrying bending coefficient page [nq1 x nq2 x 9]
Dbf = YOUNG / (12 * (1 - POISSON^2));
switch lower(bend_model)
    case 't1',    with_nu = false;
    case 'plate', with_nu = true;
    otherwise, error('harmonic_kd_setup: unknown bend_model ''%s''.', bend_model);
end
% The map Jacobian is evaluated directly from the geometry at all WQ points
% (msh.map_der returns a flat point list, not a tensor grid). Point ordering
% matches the LIN gather: index = i1 + (i2-1)*nqtot(1).
[Q1g, Q2g] = meshgrid(qn{1}, qn{2});
P_all = [Q1g(:).'; Q2g(:).'];
[~, Jall] = geometry.map_der(P_all);
Jall = reshape(Jall, 2, 2, nqtot(1), nqtot(2));
bendK = zeros(nqtot(1)*nqtot(2), 9);      % flat, in the LIN point ordering
for i2 = 1:nqtot(2)
    for i1 = 1:nqtot(1)
        Jm = Jall(:, :, i1, i2);
        dJ = Jm(1,1)*Jm(2,2) - Jm(1,2)*Jm(2,1);
        Ji = inv(Jm);
        av = [Ji(1,1)^2; 2*Ji(1,1)*Ji(2,1); Ji(2,1)^2];
        bv = [Ji(1,2)^2; 2*Ji(1,2)*Ji(2,2); Ji(2,2)^2];
        cv = [Ji(1,1)*Ji(1,2); Ji(1,1)*Ji(2,2) + Ji(1,2)*Ji(2,1); Ji(2,1)*Ji(2,2)];
        Mm = av*av.' + bv*bv.' + 2*cv*cv.';
        if with_nu
            Mm = Mm + POISSON*(av*bv.' + bv*av.') - 2*POISSON*(cv*cv.');
        end
        p_flat = i1 + (i2-1)*nqtot(1);
        bendK(p_flat, :) = reshape(Dbf * dJ * Mm, 1, 9);
    end
end
bend = struct('model', lower(bend_model), 'Dbf', Dbf);

% 3c. Mass coefficient page (scalar, RHO*|J| at WQ points)
m1 = RHO * (squeeze(jac(1,1,:,:)) .* squeeze(jac(2,2,:,:)) - ...
            squeeze(jac(1,2,:,:)) .* squeeze(jac(2,1,:,:)));
m1 = m1(:);

%% 4. Gather index and pre-gathered pages
LIN = reshape(P{1}(:, I1), nqmax(1), 1, N) + ...
      (reshape(P{2}(:, I2), 1, nqmax(2), N) - 1) * nqtot(1);
M0p = m1(LIN);   % [nq1 x nq2 x N]

% Bending coefficient pages for the contraction: [nq1 x nq2 x N] per slot.
% For an AFFINE map the coefficient is spatially constant, so it is kept as a
% 1x1 page and broadcast against the thickness pages (cheap); otherwise it is
% gathered over the sparsity grid (general, memory heavier).
benchk = max(abs(bendK(:)));
isconst = all(max(abs(bendK - bendK(1,:)), [], 1) <= 1e-12 * max(benchk, eps));
bendKp = cell(1, 9);
if isconst
    for k = 1:9, bendKp{k} = reshape(bendK(1,k), 1, 1); end
else
    for k = 1:9, bendKp{k} = reshape(bendK(LIN, k), nqmax(1), nqmax(2), N); end
end

C0p = cell(nsd, nsd);
for k1 = 1:nsd
    for k2 = 1:nsd
        blk = zeros(nqmax(1), nqmax(2), N, 4);
        c = 0;
        for j = 1:nsd
            for i = 1:nsd
                c = c + 1;
                Cij = Cm{i, j}(:, :, k1, k2);
                blk(:, :, :, c) = Cij(LIN);
            end
        end
        C0p{k1, k2} = blk;
    end
end

% Mass rule (global_mid_nonzero) shares the stiffness point set: its weights
% solve the identical least-squares system as quad_weights_00.
S_mass_rule = struct('QM', {cell(1, nsd)});
for d = 1:nsd
    S_mass_rule.QM{d} = global_mid_nonzero(sp_s.sp_univ(d));
end

% Row point/neighbour bookkeeping for element-local quadratic forms
for d = 1:nsd
    Pd{d} = P{d}; Jd{d} = J{d}; lend{d} = len{d};
    nd = n(d); nnm = nnmax(d);
    posm = sparse(nd, nd);
    for a = 1:nd
        posm(Jd{d}(1:lend{d}(a), a), a) = (1:lend{d}(a)).';
    end
    POSd{d} = posm;
end

% Self-neighbour indices per row (for the exact diagonal)
LS1 = lself{1}(I1); LS2 = lself{2}(I2);

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

% Official weight-type table (membrane family)
typ = zeros(nsd, nsd, nsd);
for k1 = 1:nsd
    for k2 = 1:nsd
        for ll = 1:nsd
            if k1 == k2 && ll == k1
                typ(k1, k2, ll) = 4;
            elseif k1 ~= k2 && ll == k1
                typ(k1, k2, ll) = 2;
            elseif k1 ~= k2 && ll == k2
                typ(k1, k2, ll) = 3;
            else
                typ(k1, k2, ll) = 1;
            end
        end
    end
end

%% 5. Sparsity (identical to wq_setup.m)
mask = (reshape((1:nnmax(1))', nnmax(1), 1, 1) <= reshape(len{1}(I1), 1, 1, N)) & ...
       (reshape(1:nnmax(2), 1, nnmax(2), 1) <= reshape(len{2}(I2), 1, 1, N));
ROWp = repmat(reshape(1:N, 1, 1, N), nnmax(1), nnmax(2), 1);
COLp = reshape(J{1}(:, I1), nnmax(1), 1, N) + (reshape(J{2}(:, I2), 1, nnmax(2), N) - 1) * n(1);
rows = ROWp(mask);
cols = COLp(mask);

%% 6. Device / precision
cast_fn = @(x) cast(x, precision);
if strcmpi(device, 'gpu')
    mv = @(x) gpuArray(cast_fn(x));
else
    mv = cast_fn;
end
for d = 1:nsd
    for k1 = 0:2
        for k2 = 0:2
            WBe{d}{k1+1, k2+1} = mv(WBe{d}{k1+1, k2+1});
        end
    end
end
for k1 = 1:nsd
    for k2 = 1:nsd
        C0p{k1, k2} = mv(C0p{k1, k2});
    end
end
M0p = mv(M0p);
for k = 1:9, bendKp{k} = mv(bendKp{k}); end
if strcmpi(device, 'gpu'), LINd = gpuArray(uint32(LIN)); else, LINd = uint32(LIN); end
maskd = mask;
if strcmpi(device, 'gpu'), maskd = gpuArray(mask); end
if strcmpi(device, 'gpu'), wait(gpuDevice); end

%% 7. Gauss-grid evaluation bundle (objective + sensitivity post-processing)
% The solution's value / gradient / Hessian fields are needed on the standard
% per-element Gauss grid to read the element energies off the state, exactly as
% fast_sensitivities.m reads the strain energy for elastostatics. The physical
% Hessian is recovered with the metric coefficients evaluated on that grid.
sp1g = sp_s.sp_univ(1); sp2g = sp_s.sp_univ(2);
sp1h = sp_bspline_1d_param(sp1g.knots, sp1g.degree, msh.qn{1}, 'hessian', true);
sp2h = sp_bspline_1d_param(sp2g.knots, sp2g.degree, msh.qn{2}, 'hessian', true);
Ev.B1   = sp1g.shape_functions;          Ev.B2   = sp2g.shape_functions;
Ev.Bx1  = sp1g.shape_function_gradients; Ev.Bx2  = sp2g.shape_function_gradients;
Ev.Bxx1 = sp1h.shape_function_hessians;  Ev.Bxx2 = sp2h.shape_function_hessians;
Ev.conn1 = sp1g.connectivity; Ev.conn2 = sp2g.connectivity;
Ev.qw1 = msh.qw{1}; Ev.qw2 = msh.qw{2};
Ev.nqn1 = msh.nqn_dir(1); Ev.nqn2 = msh.nqn_dir(2);
Ev.nel1 = msh.nel_dir(1); Ev.nel2 = msh.nel_dir(2);
ng1 = Ev.nqn1; ng2 = Ev.nqn2;
Ev.detJ = zeros(ng1, ng2, Ev.nel1, Ev.nel2);
Ev.CXX  = zeros(ng1, ng2, Ev.nel1, Ev.nel2, 3);
Ev.CYY  = Ev.CXX; Ev.CXY = Ev.CXX;
for el1 = 1:Ev.nel1
    for el2 = 1:Ev.nel2
        for g1 = 1:ng1
            for g2 = 1:ng2
                P = [msh.qn{1}(g1,el1); msh.qn{2}(g2,el2)];
                [~, jp] = geometry.map_der(P);
                Ev.detJ(g1,g2,el1,el2) = jp(1,1)*jp(2,2) - jp(1,2)*jp(2,1);
                Ji = inv(jp);
                Ev.CXX(g1,g2,el1,el2,:) = [Ji(1,1)^2, 2*Ji(1,1)*Ji(2,1), Ji(2,1)^2];
                Ev.CYY(g1,g2,el1,el2,:) = [Ji(1,2)^2, 2*Ji(1,2)*Ji(2,2), Ji(2,2)^2];
                Ev.CXY(g1,g2,el1,el2,:) = [Ji(1,1)*Ji(1,2), Ji(1,1)*Ji(2,2)+Ji(1,2)*Ji(2,1), Ji(2,1)*Ji(2,2)];
            end
        end
    end
end
Ev.model = lower(bend_model);

% Density-to-WQ-point element maps
knots_u1 = unique(sp_s.knots{1});
knots_u2 = unique(sp_s.knots{2});
e1 = discretize(qn{1}, knots_u1);
e2 = discretize(qn{2}, knots_u2);

S = struct();
S.device = lower(device); S.precision = precision;
S.N = N; S.n = n; S.nqtot = nqtot; S.nnmax = nnmax; S.nqmax = nqmax;
S.WBe = WBe; S.C0p = C0p; S.M0p = M0p; S.bendKp = bendKp; S.bend = bend;
S.bendK_const = isconst;
S.lself = {LS1, LS2};
% Element connectivity (for element-wise quadratic forms in the sensitivities)
S.P{1} = Pd{1}; S.P{2} = Pd{2};
S.J{1} = Jd{1}; S.J{2} = Jd{2};
S.length{1} = lend{1}; S.length{2} = lend{2};
S.POS{1} = POSd{1}; S.POS{2} = POSd{2};
S.conn = {sp_s.sp_univ(1).connectivity, sp_s.sp_univ(2).connectivity};
S.nel_dir = [size(sp_s.sp_univ(1).connectivity,2), size(sp_s.sp_univ(2).connectivity,2)];
S.typ = typ;
S.LIN = LINd; S.mask = maskd; S.rows = rows; S.cols = cols;
S.YOUNG = YOUNG; S.POISSON = POISSON; S.RHO = RHO;
S.nnz_block = numel(rows);
S.e1 = e1; S.e2 = e2;
S.Ev = Ev;
S.mass_rule = S_mass_rule;
S.t_setup = toc(t_setup);
end
