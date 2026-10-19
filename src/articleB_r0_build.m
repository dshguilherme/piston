function R0op = articleB_r0_build(ctx, rho0, c0, mrank, fquery)
% ARTICLEB_R0_BUILD  Build the thesis radiation-resistance weighting operator
% R0 (thesis eqs. 3.63-3.66) for the 'surface-averaged quadratic velocity' /
% radiated-sound-power objective
%
%       W_rad = V^H * R0 * V                                          (3.65)
%       R_ef  = (omega^2*rho0*S^2)/(4*pi*c0) * sin(k*r_ef)/(k*r_ef)   (3.66)
%
% on the Article B curved-shell context (Scordelis-Lo roof).
%
% OBSERVATION POINTS (the 'spatial data' in R0). One equivalent point per
% knot-span (finite element), placed at the element CENTRE - the
% quadrature-weighted centroid of the mapped element. The system velocity
% vector V is the set of surface-NORMAL velocities at those centres,
%
%       V_e = 1i*omega * (C*U)_e ,   (C*U)_e = element-average of (u . n_e) ,
%
% where U is the FE displacement (3 components, N dofs), n_e the element-average
% surface normal, and C the (ne x N) centre-averaging operator. The element
% average uses the B-spline partition of unity (sum over local basis of the mass
% matrix row) so C is exact for degree-2 elements.
%
%   A_e   = element area  = sum_q jacdet*qw
%   x_e   = element centre (quadrature-weighted centroid)
%   r_ef  = |x_e - x_f|   ,  k = omega/c0
%   sinc  = sin(k r)/(k r)  ( = 1 on the diagonal r=0 )
% With uniform A_e = S this is exactly eq. (3.66); we carry the per-element
% areas A_e in place of S^2 (S_i S_j weighting), which reduces to 3.66 when the
% areas are uniform. This is the ONE interpretation we adopt; the thesis states
% S as the surface area but does not fix whether it is a single reference area
% or the per-span area - we use the per-span area (the standard discretisation
% of the Rayleigh integral W = (rho0 c0 k^2/4pi) int int v* v sin(kr)/(kr)).
%
% DESIGN-INDEPENDENCE. Centres, areas, the normal directions and omega do NOT
% depend on the thickness field t, so R0 is assembled ONCE per frequency and
% reused for every design. (A thickness-dependent R0 would add a product-rule
% term that the sensitivity routine does not compute.)
%
% R0op also carries the two CONTROL functionals of the phase-cancellation
% ladder (see articleB_r0_eval):
%   'I'    R0 = I         ->  W = sum_e |V_e|^2        (phase-BLIND)
%   'ones' R0 = C*ones    ->  W = C*|sum_e V_e|^2      (net volume velocity,
%                            the true r_ef -> 0 point-source / monopole limit of
%                            the FULL R0; MAXIMALLY cancellation-sensitive)
%   'rad'  R0 = eq.(3.66) ->  W = V^H R0 V             (physical radiated power)
% NOTE the rank-1 'ones' operator is NOT the identity and NOT its point-source
% limit: the identity is a different operator with opposite (phase-blind)
% behaviour. They must never be conflated.
%
% Inputs: ctx from articleB_build_ctx (must carry ctx.w); rho0, c0 fluid.

if nargin < 2 || isempty(rho0),  rho0  = 1.21; end   % air
if nargin < 3 || isempty(c0),    c0    = 343;  end   % air
if nargin < 4 || isempty(mrank), mrank = 96;   end
if nargin < 5, fquery = []; end

m = ctx.msh; ne = ctx.ne; N = ctx.ndof; n = size(ctx.lm,2);
ncomp = 3; q = n/ncomp;    % local basis per component (degree+1)^2

% ---- element area weights, centres, averaged normals ----------------------
Jw   = m.jacdet .* m.quad_weights;                 % [nqn ne]
area = sum(Jw,1);                                  % [1 ne]
Jw3  = reshape(Jw, 1, [], ne);
gm   = reshape(m.geo_map, 3, [], ne);
cen  = reshape(squeeze(sum(gm.*Jw3,2)) ./ area, 3, ne);    % [3 ne] centres
nm   = reshape(m.normal, 3, [], ne);
navg = reshape(squeeze(sum(nm.*Jw3,2)) ./ area, 3, ne);    % [3 ne] avg normals

% ---- centre-averaging operator C  (ne x N) --------------------------------
% (C*U)_e = (1/A_e) sum_c navg(c,e) * sum_a (int_e phi_a) U_{a,c}
% int_e phi_a  =  row-sum of the element mass block / RHO  (partition of unity)
RHO = ctx.RHO;
ii = zeros(ne*n,1); jj = zeros(ne*n,1); vv = zeros(ne*n,1); p = 0;
for e = 1:ne
    Me = reshape(ctx.Me(e,:), n, n);       % n x n element mass (carries RHO)
    for c = 1:ncomp
        blk = Me((c-1)*q+(1:q), (c-1)*q+(1:q));
        intphi = sum(blk,2)/RHO;           % q-vector: int_e phi_a
        g = intphi * navg(c,e) / area(e);  % q-vector of C entries
        cidx = (c-1)*q + (1:q);
        for a = 1:q
            p = p+1; ii(p) = e; jj(p) = ctx.lm(e, cidx(a)); vv(p) = g(a);
        end
    end
end
C = sparse(ii(1:p), jj(1:p), vv(1:p), ne, N);

% ---- frequency, wave number, prefactor, Nystrom factors of K ---------------
if ~isempty(fquery), omega = 2*pi*fquery; else, omega = ctx.w; end
if isempty(omega) || isnan(omega)
    error('articleB_r0_build: ctx.w is empty/NaN - build the context WITH a frequency.');
end
k    = omega/c0;
pref = (omega^2 * rho0)/(4*pi*c0);
% sinc kernel K_ef = sin(k r_ef)/(k r_ef) in LOW-RANK (Nystrom) form so R0
% scales to fine meshes (no ne x ne dense matrix).
rng(12345); S = randperm(ne, min(mrank,ne)).';
Ds = distmatrix(cen, cen(:,S));            % ne x m
KS = sincmat(k, Ds.');                     % m x ne
KSS = KS(:,S);                             % m x m
KSSinv = pinv((KSS+KSS.')/2);
if ne <= 6000
    Kfull = sincmat(k, distmatrix(cen,cen));   % exact, small meshes only
else
    Kfull = [];
end

R0op = struct();
R0op.C = C; R0op.area = area(:); R0op.centres = cen; R0op.navg = navg;
R0op.k = k; R0op.omega = omega; R0op.pref = pref;
R0op.rho0 = rho0; R0op.c0 = c0; R0op.ne = ne; R0op.freq = omega/(2*pi);
R0op.S = S; R0op.KS = KS; R0op.KSS = KSS; R0op.KSSinv = KSSinv; R0op.Kfull = Kfull; R0op.mrank = numel(S);
end

% ---------------------------------------------------------------------------
function D = distmatrix(A,B)
DX = A(1,:).' - B(1,:); DY = A(2,:).' - B(2,:); DZ = A(3,:).' - B(3,:);
D = sqrt(DX.^2 + DY.^2 + DZ.^2);
end
function KM = sincmat(k,D)
KM = ones(size(D));
idx = D ~= 0;
KM(idx) = sin(k*D(idx))./(k*D(idx));
end
