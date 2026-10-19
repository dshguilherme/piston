function cts = articleB_band_ctxs(ctx0, band, zeta)
% ARTICLEB_BAND_CTXS  Build the per-frequency context cell for a band, reusing
% ONE precomputed shell context (the element matrices Bke/Ske/Me are
% frequency-INDEPENDENT, so only w, alpha and the damping bookkeeping change).
%
%   ctx0  context from articleB_build_ctx(...,NaN,0) (no frequency)
%   band  struct with band.Nf and band.f (Hz)
%   zeta  modal damping ratio; alpha = 2*zeta*w (Rayleigh stiffness part), beta=0
if nargin < 3 || isempty(zeta), zeta = 0.05; end
if ~isnan(ctx0.w)
    error('articleB_band_ctxs: ctx0 must be built with freq=NaN (got w=%g).', ctx0.w);
end
Nf = band.Nf;
cts = cell(Nf,1);
for k = 1:Nf
    ck = ctx0;
    ck.w     = 2*pi*band.f(k);
    ck.alpha = 2*zeta*ck.w;
    ck.beta  = 0;
    ck.zeta  = zeta;
    ck.R0    = [];
    cts{k}   = ck;
end
end
