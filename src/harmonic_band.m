function [Jb, dJb, info] = harmonic_band(band, t, eval_fun, mode, rho)
% HARMONIC_BAND  Frequency-RANGE (band) aggregation of a harmonic objective and
% its element-thickness gradient. This is the de-aliasing objective: every
% per-frequency metric in this project is evaluated at a SINGLE frequency, so an
% optimizer can "win" by shifting a resonance peak away from the excitation
% frequency -- a trick, not a design. Aggregating a metric over a BAND of
% frequencies removes that loophole: the only way to reduce the aggregate is to
% lower the metric across the whole range (equivalently, to suppress -- not
% relocate -- the resonant peaks inside the band).
%
% USAGE
%   band      struct: band.Nf (number of samples), band.f (1 x Nf frequencies,
%             Hz, ascending). band.f_lo/band.f_hi are carried for bookkeeping.
%   t         design (thickness) field, any shape; only its SIZE matters here.
%   eval_fun  per-frequency evaluator  [Jk, dJk] = eval_fun(t, k)  returning the
%             REAL per-frequency objective Jk (scalar) and its element-wise
%             gradient dJk (same size as t). This is the ONLY coupling to the
%             underlying per-frequency machinery: pass
%             @(tt,k) pair(articleB_shell_eval(cts{k}, tt, metric)) to reuse the
%             verified shell metrics (dyncomp / inputpower / quadvel) -- the
%             metrics themselves are NEVER re-derived here.
%   mode      'mean'  -> band AVERAGE        Jb = mean_k Jk
%             'ks'    -> worst-case / KS soft-max approximating the band MAXIMUM.
%   rho       KS parameter (only for 'ks'); default 60.
%
% AGGREGATED QUANTITY (explicit): every metric here is a REAL scalar once its
% (complex) state is contracted, and that real scalar is what is minimised:
%     dyncomp     Jk = |F.'U|                      (radius of the bilinear pairing)
%     inputpower  Jk = -0.5*w*Im(conj(F).'U)       (active input power, real)
%     quadvel     Jk = w^2 * Re(U^H R0 U)          (real part of the quadratic form)
% Jk is real for every k, so the band aggregation acts on a REAL sequence and the
% complex state only enters through each dJk. A complex aggregator (e.g. KS on a
% complex |g| before taking the modulus) would be WRONG -- aggregate the REAL
% objective actually being minimised, which is what this routine does.
%
% KS (Kreisselmeier-Steinhauser) AGGREGATION -- the band maximum, smoothed:
%     m  = max_k Jk
%     e  = exp(rho*(Jk - m))                 (shifted => numerically safe)
%     Jb = m + log(sum_k e)/rho              -> max_k Jk  as rho -> Inf
%     wk = e / sum_l e                        (softmax weights, sum wk = 1)
%     dJb = sum_k wk * dJk                    (chain rule; the weights are the
%                                              correct ones because KS is a
%                                              function of the REAL Jk alone)
% The shift by m makes the weights invariant to it. As rho grows Jb -> max_k Jk
% from ABOVE and wk concentrate on the worst frequency; Jb <= max_k Jk +
% log(Nf)/rho, so rho >= ~20 already pins Jb to the band max to well under a
% percent for the sample counts used here (the FD test sweeps rho to confirm the
% optimum and the gradient DIRECTION are insensitive to it).
%
% info carries Jk, dJk, the aggregation weights wk, the frequencies, mode/rho and
% the band spec, so the caller can report the per-frequency breakdown.
if nargin < 4 || isempty(mode), mode = 'mean'; end
if nargin < 5 || isempty(rho),  rho  = 60;   end

Nf = band.Nf;
% evaluate the metric at every frequency sample, reusing the per-frequency code
Jk  = zeros(Nf,1);
dJk = [];
for k = 1:Nf
    [Jk(k), dk] = eval_fun(t, k);
    dk = dk(:);
    if isempty(dJk), dJk = zeros(numel(dk), Nf); end
    dJk(:,k) = dk;
end

switch lower(mode)
    case 'mean'
        wk  = ones(Nf,1) / Nf;
        Jb  = mean(Jk);
        dJb = dJk * wk;                       % average of the per-frequency gradients
    case 'ks'
        m   = max(Jk);
        e   = exp(rho*(Jk - m));
        S   = sum(e);
        Jb  = m + log(S)/rho;
        wk  = e / S;                          % softmax -> focuses on the worst f
        dJb = dJk * wk;                       % sum_k wk * dJk
    otherwise
        error('harmonic_band: unknown mode ''%s'' (use ''mean'' or ''ks'').', mode);
end

info = struct('Jk', Jk, 'dJk', dJk, 'wk', wk, 'fk', band.f(:).', ...
              'mode', lower(mode), 'rho', rho, 'Jb', Jb, 'band', band, ...
              'argmax_f', band.f(find(Jk==max(Jk),1)), 'maxJk', max(Jk));
end
