function [J, dJ] = articleB_band_eval(cts, k, t, metric)
% ARTICLEB_BAND_EVAL  Thin per-frequency wrapper for harmonic_band: returns the
% first two outputs (objective, gradient) of the verified shell evaluator at the
% k-th band frequency.
[J, dJ] = articleB_shell_eval(cts{k}, t(:), metric);
end
