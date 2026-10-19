function [Kb, Vb] = shell_kd_form(H, t_e, output)
% SHELL_KD_FORM  Assembly of the thin-shell bending operator from the WQ pages
% of shell_kd_setup.m. It is the SAME 5x5-page contraction as
% harmonic_kd_form_curved.m (which is reused verbatim), documented here so the
% shell module is self-contained:
%
%   Kb = int t(x)^3 * D * sqrt(a) * Phi' G_shell Phi dxi deta ,
%   Phi = [w_xixi; w_xieta; w_etaeta; w_xi; w_eta]
%
% The five Phi components have the map 1:(2,0) 2:(1,1) 3:(0,2) 4:(1,0) 5:(0,1),
% the coefficient page slot is (ia-1)*5 + ib, and each term is the two batched
% page products
%   R = pagemtimes(T,'none',P2,'transpose');  R = pagemtimes(P1,R).
%
% OUTPUT
%   Kb : sparse NxN real bending operator (bending only).
%   Vb : raw nonzero values if requested.

if nargin < 2 || isempty(t_e), t_e = 1; end
if nargin < 3 || isempty(output), output = 'cpu_sparse'; end
[Kb, Vb] = harmonic_kd_form_curved(H, t_e, output);
end
