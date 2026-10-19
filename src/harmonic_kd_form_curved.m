function [Kb, Vb] = harmonic_kd_form_curved(H, t_e, output)
% HARMONIC_KD_FORM_CURVED Assembly of the curved (5-component) bending operator
% from the WQ pages of harmonic_kd_setup_curved.m.
%
%   Kb = int t(x)^3 * D * det(J) * Phi' G Phi dx ,
%   Phi = [w_xixi; w_xieta; w_etaeta; w_xi; w_eta]  (parametric second/first derivs)
%
% The five Phi components have the map (see quadrule_bend_curved.m)
%     1: xi xi = (2,0)  2: xi eta = (1,1)  3: eta eta = (0,2)
%     4: xi    = (1,0)  5: eta    = (0,1)
% and, for test component alpha=(i,j) and trial component beta=(k,l), direction 1
% uses WQ rule (i,k) and direction 2 uses rule (j,l); the coefficient page slot
% is (ia-1)*5 + ib. Each term is the two batched page products of wq_form.m,
%     R = pagemtimes(T,'none',P2,'transpose');  R = pagemtimes(P1,R)
% exactly as harmonic_kd_form.m, just over 5x5 = 25 slots instead of 3x3 = 9.
%
% OUTPUT
%   Kb : sparse NxN real bending operator (bending ONLY; mass/membrane are the
%        affine-safe pages of harmonic_kd_setup.m and are not part of this file).
%   Vb : the raw nonzero values (before sparse assembly), if requested.
%
% t_e : thickness per element ([nel1 x nel2]); default = 1 everywhere.

if nargin < 2 || isempty(t_e), t_e = 1; end
if nargin < 3 || isempty(output), output = 'cpu_sparse'; end

N = H.N;
tq = t_e(H.e1, H.e2);          % thickness at WQ points [nq1 x nq2] (scalar -> all)
s3 = tq(H.LIN);                % t^3 pages [nq1 x nq2 x N]

comps = H.comps;
V = [];
for ia = 1:5
    for ib = 1:5
        C = H.bendKp{(ia-1)*5 + ib};
        if ~any(C(:)), continue; end
        i = comps(ia,1); j = comps(ia,2);
        k = comps(ib,1); l = comps(ib,2);
        P1 = H.WBe{1}{i+1, k+1};
        P2 = H.WBe{2}{j+1, l+1};
        T = s3 .* C;
        R = pagemtimes(T, 'none', P2, 'transpose');
        R = pagemtimes(P1, R);
        if isempty(V), V = R; else, V = V + R; end
    end
end
Vb = V(H.mask);

if strcmpi(output, 'values')
    Kb = Vb;
    return;
end
Kb = sparse(H.rows, H.cols, double(Vb), N, N);
Kb = 0.5 * (Kb + Kb.');
end
