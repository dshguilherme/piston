function y = harmonic_kd_matvec_curved(H, t_e, x)
% HARMONIC_KD_MATVEC_CURVED Matrix-free application of the curved bending
% operator Kb of harmonic_kd_form_curved.m:
%   y = Kb * x ,  Kb = int t^3 D det(J) Phi' G Phi dx  (25-slot curved page).
% The trial vector is gathered once over the sparsity block and each of the 25
% terms is applied as two batched page products in the MATRIX-VECTOR direction
% (the transpose of the forward formation), exactly as harmonic_kd_matvec.m.
% Bilinear pairing only (x.' not x').

if nargin < 2 || isempty(t_e), t_e = 1; end
N = H.N;
tq = t_e(H.e1, H.e2);
s3 = tq(H.LIN);

Xb = zeros(H.nnmax(1), H.nnmax(2), N);
Xb(H.mask) = x(H.cols);

acc = zeros(1, 1, N);
comps = H.comps;
for ia = 1:5
    for ib = 1:5
        C = H.bendKp{(ia-1)*5 + ib};
        if ~any(C(:)), continue; end
        i = comps(ia,1); j = comps(ia,2);
        k = comps(ib,1); l = comps(ib,2);
        P1 = H.WBe{1}{i+1, k+1};
        P2 = H.WBe{2}{j+1, l+1};
        T = s3 .* C;
        acc = acc + contract(T, P1, P2, Xb);
    end
end
y = acc(:);
end

function out = contract(T, P1, P2, Xb)
% S(q1,q2) = sum_{j1,j2} P1(j1,q1) Xb(j1,j2) P2(j2,q2); out = sum T.*S
M1  = pagemtimes(P1, 'transpose', Xb, 'none');
SS  = pagemtimes(M1, P2);
out = sum(sum(T .* SS, 1), 2);
end
