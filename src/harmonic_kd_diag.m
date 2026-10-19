function [d, t] = harmonic_kd_diag(H, t_e, w, alpha, beta)
% HARMONIC_KD_DIAG Exact diagonal of the harmonic K_D operator (matrix-free).
%
%   d = diag(K_D) ,  K_D = cK*(t^3 Kb) + cM*M(t)
%
% Row n contributes only through the local entry in which the trial dof equals
% the test dof, i.e. at the local neighbour index H.lself{d}(n_d); therefore
%     d(n) = sum_{alpha,beta} sum_{q1,q2} P1(ls1,q1,n) T(q1,q2,n) P2(ls2,q2,n)
% with the same page pairs and coefficient pages as harmonic_kd_form. The
% diagonal is used as the Jacobi preconditioner; it is recomputed per iteration
% because it depends on the thickness.

t0 = tic;
tq = cast(t_e, H.precision);
tqp = tq(H.e1, H.e2);
s1 = tqp(H.LIN);
s3 = s1 .^ 3;
cK = 1 + 1i * w * beta;
cM = 1i * w * alpha - w^2;
N = H.N;
ls1 = H.lself{1}; ls2 = H.lself{2};          % [N x 1]

acc = zeros(1, 1, N, 'like', s1);
comps = [2 0; 1 1; 0 2];
for ia = 1:3
    for ib = 1:3
        C = H.bendKp{(ia-1)*3 + ib};
        if ~any(C(:)), continue; end
        i = comps(ia,1); j = comps(ia,2);
        k = comps(ib,1); l = comps(ib,2);
        T = cK * (s3 .* C);
        [A1, A2] = self_entries(H, i+1, k+1, j+1, l+1, ls1, ls2);
        Y = pagemtimes(T, 'none', A2, 'transpose');   % [nq1 x 1 x N]
        acc = acc + pagemtimes(A1, Y);                % [1 x 1 x N]
    end
end
Tm = cM * (s1 .* H.M0p);
[A1, A2] = self_entries(H, 1, 1, 1, 1, ls1, ls2);
Y = pagemtimes(Tm, 'none', A2, 'transpose');
acc = acc + pagemtimes(A1, Y);

d = acc(:);
t = toc(t0);
end

function [A1, A2] = self_entries(H, k1, k2, l1, l2, ls1, ls2)
% extract page entries at the self neighbour index: [1 x nq x N]
P1 = H.WBe{1}{k1, k2}; P2 = H.WBe{2}{l1, l2};
nq1 = size(P1,2); nq2 = size(P2,2); N = H.N;
IL1 = reshape(ls1, 1, 1, N) + reshape(0:nq1-1, 1, nq1, 1) * H.nnmax(1) + ...
      reshape((0:N-1) * H.nnmax(1) * nq1, 1, 1, N);
IL2 = reshape(ls2, 1, 1, N) + reshape(0:nq2-1, 1, nq2, 1) * H.nnmax(2) + ...
      reshape((0:N-1) * H.nnmax(2) * nq2, 1, 1, N);
A1 = reshape(P1(IL1), 1, nq1, N);
A2 = reshape(P2(IL2), 1, nq2, N);
end
