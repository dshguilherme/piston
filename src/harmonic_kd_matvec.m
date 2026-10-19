function [y, t] = harmonic_kd_matvec(H, t_e, x, w, alpha, beta)
% HARMONIC_KD_MATVEC Matrix-free application of the harmonic K_D operator.
%
%   y = K_D(t, w) * x ,  K_D = cK*(t^3 Kb) + cM*M(t)
%   cK = 1 + j*w*beta ,  cM = j*w*alpha - w^2
%
% The trial vector is gathered once over the sparsity block
% (Xb(H.mask) = x(H.cols)) and each of the nine metric bending terms plus the
% mass term is applied as two batched page products in the matrix-vector
% direction (the transpose of the forward formation):
%   S = pagemtimes(pagemtimes(P1,'transpose',Xb,'none'), P2);  y += sum(T.*S)
% No sparse matrix is ever formed.

on_gpu = strcmp(H.device, 'gpu');
sync = @() [];
if on_gpu
    dev = gpuDevice; sync = @() wait(dev);
end

t0 = tic;
x = cast(x, H.precision);
if on_gpu, x = gpuArray(x); end
sync(); t.upload = toc(t0);

t0 = tic;
Xb = zeros(H.nnmax(1), H.nnmax(2), H.N, 'like', x);
Xb(H.mask) = x(H.cols);

tq = cast(t_e, H.precision);
if on_gpu, tq = gpuArray(tq); end
tqp = tq(H.e1, H.e2);
s1 = tqp(H.LIN);
s3 = s1 .^ 3;
cK = 1 + 1i * w * beta;
cM = 1i * w * alpha - w^2;

acc = zeros(1, 1, H.N, 'like', x);
comps = [2 0; 1 1; 0 2];
for ia = 1:3
    for ib = 1:3
        C = H.bendKp{(ia-1)*3 + ib};
        if ~any(C(:)), continue; end
        i = comps(ia,1); j = comps(ia,2);
        k = comps(ib,1); l = comps(ib,2);
        P1 = H.WBe{1}{i+1, k+1};
        P2 = H.WBe{2}{j+1, l+1};
        T = cK * (s3 .* C);
        acc = acc + contract(T, P1, P2, Xb);
    end
end
Tm = cM * (s1 .* H.M0p);
acc = acc + contract(Tm, H.WBe{1}{1,1}, H.WBe{2}{1,1}, Xb);

y = acc(:);
sync(); t.kernel = toc(t0);
t.total = t.upload + t.kernel;
end

function out = contract(T, P1, P2, Xb)
% One sum-factorized term in the MATRIX-VECTOR direction:
%   S(q1,q2) = sum_{j1,j2} P1(j1,q1) Xb(j1,j2) P2(j2,q2)
%   out      = sum_{q1,q2} T(q1,q2) S(q1,q2)
M1  = pagemtimes(P1, 'transpose', Xb, 'none');   % [nq1 x nn1] x [nn1 x nn2]
S   = pagemtimes(M1, P2);                        % [nq1 x nn2] x [nn2 x nq2]
out = sum(sum(T .* S, 1), 2);                    % [1 x 1 x N]
end
