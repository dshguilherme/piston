function [KD, t] = harmonic_kd_form(H, t_e, w, alpha, beta, output)
% HARMONIC_KD_FORM Per-iteration formation of the complex-symmetric dynamic
% stiffness operator (flat Kirchhoff-Love plate in physical coordinates,
% decoupled transverse field, N-DOF scalar operator used by the T1 suite).
%
%   KD = cK * (t^3 Kb) + cM * M(t)
%   cK = 1 + j*w*beta ,  cM = j*w*alpha - w^2
%
% Bending runs on the metric-carrying weighted-quadrature pages: the physical
% Hessian is expanded in the parametric components h = [w_xi xi; w_xi eta;
% w_eta eta], giving nine (alpha,beta) terms with component map
%   xi xi = (2,0) , xi eta = (1,1) , eta eta = (0,2)
% and, for alpha = (i,j), beta = (k,l),
%   direction 1 uses rule (i,k) , direction 2 uses rule (j,l)
% with coefficient page H.bendK(:,:,(ia-1)*3+ib). Each term is the two batched
% page products of wq_form.m:
%   R = pagemtimes(T,'none',P2,'transpose'); R = pagemtimes(P1,R)
%
% output: 'cpu_sparse' (default) -> sparse KD on the host
%         'values'               -> raw nonzero values (Kb, Mm)

if nargin < 6 || isempty(output), output = 'cpu_sparse'; end

on_gpu = strcmp(H.device, 'gpu');
sync = @() [];
if on_gpu
    dev = gpuDevice; sync = @() wait(dev);
end

t0 = tic;
x = cast(t_e, H.precision);
if on_gpu, x = gpuArray(x); end
sync(); t.upload = toc(t0);

t0 = tic;
tq = x(H.e1, H.e2);                 % thickness at WQ points  [nq1 x nq2]
s1 = tq(H.LIN);                     % t   pages [nq1 x nq2 x N]
s3 = s1 .^ 3;                       % t^3 pages
cK = 1 + 1i * w * beta;
cM = 1i * w * alpha - w^2;

% ---- bending family (t^3), nine metric terms ----------------------------
comps = [2 0; 1 1; 0 2];
V = [];
for ia = 1:3
    for ib = 1:3
        C = H.bendKp{(ia-1)*3 + ib};
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

% ---- mass family (t), WQ type-00 pages ----------------------------------
Tm = s1 .* H.M0p;
Rm = pagemtimes(Tm, 'none', H.WBe{2}{1,1}, 'transpose');
Rm = pagemtimes(H.WBe{1}{1,1}, Rm);
Vm = Rm(H.mask);
sync(); t.kernel = toc(t0);

switch output
    case 'values'
        if on_gpu, Vb = gather(Vb); Vm = gather(Vm); end
        KD = struct('Kb', Vb, 'Mm', Vm);
        t.total = t.upload + t.kernel;
        return;
    otherwise
        t0 = tic;
        if on_gpu, Vb = gather(Vb); Vm = gather(Vm); end
        t.download = toc(t0);
        t0 = tic;
        K = sparse(H.rows, H.cols, cK * double(Vb) + cM * double(Vm), H.N, H.N);
        KD = 0.5 * (K + K.');
        t.sparse = toc(t0);
end
t.extract = 0;
t.total = t.upload + t.kernel + t.download + t.sparse;
end
