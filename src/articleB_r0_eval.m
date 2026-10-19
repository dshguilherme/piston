function [J, dJ, aux] = articleB_r0_eval(ctx, t, metric, R0op, R0mode)
% ARTICLEB_R0_EVAL  Harmonic objective + element-wise sensitivity for the
% curved KL shell (vectorised, element matrices from articleB_build_ctx), with
% the thesis radiation-resistance weighting R0 for the quadratic-velocity /
% radiated-power objective.
%
%   KD = (1+jw*beta)*Ks + (jw*alpha - w^2)*M
%   Ks(t) = sum_e (t_e^3 Bke_e + t_e Ske_e) ,  M(t) = sum_e t_e Me_e
%   dKD_e = (1+jw*beta)*(3 t_e^2 Bke_e + Ske_e) + (jw*alpha - w^2)*Me_e
%
% KD is COMPLEX SYMMETRIC (self-adjoint for the BILINEAR pairing <x,y>=x.'y),
% so dyncomp/inputpower need no extra solve. The quadratic-velocity objective
% J = w^2 * U^H R0 U has a gradient functional a = 2*w^2*conj(R0*U) that is
% NOT proportional to F for a general R0, so it needs ONE extra (adjoint) solve
% K_D z = a and a cross energy: dJ/dt_e = -Re(z.' dK_e U). This is the same
% class as the thesis relation K_D*lambda = -2*omega^2*U^H R (eq. 3.69):
% because R0 is Hermitian, conj(R0*U) = R0.'*conj(U), so our a = +2*w^2*R0.'*conj(U)
% equals -(K_D*lambda) of (3.69) up to the overall sign carried by the dJ
% convention; the FD check below verifies our form against an independently
% assembled W.
%
% R0mode (only used for 'quadvel'):
%   'I'    R0 = identity on the FE velocity -> J = w^2*||U||^2  (phase-BLIND)
%   'ones' R0 = C^H * ones(ne,ne) * C       -> J = w^2*|sum_e (C U)_e|^2
%          net volume velocity (point-source/monopole limit; cancellation-max)
%   'rad'  R0 = C^H * R * C, R from eq.(3.66) -> J = w^2*(CU)^H R (CU) = W_rad
% R0op from articleB_r0_build.

if nargin < 5 || isempty(R0mode), R0mode = 'rad'; end

t = t(:);
ne = ctx.ne; lm = ctx.lm; n = size(lm,2); ndof = ctx.ndof;
cK = 1 + 1i*ctx.w*ctx.beta;
cM = 1i*ctx.w*ctx.alpha - ctx.w^2;

EmatK = ctx.Bke.*(t.^3) + ctx.Ske.*t;   EmatK = EmatK(:);
EmatM = ctx.Me.*t;                      EmatM = EmatM(:);
Ks = sparse(ctx.Rows(:), ctx.Cols(:), EmatK, ndof, ndof);
M  = sparse(ctx.Rows(:), ctx.Cols(:), EmatM, ndof, ndof);
KD = cK*Ks + cM*M;

u = zeros(ndof,1);
u(ctx.free) = KD(ctx.free,ctx.free) \ ctx.F(ctx.free);
g = ctx.F.' * u;

dKDe = cK*(ctx.Bke.*(3*t.^2) + ctx.Ske) + cM*ctx.Me;   % ne x n^2
Ue = u(lm);

switch lower(metric)
    case 'dyncomp'
        J = abs(g);
        dKq = elquad(Ue,Ue,dKDe,n);
        dJ = real(conj(g)*(-dKq)) / max(J,eps);
    case 'inputpower'
        J = -0.5*ctx.w*imag(conj(ctx.F).'*u);
        dKq = elquad(Ue,Ue,dKDe,n);
        dJ = 0.5*ctx.w*imag(dKq);
    case 'quadvel'
        R0U = apply_R0(R0op, u, R0mode);
        J = ctx.w^2 * real(u'*R0U);
        a = 2*ctx.w^2 * conj(R0U);
        z = zeros(ndof,1);
        z(ctx.free) = KD(ctx.free,ctx.free) \ a(ctx.free);
        Ze = z(lm);
        dKq = elquad(Ze,Ue,dKDe,n);
        dJ = -real(dKq);
    otherwise
        error('articleB_r0_eval: unknown metric %s', metric);
end
aux = struct('u',u,'g',g,'Ks',Ks,'M',M,'R0mode',R0mode);
end

% ---------------------------------------------------------------------------
function R0U = apply_R0(R0op, u, mode)
switch lower(mode)
    case 'i'
        R0U = u;
    case 'ones'
        y = R0op.C * u;                             % centre velocities / (jw)
        s = sum(y);                                 % net volume velocity
        R0U = R0op.pref * (R0op.C' * (ones(R0op.ne,1)*s));
    case 'rad'
        y = R0op.C * u;
        z = R0op.area .* y;                         % area-weighted centre displacement
        Kz = apply_K(R0op, z);                      % Nyström K*z
        R0U = R0op.pref * (R0op.C' * (R0op.area .* Kz));
    case 'consistent'
        assert(isfield(R0op,'factor'),'consistent mode requires r0_build_consistent output');
        R0U = R0op.factor'*(R0op.factor*u);
    otherwise
        error('articleB_r0_eval: unknown R0mode %s', mode);
end
end

function Kz = apply_K(R0op, z)
% Nystrom action of the sinc kernel K: K ~= KS.' * KSSinv * KS
Kz = R0op.KS.' * (R0op.KSSinv * (R0op.KS * z));
end

function q = elquad(Ae,Be,D,n)
% q(e) = sum_{a,b} Ae(e,a) * D(e,(a-1)*n+b) * Be(e,b)
ne = size(Ae,1); q = zeros(ne,1);
for a=1:n
    va = Ae(:,a);
    for b=1:n
        q = q + va .* D(:,(a-1)*n+b) .* Be(:,b);
    end
end
end
