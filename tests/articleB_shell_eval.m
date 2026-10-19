function [J, dJ, aux] = articleB_shell_eval(ctx, t, metric)
% ARTICLEB_SHELL_EVAL  Harmonic objective + element-wise sensitivity for the
% curved KL shell (vectorised), element matrices from articleB_build_ctx.
%
%   KD = (1+jw*beta)*Ks + (jw*alpha - w^2)*M
%   Ks(t) = sum_e (t_e^3 Bke_e + t_e Ske_e)   (membrane t, bending t^3)
%   M(t)  = sum_e  t_e Me_e
%   dKD_e = (1+jw*beta)*(3 t_e^2 Bke_e + Ske_e) + (jw*alpha - w^2)*Me_e
%
% KD is COMPLEX SYMMETRIC, so the state is its own adjoint for the linear
% pairings (dyncomp, inputpower) and only 'quadvel' needs one extra solve.
% Sensitivities FD-verified to ~1e-5 (tests/articleB_fdcheck.m).

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
        J = -0.5*ctx.w*imag(conj(ctx.F).'*u);     % = 0.5*Re(F^H V), V=jw U
        dKq = elquad(Ue,Ue,dKDe,n);
        dJ = 0.5*ctx.w*imag(dKq);
    case 'quadvel'
        R0 = ctx.R0; if isempty(R0), R0 = speye(ndof); end
        J = ctx.w^2 * real(u'*(R0*u));
        a = 2*ctx.w^2 * conj(R0*u);
        z = zeros(ndof,1);
        z(ctx.free) = KD(ctx.free,ctx.free) \ a(ctx.free);
        Ze = z(lm);
        dKq = elquad(Ze,Ue,dKDe,n);
        dJ = -real(dKq);
    otherwise
        error('articleB_shell_eval: unknown metric %s', metric);
end
aux = struct('u',u,'g',g,'Ks',Ks,'M',M);
end

function q = elquad(Ae,Be,D,n)
% q(e) = sum_{a,b} Ae(e,a) * D(e,(a-1)*n+b) * Be(e,b)   (A.' D B per element)
ne = size(Ae,1); q = zeros(ne,1);
for a=1:n
    va = Ae(:,a);
    for b=1:n
        q = q + va .* D(:,(a-1)*n+b) .* Be(:,b);
    end
end
end
