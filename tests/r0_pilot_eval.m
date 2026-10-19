function [J, dJ] = r0_pilot_eval(ctx, t, R0op)
% R0_PILOT_EVAL  articleB_r0_eval(ctx,t,'quadvel',R0op,'consistent') with the same
% mathematics (state solve, J = w^2 ||factor U||^2, adjoint K_D z = a, element
% sensitivities) and two purely computational changes for the regime pilot, where the
% machine is CPU-contended:
%   (1) one sparse LU of K_D(free,free) is reused for the state and the adjoint solve
%       (K_D is complex symmetric, so the adjoint uses the same matrix);
%   (2) the adjoint and sensitivities are skipped when only J is requested (nargout == 1),
%       which is what the GCMMA inner loop needs.
% Agreement with articleB_r0_eval is checked in r0_pilot_runs.m stage 'ops'.
t = t(:);
ne = ctx.ne; lm = ctx.lm; n = size(lm,2); ndof = ctx.ndof;
cK = 1 + 1i*ctx.w*ctx.beta;
cM = 1i*ctx.w*ctx.alpha - ctx.w^2;
EmatK = ctx.Bke.*(t.^3) + ctx.Ske.*t;
EmatM = ctx.Me.*t;
KD = sparse(ctx.Rows(:), ctx.Cols(:), cK*EmatK(:) + cM*EmatM(:), ndof, ndof);
A = KD(ctx.free,ctx.free);
[L,U,P,Q] = lu(A);
u = zeros(ndof,1);
u(ctx.free) = Q*(U\(L\(P*ctx.F(ctx.free))));
R0U = R0op.factor'*(R0op.factor*u);
J = ctx.w^2 * real(u'*R0U);
if nargout < 2, return; end
a = 2*ctx.w^2 * conj(R0U);
z = zeros(ndof,1);
z(ctx.free) = Q*(U\(L\(P*a(ctx.free))));
dKDe = cK*(ctx.Bke.*(3*t.^2) + ctx.Ske) + cM*ctx.Me;
Ue = u(lm); Ze = z(lm);
q = zeros(ne,1);
for ia=1:n
    va = Ze(:,ia);
    for ib=1:n
        q = q + va .* dKDe(:,(ia-1)*n+ib) .* Ue(:,ib);
    end
end
dJ = -real(q);
end
