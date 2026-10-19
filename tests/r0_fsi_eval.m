function [J, dJ] = r0_fsi_eval(ctx, t, R0op)
% R0_FSI_EVAL  r0_pilot_eval with fluid loading: the state solve uses
%   A = K_D(free,free) + R0op.iwZ,   R0op.iwZ = i w s Z_plus (dense, complex symmetric, free dofs;
%   tests/r0_fsi_export_z.py), and J = w^2 ||factor U||^2 as in the pilot.
% Z does not depend on the thickness, so dA/dt_e = dK_D/dt_e as in the pilot and, because A is
% complex symmetric, the adjoint uses the same LU (as in r0_pilot_eval).  Dense LU (2952 dofs).
t = t(:);
ne = ctx.ne; lm = ctx.lm; n = size(lm,2); ndof = ctx.ndof;
cK = 1 + 1i*ctx.w*ctx.beta;
cM = 1i*ctx.w*ctx.alpha - ctx.w^2;
EmatK = ctx.Bke.*(t.^3) + ctx.Ske.*t;
EmatM = ctx.Me.*t;
KD = sparse(ctx.Rows(:), ctx.Cols(:), cK*EmatK(:) + cM*EmatM(:), ndof, ndof);
A = full(KD(ctx.free,ctx.free)) + R0op.iwZ;
[L,U,P] = lu(A);
u = zeros(ndof,1);
u(ctx.free) = U\(L\(P*ctx.F(ctx.free)));
R0U = R0op.factor'*(R0op.factor*u);
J = ctx.w^2 * real(u'*R0U);
if nargout < 2, return; end
a = 2*ctx.w^2 * conj(R0U);
z = zeros(ndof,1);
z(ctx.free) = U\(L\(P*a(ctx.free)));
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
