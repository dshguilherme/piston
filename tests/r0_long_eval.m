function [J, dJ, u] = r0_long_eval(ctx, t, R0op)
% R0_LONG_EVAL  r0_pilot_eval (= articleB_r0_eval(ctx,t,'quadvel',R0op,'consistent'), one LU per
% evaluation, no adjoint when only J is requested) plus a one-entry cache of the sparse LU and
% the state u keyed on the design t (exact equality) and on ctx.cachekey.
% Why: in GCMMA the accepted trial point of outer iteration k (already evaluated, J only, in the
% inner loop) is the point where outer iteration k+1 needs J and the gradient; the cache saves
% that repeated LU. It also lets several factors (exact, sinc, ERP) be evaluated, with gradients,
% at the same design for the cost of one LU. The arithmetic is the same as r0_pilot_eval line by
% line, so J and dJ are bitwise identical (checked in r0_long_runs('validate')).
persistent C
t = t(:);
if isfield(ctx,'cachekey'), key = ctx.cachekey; else, key = sprintf('%d_%.17g',ctx.ndof,ctx.w); end
ndof = ctx.ndof;
cK = 1 + 1i*ctx.w*ctx.beta;
cM = 1i*ctx.w*ctx.alpha - ctx.w^2;
if isempty(C) || ~strcmp(C.key,key) || numel(C.t)~=numel(t) || ~isequal(C.t,t)
    EmatK = ctx.Bke.*(t.^3) + ctx.Ske.*t;
    EmatM = ctx.Me.*t;
    KD = sparse(ctx.Rows(:), ctx.Cols(:), cK*EmatK(:) + cM*EmatM(:), ndof, ndof);
    A = KD(ctx.free,ctx.free);
    [L,U,P,Q] = lu(A);
    uu = zeros(ndof,1);
    uu(ctx.free) = Q*(U\(L\(P*ctx.F(ctx.free))));
    C = struct('key',key,'t',t,'L',L,'U',U,'P',P,'Q',Q,'u',uu);
end
u = C.u;
R0U = R0op.factor'*(R0op.factor*u);
J = ctx.w^2 * real(u'*R0U);
if nargout < 2, return; end
ne = ctx.ne; lm = ctx.lm; n = size(lm,2);
a = 2*ctx.w^2 * conj(R0U);
z = zeros(ndof,1);
z(ctx.free) = C.Q*(C.U\(C.L\(C.P*a(ctx.free))));
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
