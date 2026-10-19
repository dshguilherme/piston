% FD gate for the two NEW harmonic sensitivities in src/harmonic_kd_sens.m:
%   'inputpower'  active input power   J2 = -0.5*w*Im(conj(F).'U)  (TRICK, no extra solve)
%   'quadvel'     quadratic velocity   J3 = w^2 * U^H R0 U        (ONE extra solve, cross energy)
%
% VERIFIED HERE (raw numbers printed):
%  (1) energy partition: sum_e Eb  == U.'*Kb*U ,  sum_e Em  == U.'*M*U   (single field)
%      and, for quadvel,    sum_e Ebz == z.'*Kb*U ,  sum_e Emz == z.'*M*U (CROSS term)
%  (2) UNIFORM-perturbation central FD of the TOTAL gradient (all elements +h/-h)
%  (3) single-element central FD on random elements: max & median RELATIVE error
%  (4) element-wise cosine diagnostic (documented, not asserted) -- the per-element
%      split does NOT match element-wise FD (knot-boundary WQ points vs per-element
%      thickness); this is the INHERITED limitation of the dyncomp mode.
%  (5) quadvel with a mass-weighted R0 ('mass') checked the same way.
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');

E = 210e9; RHO = 7850; NU = 0.3;
pd = cantilever_beam(1.0, 1.0);
md = struct('degree',[3 3],'regularity',[2 2],'nsub',[8 4],'nquad',[4 4]);
[gm, mh, sp] = buildSpaces(pd, md);
H = harmonic_kd_setup(mh, sp, gm, E, NU, RHO, 'cpu', 'double', 'plate');
t0 = 2e-3 * ones(8, 4);
rng(1);
F  = randn(H.N, 1);                 % REAL loading
al = 25; bet = 0;
h  = 1e-8 * mean(t0(:));
freqs = [30 60 120];

Vv = harmonic_kd_form(H, ones(8,4), 0,0,0, 'values');
KB = sparse(H.rows, H.cols, Vv.Kb, H.N, H.N);   % unit-thickness bending
MB = sparse(H.rows, H.cols, Vv.Mm, H.N, H.N);   % unit-thickness mass

allok = true;
modes = {'inputpower','quadvel'};
for imode = 1:numel(modes)
  mode = modes{imode};
  fprintf('\n==== METRIC %s ====\n', mode);
  for ifr = 1:numel(freqs)
    w = 2*pi*freqs(ifr);
    U = harmonic_kd_form(H, t0, w, al, bet, 'cpu_sparse') \ F;

    if strcmpi(mode,'inputpower')
      [dJ, aux] = harmonic_kd_sens(H, t0, F, U, w, al, bet, 'inputpower');
      Jfun   = @(tt) ipJ(H, tt, w, al, bet, F);
      rEb = abs(sum(aux.Eb(:)) - U.'*KB*U) / abs(U.'*KB*U);
      rEm = abs(sum(aux.Em(:)) - U.'*MB*U) / abs(U.'*MB*U);
      extra = '';
    else
      [dJ, aux] = harmonic_kd_sens(H, t0, F, U, w, al, bet, 'quadvel');   % default R0 = I
      Jfun   = @(tt) qvJ(H, tt, w, al, bet, F, []);
      zt = aux.z.';
      rEb = abs(sum(aux.Ebz(:)) - zt*KB*U) / abs(zt*KB*U);
      rEm = abs(sum(aux.Emz(:)) - zt*MB*U) / abs(zt*MB*U);
      extra = sprintf(' | Im(z.dK.U) rel %+.1e', ...
                      norm(imag(aux.dKq_cross))/max(norm(aux.dKq_cross),eps));
    end

    % (2) uniform FD on the TOTAL
    dJfd_u = (Jfun(t0+h) - Jfun(t0-h)) / (2*h);
    eu = abs(sum(dJ(:)) - dJfd_u) / abs(dJfd_u);

    % (3) random single-element FD
    nsel = 6;
    idx  = randperm(numel(t0), nsel);
    err  = zeros(nsel,1);
    for k = 1:nsel
      tp = t0; tp(idx(k)) = tp(idx(k)) + h;
      tm = t0; tm(idx(k)) = tm(idx(k)) - h;
      dfd = (Jfun(tp) - Jfun(tm)) / (2*h);
      err(k) = abs(dJ(idx(k)) - dfd) / abs(dfd);
    end

    % (4) full element-wise FD for the cosine / mean-ratio diagnostic
    dJfd_full = zeros(size(t0));
    for a1 = 1:8
      for a2 = 1:4
        tp = t0; tp(a1,a2) = tp(a1,a2) + h;
        tm = t0; tm(a1,a2) = tm(a1,a2) - h;
        dJfd_full(a1,a2) = (Jfun(tp) - Jfun(tm)) / (2*h);
      end
    end
    cs = abs(dJ(:).'*dJfd_full(:)) / (norm(dJ(:))*norm(dJfd_full(:)));
    ratio = mean(dJfd_full(:) ./ dJ(:));
    relfull = abs(dJ(:) - dJfd_full(:)) ./ abs(dJfd_full(:));

    fprintf(['f=%3.0f | partition Eb %.1e Em %.1e%s | uniform-FD rel %.2e | ' ...
             'single-elem: max %.2e median %.2e | elemwise-FD rel: max %.2e median %.2e | ' ...
             'cos %.4f | mean FD/an %.4f\n'], ...
            freqs(ifr), rEb, rEm, extra, eu, max(err), median(err), ...
            max(relfull), median(relfull), cs, ratio);
    if eu < 1e-5 && max([rEb rEm]) < 1e-10
      fprintf('   OK\n');
    else
      fprintf('   FAIL\n'); allok = false;
    end
  end
end

% (5) mass-weighted R0 for quadvel
fprintf('\n==== quadvel, R0 = ''mass'' (unit-thickness M0) ====\n');
w = 2*pi*60;
U = harmonic_kd_form(H, t0, w, al, bet, 'cpu_sparse') \ F;
[dJ, aux] = harmonic_kd_sens(H, t0, F, U, w, al, bet, 'quadvel', 'mass');
M0 = aux.R0;
Jfunm = @(tt) qvJ(H, tt, w, al, bet, F, M0);
zt = aux.z.';
rEb = abs(sum(aux.Ebz(:)) - zt*KB*U) / abs(zt*KB*U);
rEm = abs(sum(aux.Emz(:)) - zt*MB*U) / abs(zt*MB*U);
dJfd_u = (Jfunm(t0+h) - Jfunm(t0-h)) / (2*h);
eu = abs(sum(dJ(:)) - dJfd_u) / abs(dJfd_u);
nsel = 6; idx = randperm(numel(t0), nsel); err = zeros(nsel,1);
for k = 1:nsel
  tp = t0; tp(idx(k)) = tp(idx(k)) + h;
  tm = t0; tm(idx(k)) = tm(idx(k)) - h;
  dfd = (Jfunm(tp) - Jfunm(tm)) / (2*h);
  err(k) = abs(dJ(idx(k)) - dfd) / abs(dfd);
end
fprintf('partition Eb %.1e Em %.1e | uniform-FD rel %.2e | single-elem max %.2e median %.2e\n', ...
        rEb, rEm, eu, max(err), median(err));
if eu < 1e-5 && max([rEb rEm]) < 1e-10
  fprintf('   OK\n');
else
  fprintf('   FAIL\n'); allok = false;
end

if allok, disp('SENS_METRICS_ALL_OK'); else, disp('SENS_METRICS_FAIL'); end

% ---------------------------------------------------------------------------
function J = ipJ(H, tt, w, al, bet, F)
U = harmonic_kd_form(H, tt, w, al, bet, 'cpu_sparse') \ F;
J = -0.5 * w * imag(conj(F).' * U);
end
function J = qvJ(H, tt, w, al, bet, F, R0)
U = harmonic_kd_form(H, tt, w, al, bet, 'cpu_sparse') \ F;
if isempty(R0), J = w^2 * (U'*U); else, J = w^2 * (U'*(R0*U)); end
J = real(J);
end
