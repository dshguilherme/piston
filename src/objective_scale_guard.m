function [f0s, df0s, g] = objective_scale_guard(J, dJ, p0, opt)
% OBJECTIVE_SCALE_GUARD  Monotone, positivity-preserving scaling of a harmonic
% objective for GCMMA.  Repairs the dB-scaling DESCENT INVERSION that occurs
% when the raw objective J falls below the dB threshold Jthr = 1e-10.
%
% ---------------------------------------------------------------------------
% THE DEFECT (legacy driver scaling, tests/articleB_obstacle_run.m)
%   Jdb  = Abar + dbfac*log10(J)          % Abar=100, dbfac=10
%   f0s  = 100*Jdb/p0                     % normalised objective for GCMMA
%   df0s = 100*(dbfac/ln10)*dJ/(J*p0)     % its gradient
%   with the reference normaliser p0 = Jdb(J0) frozen at the first iteration.
% Jdb is strictly increasing in J, but it crosses ZERO at
%   Jthr = 10^(-Abar/dbfac) = 1e-10.
% For J0 < Jthr the reference p0 = Jdb(J0) is NEGATIVE, so BOTH f0s and df0s
% carry the factor 1/p0 < 0.  The pair (f0s, df0s) is then monotonically
% DECREASING in J while its gradient points the wrong way: GCMMA, which only
% sees df0s, descends f0s and therefore CLIMBS the true objective J.  Below the
% threshold the scaling is not merely badly conditioned, it is INVERTED.
% (Observed on the cylindrical panel 'quadvel' objective, J0 ~ 8.4e-13 =
%  -120.7 dB < 0: the design stayed uniform, volume drifted to 0.44.)
%
% ---------------------------------------------------------------------------
% THE FIX
%   psi(J) = Abar + dbfac*log10(J)                             , J >= Jthr
%   psi(J) = psi(Jthr) + s0*(J - Jthr)/Jthr ,  s0 = dbfac/ln10 , J <  Jthr
% This map is C1 (both value and slope agree at Jthr), STRICTLY INCREASING for
% every J>0 because psi'(J) = s0/J >= s0/Jthr > 0, and it replaces the unbounded
% 1/J gradient blow-up below the threshold by the constant tangent slope
% s0/Jthr.  Above the threshold psi is the legacy DdB map VERBATIM, so cases with
% J >> 1e-10 are untouched.
%
% The objective is normalised by the POSITIVE reference p0 = |psi(J0)|.  Since
% psi(J0) = Jdb(J0) > 0 whenever J0 >= Jthr, this is IDENTICAL to the legacy p0
% on the well-scaled branch; it only changes (to a positive number) on the
% inverted branch, which is exactly the defect.  With p0 > 0 the factor 1/p0 is
% positive, so
%   f0s  = 100*psi(J)/p0
%   df0s = 100*psi'(J)*dJ/p0
% are both monotone increasing in J and the gradient sign matches descent on J.
% (A constant shift of psi leaves GCMMA's subproblem unchanged -- the model
%  offset r0 absorbs it -- so the C1 continuation is free of side effects.)
%
% ---------------------------------------------------------------------------
% SIGNATURE
%   [f0s, df0s, g] = objective_scale_guard(J, dJ, p0, opt)
%     J    scalar      raw harmonic objective at the current design (J > 0)
%     dJ   (ne x 1)    dJ/dx at the current design
%     p0   scalar or []  frozen normaliser: pass [] on the FIRST iteration to
%                      have it computed here as |psi(J0)|; pass g.p0 afterwards
%                      so the reference stays frozen across iterations.
%     opt  (optional)  .Abar (100) .dbfac (10) .psimin (1e-6)
%   f0s  scalar        scaled objective to hand to gcmma_solver
%   df0s (ne x 1)      scaled gradient to hand to gcmma_solver
%   g    struct        .psi .dpsi .p0 .branch('log'|'linear') .Jthr .Jdb_legacy
%
% The companion opt.eval given to gcmma_solver must return the SAME scaled
% objective, i.e.
%   opt.eval = @(tq) deal(scaled_obj(tq), vol(tq))
% where scaled_obj calls objective_scale_guard with the SAME frozen p0 (see
% tests/articleB_obstacle_hi_run.m, which wires this through).
%
% C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src/

if nargin < 4 || isempty(opt), opt = struct(); end
if ~isfield(opt,'Abar'),   opt.Abar   = 100;  end
if ~isfield(opt,'dbfac'),  opt.dbfac  = 10;   end
if ~isfield(opt,'psimin'), opt.psimin = 1e-6; end

Abar = opt.Abar; dbfac = opt.dbfac;
s0   = dbfac/log(10);                 % slope of the dB map wrt J (constant factor /J)
Jthr = 10^(-Abar/dbfac);              % = 1e-10 : where the legacy dB map crosses 0

J = J(:);

% --- strictly increasing, C1 map psi(J) -----------------------------------
psi  = Abar + dbfac*log10(J);
dpsi = s0 ./ J;                       % > 0
below = J < Jthr;
if any(below)
    psi(below)  = (Abar + dbfac*log10(Jthr)) + s0*(J(below)-Jthr)/Jthr;  % = s0*(J-Jthr)/Jthr
    dpsi(below) = s0/Jthr;            % tangent slope: monotone, no 1/J blow-up
end
Jdb_legacy = Abar + dbfac*log10(J);   % for reporting/onset comparison

% --- positivity-preserving frozen reference --------------------------------
if nargin < 3 || isempty(p0)
    p0 = abs(psi(1));                 % |psi(J0)| : equals legacy p0 above threshold
end
p0 = max(p0, opt.psimin);

% --- scaled objective + gradient (factor 1/p0 now provably > 0) ------------
f0s  = 100*psi(1)/p0;
df0s = 100*(dpsi(1)/p0) * dJ(:);

branch = 'log'; if J(1) < Jthr, branch = 'linear'; end
g = struct('psi',psi(1),'dpsi',dpsi(1),'p0',p0,'branch',branch, ...
           'Jthr',Jthr,'Jdb_legacy',Jdb_legacy(1));
end
