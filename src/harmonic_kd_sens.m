function [dJ, aux] = harmonic_kd_sens(H, t_e, F, U, w, alpha, beta, metric, R0)
% HARMONIC_KD_SENS Element-wise thickness sensitivity of the harmonic
% objectives -- the "strain-energy trick" of fast_sensitivities.m, adapted to
% the harmonic complex-symmetric operator.
%
% WHY IT WORKS. K_D = cK(t^3 Kb + t Km) + cM M(t) is COMPLEX SYMMETRIC, so it is
% self-adjoint for the BILINEAR pairing <x,y> = x.'*y. With K_D U = F,
%     dU = -K_D^{-1} dK U  =>  F.' dU = -(K_D^{-1} F).' dK U = -U.' dK U ,
% i.e. THE STATE SOLUTION IS ALSO THE ADJOINT: no second (adjoint) solve.
% Because dK/dt_e acts only on element e, U.' dK_e U is the ELEMENT ENERGY of the
% state, so the whole sensitivity is a POST-PROCESSING of U.
%
% GENERAL (bilinear-adjoint) FORM. For an objective J(U) with gradient
% functional a := dJ/dU (so that dJ = Re(a.' dU)) and K_D z = a,
%     dJ/dt_e = -Re( z.' dK_e U ) ,   dK_e = cK*3 t_e^2 Kb_e + cM M_e .
% When a is proportional to F (a = c*F) then z = c*U and NO extra solve is
% needed -- the "trademark trick".
%
% METRICS IMPLEMENTED (metric argument, case-insensitive):
%
%  'dyncomp'  dynamic compliance  J = |F.'U| = |g|,  g = F.'U.
%             dJ/dt_e = Re( conj(g) * (-U.' dK_e U) ) / |g|   (a = conj(g)/|g|*F).
%
%  'inputpower'  ACTIVE INPUT POWER. With V = j*w*U,
%             J2 = 0.5*Re(F^H V) = -0.5*w*Im(conj(F).' U)  (= -0.5*w*Im(g) for
%             REAL loading F). For REAL F, a is proportional to F, so THE TRICK
%             APPLIES:
%             dJ2/dt_e = +0.5*w*Im( U.' dK_e U ).
%             This is the SAME curvature/mass field post-processing the dyncomp
%             mode already does, with Im in place of Re(conj(g)*(.)). NO extra
%             solve. Note J1 = |g| and J2 = -0.5*w*Im(g) are the real and
%             imaginary parts of one bilinear pairing g.
%             LIMITATION: the closed form assumes REAL F (conj(F)=F), the case in
%             every test harness here.
%
%  'quadvel'  QUADRATIC VELOCITY. For a design-INDEPENDENT Hermitian weighting R0,
%             J3 = V^H R0 V = w^2 * U^H R0 U .
%             Here a = 2*w^2*conj(R0*U) contains NO F, so THE TRICK DOES NOT
%             APPLY: one EXTRA solve K_D z = a is required (SAME operator), then
%             dJ3/dt_e = -Re( z.' dK_e U ) -- a CROSS energy between the two
%             fields, accumulated over the weighted-quadrature points exactly as
%             the single-field energy is.
%             R0 (9th arg): [] -> identity (J3 = w^2 ||U||^2, the default), a
%             Hermitian MATRIX -> used as-is, or the string 'mass' -> the
%             unit-thickness mass matrix M0 = M(t)/t (mass-weighted velocity).
%             R0 MUST be design-independent: a thickness-dependent R0 would add a
%             product-rule term to dJ3 that this routine does not compute.
%
% The element integrals (bending + mass, single- or cross-field) are read off
% the standard per-element Gauss grid, so no element matrices are formed and no
% page-level assembly is needed. Both metrics share the SAME weighted-quadrature
% / curvature-field machinery (local function elem_bilinear).
%
% aux carries: g, J, Eb, Em (single-field element energies), and, for 'quadvel',
% the adjoint z, the cross energies Ebz/Emz, the gradient functional a and R0.

if nargin < 8 || isempty(metric), metric = 'dyncomp'; end
if nargin < 9, R0 = []; end

cK = 1 + 1i*w*beta;
cM = 1i*w*alpha - w^2;
F = F(:); U = U(:);
Umat = reshape(U, H.n(1), H.n(2));
Ev = H.Ev;

% Single-field element energies of the state (shared machinery).
[Eb, Em] = elem_bilinear(H, Ev, Umat, Umat);
g = F.' * U;  Jg = abs(g);
dKq = cK * 3 * (t_e .^ 2) .* Eb + cM * Em;     % = U.' dK_e U per element

switch lower(metric)
    case 'dyncomp'
        dJ = real(conj(g) * (-dKq)) / Jg;

    case 'inputpower'
        % J2 = 0.5*Re(F^H V) = -0.5*w*Im(conj(F).'U); real F => -0.5*w*Im(g)
        dJ = 0.5 * w * imag(dKq);               % the trick: no extra solve

    case 'quadvel'
        R0 = resolve_R0(H, R0, t_e);
        J3 = w^2 * (U' * (R0*U));               % U^H R0 U
        a  = 2 * w^2 * conj(R0*U);
        KD = harmonic_kd_form(H, t_e, w, alpha, beta, 'cpu_sparse');
        z  = KD \ a;                            % the ONE extra solve
        Zmat = reshape(z, H.n(1), H.n(2));
        [Ebz, Emz] = elem_bilinear(H, Ev, Zmat, Umat);   % cross z.' Kb_e U
        dKq_cross = cK * 3 * (t_e .^ 2) .* Ebz + cM * Emz;
        dJ = -real(dKq_cross);                  % objective is real

    case 'energy'
        % return the raw element energies (single field)
        dJ = struct('Eb', Eb, 'Em', Em, 'g', g);

    otherwise
        error('harmonic_kd_sens: unknown metric ''%s''.', metric);
end

J2d = -0.5 * w * imag(conj(F).' * U);
if strcmpi(metric, 'quadvel')
    aux = struct('g', g, 'J', real(J3), 'Eb', Eb, 'Em', Em, 'mode', lower(metric), ...
                 'z', z, 'Ebz', Ebz, 'Emz', Emz, 'a', a, 'R0', R0, 'dKq_cross', dKq_cross);
elseif strcmpi(metric, 'inputpower')
    aux = struct('g', g, 'J', J2d, 'Eb', Eb, 'Em', Em, 'mode', lower(metric), 'dKq', dKq);
else
    aux = struct('g', g, 'J', Jg, 'Eb', Eb, 'Em', Em, 'mode', lower(metric));
end
end

% ---------------------------------------------------------------------------
function [Eb, Em] = elem_bilinear(H, Ev, Amat, Bmat)
% BILINEAR element energies between two state fields A and B:
%   Eb_e = int_e Dbf*( Axx Bxx + Ayy Byy + nu (Axx Byy + Ayy Bxx)
%                      + 2(1-nu) Axy Bxy ) dO     ('plate')
%   Eb_e = int_e Dbf*( Axx Bxx + Ayy Byy + 2 Axy Bxy ) dO          ('T1')
%   Em_e = RHO * int_e ( A B ) dO
% with the PHYSICAL (metric-mapped) Hessian components, on the per-element Gauss
% grid. A = B = state gives the standard (single-field) energy partition; A = z
% (adjoint), B = U gives the CROSS energy z.' Kb_e U, z.' M_e U.
Dbf = H.bend.Dbf; with_nu = strcmpi(Ev.model, 'plate');
Eb = zeros(Ev.nel1, Ev.nel2);
Em = zeros(Ev.nel1, Ev.nel2);
for el1 = 1:Ev.nel1
    for el2 = 1:Ev.nel2
        Al = Amat(Ev.conn1(:,el1), Ev.conn2(:,el2));        % [nsh1 x nsh2]
        Bl = Bmat(Ev.conn1(:,el1), Ev.conn2(:,el2));
        % parametric Hessian / value fields of A, [nqn1 x nqn2]
        a11 = Ev.Bxx1(:,:,el1) * Al * Ev.B2(:,:,el2).';
        a12 = Ev.Bx1(:,:,el1)  * Al * Ev.Bx2(:,:,el2).';
        a22 = Ev.B1(:,:,el1)   * Al * Ev.Bxx2(:,:,el2).';
        av  = Ev.B1(:,:,el1)   * Al * Ev.B2(:,:,el2).';
        % ... and of B
        b11 = Ev.Bxx1(:,:,el1) * Bl * Ev.B2(:,:,el2).';
        b12 = Ev.Bx1(:,:,el1)  * Bl * Ev.Bx2(:,:,el2).';
        b22 = Ev.B1(:,:,el1)   * Bl * Ev.Bxx2(:,:,el2).';
        bv  = Ev.B1(:,:,el1)   * Bl * Ev.B2(:,:,el2).';
        CXX = Ev.CXX(:,:,el1,el2,:); CYY = Ev.CYY(:,:,el1,el2,:); CXY = Ev.CXY(:,:,el1,el2,:);
        axx = CXX(:,:,1).*a11 + CXX(:,:,2).*a12 + CXX(:,:,3).*a22;
        ayy = CYY(:,:,1).*a11 + CYY(:,:,2).*a12 + CYY(:,:,3).*a22;
        axy = CXY(:,:,1).*a11 + CXY(:,:,2).*a12 + CXY(:,:,3).*a22;
        bxx = CXX(:,:,1).*b11 + CXX(:,:,2).*b12 + CXX(:,:,3).*b22;
        byy = CYY(:,:,1).*b11 + CYY(:,:,2).*b12 + CYY(:,:,3).*b22;
        bxy = CXY(:,:,1).*b11 + CXY(:,:,2).*b12 + CXY(:,:,3).*b22;
        if with_nu
            dens = axx.*bxx + ayy.*byy + H.POISSON*(axx.*byy + ayy.*bxx) ...
                   + 2*(1-H.POISSON)*axy.*bxy;
        else
            dens = axx.*bxx + ayy.*byy + 2*axy.*bxy;
        end
        wgt = Ev.qw1(:,el1) * Ev.qw2(:,el2).' .* Ev.detJ(:,:,el1,el2);
        Eb(el1,el2) = Dbf * sum(sum(wgt .* dens));
        Em(el1,el2) = H.RHO * sum(sum(wgt .* (av.*bv)));
    end
end
end

% ---------------------------------------------------------------------------
function R0 = resolve_R0(H, R0, t_e)
% Resolve the quadratic-velocity weighting operator R0 (design-independent).
if isempty(R0)
    R0 = speye(H.N);                       % identity: J3 = w^2 * sum|U|^2
elseif ischar(R0)
    switch lower(R0)
        case 'mass'
            % unit-thickness mass matrix M0 = M(t)/t (mass-weighted velocity)
            Vv = harmonic_kd_form(H, ones(size(t_e)), 0, 0, 0, 'values');
            R0 = sparse(H.rows, H.cols, Vv.Mm, H.N, H.N);
        otherwise
            error('harmonic_kd_sens: unknown R0 spec ''%s''.', R0);
    end
elseif isscalar(R0) && isnumeric(R0)
    R0 = R0 * speye(H.N);
end
end
