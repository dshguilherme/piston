function op = r0_angular_gram(x, weights, B, omega, rho, c, nmu, nphi, normalization)
%R0_ANGULAR_GRAM Consistent quadrature of the sinc source-layer power form.
% x: 3 x nq physical surface points; weights: nq x 1 surface measures.
% B: nq x ndof map from velocity coefficients to local normal velocity.
% This is a source-layer model, NOT an exterior Neumann solver for a shell.
% With normalization='rayleigh', W=v'*R*v uses peak complex velocities and
% R=(rho*omega^2/(4*pi*c))*B'*diag(w)*sinc(k*r)*diag(w)*B.
% That power normalization is physical for a planar baffled radiator only.
% 'monopole-fullspace' gives half this power for prescribed free monopoles;
% it does NOT turn B*v into a physical unbaffled-shell boundary solution.
% The angular identity sinc(k*r)=1/(4*pi)*int_S2 exp(i*k*s.r) dOmega
% gives a positive-semidefinite Gram factor, without forming nq x nq K.
if nargin < 9, normalization = 'rayleigh'; end
validateattributes(x, {'double'}, {'2d','nrows',3,'finite','real'});
validateattributes(weights, {'double'}, {'vector','positive','finite','real'});
validateattributes(omega, {'double'}, {'scalar','positive','finite'});
validateattributes(rho, {'double'}, {'scalar','positive','finite'});
validateattributes(c, {'double'}, {'scalar','positive','finite'});
validateattributes(nmu, {'double'}, {'scalar','integer','>=',2});
validateattributes(nphi, {'double'}, {'scalar','integer','>=',4});
w = weights(:); nq = numel(w);
assert(size(x,2)==nq && size(B,1)==nq, 'Surface quadrature dimensions disagree.');
% Gauss-Legendre in mu=cos(theta), trapezoidal in azimuth.
i=(1:nmu-1)'; beta=i./sqrt(4*i.^2-1);
[V,D]=eig(diag(beta,1)+diag(beta,-1));
[mu,order]=sort(diag(D)); wm=2*(V(1,order).^2)';
phi=2*pi*(0:nphi-1)/nphi;
[MU,PHI]=ndgrid(mu,phi);
s=[sqrt(1-MU(:).^2).*cos(PHI(:)), ...
   sqrt(1-MU(:).^2).*sin(PHI(:)), MU(:)];
aw=repmat(wm,nphi,1)/(2*nphi); % angular weights / (4*pi)
pref=rho*omega^2/(4*pi*c);
switch lower(normalization)
    case 'rayleigh'
    case 'monopole-fullspace', pref=pref/2;
    otherwise, error('Unknown power normalization: %s',normalization);
end
WB=spdiags(w,0,nq,nq)*B;
H=complex(zeros(size(s,1),size(B,2)));
for first=1:128:size(s,1)
    rows=first:min(first+127,size(s,1));
    H(rows,:)=exp(-1i*(omega/c)*(s(rows,:)*x))*WB;
end
L=bsxfun(@times,sqrt(pref*aw),H);
op=struct('factor',L,'directions',s,'angular_weights',aw, ...
    'omega',omega,'rho',rho,'c',c,'normalization',normalization, ...
    'nmu',nmu,'nphi',nphi,'nq',nq);
end
