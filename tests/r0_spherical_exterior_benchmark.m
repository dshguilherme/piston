function r0_spherical_exterior_benchmark()
% Exact exterior-Neumann benchmark for a uniformly breathing sphere.
% Compare the paper's prescribed-source sinc quadratic with the physical
% outgoing-wave impedance; this is not a cap or shell prediction.
root=fileparts(fileparts(mfilename('fullpath'))); out=fullfile(root,'results','r0_revision');
if ~exist(out,'dir'), mkdir(out); end
R=1; rho0=1.21; c0=343; V0=1; xList=[0.1 0.5 1 2 3 4 5 6 8];
% Gauss-Legendre rule in mu=cos(theta), periodic trapezoid in azimuth.
Nmu=20; Nphi=40; b=(1:Nmu-1)./sqrt(4*(1:Nmu-1).^2-1);
[Q,D]=eig(diag(b,1)+diag(b,-1)); [mu,ord]=sort(diag(D)); wmu=2*(Q(1,ord).^2).';
phi=(0:Nphi-1)*2*pi/Nphi; [MU,PHI]=ndgrid(mu,phi);
muPts=MU(:); phiPts=PHI(:);
X=R*[(sqrt(1-muPts.^2).*cos(phiPts)).'; (sqrt(1-muPts.^2).*sin(phiPts)).'; muPts.'];
w=R^2*repmat(wmu(:),Nphi,1)*(2*pi/Nphi); S=4*pi*R^2;
r2=max(0,sum(X.^2,1).'+sum(X.^2,1)-2*(X.'*X)); r=sqrt(r2);
qnum=zeros(numel(xList),1); qexact=qnum; ratio_quad=qnum; ratio_formula=qnum; relerr=qnum; zerr=qnum;
for n=1:numel(xList)
  x=xList(n); k=x/R; kr=k*r; K=ones(size(kr)); nz=kr>1e-14; K(nz)=sin(kr(nz))./kr(nz);
  qnum(n)=real(w.'*K*w);
  j0=sin(x)/x; qexact(n)=S^2*j0^2*abs(V0)^2;
  omega=c0*k; alpha=rho0*omega^2/(4*pi*c0);
  Psrc=alpha*qnum(n);
  Z0=rho0*c0*(x^2+1i*x)/(1+x^2);
  h0=1i*exp(-1i*x)/x; h0prime=h0*(-1i-1/x);
  Zhankel=-1i*rho0*c0*h0/h0prime; zerr(n)=abs(Zhankel-Z0)/abs(Z0);
  Pext=0.5*S*real(Z0)*abs(V0)^2;
  ratio_quad(n)=Psrc/Pext;
  ratio_formula(n)=2*sin(x)^2*(1+x^2)/x^2;
  relerr(n)=abs(qnum(n)-qexact(n))/abs(qexact(n));
end
T=table(xList(:),qnum,qexact,relerr,ratio_quad,ratio_formula,zerr, ...
 'VariableNames',{'ka','quadrature_sinc_integral','analytic_sinc_integral', ...
 'relative_quadrature_error','prescribed_sinc_to_exterior_power_ratio','analytic_ratio','relative_impedance_formula_error'});
T.relative_impedance_formula_error=zerr;
writetable(T,fullfile(out,'spherical-exterior-benchmark.csv'));
fprintf('max direct-to-analytic double-integral error %.3e\n',max(relerr));
fprintf('max ratio formula discrepancy %.3e\n',max(abs(ratio_quad-ratio_formula)));
fprintf('max impedance formula discrepancy %.3e\n',max(zerr));
disp(T); fprintf('R0_SPHERICAL_EXTERIOR_BENCHMARK_PASS\n');
end
