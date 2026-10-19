function r0_sphere_modal_spectrum()
% Compare prescribed-source sinc and exact exterior operators by sphere mode.
root=fileparts(fileparts(mfilename('fullpath'))); out=fullfile(root,'results','r0_revision');
if ~exist(out,'dir'), mkdir(out); end
R=1; rho0=1.21; c0=343; xList=[0.1 0.25 0.5 1 1.5 2 2.5 3 3.5 4 5 6 8]; Lmax=10;
rows=cell(numel(xList)*(Lmax+1),5); it=0; ratioGrid=zeros(Lmax+1,numel(xList));
for ix=1:numel(xList)
 x=xList(ix);
 for ell=0:Lmax
   jl=sqrt(pi/(2*x))*besselj(ell+0.5,x);
   jl1=sqrt(pi/(2*x))*besselj(ell+1.5,x);
   yl=sqrt(pi/(2*x))*bessely(ell+0.5,x);
   yl1=sqrt(pi/(2*x))*bessely(ell+1.5,x);
   hl=jl-1i*yl; hl1=jl1-1i*yl1;
   hprime=ell/x*hl-hl1;
   Zhat=-1i*hl/hprime; F=real(Zhat);
   assert(isfinite(F) && F>0,'Nonpositive exact exterior resistance at l=%d, ka=%g.',ell,x);
   ratio=2*x^2*jl^2/F;
   ratioGrid(ell+1,ix)=ratio;
   it=it+1; rows(it,:)={x,ell,4*pi*R^2*jl^2,F,ratio};
 end
end
T=cell2table(rows,'VariableNames',{'ka','degree_l','sinc_integral_operator_eigenvalue_R1', ...
 'normalized_exterior_radiation_resistance','sinc_to_exterior_power_ratio'});
writetable(T,fullfile(out,'sphere-modal-spectrum.csv'));
% Independent direct surface quadrature for normalized Y_l0 modes at ka=1.
Nmu=20; Nphi=40; b=(1:Nmu-1)./sqrt(4*(1:Nmu-1).^2-1);
[Q,D]=eig(diag(b,1)+diag(b,-1)); [mu,ord]=sort(diag(D)); wmu=2*(Q(1,ord).^2).';
phi=(0:Nphi-1)*2*pi/Nphi; [MU,PHI]=ndgrid(mu,phi); muv=MU(:); phiv=PHI(:);
X=R*[(sqrt(1-muv.^2).*cos(phiv)).'; (sqrt(1-muv.^2).*sin(phiv)).'; muv.'];
w=R^2*repmat(wmu(:),Nphi,1)*(2*pi/Nphi); rr=sqrt(max(0,sum(X.^2,1).'+sum(X.^2,1)-2*(X.'*X)));
x=1; kr=x/R*rr; K=ones(size(kr)); nz=kr>1e-14; K(nz)=sin(kr(nz))./kr(nz);
checks=cell(6,5);
for ell=0:5
  P=legendre(ell,muv.'); Y=sqrt((2*ell+1)/(4*pi))*reshape(P(1,:),[],1);
  q=real((w.*Y)'*K*(w.*Y)); jl=sqrt(pi/2)*besselj(ell+0.5,1);
  qExact=4*pi*R^4*jl^2;
  rowsel=find(T.ka==1 & T.degree_l==ell,1); ratioFormula=T.sinc_to_exterior_power_ratio(rowsel);
  ratioDirect=q/(4*pi*R^4*jl^2)*ratioFormula;
  checks(ell+1,:)={ell,q,qExact,abs(q-qExact)/abs(qExact),ratioDirect};
end
C=cell2table(checks,'VariableNames',{'degree_l','direct_quadratic','analytic_quadratic','relative_quadrature_error','direct_ratio'});
writetable(C,fullfile(out,'sphere-modal-direct-check.csv'));
fprintf('max direct modal integral relative error %.3e\n',max(C.relative_quadrature_error));
fprintf('ka=1, l=0 ratio %.8g; max ratio difference from direct %.3e\n', ...
 ratioGrid(1,find(xList==1,1)),max(abs(C.direct_ratio-T.sinc_to_exterior_power_ratio(ismember(T.degree_l,0:5)&T.ka==1))));
disp(C); fprintf('R0_SPHERE_MODAL_SPECTRUM_PASS\n');
end
