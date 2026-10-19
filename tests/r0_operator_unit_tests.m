function r0_operator_unit_tests(output_dir)
% Independent analytic and algebraic tests; no legacy solver modification.
root=fileparts(fileparts(mfilename('fullpath'))); addpath(fullfile(root,'src'));
if nargin<1, output_dir=fullfile(root,'results','r0_revision'); end
if ~exist(output_dir,'dir'), mkdir(output_dir); end
rng(921,'twister'); rho=1.21; c=343; omega=2*pi*100;
x=randn(3,13); w=0.1+rand(13,1); B=randn(13,7);
op=r0_angular_gram(x,w,B,omega,rho,c,32,64);
D=sqrt(max(0,sum(x.^2,1)'+sum(x.^2,1)-2*(x'*x)));
z=omega/c*D; K=ones(size(z)); mask=z~=0; K(mask)=sin(z(mask))./z(mask);
WB=w.*B; R=rho*omega^2/(4*pi*c)*(WB'*K*WB);
Rg=op.factor'*op.factor;
err_identity=norm(Rg-R,'fro')/norm(R,'fro');
assert(err_identity<1e-11,'Angular identity failed: %.3e',err_identity);
hermitian=norm(Rg-Rg','fro')/norm(Rg,'fro');
assert(hermitian<1e-13); assert(min(eig((Rg+Rg')/2))>-1e-11*norm(Rg));
v=randn(7,1)+1i*randn(7,1); du=randn(7,1)+1i*randn(7,1);
h=1e-5; fd=(norm(op.factor*(v+h*du))^2-norm(op.factor*(v-h*du))^2)/(2*h);
exact=2*real((op.factor*v)'*(op.factor*du));
err_derivative=abs(fd-exact)/max(abs(exact),eps);
assert(err_derivative<1e-8);
% Uniform circular piston in infinite rigid baffle: W=0.5*rho*c*A*
% (1-J1(2ka)/(ka))*|v|^2. Surface polar quadrature converged independently.
a=0.7; nr=32; nt=96;
[r,wr]=gauss01(nr); r=a*r; wr=a*wr;
[RR,TT]=ndgrid(r,2*pi*(0:nt-1)/nt);
xp=[RR(:).*cos(TT(:)),RR(:).*sin(TT(:)),zeros(numel(RR),1)]';
wp=repmat(wr.*r,nt,1)*(2*pi/nt);
ka_list=[0.1,1,3,6]; numerical=zeros(size(ka_list)); analytical=numerical;
for j=1:numel(ka_list)
    ka=ka_list(j); om=ka*c/a;
    P=r0_angular_gram(xp,wp,ones(numel(wp),1),om,rho,c,48,96);
    numerical(j)=norm(P.factor)^2;
    analytical(j)=0.5*rho*c*pi*a^2*(1-besselj(1,2*ka)/ka);
end
piston_errors=abs(numerical-analytical)./analytical;
assert(max(piston_errors)<1e-9,'Piston analytic comparison failed.');
T=table(ka_list',analytical',numerical',piston_errors', ...
    'VariableNames',{'ka','analytic_power_W','quadrature_power_W','relative_error'});
writetable(T,fullfile(output_dir,'unit_piston.csv'));
save(fullfile(output_dir,'unit_tests.mat'),'T','err_identity','hermitian','err_derivative');
fprintf('Angular vs dense sinc relative error: %.6e\n',err_identity);
fprintf('Hermitian residual: %.6e\n',hermitian);
fprintf('Power directional derivative error: %.6e\n',err_derivative);
disp(T); fprintf('R0_OPERATOR_UNIT_TESTS_PASS\n');
end
function [x,w]=gauss01(n)
i=(1:n-1)'; b=i./sqrt(4*i.^2-1); [V,D]=eig(diag(b,1)+diag(b,-1));
[x,idx]=sort(diag(D)); w=V(1,idx)'.^2; x=(x+1)/2;
end
