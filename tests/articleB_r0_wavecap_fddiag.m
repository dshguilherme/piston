% ARTICLEB_R0_WAVECAP_FDDIAG  Diagnose the 'rad' gradient failure in the wave
% regime: is it (a) the Nystrom rank m, (b) kernel conditioning, or (c) the
% adjoint formula?  Compares, at t and t+-h, the eval forward J against an
% independent Wref using (i) the same Nystrom kernel, (ii) the EXACT dense
% kernel, and (iii) higher Nystrom rank; and checks the analytic gradient
% against dense-kernel finite differences.
function articleB_r0_wavecap_fddiag(nsub, freq)
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if nargin<1||isempty(nsub), nsub=[24 12]; end
if nargin<2||isempty(freq), freq=100; end
rho0=1.21; c0=343; zeta=0.05;
spec=struct('R',10,'alpha',deg2rad(30),'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,...
  'deg',[2 2],'reg',[1 1],'nsub',nsub,'nquad',[4 4],'bc',{2},...
  'loadtype','pole','loaddir',[0;0;1],'loadmag',1,'freq',freq,'zeta',zeta);
ctx=obstacle_build_ctx('spherecap',spec);
ne=ctx.ne; t0=ctx.t0; tu=t0*ones(ne,1);
rng(3); pe=randperm(ne,6); pert=zeros(ne,1); pert(pe)=1;

fid=fopen(fullfile(RES,'articleB_r0_wavecap_fddiag.txt'),'w');
lg=@(varargin) both(fid,varargin{:});
lg('FD DIAG cap mesh %s ne=%d f=%g\n',mat2str(nsub),ne,freq);
mranks=[128 256 512];
for mr=mranks
  R0op=articleB_r0_build(ctx,rho0,c0,mr,freq);
  sv=svd(R0op.KSS); 
  lg('\n--- Nystrom m=%d : KSS cond=%.3e  s2/s1=%.4f ---\n',mr,cond(R0op.KSS),sv(2)/sv(1));
  [J,dJ]=articleB_r0_eval(ctx,tu,'quadvel',R0op,'rad');
  lg('  eval J_rad=%.10e\n',J);
  for h=[1e-3 1e-4 1e-5]
    % (i) Wref with SAME nystrom kernel (matches eval's operator)
    Jp=refn(ctx,tu+h*pert,R0op,'nys'); Jm=refn(ctx,tu-h*pert,R0op,'nys');
    fd1=(Jp-Jm)/(2*h);
    % (ii) Wref with EXACT dense kernel
    Jpd=refn(ctx,tu+h*pert,R0op,'dense'); Jmd=refn(ctx,tu-h*pert,R0op,'dense');
    fd2=(Jpd-Jmd)/(2*h);
    lg('  h=%.0e | eval dnys(2h)=%.3e fd_nys=%.6e (rel %.3e) | fd_dense=%.6e (rel %.3e) | Jdense-Jnys=%.3e\n',...
      h, Jp-Jm, fd1, abs(sum(dJ(pe))-fd1)/max(abs(fd1),eps), fd2, abs(sum(dJ(pe))-fd2)/max(abs(fd2),eps), Jpd-Jp);
  end
  lg('  analytic sum(dJ(pe))=%.6e\n',sum(dJ(pe)));
end
lg('\nFD_DIAG_DONE\n');
fclose(fid);
fprintf('wrote %s\n',fullfile(RES,'articleB_r0_wavecap_fddiag.txt'));
end

function W=refn(ctx,t,R0op,mode)
t=t(:); ne=ctx.ne; lm=ctx.lm; ndof=ctx.ndof;
cK=1+1i*ctx.w*ctx.beta; cM=1i*ctx.w*ctx.alpha-ctx.w^2;
EK=ctx.Bke.*(t.^3)+ctx.Ske.*t; EK=EK(:); EM=ctx.Me.*t; EM=EM(:);
Ks=sparse(ctx.Rows(:),ctx.Cols(:),EK,ndof,ndof); M=sparse(ctx.Rows(:),ctx.Cols(:),EM,ndof,ndof);
KD=cK*Ks+cM*M; u=zeros(ndof,1); u(ctx.free)=KD(ctx.free,ctx.free)\ctx.F(ctx.free);
Cu=R0op.C*u; z=R0op.area.*Cu;
switch lower(mode)
  case 'nys',   Kz=R0op.KS.'*(R0op.KSSinv*(R0op.KS*z));
  case 'dense', Kz=R0op.Kfull*z;
end
W=ctx.w^2*R0op.pref*real(z'*Kz);
end
function both(fid,varargin), fprintf(fid,varargin{:}); fprintf(1,varargin{:}); end
