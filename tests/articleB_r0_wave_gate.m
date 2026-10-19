% ARTICLEB_R0_WAVE_GATE  Regime gate for the phase-cancellation ladder in the
% acoustic WAVE regime.  VERIFICATION FIRST, before any optimisation.
%
% For a sweep of assembly/excitation frequencies it reports, for the UNIFORM
% (initial) design:
%   - k*Lmax and sinc(k*Lmax)         (acoustic compactness)
%   - the rank-1 indicator of the sinc kernel K  (s2/s1 of KSS; and
%     ||K-ones||/||K|| when the dense kernel is available)
%   - rho_coh = |sum_e V_e| / sum_e |V_e|  (coherence of the surface velocity)
%   - the DIAGONAL vs OFF-DIAGONAL split of z^H K z (interference term)
% A frequency is a "genuine wave regime" only when the kernel has departed from
% rank-1 AND rho_coh is clearly below 1.
%
%   articleB_r0_wave_gate([128 64])
%   articleB_r0_wave_gate([128 64], [2 3 5 8 10])
%
% Writes results/articleB_r0_wave_gate.txt (new file only).
function articleB_r0_wave_gate(nsub, fscan)
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(RES,'dir'), mkdir(RES); end
if nargin<1 || isempty(nsub),  nsub=[128 64]; end
if nargin<2 || isempty(fscan), fscan=[2 3 4 5 6 8 10 14]; end
tag=sprintf('%dx%d',nsub(1),nsub(2));

E=4.32e8; NU=0.3; RHO=2500; t0=0.25; rho0=1.21; c0=343; zeta=0.05;

ts=tic;
ctxc=articleB_build_ctx([32 16],E,NU,RHO,t0,1,0);
Kc=sparse(ctxc.Rows(:),ctxc.Cols(:),reshape(ctxc.Bke.*t0^3+ctxc.Ske.*t0,[],1),ctxc.ndof,ctxc.ndof);
Mc=sparse(ctxc.Rows(:),ctxc.Cols(:),reshape(ctxc.Me(:)*t0,[],1),ctxc.ndof,ctxc.ndof);
[~,D]=eigs(Kc(ctxc.free,ctxc.free),Mc(ctxc.free,ctxc.free),6,'sm');
f1=min(sqrt(abs(diag(D))))/(2*pi); clear ctxc Kc Mc D
fprintf('f1 (coarse) = %.6g Hz (first structural mode)\n',f1);

ctx=articleB_build_ctx(nsub,E,NU,RHO,t0,fscan(1),zeta);
fprintf('ctx %s built in %.1f s  N=%d ne=%d\n',tag,toc(ts),ctx.ndof,ctx.ne);

Lmax=blockmax(ctx);
fid=fopen(fullfile(RES,'articleB_r0_wave_gate.txt'),'w');
lg=@(varargin) both(fid,varargin{:});
lg('# ARTICLEB_R0_WAVE_GATE  regime gate, mesh %s, N=%d, ne=%d, f1=%.6g Hz\n',tag,ctx.ndof,ctx.ne,f1);
lg('# fluid rho0=%.3g c0=%.3g ; zeta=%.2f ; uniform design t=%.3g ; E=%.3g RHO=%.3g\n',rho0,c0,zeta,t0,E,RHO);
lg('# Lmax(centres) = %.4f m\n',Lmax);
lg('# f_Hz           k          kLmax     sinc     lambda      s2/s1(KSS)  ||K-ones||/||K||   rho_coh      diag         offdiag      od/tot\n');
fprintf('Lmax(centres)=%.4f m\n',Lmax);

for i=1:numel(fscan)
  f=fscan(i);
  ctx.w=2*pi*f; ctx.alpha=2*zeta*ctx.w; ctx.beta=0;
  R0op=articleB_r0_build(ctx,rho0,c0,128,f);
  sv=svd(R0op.KSS);
  relK1=sv(2)/sv(1);
  if ~isempty(R0op.Kfull)
    relKfro=norm(R0op.Kfull-ones(ctx.ne),'fro')/norm(R0op.Kfull,'fro');
  else
    Kn=R0op.KS.'*pinv((R0op.KSS+R0op.KSS.')/2)*R0op.KS;
    relKfro=norm(Kn-ones(ctx.ne),'fro')/norm(Kn,'fro');
  end
  k=R0op.k; kl=k*Lmax; s=mysin(kl); lam=2*pi/k;
  t=t0*ones(ctx.ne,1);
  [~,~,aux]=articleB_r0_eval(ctx,t,'quadvel',R0op,'I');
  u=aux.u; Cu=R0op.C*u; z=R0op.area.*Cu;
  sumV=abs(sum(Cu)); sumabs=sum(abs(Cu)); rho=sumV/max(sumabs,eps);
  Kz=applyK(R0op,z); tot=real(z'*Kz); dg=sum(abs(z).^2); od=tot-dg;
  lg('%12.5g %12.5g %8.4f %9.4f %10.3f %12.5g %14.5g %11.6f %13.5e %13.5e %10.4f\n',...
     f,k,kl,s,lam,relK1,relKfro,rho,dg,od,od/max(abs(tot),eps));
  fprintf('f=%7.4g  kLmax=%7.4f  s2/s1=%.4g  ||K-1||=%.4g  rho_coh=%.6f  od/tot=%+.4f\n',...
     f,kl,relK1,relKfro,rho,od/max(abs(tot),eps));
end
fclose(fid);
fprintf('gate written: %s\n',fullfile(RES,'articleB_r0_wave_gate.txt'));
fprintf('ARTICLEB_R0_WAVE_GATE_DONE\n');
end

function Kz=applyK(R0op,z)
Kz=R0op.KS.'*(R0op.KSSinv*(R0op.KS*z));
end
function d=blockmax(ctx)
m=ctx.msh; Jw=m.jacdet.*m.quad_weights; area=sum(Jw,1);
Jw3=reshape(Jw,1,[],ctx.ne); gm=reshape(m.geo_map,3,[],ctx.ne);
X=reshape(squeeze(sum(gm.*Jw3,2))./area,3,ctx.ne);
n=size(X,2); dmax=0; blk=512;
for i=1:blk:n
  ii=i:min(i+blk-1,n);
  DX=X(1,ii).'-X(1,:); DY=X(2,ii).'-X(2,:); DZ=X(3,ii).'-X(3,:);
  DD=DX.^2+DY.^2+DZ.^2; dmax=max(dmax,max(DD(:)));
end
d=sqrt(dmax);
end
function s=mysin(z), s=ones(size(z)); idx=z~=0; s(idx)=sin(z(idx))./z(idx); end
function both(fid,varargin), fprintf(fid,varargin{:}); fprintf(1,varargin{:}); end
