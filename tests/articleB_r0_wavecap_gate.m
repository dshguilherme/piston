% ARTICLEB_R0_WAVECAP_GATE  Regime gate for the R0 phase-cancellation ladder on
% the MODALLY DENSE clamped spherical cap (obstacle-course benchmark).
%
% Two conditions must BOTH hold for phase cancellation to be exploitable:
%   (A) STRUCTURAL / modal density: many modes with different spatial characters
%       overlap at the excitation frequency (short-wavelength, different nodal
%       patterns) -> surface velocity field is INCOHERENT (rho_coh < 1).
%   (B) ACOUSTIC / non-compactness: k*Lmax >~ 1 so the sinc kernel R0 has real
%       off-diagonal structure (not rank-1).
%
% This gate reports BOTH, for the UNIFORM (initial) design:
%   - modal count / density in the design band [80 130] Hz
%   - k*Lmax, wavelength, sinc(k*Lmax)
%   - kernel rank-1 indicator s2/s1(KSS) and ||K-ones||/||K||
%   - rho_coh = |sum_e V_e| / sum_e |V_e|
%   - DIAGONAL vs OFF-DIAGONAL split of z^H K z (interference term)
% BEFORE any optimisation.  If rho_coh ~ 1 the field is still coherent: say so.
%
%   articleB_r0_wavecap_gate([128 64], [85 90 100 110 120 128])
%
% Writes results/articleB_r0_wavecap_gate.txt (new file only).
function articleB_r0_wavecap_gate(nsub, fscan)
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(RES,'dir'), mkdir(RES); end
if nargin<1 || isempty(nsub),  nsub=[128 64]; end
if nargin<2 || isempty(fscan), fscan=[85 90 100 110 120 128]; end
tag=sprintf('%dx%d',nsub(1),nsub(2));

rho0=1.21; c0=343; zeta=0.05;
spec=struct('R',10,'alpha',deg2rad(30),'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,...
  'deg',[2 2],'reg',[1 1],'nsub',nsub,'nquad',[4 4],'bc',{2},...
  'loadtype','pole','loaddir',[0;0;1],'loadmag',1,'freq',100,'zeta',zeta);

ts=tic;
ctx=obstacle_build_ctx('spherecap',spec);
fprintf('spherecap ctx %s built %.1f s  N=%d ne=%d\n',tag,toc(ts),ctx.ndof,ctx.ne);

% modal spectrum
Ks=sparse(ctx.Rows(:),ctx.Cols(:),reshape(ctx.Bke.*ctx.t0^3+ctx.Ske.*ctx.t0,[],1),ctx.ndof,ctx.ndof);
M =sparse(ctx.Rows(:),ctx.Cols(:),reshape(ctx.Me.*ctx.t0,[],1),ctx.ndof,ctx.ndof);
[~,D]=eigs(Ks(ctx.free,ctx.free),M(ctx.free,ctx.free),60,(2*pi*105)^2);
freqs=sort(sqrt(abs(real(diag(D)))))/(2*pi);
band=[80 130]; nband=sum(freqs>=band(1)&freqs<=band(2));
Lmax=blockmax(ctx);
fprintf('modal spectrum: %d modes in [%g %g] Hz  density %.4f modes/Hz\n',nband,band(1),band(2),nband/(band(2)-band(1)));
fprintf('band freqs: %s\n',mat2str(freqs(freqs<=band(2)+1e-9).',4));

fid=fopen(fullfile(RES,'articleB_r0_wavecap_gate.txt'),'w');
lg=@(varargin) both(fid,varargin{:});
lg('# ARTICLEB_R0_WAVECAP_GATE  modally dense clamped spherical cap, mesh %s\n',tag);
lg('# geometry: R=10 m, half-angle 30 deg, cap extent (chord) %.4f m, clamped rim, crown point load\n',2*10*sin(deg2rad(30)));
lg('# shell: E=%.3g NU=%.2f RHO=%.3g t0=%.3g ; fluid rho0=%.3g c0=%.3g ; zeta=%.2f\n',spec.E,spec.NU,spec.RHO,spec.t0,rho0,c0,zeta);
lg('# N=%d ne=%d ; Lmax(centres)=%.4f m\n',ctx.ndof,ctx.ne,Lmax);
lg('# modal band [%g %g] Hz : %d modes, density %.4f modes/Hz\n',band(1),band(2),nband,nband/(band(2)-band(1)));
lg('# band freqs: %s\n',mat2str(freqs(freqs<=band(2)+1e-9).',5));
lg('# f_Hz        k           kLmax     sinc     lambda    s2/s1(KSS)   ||K-ones||/||K||   rho_coh     diag         offdiag      od/tot\n');

for i=1:numel(fscan)
  f=fscan(i);
  ctx.w=2*pi*f; ctx.alpha=2*zeta*ctx.w; ctx.beta=0;
  R0op=articleB_r0_build(ctx,rho0,c0,128,f);
  sv=svd(R0op.KSS); relK1=sv(2)/sv(1);
  Kn=R0op.KS.'*pinv((R0op.KSS+R0op.KSS.')/2)*R0op.KS;
  relKfro=norm(Kn-ones(ctx.ne),'fro')/norm(Kn,'fro');
  k=R0op.k; kl=k*Lmax; lam=2*pi/k;
  t=spec.t0*ones(ctx.ne,1);
  [~,~,aux]=articleB_r0_eval(ctx,t,'quadvel',R0op,'I');
  u=aux.u; Cu=R0op.C*u; z=R0op.area.*Cu;
  sumV=abs(sum(Cu)); sumabs=sum(abs(Cu)); rho=sumV/max(sumabs,eps);
  Kz=applyK(R0op,z); tot=real(z'*Kz); dg=sum(abs(z).^2); od=tot-dg;
  lg('%11.4f %12.5g %8.4f %9.4f %9.3f %12.5g %14.5g %11.6f %13.5e %13.5e %10.4f\n',...
     f,k,kl,mysin(kl),lam,relK1,relKfro,rho,dg,od,od/max(abs(tot),eps));
  fprintf('f=%6.3g  kLmax=%7.4f  s2/s1=%.4g  ||K-1||=%.4g  rho_coh=%.6f  od/tot=%+.4f\n',...
     f,kl,relK1,relKfro,rho,od/max(abs(tot),eps));
end
fclose(fid);
fprintf('gate written: %s\n',fullfile(RES,'articleB_r0_wavecap_gate.txt'));
fprintf('ARTICLEB_R0_WAVECAP_GATE_DONE\n');
end

function Kz=applyK(R0op,z), Kz=R0op.KS.'*(R0op.KSSinv*(R0op.KS*z)); end
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
