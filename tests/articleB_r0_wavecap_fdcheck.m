% ARTICLEB_R0_WAVECAP_FDCHECK  Careful central-difference verification of the
% quadvel (R0) gradient on the clamped spherical cap, at the SAME frequency and
% design family as the wave-regime ladder, across several step sizes h and for
% all three operators.  Independent reference W built from a fresh solve + R0
% applied explicitly (never through articleB_r0_eval's gradient path).
%
%   articleB_r0_wavecap_fdcheck([24 12], 100)
%
% Writes results/articleB_r0_wavecap_fdcheck.txt (new file only).
function articleB_r0_wavecap_fdcheck(nsub, freq)
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if nargin<1 || isempty(nsub),  nsub=[24 12]; end
if nargin<2 || isempty(freq),  freq=100; end
rho0=1.21; c0=343; zeta=0.05;
spec=struct('R',10,'alpha',deg2rad(30),'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,...
  'deg',[2 2],'reg',[1 1],'nsub',nsub,'nquad',[4 4],'bc',{2},...
  'loadtype','pole','loaddir',[0;0;1],'loadmag',1,'freq',freq,'zeta',zeta);
ctx=obstacle_build_ctx('spherecap',spec);
R0op=articleB_r0_build(ctx,rho0,c0,128,freq);
ne=ctx.ne; t0=ctx.t0;
fid=fopen(fullfile(RES,'articleB_r0_wavecap_fdcheck.txt'),'w');
lg=@(varargin) both(fid,varargin{:});
lg('# ARTICLEB_R0_WAVECAP_FDCHECK  mesh %s  N=%d ne=%d  f=%.6g Hz  zeta=%.2f  t0=%.3g\n',mat2str(nsub),ctx.ndof,ctx.ne,freq,zeta,t0);
lg('# k=%.4e k*Lmax=%.4f\n',R0op.k, R0op.k*blockmax(ctx));

tu=t0*ones(ne,1);
rng(3); pe=randperm(ne,6); pert=zeros(ne,1); pert(pe)=1;
modes={'I','ones','rad'};
for im=1:numel(modes)
  mode=modes{im};
  [J,dJ]=articleB_r0_eval(ctx,tu,'quadvel',R0op,mode);
  Jf=@(tt) Wref(ctx,tt,R0op,mode);
  lg('\n== R0mode=%s : J=%.8e ==\n',mode,J);
  hs=[1e-3 1e-4 1e-5 1e-6 1e-7];
  for h=hs
    d_dir=sum(dJ(pe));
    dfd=(Jf(tu+h*pert)-Jf(tu-h*pert))/(2*h);
    lg('  h=%.0e  analytic dir-deriv=%.8e  FD=%.8e  rel=%.3e\n',h,d_dir,dfd,abs(d_dir-dfd)/max(abs(dfd),eps));
  end
  % per-element native-scale FD (perturb one element by h, divide by h)
  h=1e-5; err=zeros(numel(pe),1);
  for kk=1:numel(pe)
    tp=tu; tp(pe(kk))=tp(pe(kk))+h; tm=tu; tm(pe(kk))=tm(pe(kk))-h;
    dfd=(Jf(tp)-Jf(tm))/(2*h);
    err(kk)=abs(dJ(pe(kk))-dfd)/max(abs(dfd),eps);
  end
  lg('  per-element (h=1e-5): max=%.3e median=%.3e\n',max(err),median(err));
end
lg('\nARTICLEB_R0_WAVECAP_FDCHECK_DONE\n');
fclose(fid);
fprintf('wrote %s\n',fullfile(RES,'articleB_r0_wavecap_fdcheck.txt'));
end

function W=Wref(ctx,t,R0op,mode)
t=t(:); ne=ctx.ne; lm=ctx.lm; n=size(lm,2); ndof=ctx.ndof;
cK=1+1i*ctx.w*ctx.beta; cM=1i*ctx.w*ctx.alpha-ctx.w^2;
EK=ctx.Bke.*(t.^3)+ctx.Ske.*t; EK=EK(:); EM=ctx.Me.*t; EM=EM(:);
Ks=sparse(ctx.Rows(:),ctx.Cols(:),EK,ndof,ndof); M=sparse(ctx.Rows(:),ctx.Cols(:),EM,ndof,ndof);
KD=cK*Ks+cM*M; u=zeros(ndof,1); u(ctx.free)=KD(ctx.free,ctx.free)\ctx.F(ctx.free);
Cu=R0op.C*u;
switch lower(mode)
  case 'i',    W=ctx.w^2*real(u'*u);
  case 'ones', W=ctx.w^2*R0op.pref*abs(sum(Cu))^2;
  case 'rad',  z=R0op.area.*Cu; Kz=R0op.KS.'*(R0op.KSSinv*(R0op.KS*z));
               W=ctx.w^2*R0op.pref*real(z'*Kz);
end
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
function both(fid,varargin), fprintf(fid,varargin{:}); fprintf(1,varargin{:}); end
