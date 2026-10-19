% ARTICLEB_R0_FDCHECK  FD verification of the quadratic-velocity gradient for
% the three-operator phase-cancellation ladder (thesis eq. 3.66 radiation R0).
%
% For EACH operator (R0 = I / R0 = C*ones / R0 = eq.3.66) we check the analytic
% gradient of articleB_r0_eval against a central finite difference of W
% evaluated INDEPENDENTLY from V and R (fresh solve + separately assembled
% C, R). The gradient code path (adjoint cross-energy) is never used to build
% the reference W.
%
% Also cross-checks the vectorised R (=3.66) and C operators against explicit
% O(ne^2)/O(ne) loop recomputation.
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(RES,'dir'), mkdir(RES); end
fid=fopen(fullfile(RES,'articleB_r0_fdcheck.txt'),'w');
lg=@(varargin) both(fid,varargin{:});

E=4.32e8; NU=0.3; RHO=2500; t0=0.25; rho0=1.21; c0=343;
nsub=[16 8];
ctx0=articleB_build_ctx(nsub,E,NU,RHO,t0,NaN,0);
Ks0=sparse(ctx0.Rows(:),ctx0.Cols(:),reshape(ctx0.Bke.*t0^3+ctx0.Ske.*t0,[],1),ctx0.ndof,ctx0.ndof);
M0 =sparse(ctx0.Rows(:),ctx0.Cols(:),reshape(ctx0.Me(:)*t0,[],1),ctx0.ndof,ctx0.ndof);
[~,D]=eigs(Ks0(ctx0.free,ctx0.free),M0(ctx0.free,ctx0.free),3,'sm');
f1=min(sqrt(abs(diag(D))))/(2*pi); freq=0.97*f1; zeta=0.05;
ctx=articleB_build_ctx(nsub,E,NU,RHO,t0,freq,zeta);
R0op=articleB_r0_build(ctx,rho0,c0,128);
ne=ctx.ne; N=ctx.ndof;
lg('==== ARTICLE B radiation-R0 FD verification (mesh %s, N=%d, f=%.6g Hz, zeta=%.2f) ====\n',...
   mat2str(nsub),N,freq,zeta);
lg('fluid: rho0=%.3g kg/m^3  c0=%.3g m/s  k=%.4e 1/m  k*Lmax=%.4e\n',...
   rho0,c0,R0op.k, R0op.k*maxdist(R0op.centres));
lg('excitation omega=%.6g rad/s (R0 assembled at this frequency)\n\n',R0op.omega);

%% (0) operator cross-checks (Nystrom kernel vs exact dense; C vs loop)
if ~isempty(R0op.Kfull)
  Knys = R0op.KS.'*(R0op.KSSinv*R0op.KS);
  relK = norm(Knys-R0op.Kfull,'fro')/norm(R0op.Kfull,'fro');
  relRank1 = norm(R0op.Kfull-ones(ne,ne),'fro')/norm(R0op.Kfull,'fro');
else
  relK=NaN; relRank1=NaN;
end
lg('(0a) Nystrom sinc kernel vs exact dense (m=%d): rel err = %.2e ; ||K-ones||/||K|| = %.4f (rank-1 nature)\n',...
   R0op.mrank, relK, relRank1);
% C vs loop: recompute (C u)_e for element e by explicit sum over local basis
rng(7); ut=randn(N,1)+1i*randn(N,1);
Cu=R0op.C*ut; Cmaxrel=0; nloc=size(ctx.lm,2); q=nloc/3;
for e=[1 3 11 ne]
  Me=reshape(ctx.Me(e,:),nloc,nloc); val=0;
  for c=1:3
    blk=Me((c-1)*q+(1:q),(c-1)*q+(1:q)); intphi=sum(blk,2)/RHO;
    g=intphi*R0op.navg(c,e)/R0op.area(e);
    for a=1:q, val=val+g(a)*ut(ctx.lm(e,(c-1)*q+a)); end
  end
  Cmaxrel=max(Cmaxrel,abs(val-Cu(e))/max(abs(val),eps));
end
lg('(0b) C vectorised vs loop, 4 elements: max rel err = %.2e\n\n',Cmaxrel);

%% (1) gradient FD for each operator
t=t0*ones(ne,1); rng(3); pe=randperm(ne,10); pert=zeros(ne,1); pert(pe)=1; h=1e-7;
modes={'I','ones','rad'};
for im=1:numel(modes)
  mode=modes{im};
  [J,dJ]=articleB_r0_eval(ctx,t,'quadvel',R0op,mode);
  Jfun=@(tt) Wref(ctx,tt,R0op,mode);
  % uniform FD of the total
  dJfd_u=(Jfun(t+h*pert)-Jfun(t-h*pert))/(2*h);
  eu=abs(sum(dJ(pe))-dJfd_u)/abs(dJfd_u);
  % single-element FD
  err=zeros(numel(pe),1);
  for kk=1:numel(pe)
    tp=t; tp(pe(kk))=tp(pe(kk))+h; tm=t; tm(pe(kk))=-h+tm(pe(kk));
    dfd=(Jfun(tp)-Jfun(tm))/(2*h);
    err(kk)=abs(dJ(pe(kk))-dfd)/max(abs(dfd),eps);
  end
  lg('(1) R0mode=%-4s : J=%.6e | uniform-FD rel %.2e | single-elem max %.2e median %.2e\n',...
     mode,J,eu,max(err),median(err));
end
lg('\nARTICLEB_R0_FDCHECK_DONE\n');
fclose(fid);
fprintf('wrote %s\n',fullfile(RES,'articleB_r0_fdcheck.txt'));

% ---------------------------------------------------------------------------
function W = Wref(ctx,t,R0op,mode)
% INDEPENDENT evaluation of W for the given operator: fresh solve, then build W
% directly from V and R (never via articleB_r0_eval's gradient path).
t=t(:); ne=ctx.ne; lm=ctx.lm; n=size(lm,2); ndof=ctx.ndof;
cK=1+1i*ctx.w*ctx.beta; cM=1i*ctx.w*ctx.alpha-ctx.w^2;
EK=ctx.Bke.*(t.^3)+ctx.Ske.*t; EK=EK(:);
EM=ctx.Me.*t; EM=EM(:);
Ks=sparse(ctx.Rows(:),ctx.Cols(:),EK,ndof,ndof);
M =sparse(ctx.Rows(:),ctx.Cols(:),EM,ndof,ndof);
KD=cK*Ks+cM*M;
u=zeros(ndof,1); u(ctx.free)=KD(ctx.free,ctx.free)\ctx.F(ctx.free);
Cu=R0op.C*u;
switch lower(mode)
  case 'i'
    W=ctx.w^2*real(u'*u);
  case 'ones'
    W=ctx.w^2*R0op.pref*abs(sum(Cu))^2;
  case 'rad'
    z=R0op.area.*Cu;
    if ~isempty(R0op.Kfull)
      W=ctx.w^2*R0op.pref*real(z'*(R0op.Kfull*z));   % EXACT dense kernel
    else
      Kz=R0op.KS.'*(R0op.KSSinv*(R0op.KS*z));
      W=ctx.w^2*R0op.pref*real(z'*Kz);
    end
end
end
function both(fid,varargin), fprintf(fid,varargin{:}); fprintf(1,varargin{:}); end
function d = maxdist(X)
DX = X(1,:)-X(1,:).'; DY = X(2,:)-X(2,:).'; DZ = X(3,:)-X(3,:).';
d = sqrt(max(max(DX.^2+DY.^2+DZ.^2)));
end
