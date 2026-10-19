% ARTICLEB_VALIDATE  Independent checks of the Article B curved-shell operator.
%   (1) ELEMENT-MATRIX ASSEMBLY vs the monolithic geopdes operators
%       op_KL_shells_tp (bending+membrane) and op_u_v (mass)
%   (2) sensitivity FD check of all three objectives
% Writes results/articleB_validation.txt
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
fid=fopen(fullfile(RES,'articleB_validation.txt'),'w');
lg=@(varargin) both(fid,varargin{:});
E=4.32e8; NU=0.3; RHO=2500; t0=0.25;

lg('==== ARTICLE B operator validation (curved Scordelis-Lo roof shell) ====\n');
%% (1) element-matrix assembly vs monolithic assembly (matched coefficients)
nsub=[8 8]; md=struct('degree',[2 2],'regularity',[1 1],'nsub',nsub,'nquad',[3 3]);
ctx=articleB_build_ctx(nsub,E,NU,RHO,t0,NaN,0);
pd=scordelis_problem;
[~,msh_raw,sp_raw]=buildSpaces(pd,md);
Eone=@(x,y,z) ones(size(x)); nuf=@(x,y,z) NU*ones(size(x));
sp1=sp_precompute_param(ctx.space,ctx.msh,'value',true,'gradient',true,'hessian',true);
for idim=1:ctx.msh.rdim, x{idim}=reshape(ctx.msh.geo_map(idim,:,:),ctx.msh.nqn,ctx.msh.nel); end
[Bke,Ske,~]=op_KL_shells_elements(sp1,sp1,ctx.msh,Eone(x{:}),nuf(x{:}));
Me1=op_u_v_elements(sp1,sp1,ctx.msh,1);
Ks_el=sparse(ctx.Rows(:),ctx.Cols(:),reshape(Bke.*t0^3+Ske.*t0,[],1),ctx.ndof,ctx.ndof);
M_el =sparse(ctx.Rows(:),ctx.Cols(:),reshape(Me1*t0,[],1),ctx.ndof,ctx.ndof);
Ks_tp=op_KL_shells_tp(sp_raw,sp_raw,msh_raw,Eone,nuf,t0);
sp_p=sp_precompute_param(sp_raw,ctx.msh,'value',true);
[r,c,v]=op_u_v(sp_p,sp_p,ctx.msh,1); M_tp=t0*sparse(r,c,v,ctx.ndof,ctx.ndof);
eK=norm(Ks_el-Ks_tp,'fro')/norm(Ks_tp,'fro');
eM=norm(M_el-M_tp,'fro')/norm(M_tp,'fro');
area=sum(sum(ctx.msh.jacdet.*ctx.msh.quad_weights));
lg('(1) element-matrix vs monolithic assembly (nsub=%s, t=%.3f):\n',mat2str(nsub),t0);
lg('    Ks rel-Frobenius err = %.3e   (op_KL_shells_tp)\n',eK);
lg('    M  rel-Frobenius err = %.3e   (op_u_v)\n',eM);
lg('    sum(M) = %.6e   (3*rho*area = %.6e)\n',full(sum(M_el(:))),3*1*area);
lg('    %s\n\n', ternary(eK<1e-12 && eM<1e-12,'OK','FAIL'));

%% (2) FD sensitivity check on a 16x8 mesh near the first resonance
nsub=[16 8];
ctx0=articleB_build_ctx(nsub,E,NU,RHO,t0,NaN,0);
Ks0=sparse(ctx0.Rows(:),ctx0.Cols(:),reshape(ctx0.Bke.*t0^3+ctx0.Ske.*t0,[],1),ctx0.ndof,ctx0.ndof);
M0 =sparse(ctx0.Rows(:),ctx0.Cols(:),reshape(ctx0.Me(:)*t0,[],1),ctx0.ndof,ctx0.ndof);
[~,D]=eigs(Ks0(ctx0.free,ctx0.free),M0(ctx0.free,ctx0.free),3,'sm');
f1=min(sqrt(abs(diag(D))))/(2*pi); freq=0.97*f1; zeta=0.05;
ctx=articleB_build_ctx(nsub,E,NU,RHO,t0,freq,zeta);
lg('(2) sensitivity FD check (nsub=%s, f=%.4f Hz, zeta=%.2f):\n',mat2str(nsub),freq,zeta);
t=t0*ones(ctx.ne,1); rng(3); pe=randperm(ctx.ne,10); pert=zeros(ctx.ne,1); pert(pe)=1; h=1e-7;
for m={'dyncomp','inputpower','quadvel'}
  [J,dJ]=articleB_shell_eval(ctx,t,m{1});
  Jp=articleB_shell_eval(ctx,t+h*pert,m{1});
  Jm=articleB_shell_eval(ctx,t-h*pert,m{1});
  fd=(Jp-Jm)/(2*h); an=sum(dJ(pe));
  lg('    %-11s J=%.5e  FD=%+.6e  analytic=%+.6e  rel err=%.2e  %s\n',...
    m{1},J,fd,an,abs(fd-an)/max(abs(fd),eps), ternary(abs(fd-an)/max(abs(fd),eps)<1e-4,'OK','CHECK'));
end
lg('\nARTICLEB_VALIDATE_DONE\n');
fclose(fid); fprintf('wrote %s\n',fullfile(RES,'articleB_validation.txt'));

function both(fid,varargin), fprintf(fid,varargin{:}); fprintf(1,varargin{:}); end
function s=ternary(c,a,b), if c, s=a; else, s=b; end, end
