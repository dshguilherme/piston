function samples = r0_surface_samples(ctx, u, nquad)
%R0_SURFACE_SAMPLES Evaluate the actual spline displacement dot local normal.
% Acoustic quadrature can be refined without rebuilding the structural state.
% u may contain multiple fixed displacement fields as columns.
if nargin<3, nquad=ctx.spec.nquad; end
spec=ctx.spec;
knots0={[0 0 0 0 1 1 1 1],[0 0 0 0 1 1 1 1]};
[knots,breaks]=kntrefine(knots0,spec.nsub-1,spec.deg,spec.reg);
rule=msh_gauss_nodes(nquad); [qn,qw]=msh_set_quad_nodes(breaks,rule);
mesh=msh_cartesian(breaks,qn,qw,ctx.geometry);
% GeoPDEs spaces cache basis samples for the constructor's quadrature grid.
% Recreate the SAME basis with new sample nodes, not a new structural solve.
periodic=ctx.space.scalar_spaces{1}.periodic_dir;
if ~isempty(periodic), knots=kntunclamp(knots,spec.deg,spec.reg,periodic); end
scalar=sp_bspline(knots,spec.deg,mesh,[],periodic);
space=sp_vector(repmat({scalar},1,3),mesh);
assert(space.ndof==ctx.ndof,'Quadrature refresh changed the coefficient space.');
m=msh_precompute(mesh);
scalar_samples=sp_precompute_param(space.scalar_spaces{1},m,'value',true);
nsmax=scalar_samples.nsh_max; Ns=space.scalar_spaces{1}.ndof;
nq=m.nqn; ne=m.nel; nf=size(u,2);
y=complex(zeros(nq,ne,nf));
bi=zeros(ne*nq*nsmax*3,1); bj=bi; bv=bi; p=0;
for e=1:ne
    ns=scalar_samples.nsh(e); ids=scalar_samples.connectivity(1:ns,e);
    normal=reshape(m.normal(:,:,e),3,nq);
    for component=1:3
        vals=reshape(scalar_samples.shape_functions(:,1:ns,e),nq,ns);
        local=normal(component,:)'.*vals;
        gdofs=ids+(component-1)*Ns;
        y(:,e,:)=y(:,e,:)+reshape(local*u(gdofs,:),nq,1,nf);
        take=p+(1:nq*ns); bi(take)=repmat((e-1)*nq+(1:nq)',ns,1);
        bj(take)=reshape(repmat(gdofs(:)',nq,1),[],1);
        bv(take)=local(:); p=p+nq*ns;
    end
end
weights=m.jacdet.*m.quad_weights;
B=sparse(bi(1:p),bj(1:p),bv(1:p),nq*ne,ctx.ndof);
samples=struct('x',reshape(m.geo_map,3,[]),'weights',weights(:), ...
    'normal_displacement',reshape(y,nq*ne,nf),'basis',B,'nq',nq,'ne',ne, ...
    'area',sum(weights,1)', 'nquad',nquad);
end
