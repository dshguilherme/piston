function r0_bem_export_fields(levels)
% Export frozen shell normal-velocity fields to triangle meshes for BEM and
% verify each export with the sinc-model power against the consistent spline
% value from frozen_field_quadrature.csv.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fileparts(root)); addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
res=fullfile(root,'results','r0_revision'); out=fullfile(res,'bem_mesh');
if ~exist(out,'dir'), mkdir(out); end
logf=fullfile(root,'paper_jsv','supplement','bem-export.log');
if exist(logf,'file'), delete(logf); end
diary(logf); cleanupObj=onCleanup(@() diary('off')); %#ok<NASGU>
if nargin<1, levels=[2 3 4 6 8]; end
rho0=1.21; c0=343;
ref=readtable(fullfile(res,'frozen_field_quadrature.csv'));
ref=ref(ref.surface_order==8 & ref.angular_order==64,:);
rows={};
for geometry={'cylpanel','spherecap'}
    kind=geometry{1};
    spec=struct('R',25,'L',50,'th0',deg2rad(40),'alpha',deg2rad(30), ...
        'E',2e11,'NU',0.3,'RHO',7850,'t0',0.1,'deg',[2 2],'reg',[1 1], ...
        'nsub',[40 24],'nquad',[4 4],'bc',[1 2 3 4],'loadtype','point', ...
        'loaduv',[0.5 0.5],'loaddir',[0;1;0],'loadmag',1,'freq',6.53547,'zeta',0.05);
    if strcmp(kind,'spherecap')
        spec.R=10; spec.bc=2; spec.loadtype='pole'; spec.loaddir=[0;0;1]; spec.freq=100;
    end
    ctx=obstacle_build_ctx(kind,spec); ctx.alpha=0; ctx.beta=0.1/ctx.w;
    S=load(fullfile(res,[kind '_frozen_states.mat'])); U=S.U; labels=S.labels(:)'; designs=S.designs;
    nd=size(U,2);
    fprintf('\n===== %s: designs in mat: %s\n',kind,strjoin(labels,', '));
    % cheap re-solve assertion
    for j=1:nd
        t=designs(:,j); KE=ctx.Bke.*(t.^3)+ctx.Ske.*t; ME=ctx.Me.*t;
        A=sparse(ctx.Rows(:),ctx.Cols(:),(1+1i*ctx.w*ctx.beta)*KE(:)-ctx.w^2*ME(:),ctx.ndof,ctx.ndof);
        Uj=complex(zeros(ctx.ndof,1)); Uj(ctx.free)=A(ctx.free,ctx.free)\ctx.F(ctx.free);
        fprintf('re-solve check %s: rel diff %.2e\n',labels{j},norm(Uj-U(:,j))/norm(U(:,j)));
    end
    om=ctx.w; k=om/c0; lambda=c0/(om/2/pi);
    sm=r0_surface_samples(ctx,U,[4 4]); spline_area=sum(sm.area);
    fprintf('freq %.5f Hz omega %.5f lambda %.4f m k %.5f; spline area (sum of weights) %.8f m^2\n', ...
        om/2/pi,om,lambda,k,spline_area);
    for m=levels
        [V,T,vn,N,bnd,info]=build_mesh(ctx,U,m,kind);
        vn=1i*om*vn;
        ext=max(max(V)-min(V)); ed=[T(:,[1 2]);T(:,[2 3]);T(:,[3 1])];
        el=sqrt(sum((V(ed(:,1),:)-V(ed(:,2),:)).^2,2));
        area=tri_area(V,T);
        fprintf('\n-- %s m=%d: Nv=%d Nt=%d (removed degenerate %d) area %.8f (spline %.8f, rel %.2e)\n', ...
            kind,m,size(V,1),size(T,1),info.ndegenerate,sum(area),spline_area,sum(area)/spline_area-1);
        fprintf('   max extent %.4f m, k*extent %.4f; edge mean %.4f max %.4f m; h_max/lambda %.4f; elements/lambda (max edge) %.2f\n', ...
            ext,k*ext,mean(el),max(el),max(el)/lambda,lambda/max(el));
        fprintf('   orientation: %d of %d triangles initially disagree with spline normal (flipped to fix); after fix %d\n', ...
            info.nflip,size(T,1),info.nflip_after);
        fprintf('   duplicate-sample consistency: max rel |dvn| %.2e, max |dx| %.2e, max |dn| %.2e; analytic-position check %.2e\n', ...
            info.dup_vn,info.dup_x,info.dup_n,info.xcheck);
        if ~isempty(info.pole)
            fprintf('   pole (merged %d samples): relative spread of v_n per design (max|v-mean|/max(|mean|,eps)), and absolute (m/s):\n',info.pole.n);
            for j=1:nd
                fprintf('      %-14s rel %.3e  abs %.3e  |mean| %.3e\n',labels{j},info.pole.rel(j),info.pole.abs(j),info.pole.mean(j));
            end
        end
        fprintf('   boundary vertices: %d\n',nnz(bnd));
        P=sinc_power(V,T,vn,om,rho0,c0,32);
        for j=1:nd
            iref=strcmp(ref.geometry,kind)&strcmp(ref.design,labels{j});
            Pc=ref.consistent_spline_W(iref);
            fprintf('   %-20s P_mesh=%.9e P_ref=%.9e rel=%+.3e\n',labels{j},P(j),Pc,P(j)/Pc-1);
            rows(end+1,:)={kind,labels{j},m,size(V,1),size(T,1),sum(area),P(j),Pc,abs(P(j)-Pc)/Pc}; %#ok<AGROW>
        end
        vertices=V; triangles=int32(T); vn_out=vn; normals=N; boundary_vertex=bnd; freq=om/2/pi; omega=om; %#ok<NASGU>
        rho0_=rho0; c0_=c0; %#ok<NASGU>
        S2=struct('vertices',V,'triangles',double(T),'vn',vn,'normals',N,'labels',{labels}, ...
            'freq',om/2/pi,'omega',om,'rho0',rho0,'c0',c0,'boundary_vertex',bnd);
        save(fullfile(out,sprintf('%s_m%d.mat',kind,m)),'-struct','S2','-v7');
    end
end
Tab=cell2table(rows,'VariableNames',{'geometry','design','m','n_vertices','n_triangles','area', ...
    'P_sinc_mesh','P_consistent_ref','rel_diff'});
writetable(Tab,fullfile(out,'export_sinc_check.csv'));
fprintf('\n===== convergence table (rel_diff) =====\n');
for geometry={'cylpanel','spherecap'}
    g=geometry{1}; ds=unique(Tab.design(strcmp(Tab.geometry,g)),'stable');
    for i=1:numel(ds)
        s=Tab(strcmp(Tab.geometry,g)&strcmp(Tab.design,ds{i}),:);
        fprintf('%-9s %-20s',g,ds{i}); fprintf('  m=%d:%.3e',[s.m s.rel_diff]'); fprintf('\n');
    end
end
fprintf('R0_BEM_EXPORT_DONE\n');
end

function a=tri_area(V,T)
c=cross(V(T(:,2),:)-V(T(:,1),:),V(T(:,3),:)-V(T(:,1),:),2);
a=0.5*sqrt(sum(c.^2,2));
end

function [V,T,vn,N,bnd,info]=build_mesh(ctx,U,m,kind)
spec=ctx.spec; nf=size(U,2);
knots0={[0 0 0 0 1 1 1 1],[0 0 0 0 1 1 1 1]};
[knots,breaks]=kntrefine(knots0,spec.nsub-1,spec.deg,spec.reg);
nel=[numel(breaks{1})-1, numel(breaks{2})-1];
qn=cell(1,2); qw=cell(1,2);
for d=1:2
    qn{d}=zeros(m+1,nel(d)); qw{d}=zeros(m+1,nel(d));
    for e=1:nel(d), qn{d}(:,e)=linspace(breaks{d}(e),breaks{d}(e+1),m+1)'; qw{d}(:,e)=(breaks{d}(e+1)-breaks{d}(e))/(m+1);
        % GeoPDEs findspan puts a right-end node in the NEXT span; shift it inside (position error ~1e-11)
        qn{d}(end,e)=qn{d}(end,e)-1e-11*(breaks{d}(e+1)-breaks{d}(e)); end
end
mesh=msh_cartesian(breaks,qn,qw,ctx.geometry);
periodic=ctx.space.scalar_spaces{1}.periodic_dir;
if ~isempty(periodic), knots=kntunclamp(knots,spec.deg,spec.reg,periodic); end
scalar=sp_bspline(knots,spec.deg,mesh,[],periodic);
space=sp_vector(repmat({scalar},1,3),mesh);
assert(space.ndof==ctx.ndof,'coefficient space changed');
mp=msh_precompute(mesh);
ss=sp_precompute_param(space.scalar_spaces{1},mp,'value',true);
Ns=space.scalar_spaces{1}.ndof; nq=mp.nqn; ne=mp.nel;
y=complex(zeros(nq,ne,nf)); X=reshape(mp.geo_map,3,nq,ne); Nrm=zeros(3,nq,ne);
for e=1:ne
    ns=ss.nsh(e); ids=ss.connectivity(1:ns,e);
    nrm=reshape(mp.normal(:,:,e),3,nq);
    bad=any(isnan(nrm),1)|vecnorm(nrm)<0.5;
    if any(bad)
        assert(strcmp(kind,'spherecap'),'degenerate normal on non-cap');
        xb=X(:,bad,e); nrm(:,bad)=xb./vecnorm(xb);   % outward radial = spline normal direction (checked below)
    end
    Nrm(:,:,e)=nrm;
    vals=reshape(ss.shape_functions(:,1:ns,e),nq,ns);
    for c=1:3
        gd=ids+(c-1)*Ns;      % corrected component-aware indexing
        y(:,e,:)=y(:,e,:)+reshape((nrm(c,:)'.*vals)*U(gd,:),nq,1,nf);
    end
end
% global grid indices
NI=nel(1)*m+1; NJ=nel(2)*m+1; NG=NI*NJ;
gidx=zeros(nq,ne);
for e=1:ne
    [e1,e2]=ind2sub(nel,e);
    [i,j]=ndgrid(0:m,0:m);
    gidx(:,e)=((e1-1)*m+i(:)+1)+(((e2-1)*m+j(:))*NI);
end
gi=gidx(:); Yall=reshape(y,nq*ne,nf); Xall=reshape(X,3,[]).'; Nall=reshape(Nrm,3,[]).';
Yg=zeros(NG,nf); Xg=zeros(NG,3); Ng=zeros(NG,3);
Yg(gi,:)=Yall; Xg(gi,:)=Xall; Ng(gi,:)=Nall;
info.dup_vn=max(abs(Yall-Yg(gi,:)),[],'all')/max(abs(Yall),[],'all');
info.dup_x=max(abs(Xall-Xg(gi,:)),[],'all'); info.dup_n=max(abs(Nall-Ng(gi,:)),[],'all');
% analytic coordinate check of grid ordering
ug=[];vg=[];
for d=1:2
    pts=zeros(1,nel(d)*m+1);
    for e=1:nel(d), pts((e-1)*m+(1:m+1))=linspace(breaks{d}(e),breaks{d}(e+1),m+1); end
    if d==1, ug=pts; else, vg=pts; end
end
Xa=ctx.geometry.map({ug,vg}).';
info.xcheck=max(abs(Xa-Xg),[],'all');
assert(info.xcheck<1e-7,'grid ordering mismatch');
ext=max(max(Xg)-min(Xg)); tol=1e-6*ext;
% merge coincident vertices
key=round(Xg/tol);
[~,ia,ic]=unique(key,'rows','stable');
nv=numel(ia);
if strcmp(kind,'spherecap'), assert(nv==1+nel(1)*m*nel(2)*m,'cap vertex count');
else, assert(nv==NG,'panel vertex count'); end
cnt=accumarray(ic,1);
V=zeros(nv,3); N=zeros(nv,3); vn=complex(zeros(nv,nf));
for d=1:3, V(:,d)=accumarray(ic,Xg(:,d))./cnt; N(:,d)=accumarray(ic,Ng(:,d)); end
for j=1:nf
    vn(:,j)=accumarray(ic,real(Yg(:,j)))./cnt+1i*accumarray(ic,imag(Yg(:,j)))./cnt;
end
N=N./vecnorm(N,2,2);
info.pole=[];
if strcmp(kind,'spherecap')
    pg=find(Xg(:,3)>max(Xg(:,3))-tol); % pole: all samples at the apex
    assert(all(ic(pg)==ic(pg(1))),'pole not single merged vertex');
    yp=Yg(pg,:); mu=mean(yp,1);
    info.pole=struct('n',numel(pg),'abs',max(abs(yp-mu),[],1), ...
        'mean',abs(mu),'rel',max(abs(yp-mu),[],1)./max(abs(mu),eps));
end
% triangulation
vid=reshape(ic,NI,NJ); T=zeros(2*(NI-1)*(NJ-1),3); p=0;
for J=1:NJ-1
    for I=1:NI-1
        a=vid(I,J); b=vid(I+1,J); c=vid(I+1,J+1); d=vid(I,J+1);
        p=p+1; T(p,:)=[a b c]; p=p+1; T(p,:)=[a c d];
    end
end
deg=T(:,1)==T(:,2)|T(:,2)==T(:,3)|T(:,1)==T(:,3);
info.ndegenerate=nnz(deg); T=T(~deg,:);
% orientation vs spline normal
cr=cross(V(T(:,2),:)-V(T(:,1),:),V(T(:,3),:)-V(T(:,1),:),2);
nm=N(T(:,1),:)+N(T(:,2),:)+N(T(:,3),:);
neg=sum(cr.*nm,2)<0; info.nflip=nnz(neg);
T(neg,[2 3])=T(neg,[3 2]);
cr=cross(V(T(:,2),:)-V(T(:,1),:),V(T(:,3),:)-V(T(:,1),:),2);
info.nflip_after=nnz(sum(cr.*nm,2)<0);
% remove unreferenced vertices? (none expected)
assert(numel(unique(T(:)))==nv,'unreferenced vertices');
% boundary vertices: edges used once
E=sort([T(:,[1 2]);T(:,[2 3]);T(:,[3 1])],2);
[ue,~,ie]=unique(E,'rows'); ec=accumarray(ie,1);
be=ue(ec==1,:); bnd=false(nv,1); bnd(be(:))=true;
assert(all(ec<=2),'non-manifold edge');
end

function P=sinc_power(V,T,vn,om,rho,c,na)
% P1 interpolation of vn, 6-point triangle rule, angular (plane-wave) form of j0(kr).
a=0.445948490915965; b=0.091576213509771; wa=0.223381589678011; wb=0.109951743655322;
L=[1-2*a a a; a 1-2*a a; a a 1-2*a; 1-2*b b b; b 1-2*b b; b b 1-2*b]; Wq=[wa wa wa wb wb wb]';
nt=size(T,1); A=V(T(:,1),:); B=V(T(:,2),:); C=V(T(:,3),:);
ar=0.5*sqrt(sum(cross(B-A,C-A,2).^2,2));
nq=6*nt; xq=zeros(nq,3); wq=zeros(nq,1); nf=size(vn,2); vq=complex(zeros(nq,nf));
for q=1:6
    r=q:6:nq;
    xq(r,:)=L(q,1)*A+L(q,2)*B+L(q,3)*C; wq(r)=Wq(q)*ar;
    vq(r,:)=L(q,1)*vn(T(:,1),:)+L(q,2)*vn(T(:,2),:)+L(q,3)*vn(T(:,3),:);
end
i=(1:na-1)'; beta=i./sqrt(4*i.^2-1);
[Vv,D]=eig(diag(beta,1)+diag(beta,-1)); [mu,o]=sort(diag(D)); wm=2*(Vv(1,o).^2)';
nphi=2*na; phi=2*pi*(0:nphi-1)/nphi; [MU,PHI]=ndgrid(mu,phi);
s=[sqrt(1-MU(:).^2).*cos(PHI(:)), sqrt(1-MU(:).^2).*sin(PHI(:)), MU(:)];
aw=repmat(wm,nphi,1)/(2*nphi); k=om/c; pref=rho*om^2/(4*pi*c);
wv=wq.*vq; P=zeros(1,nf); bs=max(1,floor(2e7/nq));
for f=1:bs:size(s,1)
    r=f:min(f+bs-1,size(s,1));
    H=exp(-1i*k*(s(r,:)*xq.'))*wv;           % ndir_block x nf
    P=P+sum(aw(r).*abs(H).^2,1);
end
P=pref*P;
end

