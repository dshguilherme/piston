function r0_cap_apex_regular_test()
% Probe displacement single-valuedness at the collapsed spherical-cap apex.
% Also verify the implemented boundary enforces displacement, not zero slope.
root=fileparts(fileparts(mfilename('fullpath')));
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath(fullfile(root,'src')); addpath(fullfile(root,'tests'));
spec=struct('R',10,'L',20,'alpha',deg2rad(30),'E',2e11,'NU',0.3,'RHO',7850, ...
 't0',0.1,'deg',[2 2],'reg',[1 1],'nsub',[40 24],'nquad',[4 4], ...
 'bc',2,'loadtype','pole','loaduv',[0 0],'loaddir',[0;0;1], ...
 'loadmag',1,'freq',100,'zeta',0.05);
ctx=obstacle_build_ctx('spherecap',spec); ctx.alpha=0; ctx.beta=0.1/ctx.w;
D=load(fullfile(root,'results','articleB_r0_thesisdamp_spherecap_40x24_f100_partial.mat'),'L');
j=find(arrayfun(@(q) strcmp(q.label,'quadvel-rad'),D.L),1);
assert(~isempty(j),'Expected frozen radiation-objective cap design was not found.');
t=D.L(j).tbest(:); ke=ctx.Bke.*t.^3+ctx.Ske.*t; me=ctx.Me.*t;
A=sparse(ctx.Rows(:),ctx.Cols(:),((1+1i*ctx.w*ctx.beta)*ke(:)-ctx.w^2*me(:)),ctx.ndof,ctx.ndof);
U=zeros(ctx.ndof,1); U(ctx.free)=A(ctx.free,ctx.free)\ctx.F(ctx.free);
% At the collapsed pole, a regular Cartesian displacement has a single limit
% independent of azimuth. The random field diagnoses what the FE space permits.
rng(7301); Urand=randn(ctx.ndof,1);
theta=(0:63)/64; phiList=[0 1e-6 1e-4 1e-2 0.1 0.3 0.6 0.9];
fields={U,Urand}; names={'frozen_response','random_coefficients'}; rows=cell(0,3);
for f=1:numel(fields)
  fprintf('%s\n',names{f});
  for p=phiList
    [ev,~]=sp_eval(fields{f},ctx.space,ctx.geometry,{p,theta},{'value'});
    if iscell(ev), vv=ev{1}; else, vv=ev; end
    V=reshape(vv,3,[]); mu=mean(V,2);
    spread=norm(V-mu,'fro')/sqrt(size(V,2)); scale=max(norm(mu),norm(V,'fro')/sqrt(size(V,2)));
    rel=spread/max(scale,realmin);
    rows(end+1,:)={names{f},p*spec.alpha*180/pi,rel}; %#ok<AGROW>
    fprintf(' u=%9.2e phi-deg=%7.3f rel-RMS-spread=%11.4e\n',p,p*spec.alpha*180/pi,rel);
  end
end
out=fullfile(root,'paper_jsv','supplement');
T=cell2table(rows,'VariableNames',{'field','polar_angle_deg','relative_rms_angular_spread'});
writetable(T,fullfile(out,'cap_apex_angular_spread.csv'));
% Boundary 2 is the open edge of the cap. Trace values are constrained, but
% a nonzero physical gradient confirms that this is not a full rotational clamp.
[ebd,~]=sp_eval(U,ctx.space,ctx.geometry,{1,theta},{'value','gradient'});
Vbd=reshape(ebd{1},3,[]); Gbd=reshape(ebd{2},3,3,[]);
maxTrace=max(vecnorm(Vbd,2,1)); gradRMS=sqrt(mean(sum(reshape(abs(Gbd).^2,9,[]),1)));
[evint,~]=sp_eval(U,ctx.space,ctx.geometry,{0.95,theta},{'gradient'});
if iscell(evint), gi=evint{1}; else, gi=evint; end
Gint=reshape(gi,3,3,[]); gradInt=sqrt(mean(sum(reshape(abs(Gint).^2,9,[]),1)));
fprintf('boundary max displacement %.6e; boundary gradient RMS %.6e; near-edge gradient RMS %.6e\n',maxTrace,gradRMS,gradInt);
fprintf('R0_CAP_APEX_REGULAR_TEST_DONE\n');
end
