% PROF_FULL Run the real driver under the MATLAB profiler and print the
% function-level total-time breakdown (works without editing any callee).
E=210e9; RHO=7850; NU=0.3;
pd = cantilever_beam(1.0,1.0);
md = struct('degree',[3 3],'regularity',[2 2],'nsub',[64 32],'nquad',[4 4]);
[gm,mh,sp] = buildSpaces(pd,md);
H = harmonic_kd_setup(mh,sp,gm,E,NU,RHO,'cpu','double','plate');
F = zeros(H.N,1); [~,ic] = min(abs((1:H.N).'-round(H.N/2))); F(ic)=1;
opt = struct('tmin',0.4,'tmax',1.0,'volfrac',0.5,'nmax',50,'stop',false,...
             'objscale','log','optimizer','gcmma');

profile on -history
[tb,hist] = harmonic_topopt(H,F,2*pi*12.09,25,0,opt);
profile off

p = profile('info');
FT = p.FunctionTable;
[~,ord] = sort([FT.TotalTime],'descend');
fprintf('\n=== FUNCTION TOTAL TIME (top 15) ===\n');
fprintf('%-40s %10s %8s\n','function','totTime[s]','calls');
for k=1:min(15,numel(FT))
    f=FT(ord(k));
    fprintf('%-40s %10.5f %8d\n', f.FunctionName, f.TotalTime, f.NumCalls);
end
fprintf('total driver hist.times sum = %.4f s over %d iters\n', sum(hist.times), numel(hist.times));
