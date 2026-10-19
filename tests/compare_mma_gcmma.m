% COMPARE_MMA_GCMMA Side-by-side comparison of the two optimizers on identical
% problems, including the discriminating case: the real harmonic plate.
%
% WHY: for vibration objectives MMA's subproblem is not conservative away from
% the current point, which can stall or oscillate. GCMMA adds an inner
% conservatism loop. This script measures the difference instead of asserting it.
%
% Reported per case: outer iterations, final objective, feasibility, and the
% number of objective-INCREASE steps (monotonicity).
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
fprintf('\n%-30s %6s %6s %12s %11s %12s %12s %6s\n', 'case','opt','iters','best J','it@best','max up-jump','tail spread','feas');
fprintf('%s\n', repmat('-',1,96));

% ---- case A: analytic KKT problem ----
for who = {'mma','gcmma'}
    x = [1;1]; xmin = zeros(2,1); xmax = 5*ones(2,1); st = []; o = struct();
    for it = 1:400
        f0 = (x(1)-3)^2 + (x(2)-2)^2; df0 = [2*(x(1)-3); 2*(x(2)-2)];
        gv = x(1)+x(2)-4; dgv = [1;1];
        if strcmp(who{1},'gcmma')
            [xn, st] = gcmma_solver(x, xmin, xmax, f0, df0, gv, dgv, st, it, o);
        else
            [xn, st] = mma_solver(x, xmin, xmax, f0, df0, gv, dgv, st, it);
        end
        if max(abs(xn-x)) < 1e-10, x = xn; break; end
        x = xn;
        itv(it) = f0;
    end
    f = (x(1)-3)^2 + (x(2)-2)^2;
    fprintf('%-30s %6s %6d %12.5e %11s %12s %12s %6d\n', 'A analytic KKT (2.5,1.5)', who{1}, it, f, '-', '-', '-', x(1)+x(2)<=4+1e-9);
end

% ---- case D: the REAL harmonic plate, both optimizers, same start/load/freq ----
E = 210e9; RHO = 7850; NU = 0.3;
pd = cantilever_beam(1.0, 1.0);
md = struct('degree',[3 3],'regularity',[2 2],'nsub',[16 8],'nquad',[4 4]);
[gm, mh, sp] = buildSpaces(pd, md);
H = harmonic_kd_setup(mh, sp, gm, E, NU, RHO, 'cpu', 'double', 'plate');
F = zeros(H.N,1); [~,ic] = min(abs((1:H.N).' - round(H.N/2))); F(ic) = 1;
for f = [12.09 30.44]
    for who = {'mma','gcmma'}
        opt = struct('tmin',0.4,'tmax',1.0,'volfrac',0.5,'nmax',300,'tol',0,'optimizer',who{1});  % tol=0: full budget for both
        [topt, hist] = harmonic_topopt(H, F, 2*pi*f, 25, 0, opt);
        gv = sum(topt(:))/(opt.volfrac*numel(topt)*opt.tmax) - 1;
        fprintf('%-30s %6s %6d %12.5e %11d %12.3e %12.3e %6d\n', ...
            sprintf('D plate f=%.2f Hz', f), who{1}, numel(hist.J), min(hist.J), ...
            hist.iters_to_best, hist.upjump, hist.tailspread, gv<=1e-9);
        tbest = hist.tbest; J = hist.J; Jbest = hist.Jbest;   %#ok<NASGU>
        save(sprintf('cmp_%s_f%05.1f.mat', who{1}, f), 'tbest', 'J', 'Jbest', 'f', 'who');
    end
end

% ---- figure: the two BEST designs at 12.09 Hz and their trajectories ----
Sm = load('cmp_mma_f012.1.mat'); Sg = load('cmp_gcmma_f012.1.mat');
fig = figure('Position',[80 80 1250 800],'Color','w');
subplot(2,2,1); imagesc(flipud(Sm.tbest.')); axis image; colorbar;
title(sprintf('MMA best design, J_{min}=%.3e', Sm.Jbest(end))); ylabel('element j');
subplot(2,2,2); imagesc(flipud(Sg.tbest.')); axis image; colorbar;
title(sprintf('GCMMA best design, J_{min}=%.3e', Sg.Jbest(end))); ylabel('element j');
subplot(2,2,3); semilogy(Sm.J,'.','MarkerSize',4); hold on; semilogy(Sm.Jbest,'-','LineWidth',1.6);
grid on; xlabel('iteration'); ylabel('J'); title('MMA trajectory');
legend('running J','best J','Location','northeast');
subplot(2,2,4); semilogy(Sg.J,'.','MarkerSize',4); hold on; semilogy(Sg.Jbest,'-','LineWidth',1.6);
grid on; xlabel('iteration'); ylabel('J'); title('GCMMA trajectory');
legend('running J','best J','Location','northeast');
figp = 'C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results/mma_vs_gcmma_plate_12Hz.png';
print(fig,'-dpng','-r140', figp);
fprintf('figure saved: %s\n', figp);
disp('COMPARE_DONE');
