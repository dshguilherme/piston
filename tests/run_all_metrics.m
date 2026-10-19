function run_all_metrics(nsub, nmax)
% RUN_ALL_METRICS End-to-end benchmark: optimise the harmonic plate for each of
% the three objectives (dynamic compliance, input power, quadratic velocity) on
% the SAME mesh/load/frequency, then render readable design + trajectory +
% volume figures, a side-by-side comparison, a table, and save a .mat.
%
%   run_all_metrics([128 64], 100)
%
% Benchmarks the OTHER two metrics through the shared building blocks
% (tests/run_metric_opt.m); src/ is untouched.

if nargin < 1 || isempty(nsub), nsub = [128 64]; end
if nargin < 2 || isempty(nmax), nmax = 100;      end

addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');

RES = 'C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
tag = sprintf('%dx%d', nsub(1), nsub(2));

% ---- benchmark definition (identical to the driver's harness) --------------
E = 210e9; RHO = 7850; NU = 0.3;
pd = cantilever_beam(1.0, 1.0);
md = struct('degree',[3 3],'regularity',[2 2],'nsub',nsub,'nquad',[4 4]);
[gm, mh, sp] = buildSpaces(pd, md);
tS = tic;
H = harmonic_kd_setup(mh, sp, gm, E, NU, RHO, 'cpu', 'double', 'plate');
tsetup = toc(tS);
F = zeros(H.N,1); [~,ic] = min(abs((1:H.N).'-round(H.N/2))); F(ic) = 1;
w = 2*pi*12.09; alpha = 25; beta = 0;
fprintf('\n==== BENCHMARK mesh %s  N=%d  f=12.09 Hz  alpha=25 beta=0  setup=%.2f s ====\n', ...
        tag, H.N, tsetup);

metrics = {'dyncomp','inputpower','quadvel'};
nM = numel(metrics);
R = struct('metric',{},'tbest',{},'hist',{},'meta',{});
opt = struct('tmin',0.4,'tmax',1.0,'volfrac',0.5,'nmax',nmax,'tol',0, ...
             'stop',true,'objscale','log','optimizer','gcmma');
for k = 1:nM
    fprintf('\n---- metric %d/%d : %s ----\n', k, nM, metrics{k});
    [tbest, hist, meta] = run_metric_opt(H, F, w, alpha, beta, metrics{k}, opt);
    R(k).metric = metrics{k}; R(k).tbest = tbest; R(k).hist = hist; R(k).meta = meta;
    fprintf('  %s DONE: iters=%d evals=%d Jbest=%.6e dB=%.3f vol=%.4f tailspread=%.3f wall=%.1f s (%s)\n', ...
        metrics{k}, meta.niter, meta.numeval, meta.Jbest, 10*log10(max(meta.Jbest,realmin)), ...
        meta.final_vol, meta.tailspread, meta.walltime, meta.stopReason);
end

% ---- save numeric histories ------------------------------------------------
matpath = fullfile(RES, sprintf('harmonic_metrics_%s.mat', tag));
save(matpath, 'R', 'metrics', 'nsub', 'nmax', 'tsetup');
fprintf('\nsaved histories: %s\n', matpath);

% ---- per-metric figures ----------------------------------------------------
for k = 1:nM
    h = R(k).hist;
    fig = figure('Position',[60 60 1700 560],'Color','w');
    % (a) design
    subplot(1,3,1);
    imagesc(R(k).tbest.'); axis image; set(gca,'YDir','normal');
    cb = colorbar; cb.Label.String = 'thickness t';
    caxis([h.tmin h.tmax]);
    xlabel('element i  (x-direction)'); ylabel('element j  (y-direction)');
    title(sprintf('%s: best design  (J=%.3e)', R(k).metric, R(k).meta.Jbest), 'Interpreter','none');
    % (b) objective trajectory
    subplot(1,3,2);
    ii = 1:numel(h.J);
    semilogy(ii, R(k).hist.J, '.', 'MarkerSize', 9); hold on;
    semilogy(ii, R(k).hist.Jbest, '-', 'LineWidth', 2); grid on;
    xlabel('outer iteration'); ylabel('objective J  (log scale)');
    title(sprintf('%s: objective trajectory', R(k).metric), 'Interpreter','none');
    legend('running J', 'best J', 'Location','northeast');
    % (c) volume fraction history
    subplot(1,3,3);
    plot(ii, R(k).hist.vol, '-o', 'LineWidth', 1.5, 'MarkerSize', 3); hold on;
    yline(0.5, 'k--', 'volfrac = 0.5'); grid on;
    xlabel('outer iteration'); ylabel('volume fraction  \Sigma t / \Sigma t_{max}');
    title(sprintf('%s: volume-fraction history', R(k).metric), 'Interpreter','none');
    ylim([0.45 0.55]);
    fp = fullfile(RES, sprintf('design_%s_%s.png', R(k).metric, tag));
    print(fig, '-dpng', '-r140', fp); close(fig);
    fprintf('figure saved: %s\n', fp);
end

% ---- side-by-side comparison of the three best designs ---------------------
fig = figure('Position',[60 60 1700 560],'Color','w');
for k = 1:nM
    subplot(1,3,k);
    imagesc(R(k).tbest.'); axis image; set(gca,'YDir','normal');
    caxis([R(k).hist.tmin R(k).hist.tmax]);
    colorbar; xlabel('element i  (x)'); ylabel('element j  (y)');
    title(sprintf('%s\nJ_{best}=%.3e  (%.2f dB)', R(k).metric, ...
          R(k).meta.Jbest, 10*log10(max(R(k).meta.Jbest,realmin))), 'Interpreter','none');
end
sgtitle(sprintf('Best designs, harmonic plate, mesh %s, N=%d, f=12.09 Hz', tag, H.N));
fp = fullfile(RES, sprintf('compare_metrics_%s.png', tag));
print(fig, '-dpng', '-r140', fp); close(fig);
fprintf('figure saved: %s\n', fp);

% ---- table -----------------------------------------------------------------
hdr = {'metric','best J (raw)','best J (dB)','iters','GCMMA evals', ...
       'wall time (s)','final volfrac','tail spread','stop'};
data = cell(nM, numel(hdr));
for k = 1:nM
    data(k,:) = { R(k).metric, sprintf('%.6e', R(k).meta.Jbest), ...
        sprintf('%.3f', 10*log10(max(R(k).meta.Jbest,realmin))), ...
        sprintf('%d', R(k).meta.niter), sprintf('%d', R(k).meta.numeval), ...
        sprintf('%.1f', R(k).meta.walltime), sprintf('%.4f', R(k).meta.final_vol), ...
        sprintf('%.3f', R(k).meta.tailspread), R(k).meta.stopReason };
end
fprintf('\n==== COMPARISON TABLE (mesh %s) ====\n', tag);
fprintf('%-12s %16s %13s %7s %13s %14s %14s %12s %10s\n', hdr{:});
for k = 1:nM
    fprintf('%-12s %16s %13s %7s %13s %14s %14s %12s %10s\n', data{k,:});
end
% text file
fid = fopen(fullfile(RES, sprintf('metrics_table_%s.txt', tag)), 'w');
fprintf(fid, '%-12s %16s %13s %7s %13s %14s %14s %12s %10s\n', hdr{:});
for k = 1:nM
    fprintf(fid, '%-12s %16s %13s %7s %13s %14s %14s %12s %10s\n', data{k,:});
end
fclose(fid);
% table as a figure (rendered as monospaced text so it survives print/export)
fig = figure('Position',[60 60 1300 300],'Color','w'); ax = axes(fig); axis(ax,'off');
colw = [0.03 0.16 0.31 0.42 0.50 0.63 0.76 0.88 1.00];
yrow = 0.88; dy = 0.16;
for c = 1:numel(hdr)
    text(ax, colw(c), yrow, hdr{c}, 'FontName','Consolas','FontSize',11, ...
         'FontWeight','bold','VerticalAlignment','middle','Interpreter','none');
end
for r = 1:nM
    yy = yrow - r*dy;
    for c = 1:numel(hdr)
        text(ax, colw(c), yy, data{r,c}, 'FontName','Consolas','FontSize',11, ...
             'VerticalAlignment','middle','Interpreter','none');
    end
end
text(ax, 0.03, yrow+dy, sprintf('Harmonic plate benchmark, mesh %s, N=%d, f=12.09 Hz', tag, H.N), ...
     'FontName','Consolas','FontSize',12,'FontWeight','bold','Interpreter','none');
print(fig, '-dpng', '-r150', fullfile(RES, sprintf('metrics_table_%s.png', tag)));
close(fig);
fprintf('figure saved: %s\n', fullfile(RES, sprintf('metrics_table_%s.png', tag)));
fprintf('table text : %s\n', fullfile(RES, sprintf('metrics_table_%s.txt', tag)));

disp('RUN_ALL_METRICS_DONE');
end
