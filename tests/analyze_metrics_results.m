% POST_PROCESS: recompute a settled-window tail spread from the saved histories,
% regenerate the table (with both spreads), add a combined-trajectory figure,
% and update the .mat. Reads results/harmonic_metrics_128x64.mat.
RES = 'C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
S = load(fullfile(RES,'harmonic_metrics_128x64.mat'));
R = S.R; nM = numel(R); tag = sprintf('%dx%d', S.nsub(1), S.nsub(2));
metrics = S.metrics; nsub = S.nsub; nmax = S.nmax; tsetup = S.tsetup;

% settled-window spread: last min(n,10) iterations (the plateau after the stop
% criterion was satisfied). Driver's "tail spread" (last 50, or all if fewer)
% is kept for reference but includes the initial transient on short runs.
for k = 1:nM
    J = R(k).hist.J(:); n = numel(J);
    w10 = J(max(1,n-9):end);
    R(k).meta.tailspread10 = max(w10)/max(min(w10),realmin);
    w50 = J(max(1,n-49):end);
    R(k).meta.tailspread50 = max(w50)/max(min(w50),realmin);
    fprintf('DIAG %-11s n=%d  last10 J: %s\n', R(k).metric, n, ...
            strjoin(arrayfun(@(v) sprintf('%.4e',v), w10(:).', 'uni', 0), ' '));
end

hdr = {'metric','best J (raw)','best J (dB)','iters','GCMMA evals', ...
       'wall time (s)','final volfrac','tail spread (last10)','tail spread (last50)','stop'};
data = cell(nM, numel(hdr));
for k = 1:nM
    m = R(k).meta;
    data(k,:) = { R(k).metric, sprintf('%.6e', m.Jbest), ...
        sprintf('%.3f', 10*log10(max(m.Jbest,realmin))), ...
        sprintf('%d', m.niter), sprintf('%d', m.numeval), ...
        sprintf('%.1f', m.walltime), sprintf('%.4f', m.final_vol), ...
        sprintf('%.6f', m.tailspread10), sprintf('%.3e', m.tailspread50), m.stopReason };
end

fprintf('\n==== FINAL COMPARISON TABLE (mesh %s, f=12.09 Hz) ====\n', tag);
fprintf('%-12s %15s %12s %6s %12s %13s %14s %20s %18s %10s\n', hdr{:});
for k = 1:nM
    fprintf('%-12s %15s %12s %6s %12s %13s %14s %20s %18s %10s\n', data{k,:});
end

fid = fopen(fullfile(RES, sprintf('metrics_table_%s.txt',tag)),'w');
fprintf(fid,'Harmonic plate benchmark, mesh %s, f=12.09 Hz, alpha=25, beta=0, volfrac=0.5\n\n', tag);
fprintf(fid,'%-12s %15s %12s %6s %12s %13s %14s %20s %18s %10s\n', hdr{:});
for k = 1:nM
    fprintf(fid,'%-12s %15s %12s %6s %12s %13s %14s %20s %18s %10s\n', data{k,:});
end
fclose(fid);

fig = figure('Position',[40 40 1900 360],'Color','w'); ax = axes(fig); axis(ax,'off');
colw = [0.005 0.10 0.235 0.345 0.42 0.53 0.635 0.76 0.90 1.00];
yrow = 0.84; dy = 0.22;
for c = 1:numel(hdr)
    text(ax, colw(c), yrow, hdr{c}, 'FontName','Consolas','FontSize',12, 'Color','k', ...
         'FontWeight','bold','VerticalAlignment','middle','Interpreter','none');
end
for r = 1:nM
    yy = yrow - r*dy;
    for c = 1:numel(hdr)
        text(ax, colw(c), yy, data{r,c}, 'FontName','Consolas','FontSize',12, 'Color','k', ...
             'VerticalAlignment','middle','Interpreter','none');
    end
end
print(fig,'-dpng','-r150', fullfile(RES, sprintf('metrics_table_%s.png',tag))); close(fig);
fprintf('table figure/ txt written.\n');

% combined best-J trajectory
fig = figure('Position',[60 60 900 560],'Color','w');
hold on; cols = lines(nM);
for k = 1:nM
    semilogy(1:numel(R(k).hist.Jbest), R(k).hist.Jbest, '-', 'LineWidth', 2, 'Color', cols(k,:));
end
grid on; set(gca,'YScale','log');
xlabel('outer iteration'); ylabel('best objective J_{best}  (log scale)');
legend(R(1).metric, R(2).metric, R(3).metric, 'Location','northeast');
title(sprintf('Best-so-far objective, three metrics (mesh %s)', tag));
print(fig,'-dpng','-r140', fullfile(RES, sprintf('trajectories_%s.png',tag))); close(fig);
fprintf('trajectories figure written.\n');

% ---- regenerate per-metric figures + comparison with TRUE dB (no eps clamp) --
for k = 1:nM
    h = R(k).hist; Jb = R(k).meta.Jbest; db = 10*log10(max(Jb,realmin));
    fig = figure('Position',[60 60 1700 560],'Color','w');
    subplot(1,3,1);
    imagesc(R(k).tbest.'); axis image; set(gca,'YDir','normal');
    cb = colorbar; cb.Label.String = 'thickness t';
    caxis([h.tmin h.tmax]);
    xlabel('element i  (x-direction)'); ylabel('element j  (y-direction)');
    title(sprintf('%s: best design  (J=%.3e = %.2f dB)', R(k).metric, Jb, db), 'Interpreter','none');
    subplot(1,3,2);
    ii = 1:numel(h.J);
    semilogy(ii, h.J, '.', 'MarkerSize', 9, 'Color',[0 0.45 0.74]); hold on;
    semilogy(ii, h.Jbest, '-', 'LineWidth', 2, 'Color',[0.85 0.33 0.10]); grid on;
    xlabel('outer iteration'); ylabel('objective J  (log scale)');
    title(sprintf('%s: objective trajectory', R(k).metric), 'Interpreter','none');
    legend('running J','best J','Location','northeast');
    subplot(1,3,3);
    plot(ii, h.vol, '-o', 'LineWidth', 1.5, 'MarkerSize', 3); hold on;
    yline(0.5, 'k--', 'volfrac = 0.5'); grid on; ylim([0.45 0.55]);
    xlabel('outer iteration'); ylabel('volume fraction  \Sigma t / \Sigma t_{max}');
    title(sprintf('%s: volume-fraction history', R(k).metric), 'Interpreter','none');
    print(fig,'-dpng','-r140', fullfile(RES, sprintf('design_%s_%s.png', R(k).metric, tag))); close(fig);
    fprintf('figure saved: %s\n', fullfile(RES, sprintf('design_%s_%s.png', R(k).metric, tag)));
end

fig = figure('Position',[60 60 1700 560],'Color','w');
for k = 1:nM
    subplot(1,3,k);
    imagesc(R(k).tbest.'); axis image; set(gca,'YDir','normal');
    caxis([R(k).hist.tmin R(k).hist.tmax]); colorbar;
    xlabel('element i  (x)'); ylabel('element j  (y)');
    title(sprintf('%s\nJ_{best}=%.3e  (%.2f dB)', R(k).metric, R(k).meta.Jbest, ...
          10*log10(max(R(k).meta.Jbest,realmin))), 'Interpreter','none');
end
sgtitle(sprintf('Best designs, harmonic plate, mesh %s, N=8777, f=12.09 Hz', tag));
print(fig,'-dpng','-r140', fullfile(RES, sprintf('compare_metrics_%s.png', tag))); close(fig);
fprintf('figure saved: %s\n', fullfile(RES, sprintf('compare_metrics_%s.png', tag)));

save(fullfile(RES,sprintf('harmonic_metrics_%s.mat',tag)), 'R', 'metrics', 'nsub', 'nmax', 'tsetup');
disp('POST_PROCESS_DONE');
