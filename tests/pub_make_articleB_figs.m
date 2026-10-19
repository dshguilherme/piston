function pub_make_articleB_figs
% PUB_MAKE_ARTICLEB_FIGS  Publication figures for Article B (harmonic topology
% optimization of a curved Scordelis-Lo roof shell).  All numbers come from
%   results/articleB_shell_metrics_128x64.mat  (variable R)
%   results/articleB_metrics_table_128x64.txt
% Writes results/pub_fig_B1_designs.{png,pdf}
%        results/pub_fig_B2_trajectories.{png,pdf}
%        results/pub_fig_B3_degeneracy.{png,pdf}
%        results/pub_fig_B4_metrics.{png,pdf}

addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
RES = 'C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
S = pub_style();
out = @(b) fullfile(RES, b);

fprintf('loading R from articleB_shell_metrics_128x64.mat ...\n');
L = load(fullfile(RES,'articleB_shell_metrics_128x64.mat'), 'R','nsub','freq','f1');
R = L.R; nsub = L.nsub; f = L.freq;
n1 = nsub(1); n2 = nsub(2);
lab = {'Dynamic compliance','Input power','Quadratic velocity'};
Mcol = {S.c.b1, S.c.b2, S.c.b3};
cmap = parula(256);   % perceptually uniform sequential map

% ============================================================ B1 designs
fig = figure('Units','inches','Position',[1 1 7.6 2.6]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
tlo = min(arrayfun(@(k) min(R(k).tbest), 1:3));
thi = max(arrayfun(@(k) max(R(k).tbest), 1:3));
for k = 1:3
    nexttile;
    tg = reshape(R(k).tbest, n1, n2);
    Ce = tg.'; Ce(end+1,:) = Ce(end,:); Ce(:,end+1) = Ce(:,end);
    pcolor(0.5:1:(n1+0.5), 0.5:1:(n2+0.5), Ce); shading flat;   % vector heatmap
    axis image; set(gca,'YDir','normal','Layer','top');
    colormap(gca, cmap); caxis([tlo thi]);
    xlabel('element $i$ (circumferential $\xi$)','Interpreter','latex','FontSize',S.fs_lab-1);
    ylabel('element $j$ (axial $\eta$)','Interpreter','latex','FontSize',S.fs_lab-1);
    ttl = sprintf('%s\n$J^{\\star}=%.3f$ dB', lab{k}, 10*log10(R(k).meta.Jbest));
    title(ttl,'Interpreter','latex','FontSize',S.fs_title,'FontWeight','normal');
    set(gca,'FontSize',S.fs-1);
end
cb = colorbar; cb.Layout.Tile = 'east';
cb.Label.String = 'thickness t (m)';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = S.fs_lab;
cb.Color = [0 0 0]; cb.Label.Color = [0 0 0]; cb.TickLabelInterpreter = 'latex';
tB1 = title(tl, sprintf('Optimal thickness fields, curved Scordelis--Lo roof shell ($%d\\times%d$, $f=%.3f$ Hz)',n1,n2,f), ...
      'Interpreter','latex','FontSize',S.fs_title,'FontWeight','normal'); set(tB1,'Color',[0 0 0]);
pub_save(fig, out('pub_fig_B1_designs')); close(fig);
fprintf('B1 written\n');

% ============================================================ B2 trajectories
fig = figure('Units','inches','Position',[1 1 7.6 2.7]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for k = 1:3
    nexttile; hold on; grid on; box on;
    h = R(k).hist; ii = 1:numel(h.J);
    plot(ii, h.J, '.', 'Color', Mcol{k}, 'MarkerSize', 8, 'DisplayName','running $J$');
    plot(ii, max(h.Jbest, realmin), '-', 'Color', [0.15 0.15 0.15], 'LineWidth', S.lw, ...
         'DisplayName','best $J$');
    set(gca,'YScale','log','FontSize',S.fs);
    xlabel('outer iteration','Interpreter','latex','FontSize',S.fs_lab);
    if k==1, ylabel('objective $J$ (log scale)','Interpreter','latex','FontSize',S.fs_lab); end
    title(sprintf('%s ($J^\\star=%.3f$ dB)', lab{k}, 10*log10(R(k).meta.Jbest)), ...
          'Interpreter','latex','FontSize',S.fs_title-1,'FontWeight','normal');
    if k==1, legend('Location','northeast','Interpreter','latex','FontSize',S.fs_leg); end
end
tB2 = title(tl,'Objective trajectories','Interpreter','latex','FontSize',S.fs_title,'FontWeight','normal'); set(tB2,'Color',[0 0 0]);
pub_save(fig, out('pub_fig_B2_trajectories')); close(fig);
fprintf('B2 written\n');

% ============================================================ B3 degeneracy
fig = figure('Units','inches','Position',[1 1 9.6 2.6]);
tl = tiledlayout(1,4,'TileSpacing','compact','Padding','compact');
for k = 1:3
    nexttile; hold on; box on; grid on;
    histogram(R(k).tbest, 20, 'FaceColor', Mcol{k}, 'EdgeColor','w');
    xline(min(R(k).tbest),'k--','LineWidth',0.9);
    xline(max(R(k).tbest),'k--','LineWidth',0.9);
    xlabel('thickness $t$ (m)','Interpreter','latex','FontSize',S.fs_lab-1);
    if k==1, ylabel('element count','Interpreter','latex','FontSize',S.fs_lab-1); end
    title(sprintf('%s\n$\\mathrm{CV}=%.2f$, $\\sigma=%.3f$', lab{k}, R(k).deg.CV, R(k).deg.std), ...
          'Interpreter','latex','FontSize',S.fs_title-1,'FontWeight','normal');
    xlim([0.24 1.02]); set(gca,'FontSize',S.fs-1);
end
% concentration metric: frac of deviation in the hottest 1% of elements
nexttile; hold on; box on; grid on;
cc = arrayfun(@(k) R(k).deg.conc1pct, 1:3);
b = bar(1:3, cc, 0.55, 'FaceColor','flat');
for k=1:3, b.CData(k,:) = Mcol{k}; end
yline(1.0,'--','Color',S.c.WQ,'LineWidth',S.lw_ref);
text(0.55, 1.04,'degenerate design $=1.00$','Color',S.c.WQ,'FontSize',8,'Interpreter','latex');
yline(0.01,'-.','Color',[0.35 0.35 0.35],'LineWidth',1.0);
text(0.55, 0.20,'uniform ideal $=0.01$','Color',[0.35 0.35 0.35],'FontSize',8,'Interpreter','latex');
for k=1:3
    text(k, cc(k)+0.045, sprintf('%.3f',cc(k)), 'HorizontalAlignment','center', ...
         'FontSize',9,'Interpreter','none');
end
set(gca,'XTick',1:3,'XTickLabel',{'dyn. comp.','input power','quad. vel.'}, ...
    'TickLabelInterpreter','none','FontSize',S.fs-1);
ylim([0 1.12]);
ylabel('deviation fraction in hottest 1\% of elements','Interpreter','latex','FontSize',S.fs_lab-1);
title({'Concentration metric','(distributed, not degenerate)'}, ...
      'Interpreter','latex','FontSize',S.fs_title-1,'FontWeight','normal');
tB3 = title(tl,'Non-degeneracy of the three optimal designs ($128\times64$)','Interpreter','latex', ...
      'FontSize',S.fs_title,'FontWeight','normal'); set(tB3,'Color',[0 0 0]);
pub_save(fig, out('pub_fig_B3_degeneracy')); close(fig);
fprintf('B3 written\n');

% ============================================================ B4 metrics
Jbest = arrayfun(@(k) R(k).meta.Jbest, 1:3);
JdB   = 10*log10(Jbest);
fig = figure('Units','inches','Position',[1 1 4.6 3.2]);
hold on; box on; grid on;
b = bar(1:3, JdB, 0.55, 'FaceColor','flat');
for k=1:3, b.CData(k,:) = Mcol{k}; end
for k=1:3
    text(k, JdB(k), sprintf('  %.3f dB\n  ($J=%.3g$)', JdB(k), Jbest(k)), ...
         'FontSize',9,'Interpreter','latex','HorizontalAlignment','center', ...
         'VerticalAlignment','bottom');
end
set(gca,'XTick',1:3,'XTickLabel',{'dyn. comp.','input power','quad. vel.'}, ...
    'TickLabelInterpreter','none','FontSize',S.fs);
ylim([min(JdB)-6, max(JdB)+12]);
ylabel('best objective $J^{\star}$ (dB)','Interpreter','latex','FontSize',S.fs_lab);
title('Best objective per metric (128$\times$64)','Interpreter','latex', ...
      'FontSize',S.fs_title,'FontWeight','normal');
pub_save(fig, out('pub_fig_B4_metrics')); close(fig);
fprintf('B4 written\n');

fprintf('ARTICLEB_PUB_FIGS_DONE\n');
end
