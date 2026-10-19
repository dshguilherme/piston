function r0_nystrom_figures()
% Compare global kernel error against errors in cancellation-sensitive power.
root=fileparts(fileparts(mfilename('fullpath'))); out=fullfile(root,'results','r0_revision');
T=readtable(fullfile(out,'nystrom_rank_sweep.csv')); T.geometry=string(T.geometry); T.field=string(T.field);
P=T(T.geometry=="cylpanel" & T.rank<=256,:); C=T(T.geometry=="spherecap" & T.field=="quadvel-rad" & T.rank<=256,:);
fig=figure('Color','w','Units','centimeters','Position',[2 2 16 10]); ax=axes(fig); hold(ax,'on');
ax.Toolbar.Visible='off'; ax.Color='w'; ax.FontName='Arial'; ax.FontSize=9; ax.LineWidth=.9;
colors=[.80 .28 .12;.05 .35 .62];
semilogy(ax,P.rank,abs(P.relative_quadratic_error),'-o','Color',colors(1,:), ...
 'LineWidth',1.8,'MarkerFaceColor','w','DisplayName','Panel: power error');
semilogy(ax,P.rank,P.relative_kernel_frobenius_error,'--s','Color',colors(1,:), ...
 'LineWidth',1.5,'MarkerFaceColor','w','DisplayName','Panel: kernel Frobenius error');
semilogy(ax,C.rank,abs(C.relative_quadratic_error),'-o','Color',colors(2,:), ...
 'LineWidth',1.8,'MarkerFaceColor','w','DisplayName','Cap: power error');
semilogy(ax,C.rank,C.relative_kernel_frobenius_error,'--s','Color',colors(2,:), ...
 'LineWidth',1.5,'MarkerFaceColor','w','DisplayName','Cap: kernel Frobenius error');
yline(ax,1e-2,':k','1% power tolerance','LabelHorizontalAlignment','left','HandleVisibility','off');
xlabel(ax,'Nystrom sample rank'); ylabel(ax,'Relative error (absolute for quadratic power)');
title(ax,'A small global matrix error can hide a large power error','FontWeight','normal','Color',[.08 .08 .08]);
legend(ax,'Location','northeast','Box','off','Color','w','TextColor',[.1 .1 .1]);
grid(ax,'on'); ax.GridAlpha=.16; ax.XLim=[8 256]; ax.YLim=[1e-7 1]; ax.XColor=[.1 .1 .1]; ax.YColor=[.1 .1 .1];
exportgraphics(fig,fullfile(out,'nystrom-rank-fidelity.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(out,'nystrom-rank-fidelity.png'),'Resolution',400); close(fig);
fprintf('R0_NYSTROM_FIGURES_DONE\n');
end
