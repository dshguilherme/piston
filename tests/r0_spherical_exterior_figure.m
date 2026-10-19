function r0_spherical_exterior_figure()
% Plot the prescribed-source sinc to exact exterior power ratio for a sphere.
root=fileparts(fileparts(mfilename('fullpath'))); out=fullfile(root,'results','r0_revision');
T=readtable(fullfile(out,'spherical-exterior-benchmark.csv'));
x=linspace(0.02,8,1600); ratio=2*sin(x).^2.*(1+x.^2)./x.^2;
fig=figure('Color','w','Units','centimeters','Position',[2 2 16 10]); ax=axes(fig); hold(ax,'on');
ax.Toolbar.Visible='off'; ax.Color='w'; ax.FontName='Arial'; ax.FontSize=9; ax.LineWidth=.9;
plot(ax,x,ratio,'-','Color',[.05 .35 .62],'LineWidth',1.8,'DisplayName','Analytical ratio');
plot(ax,T.ka,T.prescribed_sinc_to_exterior_power_ratio,'o','Color',[.80 .28 .12], ...
 'LineWidth',1.3,'MarkerFaceColor','w','MarkerSize',5,'DisplayName','Direct surface quadrature');
yline(ax,1,'--k','Equal power','LabelHorizontalAlignment','left','HandleVisibility','off');
xline(ax,pi,':','Color',[.4 .4 .4],'HandleVisibility','off');
text(ax,pi+.10,2.55,'ka=3.14: sinc-model zero','FontSize',8,'Color',[.2 .2 .2]);
xlabel(ax,'ka'); ylabel(ax,'Prescribed-source sinc power / exact exterior power');
title(ax,'A curved-surface sinc kernel is not the exterior Neumann operator','FontWeight','normal','Color',[.08 .08 .08]);
legend(ax,'Location','northeast','Box','off','Color','w','TextColor',[.1 .1 .1]);
grid(ax,'on'); ax.GridAlpha=.16; ax.XLim=[0 8]; ax.YLim=[0 3.1]; ax.XColor=[.1 .1 .1]; ax.YColor=[.1 .1 .1];
exportgraphics(fig,fullfile(out,'spherical-exterior-contrast.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(out,'spherical-exterior-contrast.png'),'Resolution',400); close(fig);
fprintf('R0_SPHERICAL_EXTERIOR_FIGURE_DONE\n');
end
