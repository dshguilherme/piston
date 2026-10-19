function r0_cap_apex_figure()
% Publication-style diagnostic of angular displacement spread on the cap.
root=fileparts(fileparts(mfilename('fullpath'))); out=fullfile(root,'results','r0_revision');
T=readtable(fullfile(root,'paper_jsv','supplement','cap_apex_angular_spread.csv'));
T.field=string(T.field); F=T(T.field=="frozen_response",:); R=T(T.field=="random_coefficients",:);
sel=[1 4 5 6 7 8]; F=F(sel,:); R=R(sel,:);
fig=figure('Color','w','Units','centimeters','Position',[2 2 16 10]); ax=axes(fig); hold(ax,'on');
ax.Toolbar.Visible='off'; ax.Color='w'; ax.FontName='Arial'; ax.FontSize=9; ax.LineWidth=.9;
semilogy(ax,F.polar_angle_deg,F.relative_rms_angular_spread,'-o','Color',[.05 .35 .62], ...
 'LineWidth',1.8,'MarkerFaceColor','w','DisplayName','Frozen cap response');
semilogy(ax,R.polar_angle_deg,R.relative_rms_angular_spread,'--s','Color',[.80 .28 .12], ...
 'LineWidth',1.6,'MarkerFaceColor','w','DisplayName','Random coefficient vector');
xlabel(ax,'Polar angle from collapsed apex (deg)');
ylabel(ax,'Relative RMS angular displacement spread');
title(ax,'A nearly regular response does not imply a regular basis','FontWeight','normal','Color',[.08 .08 .08]);
legend(ax,'Location','northoutside','Orientation','horizontal','Box','off','Color','w','TextColor',[.1 .1 .1]);
grid(ax,'on'); ax.GridAlpha=.16; ax.XMinorGrid='off'; ax.YMinorGrid='off'; ax.XLim=[-0.5 27.5]; ax.YLim=[1e-5 1.2]; ax.YScale='log';
ax.XColor=[.1 .1 .1]; ax.YColor=[.1 .1 .1];
exportgraphics(fig,fullfile(out,'cap-apex-regularity.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(out,'cap-apex-regularity.png'),'Resolution',400); close(fig);
fprintf('R0_CAP_APEX_FIGURE_DONE\n');
end
