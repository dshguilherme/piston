function r0_operator_figures()
% Publication-style plots for the radiation quadrature study.
root=fileparts(fileparts(mfilename('fullpath'))); out=fullfile(root,'results','r0_revision');
T=readtable(fullfile(out,'frozen_field_quadrature.csv'));
T.geometry=string(T.geometry); T.design=string(T.design);
T=T(T.surface_order==8 & T.angular_order==64,:);
T=T(ismember(T.design,["uniform","quadvel-rad","quadvel-I","quadvel-ones","legacy_long_horizon"]),:);
levels=["Centroid","Corrected centroid","Integrated patches","Spline trace"];
fig=figure('Color','w','Units','centimeters','Position',[2 2 17.5 10]);
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
colors=[0.05 0.35 0.62;0.80 0.28 0.12;0.20 0.55 0.38;0.48 0.34 0.65;0.15 0.62 0.68;0.70 0.52 0.12];
for panel=1:2
    if panel==1, geom="spherecap"; ax=nexttile; else, geom="cylpanel"; ax=nexttile; end
    rows=T.geometry==geom; labels=unique(T.design(rows),'stable'); hold(ax,'on');
    for j=1:numel(labels)
        r=rows & T.design==labels(j); z=T(r,:); z=sortrows(z,'design');
        vals=[z.legacy_centroid_W(1),z.joint_average_centroid_W(1), ...
              z.finite_patch_W(1),z.consistent_spline_W(1)];
        delta=100*(vals/vals(1)-1);
        sty='-'; if j==2,sty='--';elseif j==3,sty=':';elseif j==4,sty='-.';end
        pretty=labels(j);
        switch labels(j)
            case "uniform", pretty="Uniform field";
            case "quadvel-rad", pretty="Legacy R0-optimized field";
            case "quadvel-I", pretty="Legacy coefficient-norm field";
            case "quadvel-ones", pretty="Legacy coherent-sum field";
            case "legacy_long_horizon", pretty="Panel long-horizon field";
        end
        plot(ax,1:4,delta,'Color',colors(j,:),'LineStyle',sty,'LineWidth',1.9, ...
            'Marker','o','MarkerSize',5,'MarkerFaceColor','w', ...
            'DisplayName',pretty);
    end
    yline(ax,0,'k-','LineWidth',0.75,'HandleVisibility','off');
    xlim(ax,[0.85 4.15]); xticks(ax,1:4); xticklabels(ax,levels);
    xtickangle(ax,20); ylabel(ax,'Change in predicted power from centroid (%)');
    if panel==1, title(ax,'Spherical cap'); else, title(ax,'Cylindrical panel'); end
    grid(ax,'on'); ax.GridAlpha=0.15; ax.Box='on'; ax.LineWidth=0.9;
    ax.FontName='Arial'; ax.FontSize=8.5; ax.XColor=[.12 .12 .12]; ax.YColor=[.12 .12 .12];
    ax.Color='w'; ax.Toolbar.Visible='off';
    ax.Title.Color=[.08 .08 .08]; ax.Title.FontWeight='normal'; ax.Title.FontSize=9.5;
    legend(ax,'Location','best','FontSize',7.3,'Box','off','TextColor',[.1 .1 .1],'Color','w');
end
sgtitle(fig,'Effect of radiation-matrix surface representation on fixed fields', ...
    'FontWeight','normal','FontSize',11,'Color',[.08 .08 .08]);
exportgraphics(fig,fullfile(out,'operator-fidelity.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(out,'operator-fidelity.png'),'Resolution',400);
% Conceptual diagram: separate numerical representation from physical model.
f=figure('Color','w','Units','centimeters','Position',[2 2 17.5 7.5]); ax=axes(f,'Position',[.04 .08 .92 .82]); ax.Toolbar.Visible='off'; ax.Color='w'; axis(ax,[0 1 0 1]); axis(ax,'off'); hold(ax,'on');
boxes=[.03 .43 .20 .34;.29 .43 .20 .34;.55 .43 .20 .34;.81 .43 .17 .34];
labels={['Element-centroid' char(10) 'normal trace' char(10) 'area-weighted sinc'], ...
 ['Finite-patch' char(10) 'normal-average trace' char(10) 'double area integration'], ...
 ['Spline-consistent' char(10) 'normal trace' char(10) 'B^H W K W B'], ...
 ['Exterior acoustic' char(10) 'boundary solve' char(10) 'BEM / FEM-PML']};
for j=1:4
    rectangle(ax,'Position',boxes(j,:),'Curvature',.05,'FaceColor',[.94 .96 .98], ...
        'EdgeColor',[.18 .34 .48],'LineWidth',1.2);
    text(ax,boxes(j,1)+boxes(j,3)/2,boxes(j,2)+boxes(j,4)/2,labels{j}, ...
        'HorizontalAlignment','center','VerticalAlignment','middle','FontName','Arial', ...
        'FontSize',9,'Color',[.08 .12 .16],'Interpreter','none');
end
for j=1:3
    x1=boxes(j,1)+boxes(j,3); x2=boxes(j+1,1);
    annotation(f,'arrow',[.04+.92*(x1+.005) .04+.92*(x2-.008)],[.08+.82*.60 .08+.82*.60],'Color',[.25 .25 .25],'LineWidth',1.1);
end
text(ax,.5,.91,'Two questions must not be confused','HorizontalAlignment','center', ...
    'FontName','Arial','FontSize',12,'FontWeight','bold','Color',[.08 .08 .08]);
text(ax,.39,.27,'Boxes 1–3: same prescribed-source sinc model; numerical refinement', ...
    'HorizontalAlignment','center','FontName','Arial','FontSize',8.5,'Color',[.08 .24 .40]);
text(ax,.50,.12,'Box 4: solve the exterior acoustic problem; the physical model changes', ...
    'HorizontalAlignment','center','FontName','Arial','FontSize',8.5,'Color',[.55 .20 .12]);
plot(ax,[.03 .75],[.22 .22],'-','Color',[.18 .34 .48],'LineWidth',1.1);
plot(ax,[.81 .98],[.07 .07],'-','Color',[.68 .25 .14],'LineWidth',1.1);
exportgraphics(f,fullfile(out,'operator-hierarchy.pdf'),'ContentType','vector');
exportgraphics(f,fullfile(out,'operator-hierarchy.png'),'Resolution',400);
close all; fprintf('R0_OPERATOR_FIGURES_DONE\n');
end
