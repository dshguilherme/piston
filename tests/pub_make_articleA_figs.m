function pub_make_articleA_figs
% PUB_MAKE_ARTICLEA_FIGS  Publication figures for Article A (weighted
% quadrature on the shell obstacle course).  Reads NO results file to invent
% numbers: every value below is transcribed from
%   results/articleA_quadrature_obstacle.txt.
% Writes results/pub_fig_A1_convergence.{png,pdf}
%        results/pub_fig_A2_ladder.{png,pdf}
%        results/pub_fig_A3_cost_accuracy.{png,pdf}
%        results/pub_fig_A4_geometries.{png,pdf}

addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
RES = 'C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
S = pub_style();
out = @(b) fullfile(RES, b);

% ============================================================ A1 convergence
% transcribed from articleA_quadrature_obstacle.txt section B/D.
names = {'Pinched cylinder','Scordelis--Lo roof','Pinched hemisphere (quarter)'};
nsub  = {[8 16],[8 16 24],[8 16]};
eWQ   = {[2.299e+00 4.796e+00], [6.346e-02 2.523e-02 1.144e-02], [4.629e-02 1.680e-02]};
eG    = {[7.597e-13 5.647e-14], [2.511e-13 1.232e-12 7.692e-12], [2.121e-09 2.479e-12]};
eC2   = {[9.707e-03 9.699e-03], [2.486e-04 7.406e-05 3.493e-05], [1.660e-04 4.979e-05]};
eC3   = {[1.977e-09 1.330e-10], [2.995e-09 3.083e-11 6.356e-12], [7.175e-08 2.647e-09]};

fig = figure('Units','inches','Position',[1 1 7.6 3.0]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');
for k = 1:3
    nexttile; hold on; grid on; box on;
    h = 1 ./ nsub{k};
    % reference order line O(h), anchored at the cmoment-3 last point
    hx = linspace(min(h)*0.9, max(h)*1.25, 20);
    yref = eC3{k}(end) * (hx/max(h)).^1;
    hRef = plot(hx, yref, '--', 'Color',[0.45 0.45 0.45], 'LineWidth', S.lw_ref, ...
                'DisplayName','$O(h)$ reference');
    % machine-precision floor of the analytic reference
    yline(1e-13,'-.','Color',[0.6 0.6 0.6],'LineWidth',0.9, ...
          'Label','machine precision','LabelHorizontalAlignment','left', ...
          'FontSize',9,'Interpreter','tex');
    p1 = loglog(h, eWQ{k}, '-o', 'Color',S.c.WQ, 'LineWidth',S.lw, 'MarkerSize',S.ms, ...
                'MarkerFaceColor',S.c.WQ,'MarkerEdgeColor','w', 'DisplayName','reduced WQ');
    p2 = loglog(h, eG{k},  '-s', 'Color',S.c.G,  'LineWidth',S.lw, 'MarkerSize',S.ms, ...
                'MarkerFaceColor',S.c.G,'MarkerEdgeColor','w', 'DisplayName','Gauss with metric');
    p3 = loglog(h, eC2{k}, '-^', 'Color',S.c.C2, 'LineWidth',S.lw, 'MarkerSize',S.ms, ...
                'MarkerFaceColor',S.c.C2,'MarkerEdgeColor','w','DisplayName','C-weighted $\mathrm{nr}=2$');
    p4 = loglog(h, eC3{k}, '-d', 'Color',S.c.C3, 'LineWidth',S.lw, 'MarkerSize',S.ms, ...
                'MarkerFaceColor',S.c.C3,'MarkerEdgeColor','w','DisplayName','C-weighted $\mathrm{nr}=3$');
    set(gca,'XScale','log','YScale','log','FontSize',S.fs);
    xlim([min(h)*0.85 max(h)*1.35]); ylim([1e-14 2e1]);
    xlabel('element size $h$','Interpreter','latex','FontSize',S.fs_lab);
    if k==1, ylabel('relative bending-energy error','Interpreter','latex','FontSize',S.fs_lab); end
    title(names{k},'Interpreter','latex','FontSize',S.fs_title,'FontWeight','normal');
    if k==1
        legend([p1 p2 p3 p4],'Interpreter','latex','FontSize',S.fs_leg, ...
               'Location','southwest','NumColumns',1);
    end
end
pub_save(fig, out('pub_fig_A1_convergence')); close(fig);
fprintf('A1 written\n');

% ============================================================ A2 enrichment ladder
nr_all = [1 2 3 4 5];
elad   = [2.090e-02 1.090e-02 6.109e-03 3.399e-10 4.958e-13];   % txt section D
fig = figure('Units','inches','Position',[1 1 4.6 3.2]);
hold on; grid on; box on;
xA = linspace(1,5,200);
semilogy(xA, elad(4)*(xA/4).^-6, '--', 'Color',[0.55 0.55 0.55], 'LineWidth',S.lw_ref, ...
         'DisplayName','high-order reference');
yline(1e-12,'-.','Color',[0.6 0.6 0.6],'LineWidth',0.9,'Label','reference floor', ...
      'FontSize',9,'Interpreter','tex');
plot(nr_all, elad, '-o', 'Color',S.c.C3, 'LineWidth',S.lw, 'MarkerSize',S.ms+1, ...
     'MarkerFaceColor',S.c.C3,'MarkerEdgeColor','w','DisplayName','C-weighted WQ');
plot(4, elad(4), 'p', 'MarkerSize',14, 'MarkerFaceColor',[1 0.85 0.1], ...
     'MarkerEdgeColor',[0.4 0.3 0], 'DisplayName','floor reached ($\mathrm{nr}\geq4$)');
set(gca,'YScale','log','FontSize',S.fs);
xlim([0.6 5.4]); ylim([1e-13 1e-1]); xticks(1:5);
xlabel('enriched moments $\mathrm{nr}$','Interpreter','latex','FontSize',S.fs_lab);
ylabel('relative bending-energy error','Interpreter','latex','FontSize',S.fs_lab);
title('Moment-enrichment ladder (pinched cylinder, $n_{\rm sub}=8$)','Interpreter','latex', ...
      'FontSize',S.fs_title,'FontWeight','normal');
legend('Location','southwest','Interpreter','latex','FontSize',S.fs_leg);
pub_save(fig, out('pub_fig_A2_ladder')); close(fig);
fprintf('A2 written\n');

% ============================================================ A3 cost / accuracy
% transcribed from txt section C.  method row: [prob, nsub, scheme, time, err, work]
prob = {'pinched cylinder','Scordelis-Lo roof','pinched hemisphere'};
M = {
 'pinched cylinder', 8,'WQ',3.170,2.299e+00,1631;
 'pinched cylinder', 8,'G',8.790,7.597e-13,10000;
 'pinched cylinder', 8,'C2',11.246,9.707e-03,1631;
 'pinched cylinder', 8,'C3',10.105,1.977e-09,1631;
 'pinched cylinder',16,'WQ',9.662,4.796e+00,1421;
 'pinched cylinder',16,'G',55.769,5.647e-14,10000;
 'pinched cylinder',16,'C2',39.914,9.699e-03,1421;
 'pinched cylinder',16,'C3',37.929,1.330e-10,1421;
 'Scordelis-Lo roof',8,'WQ',0.888,6.346e-02,1650;
 'Scordelis-Lo roof',8,'G',2.034,2.511e-13,10000;
 'Scordelis-Lo roof',8,'C2',4.407,2.486e-04,1650;
 'Scordelis-Lo roof',8,'C3',5.198,2.995e-09,1650;
 'Scordelis-Lo roof',16,'WQ',3.556,2.523e-02,1430;
 'Scordelis-Lo roof',16,'G',12.477,1.232e-12,10000;
 'Scordelis-Lo roof',16,'C2',12.159,7.406e-05,1430;
 'Scordelis-Lo roof',16,'C3',10.432,3.083e-11,1430;
 'Scordelis-Lo roof',24,'WQ',4.504,1.144e-02,1360;
 'Scordelis-Lo roof',24,'G',25.868,7.692e-12,10000;
 'Scordelis-Lo roof',24,'C2',18.136,3.493e-05,1360;
 'Scordelis-Lo roof',24,'C3',22.443,6.356e-12,1360;
 'pinched hemisphere',8,'WQ',0.423,4.629e-02,1650;
 'pinched hemisphere',8,'G',0.733,2.121e-09,10000;
 'pinched hemisphere',8,'C2',2.114,1.660e-04,1650;
 'pinched hemisphere',8,'C3',1.734,7.175e-08,1650;
 'pinched hemisphere',16,'WQ',2.005,1.680e-02,1430;
 'pinched hemisphere',16,'G',6.591,2.479e-12,10000;
 'pinched hemisphere',16,'C2',9.219,4.979e-05,1430;
 'pinched hemisphere',16,'C3',8.567,2.647e-09,1430;
};
pcol = containers.Map({'pinched cylinder','Scordelis-Lo roof','pinched hemisphere'}, ...
                      {S.c.p1,S.c.p2,S.c.p3});
mk = containers.Map({'WQ','G','C2','C3'}, {'o','s','^','d'});
methname = containers.Map({'WQ','G','C2','C3'}, ...
   {'reduced WQ','Gauss','C-weighted $\mathrm{nr}=2$','C-weighted $\mathrm{nr}=3$'});

fig = figure('Units','inches','Position',[1 1 7.4 3.1]);
tl = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
xcol = [6 4];    % column 6 = work/el, 4 = time_s
xlab = {'stored work per element (WQ page slots / Gauss entries)','assembly wall time (s)'};
for a = 1:2
    nexttile; hold on; grid on; box on;
    hh = gobjects(0);
    for i = 1:size(M,1)
        xn = M{i,xcol(a)}; ye = M{i,5};
        h = plot(xn, ye, mk(M{i,3}), 'Color', pcol(M{i,1}), 'MarkerSize', S.ms+1, ...
                 'MarkerFaceColor', pcol(M{i,1}), 'MarkerEdgeColor','w', 'LineWidth',1.0);
        if a==1 && strcmp(M{i,3},'WQ'), hh(end+1) = h; end %#ok<AGROW>
    end
    set(gca,'XScale','log','YScale','log','FontSize',S.fs);
    xlabel(xlab{a},'Interpreter','latex','FontSize',S.fs_lab);
    if a==1, ylabel('achieved relative error','Interpreter','latex','FontSize',S.fs_lab); end
    ylim([1e-14 1e1]);
end
% shared legend: colour = problem, marker = scheme (built from proxy handles)
nexttile(1); %#ok<NASGU>
ax0 = nexttile(1);
leg1 = gobjects(3,1); leg2 = gobjects(4,1);
for i=1:3, leg1(i) = plot(ax0, NaN,NaN,'s','Color',pcol(prob{i}),'MarkerFaceColor',pcol(prob{i}),'MarkerEdgeColor','w','MarkerSize',S.ms); end
schemes = {'WQ','G','C2','C3'};
for i=1:4, leg2(i) = plot(ax0, NaN,NaN,mk(schemes{i}),'Color',0.3*[1 1 1],'MarkerFaceColor',0.4*[1 1 1],'MarkerEdgeColor','w','MarkerSize',S.ms+1); end
lgd = legend([leg1;leg2], ...
   {'pinched cylinder','Scordelis--Lo roof','pinched hemisphere', ...
    'reduced WQ','Gauss','C-weighted $\mathrm{nr}=2$','C-weighted $\mathrm{nr}=3$'}, ...
   'Interpreter','latex','FontSize',S.fs_leg-1,'Location','southoutside','NumColumns',4, ...
   'Orientation','horizontal');
lgd.Layout.Tile = 'south';
ttlA = title(tl,'Cost vs achieved accuracy of the bending-energy scheme','Interpreter','latex', ...
      'FontSize',S.fs_title,'FontWeight','normal'); set(ttlA,'Color',[0 0 0]);
pub_save(fig, out('pub_fig_A3_cost_accuracy')); close(fig);
fprintf('A3 written\n');

% ============================================================ A4 geometries
pg = obstacle_course_geom('cylinder');
Rr=10; hole=18*pi/180;
geo_hem = geo_load(nrbrevolve(nrbcirc(Rr,[0 0 0],hole,pi/2),[0 0 0],[1 0 0],pi/2));
G = { pg(1).nurbs, ...                                             % pinched cylinder
      obstacle_course_geom('roof').nurbs, ...                      % Scordelis-Lo roof
      geo_hem.nurbs };                                             % hemisphere quarter
ttl = {'Pinched cylinder','Scordelis--Lo roof','Pinched hemisphere (quarter)'};
tgt = { [300;0;300], ...                                           % Pc = [R;0;L/2]
        [25*cos(pi/2-2*pi/9); 25; 25*sin(pi/2-2*pi/9)], ...        % free-edge midpoint
        [0;Rr;0] };                                               % equator load point
axlab = {'$x$, $y$, $z$ (mm)','$x$, $y$, $z$ (units)','$x$, $y$, $z$ (units)'};

fig = figure('Units','inches','Position',[1 1 8.0 3.0]);
tl = tiledlayout(1,3,'TileSpacing','compact','Padding','loose');
for k = 1:3
    nexttile; hold on; box on;
    nr = G{k};
    u = linspace(nr.knots{1}(1), nr.knots{1}(end), 80);
    v = linspace(nr.knots{2}(1), nr.knots{2}(end), 80);
    X=zeros(numel(v),numel(u)); Y=X; Z=X;
    for j=1:numel(v)
        for i=1:numel(u)
            p = nrbeval(nr, {u(i), v(j)});
            X(j,i)=p(1); Y(j,i)=p(2); Z(j,i)=p(3);
        end
    end
    [Nx,Ny,Nz] = surfnorm(X,Y,Z);
    Ldir = [-0.4 -0.6 0.7]; Ldir = Ldir/norm(Ldir);
    ill = Nx*Ldir(1)+Ny*Ldir(2)+Nz*Ldir(3);
    ill = max(0, min(1, 0.45 + 0.55*abs(ill)));
    cmapshade = [linspace(0.55,0.95,256)', linspace(0.60,0.96,256)', linspace(0.70,0.98,256)'];
    surf(X,Y,Z,'CData',ill,'FaceColor','flat','EdgeColor','none');
    colormap(gca, cmapshade); caxis(gca,[0 1]);
    set(gca,'Color','w');
    tp = tgt{k};
    plot3(tp(1),tp(2),tp(3),'p','MarkerSize',15,'MarkerFaceColor',S.c.WQ, ...
          'MarkerEdgeColor','k','LineWidth',0.8);
    text(tp(1),tp(2),tp(3),'  target','FontSize',9,'Color','k','Interpreter','none');
    axis equal; view(-35,22); grid off;
    xlabel(axlab{k},'Interpreter','latex','FontSize',S.fs_lab-1);
    title(ttl{k},'Interpreter','latex','FontSize',S.fs_title,'FontWeight','normal','Color',[0 0 0]);
    set(gca,'FontSize',S.fs-1);
    zl = get(gca,'ZLim'); set(gca,'ZTick',[zl(1) zl(2)]);   % avoid crowded z ticks
end
pub_save(fig, out('pub_fig_A4_geometries')); close(fig);
fprintf('A4 written\n');

fprintf('ARTICLEA_PUB_FIGS_DONE\n');
end
