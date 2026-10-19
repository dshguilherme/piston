function S = pub_style()
% PUB_STYLE  Unified publication style for Article A / Article B figures.
% Applies consistent fonts, line widths, tick direction, and box for the whole
% figure set, and returns a struct of named colours/line widths so every panel
% in the set matches. Call once at the top of each figure script.

S.font      = 'Helvetica';
S.fs        = 12;   % base axes font size
S.fs_lab    = 13;   % axis label
S.fs_title  = 13;   % panel / figure title
S.fs_leg    = 11;   % legend
S.lw        = 1.8;  % line width for data lines
S.lw_ref    = 1.1;  % reference / guide lines
S.ms        = 8;    % marker size
S.lw_ax     = 0.9;  % axis line width

% --- scheme colours (Article A) : colour-blind safe, one per quadrature rule -
S.c.WQ   = [0.835 0.153 0.153];   % reduced weighted quadrature (red)
S.c.G    = [0.122 0.403 0.678];   % Gauss-with-metric          (blue)
S.c.C2   = [0.176 0.615 0.271];   % C-weighted, nr = 2         (green)
S.c.C3   = [0.549 0.235 0.706];   % C-weighted, nr = 3         (purple)

% --- problem colours (Article A cost scatter) -------------------------------
S.c.p1   = [0.000 0.310 0.780];   % pinched cylinder
S.c.p2   = [0.850 0.470 0.090];   % Scordelis-Lo roof
S.c.p3   = [0.180 0.560 0.300];   % pinched hemisphere

% --- Article B metric/panel accent colours ----------------------------------
S.c.b1   = [0.122 0.403 0.678];
S.c.b2   = [0.850 0.470 0.090];
S.c.b3   = [0.176 0.615 0.271];

% --- global defaults (WHITE backgrounds, dark-on-white text/axes) -----------
set(groot,'defaultFigureColor','w');
set(groot,'defaultAxesColor','w');
set(groot,'defaultLegendColor','w');
set(groot,'defaultColorbarColor','w');
set(groot,'defaultAxesFontName',S.font);
set(groot,'defaultTextFontName',S.font);
set(groot,'defaultLegendFontName',S.font);
set(groot,'defaultColorbarFontName',S.font);
set(groot,'defaultAxesFontSize',S.fs);
set(groot,'defaultAxesLineWidth',S.lw_ax);
set(groot,'defaultLineLineWidth',S.lw);
set(groot,'defaultAxesTickDir','out');
set(groot,'defaultAxesTickLength',[0.018 0.018]);
set(groot,'defaultAxesBox','on');
set(groot,'defaultAxesXColor',[0 0 0]);
set(groot,'defaultAxesYColor',[0 0 0]);
set(groot,'defaultAxesZColor',[0 0 0]);
set(groot,'defaultTextColor',[0 0 0]);
set(groot,'defaultAxesGridColor',[0.75 0.75 0.75]);
set(groot,'defaultAxesGridAlpha',0.5);
set(groot,'defaultAxesXGrid','off');
set(groot,'defaultAxesYGrid','off');
set(groot,'defaultAxesLayer','top');
end
