% ARTICLEB_BAND_FIGURE  Regenerate the de-aliasing SPECTRUM figures from saved
% band-optimization results (no solves). Left: metric spectrum over a range
% wider than the band, band shaded, f_s marked. Right: grouped bars of the band
% average, band maximum AND the wide-sweep maximum, so the RELOCATION (a peak
% pushed outside the band) is explicit. PNG + PDF, forced WHITE background.
%
%   articleB_band_figure
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
files={'articleB_band_32x16_z0.01_bh0.1.mat','articleB_band_32x16_z0.01_bh0.25.mat'};
S=pub_style();
for fi=1:numel(files)
  Sd=load(fullfile(RES,files{fi}));
  band=Sd.band; fsw=Sd.fsw(:); Jsw=Sd.Jsw; T=Sd.T; f1=Sd.f1; f_s=Sd.f_s;
  tag=regexprep(files{fi},'^articleB_band_',''); tag=regexprep(tag,'\.mat$','');
  fig=figure('Position',[60 60 1500 560],'Color','w');

  ax1=subplot(1,2,1); hold(ax1,'on');
  yl=[min(Jsw(:))*0.6, max(Jsw(:))*1.8];
  patch(ax1,[band.f_lo band.f_hi band.f_hi band.f_lo],[yl(1) yl(1) yl(2) yl(2)],...
        [0.85 0.90 0.97],'EdgeColor','none','FaceAlpha',0.9,'HandleVisibility','off');
  set(ax1,'YScale','log'); ylim(ax1,yl);
  h1=semilogy(ax1,fsw,Jsw(:,1),'-','Color',[0.35 0.35 0.35],'LineWidth',S.lw,'DisplayName','initial design');
  h2=semilogy(ax1,fsw,Jsw(:,2),'-','Color',S.c.b2,'LineWidth',S.lw,'DisplayName','single-freq obj (f_s)');
  h3=semilogy(ax1,fsw,Jsw(:,3),'-','Color',S.c.b1,'LineWidth',S.lw,'DisplayName','band obj');
  xline(ax1,f_s,'k--','f_s','LineWidth',S.lw_ref,'HandleVisibility','off');
  xline(ax1,band.f_lo,':','Color',[0.3 0.4 0.6],'LineWidth',S.lw_ref,'HandleVisibility','off');
  xline(ax1,band.f_hi,':','Color',[0.3 0.4 0.6],'LineWidth',S.lw_ref,'HandleVisibility','off');
  xlabel(ax1,'frequency f (Hz)','FontSize',S.fs_lab);
  ylabel(ax1,'dynamic compliance |F^T U|','FontSize',S.fs_lab);
  title(ax1,'metric spectrum (band shaded)','FontSize',S.fs_title);
  legend(ax1,[h1 h2 h3],'Location','northwest','FontSize',S.fs_leg);
  grid(ax1,'on'); box(ax1,'on');

  ax2=subplot(1,2,2);
  M=[T(1).band_avg T(2).band_avg T(3).band_avg; ...
     T(1).band_max T(2).band_max T(3).band_max; ...
     T(1).sweep_max T(2).sweep_max T(3).sweep_max].';
  b=bar(ax2,M,'grouped');
  b(1).FaceColor=[0.30 0.45 0.70]; b(2).FaceColor=[0.80 0.55 0.15]; b(3).FaceColor=[0.70 0.25 0.20];
  set(ax2,'YScale','log'); set(ax2,'XTickLabel',{'initial','single-freq','band'});
  ylabel(ax2,'objective (log scale)','FontSize',S.fs_lab);
  title(ax2,'band avg, band max, wide-sweep max','FontSize',S.fs_title);
  legend(ax2,{'band average','band maximum (in band)','wide-sweep maximum'},...
         'Location','northwest','FontSize',S.fs_leg);
  grid(ax2,'on'); box(ax2,'on'); yl2=get(ax2,'YLim');
  for di=1:3
    for bi=1:3
      xp=di+(bi-2)*0.22;
      text(ax2,xp,M(di,bi),sprintf('%.2g',M(di,bi)),'HorizontalAlignment','center',...
           'VerticalAlignment','bottom','FontSize',8,'Rotation',0);
    end
  end
  sgtitle(sprintf('Band de-aliasing: curved Scordelis-Lo roof, %s, f_1=%.4f Hz, band [%.3f,%.3f] Hz',...
          tag,f1,band.f_lo,band.f_hi),'FontSize',S.fs_title);
  base=fullfile(RES,sprintf('articleB_band_spectrum_%s',tag));
  pub_save(fig,base,300); close(fig);
  fprintf('figure: %s.png / .pdf\n',base);
end
disp('ARTICLEB_BAND_FIGURE_DONE');
