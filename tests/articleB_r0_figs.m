function articleB_r0_figs(tag)
% ARTICLEB_R0_FIGS  Figures for the phase-cancellation ladder, from the saved
% results/articleB_r0_metrics_<tag>.mat. Smooth (interpolated) design fields,
% 300 dpi PNG + raster PDF, forced WHITE background, subplot layout (robust on
% the shared machine). Safe to run on its own after articleB_r0_run.
%
%   articleB_r0_figs('256x128')
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
pub_style();
M=load(fullfile(RES,sprintf('articleB_r0_metrics_%s.mat',tag)),'L','nsub','freq');
L=M.L; nsub=M.nsub; freq=M.freq;
nJ=numel(L); ne=prod(nsub); ndof=numel(L(1).meas.u);
lab3={'dyncomp','inputpower','quadvel-rad'}; idx3=[1 2 3];
tmin=0.25; tmax=1.0; mu=L(1).meas0; cmap=parula(256);
fprintf('figs: mesh %s  N=%d ne=%d f=%.6g Hz\n',tag,ndof,ne,freq);

% F1: five designs, smooth, same colour scale
try, fig=figure('Position',[10 10 1600 300],'Color','w');
 for j=1:nJ
  subplot(1,nJ,j); smooth_field(reshape(L(j).tbest,nsub(1),nsub(2)));
  colormap(gca,cmap); caxis([tmin tmax]); colorbar;
  xlabel('i (\xi)'); ylabel('j (\eta)');
  title(sprintf('%s (%.2f dB)',L(j).label,10*log10(max(L(j).meta.Jbest,eps))),'Interpreter','none','FontWeight','normal');
 end
 sgtitle(sprintf('Optimised thickness fields, Scordelis-Lo roof, %s (%d elements, N=%d), f=%.4g Hz',tag,ne,ndof,freq),'Interpreter','none');
 save_fig(fig, fullfile(RES,sprintf('articleB_r0_designs_%s',tag))); close(fig);
catch err, warning('F1: %s',err.message); end

% F2: three primary designs, smooth, same colour scale
try, fig=figure('Position',[10 10 1100 320],'Color','w');
 for jj=1:3
  j=idx3(jj); subplot(1,3,jj); smooth_field(reshape(L(j).tbest,nsub(1),nsub(2)));
  colormap(gca,cmap); caxis([tmin tmax]); colorbar;
  xlabel('i (\xi)'); ylabel('j (\eta)'); title(L(j).label,'Interpreter','none','FontWeight','normal');
 end
 sgtitle(sprintf('Three objectives, same colour scale (%s)',tag),'Interpreter','none');
 save_fig(fig, fullfile(RES,sprintf('articleB_r0_compare_%s',tag))); close(fig);
catch err, warning('F2: %s',err.message); end

% F3: bar chart relative W per operator
try, fig=figure('Position',[10 10 900 460],'Color','w');
 Rr=zeros(3,3);
 for jj=1:3, j=idx3(jj);
   Rr(jj,:)=[L(j).meas.Wrad/mu.Wrad, L(j).meas.WI/mu.WI, L(j).meas.Wones/mu.Wones]; end
 hold on; grid on; box on;
 b=bar(1:3, Rr, 0.85,'grouped');
 cols=[0.122 0.403 0.678; 0.850 0.470 0.090; 0.176 0.615 0.271];
 for c=1:3, b(c).FaceColor=cols(c,:); end
 for c=1:3, for jj=1:3
   text(b(c).XEndPoints(jj),Rr(jj,c),sprintf('%.3f',Rr(jj,c)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',8);
 end, end
 set(gca,'XTick',1:3,'XTickLabel',lab3,'TickLabelInterpreter','none');
 ylabel('W / W_{uniform}  (lower = better)');
 legend({'W_{rad} (full R0, eq.3.66)','W with R0=I (phase-BLIND)','W with R0=rank1 (net vol. vel.)'},...
        'Location','northoutside','Orientation','horizontal');
 title(sprintf('Cancellation ladder on optimised designs (%s)',tag),'Interpreter','none');
 save_fig(fig, fullfile(RES,sprintf('articleB_r0_bars_%s',tag))); close(fig);
catch err, warning('F3: %s',err.message); end

% F4: phase / coherence field for the qv-rad design
try, Cu=L(3).meas.Cu; y=Cu;
 fig=figure('Position',[10 10 1300 380],'Color','w');
 subplot(1,3,1); smooth_field(angle(reshape(Cu,nsub(1),nsub(2))));
  colormap(gca,hsv(256)); caxis([-pi pi]); cb=colorbar; cb.Label.String='arg(V) [rad]';
  xlabel('i (\xi)'); ylabel('j (\eta)'); title('phase of centre velocity, qv-rad','Interpreter','none','FontWeight','normal');
 subplot(1,3,2); hold on; box on; grid on;
  yn=y/max(abs(y)); plot(real(yn),imag(yn),'.','Color',[0.55 0.55 0.55],'MarkerSize',10);
  s=sum(yn); plot([0 real(s)],[0 imag(s)],'-o','Color',[0.85 0.15 0.15],'LineWidth',2.4,'MarkerSize',6);
  plot(cos(linspace(0,2*pi,200)),sin(linspace(0,2*pi,200)),'k:');
  axis equal; xlabel('Re V (norm.)'); ylabel('Im V');
  title(sprintf('phasor sum |\\SigmaV|/\\Sigma|V| = %.3f',abs(s)/sum(abs(yn))),'Interpreter','tex','FontWeight','normal');
 subplot(1,3,3); hold on; box on; grid on;
  plot(cumsum(sort(abs(yn),'descend')),'LineWidth',1.8,'Color',[0.122 0.403 0.678]);
  yline(abs(s),'r--','LineWidth',1.4);
  xlabel('elements (sorted by |V|)'); ylabel('cumulative |V| / max');
  legend({'incoherent \Sigma|V|','coherent |\SigmaV|'},'Location','southeast');
  title('incoherent >> coherent','Interpreter','none','FontWeight','normal');
 sgtitle(sprintf('Phase/coherence field, quadvel-rad design (%s)',tag),'Interpreter','none');
 save_fig(fig, fullfile(RES,sprintf('articleB_r0_phase_%s',tag))); close(fig);
catch err, warning('F4: %s',err.message); end

fprintf('ARTICLEB_R0_FIGS_DONE\n');
end

function smooth_field(tg)
[n1,n2]=size(tg);
Xi=linspace(1,n1,6*n1); Yi=linspace(1,n2,6*n2);
Zi=interp2(1:n1,1:n2,tg.',Xi,Yi.','linear');
imagesc(Xi,Yi,Zi); axis image; set(gca,'YDir','normal');
end
function save_fig(fig, base)
set(fig,'Color','w');
ok=false;
for a=1:4
  try, print(fig,[base '.png'],'-dpng','-r300'); ok=true; break; catch, pause(1); end
end
if ~ok, warning('png failed %s',base); end
for a=1:4
  try, print(fig,[base '.pdf'],'-dpdf','-r200'); break; catch, pause(1); end
end
end
