function articleB_r0_wavecap_figs(tag)
% ARTICLEB_R0_WAVECAP_FIGS  Figures for the WAVE-regime R0 ladder on the clamped
% spherical cap, from results/articleB_r0_wavecap_metrics_<tag>.mat and
% ..._freqscan_<tag>.txt.  White background forced and PIXEL-VERIFIED (pub_save),
% 300 dpi PNG + vector PDF.
%
%   articleB_r0_wavecap_figs('128x64_f100')
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
RES='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
pub_style();

M=load(fullfile(RES,sprintf('articleB_r0_wavecap_metrics_%s.mat',tag)),...
  'L','nsub','freq','Lmax','kLmax','relK1','relKfro','band','nband');
L=M.L; nsub=M.nsub; freq=M.freq; kLmax=M.kLmax; nJ=numel(L);
mu=L(1).meas0; tmin=0.25; tmax=1.0; cmap=parula(256);
lab3={'dyncomp','inputpower','quadvel-rad'};
fprintf('wavecap figs: %s mesh %s f=%.5g Hz kLmax=%.3f nband=%d\n',tag,mat2str(nsub),freq,kLmax,M.nband);

% F1: five designs, smooth, same colour scale
fig=figure('Position',[10 10 1700 320],'Color','w');
for j=1:nJ
  subplot(1,nJ,j); smooth_field(reshape(L(j).tbest,nsub(1),nsub(2)));
  colormap(gca,cmap); caxis([tmin tmax]); cb=colorbar; cb.Label.String='t (m)';
  xlabel('i (\xi)'); ylabel('j (\eta)');
  title(sprintf('%s (%.2f dB)',L(j).label,10*log10(max(L(j).meta.Jbest,eps))),'Interpreter','none','FontWeight','normal');
end
sgtitle(sprintf('Clamped spherical cap: wave-regime designs, %s, f=%.4g Hz, kL_{max}=%.2f, %d modes in [80,130] Hz',tag,freq,kLmax,M.nband),'Interpreter','none');
pub_save(fig, fullfile(RES,sprintf('articleB_r0_wavecap_designs_%s',tag))); close(fig);

% F2: three primary designs side by side, same colour scale
fig=figure('Position',[10 10 1100 340],'Color','w');
for jj=1:3
  subplot(1,3,jj); smooth_field(reshape(L(jj).tbest,nsub(1),nsub(2)));
  colormap(gca,cmap); caxis([tmin tmax]); cb=colorbar; cb.Label.String='t (m)';
  xlabel('i (\xi)'); ylabel('j (\eta)'); title(lab3{jj},'Interpreter','none','FontWeight','normal');
end
sgtitle(sprintf('Three objectives, same colour scale (%s)',tag),'Interpreter','none');
pub_save(fig, fullfile(RES,sprintf('articleB_r0_wavecap_compare_%s',tag))); close(fig);

% F3: operator-comparison bars
fig=figure('Position',[10 10 950 480],'Color','w');
Rr=zeros(3,3);
for jj=1:3, Rr(jj,:)=[L(jj).meas.Wrad/mu.Wrad, L(jj).meas.WI/mu.WI, L(jj).meas.Wones/mu.Wones]; end
hold on; grid on; box on;
b=bar(1:3, Rr, 0.85,'grouped');
cols=[0.122 0.403 0.678; 0.850 0.470 0.090; 0.176 0.615 0.271];
for c=1:3, b(c).FaceColor=cols(c,:); end
for c=1:3, for jj=1:3
  text(b(c).XEndPoints(jj),Rr(jj,c),sprintf('%.3f',Rr(jj,c)),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',9);
end, end
set(gca,'XTick',1:3,'XTickLabel',lab3,'TickLabelInterpreter','none');
ylabel('W / W_{uniform}  (lower = better)');
legend({'W_{rad} (full R0, eq.3.66)','W with R0=I (phase-BLIND)','W with R0=rank-1 (net vol. vel.)'},...
       'Location','northoutside','Orientation','horizontal');
title(sprintf('Cancellation ladder, wave regime, cap (%s)',tag),'Interpreter','none');
pub_save(fig, fullfile(RES,sprintf('articleB_r0_wavecap_bars_%s',tag))); close(fig);

% F4: phase / coherence field for the qv-rad design (REAL phase structure)
jrad=3; Cu=L(jrad).meas.Cu; y=Cu;
fig=figure('Position',[10 10 1350 400],'Color','w');
subplot(1,3,1); smooth_field(angle(reshape(Cu,nsub(1),nsub(2))));
 colormap(gca,hsv(256)); caxis([-pi pi]); cb=colorbar; cb.Label.String='arg(V) [rad]';
 xlabel('i (\xi)'); ylabel('j (\eta)'); title('phase of centre velocity, qv-rad','Interpreter','none','FontWeight','normal');
subplot(1,3,2); hold on; box on; grid on;
 yn=y/max(abs(y)); plot(real(yn),imag(yn),'.','Color',[0.55 0.55 0.55],'MarkerSize',11);
 s=sum(yn); plot([0 real(s)],[0 imag(s)],'-o','Color',[0.85 0.15 0.15],'LineWidth',2.4,'MarkerSize',6);
 plot(cos(linspace(0,2*pi,200)),sin(linspace(0,2*pi,200)),'k:');
 axis equal; xlabel('Re V (norm.)'); ylabel('Im V');
 title(sprintf('phasor sum |\\SigmaV|/\\Sigma|V| = %.3f',abs(s)/sum(abs(yn))),'Interpreter','tex','FontWeight','normal');
subplot(1,3,3); hold on; box on; grid on;
 plot(cumsum(sort(abs(yn),'descend')),'LineWidth',1.8,'Color',[0.122 0.403 0.678]);
 yline(abs(s),'r--','LineWidth',1.4);
 xlabel('elements (sorted by |V|)'); ylabel('cumulative |V| / max');
 legend({'incoherent \Sigma|V|','coherent |\SigmaV|'},'Location','southeast');
 title('incoherent vs coherent sum','Interpreter','none','FontWeight','normal');
sgtitle(sprintf('Phase/coherence field, quadvel-rad design, cap (%s, f=%.4g Hz)',tag,freq),'Interpreter','none');
pub_save(fig, fullfile(RES,sprintf('articleB_r0_wavecap_phase_%s',tag))); close(fig);

% F5: W_rad vs frequency (uniform) with the kernel rank-1 indicator
fs=readfreqscan(fullfile(RES,sprintf('articleB_r0_wavecap_freqscan_%s.txt',tag)));
fig=figure('Position',[10 10 900 480],'Color','w');
yyaxis left; semilogx(fs.f, fs.Wrad/fs.Wrad(end),'-o','Color',[0.122 0.403 0.678],'LineWidth',1.8,'MarkerSize',6);
 ylabel('W_{rad} / W_{rad}(uniform, 4f_{exc})');
yyaxis right; semilogx(fs.f, fs.rank1,'-s','Color',[0.850 0.470 0.090],'LineWidth',1.8,'MarkerSize',6);
 ylabel('rank-1 indicator  s_2/s_1(K)');
xline(freq,'k--','LineWidth',1.2);
grid on; box on; xlabel('assembly frequency f (Hz)');
title(sprintf('W_{rad} and kernel rank-1 departure vs frequency, cap (%s)',tag),'Interpreter','none');
legend({'W_{rad} (normalised)','s_2/s_1(K)','operating f'},'Location','best');
pub_save(fig, fullfile(RES,sprintf('articleB_r0_wavecap_freqcurve_%s',tag))); close(fig);

fprintf('ARTICLEB_R0_WAVECAP_FIGS_DONE\n');
end

function S=readfreqscan(p)
fid=fopen(p,'r'); C=textscan(fid,'%f %f %f %f %f %f %f','CommentStyle','#','CollectOutput',true); fclose(fid);
A=C{1}; S.f=A(:,1); S.kl=A(:,2); S.rank1=A(:,3); S.Wrad=A(:,5); S.WI=A(:,6); S.odfrac=A(:,7);
end
function smooth_field(tg)
[n1,n2]=size(tg);
Xi=linspace(1,n1,6*n1); Yi=linspace(1,n2,6*n2);
Zi=interp2(1:n1,1:n2,tg.',Xi,Yi.','linear');
imagesc(Xi,Yi,Zi); axis image; set(gca,'YDir','normal');
end
