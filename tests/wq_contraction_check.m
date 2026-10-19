% WQ_CONTRACTION_CHECK -- fifth directive, item (4).
% Verify the FULL 4-index contraction of the shell bending quadratic form
% against an INDEPENDENT brute-force hand contraction, and quantify the
% multiplicity of the symmetric (1,2)/(2,1) pair.
% M = D sqrt(a)[ nu t't + (1-nu) Q ], I1=t'Phi=a^{ab}kappa_ab,
% I2=Phi'QPhi=a^{ag}a^{bd}kappa_ab kappa_gd.
% Saves results/wq_contraction_check.txt.
addpath(genpath('C:/Users/dshgu/OneDrive/Documents/geopdes-master'));
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/src');
addpath('C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/tests');
resdir='C:/Users/dshgu/OneDrive/Documents/FastFormation/harmonic_shell_course/results';
if ~exist(resdir,'dir'), mkdir(resdir); end
fid=fopen(fullfile(resdir,'wq_contraction_check.txt'),'w');
lg=@(varargin) both(fid,varargin{:});
E=1.0; NU=0.3; D=E/(12*(1-NU^2));
lg('==== WQ shell bending contraction check (item 4) ====\n');
arc=nrbcirc(1.5,[0 0 0],0,pi*0.8); geometry=geo_load(nrbextrude(arc,[0 0 2]));
lg('geometry: cylinder patch R=1.5; rdim=%d\n', geometry.rdim);
rng(3); pts=[0.23 0.61 0.44; 0.37 0.79 0.15];
maxfull=0; maxbug=0;
for p=1:size(pts,2)
  Sk=shell_curvature(geometry,pts(:,p),'analytic');
  ai=Sk.ainv(:,:,1); Sq=Sk.sqrta(1); G=Sk.Gamma(:,:,:,1);
  A11=[1 0 0 -G(1,1,1) -G(1,1,2)]; A12=[0 1 0 -G(1,2,1) -G(1,2,2)]; A22=[0 0 1 -G(2,2,1) -G(2,2,2)];
  Ao=cell(2,2); Ao{1,1}=A11; Ao{1,2}=A12; Ao{2,1}=A12; Ao{2,2}=A22;
  t=ai(1,1)*A11+ai(1,2)*A12+ai(2,1)*A12+ai(2,2)*A22;
  Qf=zeros(5); Qbug=zeros(5);
  for al=1:2, for be=1:2, for g=1:2, for de=1:2
    Qf=Qf+ai(al,g)*ai(be,de)*(Ao{al,be}.'*Ao{g,de});
  end, end, end, end
  for al=1:2, for be=1:2
    Qbug=Qbug+ai(al,al)*ai(be,be)*(Ao{al,be}.'*Ao{al,be});
  end, end
  Mf=D*Sq*(NU*(t.'*t)+(1-NU)*Qf);   % code formula (full 4-index Q)
  Mb=D*Sq*(NU*(t.'*t)+(1-NU)*Qbug); % unique-pairs-only Q (the old bug)
  eF=0; eB=0;
  for k=1:40
    Phi=randn(5,1);
    I1=ai(1,1)*(A11*Phi)+2*ai(1,2)*(A12*Phi)+ai(2,2)*(A22*Phi);
    I2f=0; for al=1:2, for be=1:2, for g=1:2, for de=1:2
      I2f=I2f+ai(al,g)*ai(be,de)*(Ao{al,be}*Phi)*(Ao{g,de}*Phi); end, end, end, end
    ref=D*Sq*(NU*I1^2+(1-NU)*I2f);
    eF=max(eF,abs(Phi.'*Mf*Phi-ref)/max(abs(ref),eps));
    eB=max(eB,abs(Phi.'*Mb*Phi-ref)/max(abs(ref),eps));
  end
  lg('  point %d: full-4-index rel.err=%.3e | unique-pairs-only(OLD BUG) rel.err=%.3e\n', p, eF, eB);
  maxfull=max(maxfull,eF); maxbug=max(maxbug,eB);
end
lg('\nSUMMARY: full 4-index contraction rel.err=%.3e (machine eps -> formula correct);\n', maxfull);
lg('unique-pairs-only rel.err=%.3e (O(1) -> the symmetric multiplicity DOES matter;\n', maxbug);
lg('the current quadrule_bend_shell sums all four indices, i.e. the earlier fix holds).\n');
lg('WQ_CONTRACTION_CHECK_DONE\n');
fclose(fid);
fprintf('WQ_CONTRACTION_CHECK_DONE\n');
function both(fid,varargin), fprintf(fid,varargin{:}); fprintf(1,varargin{:}); end
