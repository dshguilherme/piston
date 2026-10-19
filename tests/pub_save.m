function pub_save(fig, base, res, mode)
% PUB_SAVE  Write a figure as 300 dpi PNG plus a PDF and an EPS sibling, same
% basename, with a FORCED WHITE background on the saved file (not just on
% screen).  Then read the PNG back and print corner + margin pixel RGB as
% white-check evidence, and report the EPS size / PostScript header.
%
%   pub_save(fig, base)              % vector PDF + vector EPS (default)
%   pub_save(fig, base, res, mode)   % mode: 'vector' (default) or 'raster'
%
% mode='raster' is for HEAVILY SHADED 3D SURFACES: MATLAB's painters renderer
% cannot keep them vector and tries (very slowly) to vectorise thousands of
% polygons.  In raster mode the PDF and EPS are written with ContentType
% 'image' (a high-resolution raster inside the PDF/EPS wrapper) which is fast
% and of sane size -- the honest choice for shaded 3D, per the figure brief.
if nargin < 3 || isempty(res), res = 300; end
if nargin < 4 || isempty(mode), mode = 'vector'; end
isRaster = strcmpi(mode,'raster');
set(fig, 'Color', 'w');
axs = findall(fig, 'type', 'axes');
for i = 1:numel(axs)
    set(axs(i), 'Color', 'w');
    try, set(axs(i), 'Toolbar', []); catch, end
end
legs = findall(fig, 'type', 'legend');
for i = 1:numel(legs), set(legs(i), 'Color', 'w'); end
if ~isRaster
    % Force the painters renderer so the EPS/PDF stay genuinely VECTOR (MATLAB
    % otherwise auto-selects OpenGL for object-heavy figures and rasterises).
    try, set(fig, 'Renderer', 'painters'); catch, end
end
drawnow;
% PDF + EPS FIRST, while the figure is certainly alive (the exportgraphics
% Exporter has an onCleanup race on this build that can invalidate the figure,
% so it runs last).  Vector by default; raster-embedded for shaded 3D.
if isRaster
    exportgraphics(fig, [base '.pdf'], 'ContentType', 'image', 'Resolution', res, ...
                   'BackgroundColor', 'white');
    exportgraphics(fig, [base '.eps'], 'ContentType', 'image', 'Resolution', res, ...
                   'BackgroundColor', 'white');
else
    % PDF vector via the stable printers path (InvertHardcopy whitens the page).
    pub_hardcopy(fig, [base '.pdf'], '-dpdf', '-vector');
    % EPS via the same printers path (Elsevier submission vector format).  NOTE:
    % MATLAB's painters renderer cannot keep a heavily shaded 3D surface vector
    % and embeds a raster inside the EPS wrapper; we still emit it and report
    % its size so the caller can decide whether it is usable.
    pub_hardcopy(fig, [base '.eps'], '-depsc', '-vector');
end
% PNG (on-screen raster, 300 dpi, forced white).  exportgraphics is fast but can
% race on this build, so fall back to the printers path.
okPng = false;
if isvalid(fig)
    for attempt = 1:5
        try
            exportgraphics(fig, [base '.png'], 'Resolution', res, 'BackgroundColor', 'white');
            okPng = true; break;
        catch
            drawnow; pause(0.5);
        end
        if ~isvalid(fig), break; end
    end
end
if ~okPng && isvalid(fig)
    pub_hardcopy(fig, [base '.png'], '-dpng', sprintf('-r%d', res));
end
if ~exist([base '.png'],'file')
    warning('pub_save:noPng','PNG export failed for %s', base); return;
end

% --- EPS sanity: non-zero, sane size, PostScript header present -------------
epsPath = [base '.eps']; epsBytes = 0; hasPS = false; epsHdr = '';
fid = fopen(epsPath, 'rb');
if fid > 0
    raw = fread(fid, 128, '*uint8').'; fclose(fid);
    epsHdr = char(raw);
    hasPS = strncmp(epsHdr, '%!PS', 4);
end
d = dir(epsPath); if ~isempty(d), epsBytes = d.bytes; end
fprintf('  [eps-check]   %s.eps  %d bytes  PS-header=%d  %s\n', base, epsBytes, ...
        hasPS, ternary(hasPS && epsBytes > 1000, 'OK', '*** CHECK ***'));

% --- white-background verification on the saved raster ----------------------
I = imread([base '.png']);
if size(I,3) == 1, I = repmat(I, 1, 1, 3); end
H = size(I,1); W = size(I,2);
c = [ squeeze(I(1,1,:)).'; squeeze(I(1,W,:)).'; ...
      squeeze(I(H,1,:)).'; squeeze(I(H,W,:)).' ];          % TL,TR,BL,BR corners
frame = [ I(1:3,:,:); I(end-2:end,:,:) ];
bmin = min([ min(double(frame),[],'all'), ...
             min(double(I(:,1:3,:)),[],'all'), ...
             min(double(I(:,end-2:end,:)),[],'all') ]);    % darkest border pixel
whitefrac = mean(double(I(:))==255);                       % fraction pure white
fprintf('  [white-check] %s\n', base);
fprintf('    corners (TL,TR,BL,BR) = %s\n', mat2str(c));
fprintf('    outer-border darkest   = %d  |  pixels==255 : %.1f%%\n', bmin, 100*whitefrac);
if any(c(:) < 250)
    fprintf('    *** WARNING: corner background not white ***\n');
end
end

function pub_hardcopy(fig, file, varargin)
% Robust print() wrapper: retry a few times, and make sure the figure still
% exists (guards against the intermittent exporter race).
for attempt = 1:4
    try
        if ~isvalid(fig), error('pub_save:invalidFig','figure gone'); end
        print(fig, file, varargin{:});
        return;
    catch err
        if attempt == 4, rethrow(err); end
        drawnow; pause(0.5);
    end
end
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
