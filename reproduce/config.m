function c = config(addpaths)
%CONFIG  Paths and shared constants for the reproducibility package (MATLAB side).
%
%   c = config        returns the struct of paths and constants.
%   c = config(true)  also puts GeoPDEs, the NURBS toolbox, the author's GeoPDEs extension files
%                     (vendor/geopdes_topopt_functions), src/ and tests/ on the MATLAB path, in that order.
%
%   Environment variables (all optional):
%     REPRO_ROOT    folder holding results/, src/, tests/ (the repository root).
%                   Default: the parent of the folder that holds this file.
%     REPRO_OUT     where re-runs write. Default: <this folder>/out.
%     GEOPDES_ROOT  GeoPDEs distribution folder (the one that contains geopdes/ and nurbs-1.4.3/).
%                   Default: <REPRO_ROOT>/external/geopdes (see README.md).
%
%   Note. The original scripts in tests/ still contain absolute addpath calls (audit gap G04),
%   e.g. tests/obstacle_build_ctx.m:29-32. On another machine those folders do not exist and
%   addpath only warns; the paths set here then do the work. The originals were not edited.

if nargin < 1, addpaths = false; end
here = fileparts(mfilename('fullpath'));

c.PKG  = here;
c.ROOT = envor('REPRO_ROOT', fullfile(here, '..'));
c.OUT  = envor('REPRO_OUT',  fullfile(here, 'out'));
c.GEOPDES_ROOT = envor('GEOPDES_ROOT', fullfile(c.ROOT, 'external', 'geopdes'));

c.RESULTS = fullfile(c.ROOT, 'results');
c.RES  = fullfile(c.RESULTS, 'r0_revision');
c.SRC  = fullfile(c.ROOT, 'src');
c.TESTS = fullfile(c.ROOT, 'tests');

% Shared constants, copied from the code (not changed there):
c.Jref = 5.9464e-7;          % W; typed in tests/r0_pilot_runs.m:33, r0_long_runs.m:31, r0_fsi_reopt.m:16
c.rho0 = 1.21; c.c0 = 343;   % air, as in tests/r0_pilot_runs.m:34
% Thread counts used for the stored runs (audit G12): [] = MATLAB default
c.threads = struct('frozen_fields', [], 'design_point', [], 'regime_runs', 2, ...
                   'fsi_ctx', 2, 'fsi_reopt', 4);

if addpaths
    g = c.GEOPDES_ROOT;
    if ~isfolder(g)
        error('config:geopdes', 'GeoPDEs folder not found: %s (set GEOPDES_ROOT)', g);
    end
    addpath(genpath(g));
    addpath(fullfile(c.ROOT, 'vendor', 'geopdes_topopt_functions'));
    addpath(c.ROOT);              % the originals call addpath(fileparts(root))
    addpath(c.SRC);
    addpath(c.TESTS);
    if ~isfolder(c.OUT), mkdir(c.OUT); end
end
end

function v = envor(name, default)
v = getenv(name);
if isempty(v), v = default; end
v = char(java.io.File(v).getCanonicalPath());
end
