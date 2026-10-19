# REPRODUCE-PACKAGE COPY. Do not edit numerics here; edit the original and re-copy.
# Original: paper_jsv/reproduce_plan/rescued_scratch/queue.ps1  (sha256 7bb56e9fb5e4d518)
# Changed: original line 5 ($root = absolute path) -> $env:REPRO_ROOT if set, else (Resolve-Path "$PSScriptRoot\..\..").Path;
#          original line 13: MATLAB executable "C:\Program Files\MATLAB\R2025a\bin\matlab.exe" -> $env:MATLAB_EXE if set, else 'matlab' (from PATH).
#          Kept: maxNumCompThreads(2), the skip rule (a finished run_<cell>_<code>.mat is skipped unless code is 'eval'),
#          the per-run MATLAB -logfile, and the queue.log semantics (start/end/done lines appended to <out>\queue.log).
# WARNING: tests/r0_pilot_runs.m (the ORIGINAL, unedited) writes its results into results\r0_revision\pilot\ under $root, and this script's
#          $out, queue.log and log_run_*.log are in the same folder. Running this therefore OVERWRITES or adds to the stored pilot runs
#          (the skip rule only protects run_<cell>_<code>.mat files that already exist). Run it only on a scratch copy of the tree
#          (set REPRO_ROOT to that copy), never on the archived results.
# NOTE: the -List argument ("cell:code:nmax,..." items) used for the paper is not recorded in this script; see drivers\pilot_queue.txt.
param([string]$List)
# sequential pilot queue: one MATLAB session at a time; skips runs already saved.
# items: cell:code:nmax  (code 'eval' runs the evaluation stage). MATLAB writes its own -logfile,
# so it does not depend on this script's stdout.
if ($env:REPRO_ROOT) { $root = $env:REPRO_ROOT } else { $root = (Resolve-Path "$PSScriptRoot\..\..").Path }
if ($env:MATLAB_EXE) { $matlabExe = $env:MATLAB_EXE } else { $matlabExe = 'matlab' }
$out = "$root\results\r0_revision\pilot"
foreach ($item in $List.Split(',')) {
    $cell, $code, $nmax = $item.Split(':')
    if ($code -ne 'eval' -and (Test-Path "$out\run_${cell}_${code}.mat")) { continue }
    $log = "$out\log_run_${cell}_${code}.log"
    if (Test-Path $log) { $log = "$out\log_run_${cell}_${code}_$(Get-Date -Format HHmmss).log" }
    "$(Get-Date -Format s) start $cell $code $nmax" | Out-File -Append "$out\queue.log"
    $p = Start-Process $matlabExe -ArgumentList "-batch","`"maxNumCompThreads(2); cd('$root\tests'); r0_pilot_runs('$cell','$code',$nmax)`"","-logfile","`"$log`"" -WindowStyle Hidden -PassThru -Wait
    "$(Get-Date -Format s) end $cell $code exit $($p.ExitCode)" | Out-File -Append "$out\queue.log"
}
"$(Get-Date -Format s) queue done" | Out-File -Append "$out\queue.log"
