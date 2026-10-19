param([string]$List)
# sequential pilot queue: one MATLAB session at a time; skips runs already saved.
# items: cell:code:nmax  (code 'eval' runs the evaluation stage). MATLAB writes its own -logfile,
# so it does not depend on this script's stdout.
$root = "C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course"
$out = "$root\results\r0_revision\pilot"
foreach ($item in $List.Split(',')) {
    $cell, $code, $nmax = $item.Split(':')
    if ($code -ne 'eval' -and (Test-Path "$out\run_${cell}_${code}.mat")) { continue }
    $log = "$out\log_run_${cell}_${code}.log"
    if (Test-Path $log) { $log = "$out\log_run_${cell}_${code}_$(Get-Date -Format HHmmss).log" }
    "$(Get-Date -Format s) start $cell $code $nmax" | Out-File -Append "$out\queue.log"
    $p = Start-Process "C:\Program Files\MATLAB\R2025a\bin\matlab.exe" -ArgumentList "-batch","`"maxNumCompThreads(2); cd('$root\tests'); r0_pilot_runs('$cell','$code',$nmax)`"","-logfile","`"$log`"" -WindowStyle Hidden -PassThru -Wait
    "$(Get-Date -Format s) end $cell $code exit $($p.ExitCode)" | Out-File -Append "$out\queue.log"
}
"$(Get-Date -Format s) queue done" | Out-File -Append "$out\queue.log"
