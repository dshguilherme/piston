# Experiment 10 queue runner: one Experiment-10 MATLAB session at a time, never more than two
# MATLAB sessions on the machine in total (other agents' sessions are left alone, never killed).
# Each line of results/r0_revision/robust/queue.txt is a MATLAB command, e.g.
#   r0_robust_runs('c30_k18 X U 16 0.01 mean 800')
# The file is re-read before every job, so jobs can be appended while the runner works.
# Lines listed in queue_done.txt are skipped; a failed job is retried from its checkpoint up to 3
# attempts (queue_failed.txt holds one line per failed attempt). A file robust/STOP stops the
# runner (a running r0_robust_runs job also stops at its next 25-iteration mark).
# The fluid-loading operators are read from $env:R0_ROBUST_ZDIR (set below if empty).
$root = "C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course"
$dir = Join-Path $root "results\r0_revision\robust"
$qf = Join-Path $dir "queue.txt"; $donef = Join-Path $dir "queue_done.txt"; $failf = Join-Path $dir "queue_failed.txt"
$log = Join-Path $dir "queue.log"; $stop = Join-Path $dir "STOP"
if (-not $env:R0_ROBUST_ZDIR) { $env:R0_ROBUST_ZDIR = "C:\Users\dshgu\AppData\Local\Temp\claude\C--Users-dshgu-OneDrive-Documents-FastFormation-harmonic-shell-course-paper-jsv\90c1eee2-75c1-4619-af46-00827600c50f\scratchpad\z" }
function Log($m) { Add-Content -Path $log -Value ("{0} {1}" -f (Get-Date -Format "yyyy-MM-ddTHH:mm:ss"), $m) -Encoding utf8 }
Log ("runner start pid {0}" -f $PID)
$k = 0
while ($true) {
    if (Test-Path $stop) { Log "STOP file present; runner exits"; break }
    $jobs = @(Get-Content $qf | Where-Object { $_.Trim() -and -not $_.Trim().StartsWith('#') } | ForEach-Object { $_.Trim() })
    $done = @(); if (Test-Path $donef) { $done = @(Get-Content $donef) }
    $fail = @(); if (Test-Path $failf) { $fail = @(Get-Content $failf) }
    $next = $jobs | Where-Object { ($done -notcontains $_) -and (@($fail -eq $_).Count -lt 3) } | Select-Object -First 1
    if (-not $next) { Log "queue empty; runner exits"; break }
    $all = @(Get-CimInstance Win32_Process -Filter "Name = 'MATLAB.exe'" -ErrorAction SilentlyContinue | Where-Object { $_.ExecutablePath -match '\\win64\\' })
    $mine = @($all | Where-Object { $_.CommandLine -match 'r0_robust_' })
    if ($mine.Count -gt 0 -or $all.Count -ge 2) { Log ("waiting 60 s: {0} MATLAB sessions running ({1} of Experiment 10)" -f $all.Count, $mine.Count); Start-Sleep 60; continue }
    $k++
    $stamp = Get-Date -Format "yyyyMMdd_HHmmss"
    $out = Join-Path $dir ("log_job_{0}.log" -f $stamp); $err = Join-Path $dir ("log_job_{0}.err" -f $stamp)
    $p = Start-Process -FilePath "matlab" -ArgumentList @('-batch', "`"$next`"") -WorkingDirectory (Join-Path $root "tests") -RedirectStandardOutput $out -RedirectStandardError $err -NoNewWindow -PassThru
    Log ("start [{0}] launcher pid {1} log {2}" -f $next, $p.Id, (Split-Path $out -Leaf))
    $p.WaitForExit()
    $ok = Select-String -Path $out -Pattern "R0_ROBUST_.*_DONE" -Quiet
    Log ("end [{0}] exit {1} done-marker {2}" -f $next, $p.ExitCode, $ok)
    if ($ok -and -not (Test-Path $stop)) { Add-Content -Path $donef -Value $next -Encoding utf8 }
    elseif (-not $ok) { Add-Content -Path $failf -Value $next -Encoding utf8; Log ("[{0}] FAILED; see {1}" -f $next, (Split-Path $err -Leaf)); Start-Sleep 60 }
}
