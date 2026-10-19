# Experiment 8 queue runner: one MATLAB session at a time.
# Reads results/r0_revision/longruns/queue.txt before every job (so jobs can be appended while it
# runs), skips lines already listed in queue_done.txt or queue_failed.txt, runs the first remaining
# job with "matlab -batch r0_long_runs('<job>')", and logs start / PIDs / end to queue.log.
# A file longruns/STOP stops the queue (the running MATLAB job also stops at its next 25-iteration
# mark and keeps its checkpoint; re-queueing the same line resumes it).
$root = "C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course"
$dir = Join-Path $root "results\r0_revision\longruns"
$qf = Join-Path $dir "queue.txt"; $donef = Join-Path $dir "queue_done.txt"; $failf = Join-Path $dir "queue_failed.txt"
$log = Join-Path $dir "queue.log"; $stop = Join-Path $dir "STOP"
function Log($m) { Add-Content -Path $log -Value ("{0} {1}" -f (Get-Date -Format "yyyy-MM-ddTHH:mm:ss"), $m) -Encoding utf8 }
Log ("runner start pid {0}" -f $PID)
while ($true) {
    if (Test-Path $stop) { Log "STOP file present; runner exits"; break }
    $jobs = @(Get-Content $qf | Where-Object { $_.Trim() -and -not $_.Trim().StartsWith('#') } | ForEach-Object { $_.Trim() })
    $done = @(); if (Test-Path $donef) { $done = @(Get-Content $donef) }
    $fail = @(); if (Test-Path $failf) { $fail = @(Get-Content $failf) }
    # a failed job (e.g. a transient out-of-memory on the shared machine) is retried from its
    # checkpoint up to 3 attempts in total; queue_failed.txt holds one line per failed attempt
    $next = $jobs | Where-Object { ($done -notcontains $_) -and (@($fail -eq $_).Count -lt 3) } | Select-Object -First 1
    if ($next -and (@($fail -eq $next).Count -gt 0)) { Log ("retrying {0} (failed attempts so far {1}) after 120 s" -f $next, @($fail -eq $next).Count); Start-Sleep 120 }
    if (-not $next) { Log "queue empty; runner exits"; break }
    # one Experiment-8 MATLAB session at a time; sessions of other users/agents are left alone
    $mine = @(Get-CimInstance Win32_Process -Filter "Name like '%MATLAB%'" -ErrorAction SilentlyContinue | Where-Object { $_.CommandLine -match 'r0_long_' })
    if ($mine.Count -gt 0) { Log ("an Experiment-8 MATLAB is still running (pids {0}); waiting 60 s" -f (($mine | ForEach-Object { $_.ProcessId }) -join ',')); Start-Sleep 60; continue }
    $safe = $next -replace '\s+', '_'
    $out = Join-Path $dir ("log_" + $safe + ".log"); $err = Join-Path $dir ("log_" + $safe + ".err")
    if (Test-Path $out) { $k = 1; while (Test-Path (Join-Path $dir ("log_{0}_r{1}.log" -f $safe, $k))) { $k++ }; $out = Join-Path $dir ("log_{0}_r{1}.log" -f $safe, $k); $err = Join-Path $dir ("log_{0}_r{1}.err" -f $safe, $k) }
    $cmd = "r0_long_runs('$next')"
    $p = Start-Process -FilePath "matlab" -ArgumentList @('-batch', "`"$cmd`"") -WorkingDirectory (Join-Path $root "tests") -RedirectStandardOutput $out -RedirectStandardError $err -NoNewWindow -PassThru
    Log ("start {0} launcher pid {1} log {2}" -f $next, $p.Id, (Split-Path $out -Leaf))
    Start-Sleep 20
    $child = @(Get-CimInstance Win32_Process -Filter ("ParentProcessId={0}" -f $p.Id) -ErrorAction SilentlyContinue | ForEach-Object { $_.ProcessId })
    if ($child.Count -gt 0) { Log ("  {0} MATLAB child pid(s) {1}" -f $next, ($child -join ',')) }
    $p.WaitForExit()
    $code = $p.ExitCode
    $ok = Select-String -Path $out -Pattern "R0_LONG_.*_DONE" -Quiet
    Log ("end {0} exit {1} done-marker {2}" -f $next, $code, $ok)
    if ($ok) {
        if (Test-Path $stop) { Log ("{0} stopped by STOP file; not marked done" -f $next) }
        else { Add-Content -Path $donef -Value $next -Encoding utf8 }
    } else { Add-Content -Path $failf -Value $next -Encoding utf8; Log ("{0} FAILED; see {1}" -f $next, (Split-Path $err -Leaf)) }
}
