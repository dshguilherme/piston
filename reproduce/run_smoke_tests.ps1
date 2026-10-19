# Wrapper for run_smoke_tests.py. Uses $env:REPRO_PYTHON if set, else 'python' from PATH
# (activate the analysis environment first: conda activate piston-analysis).
# Extra arguments are passed through, e.g.  .\run_smoke_tests.ps1 --no-tex   or   --only factor_two,tables
$py = if ($env:REPRO_PYTHON) { $env:REPRO_PYTHON } else { 'python' }
& $py (Join-Path $PSScriptRoot 'run_smoke_tests.py') @args
exit $LASTEXITCODE
