param([string]$Src = "manuscript_proto_voice.tex", [string]$Sub = "")
$proto = "C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course\paper_jsv\prototype"
$root = "C:\Users\dshgu\OneDrive\Documents\FastFormation\harmonic_shell_course\paper_jsv"
$out = Join-Path $env:TEMP ("trim_build" + $Sub)
New-Item -ItemType Directory -Force $out | Out-Null
$env:TEXINPUTS = ".;$proto;$proto\figs;$root;"
$env:BIBINPUTS = ".;$proto;$root;"
$job = [System.IO.Path]::GetFileNameWithoutExtension($Src)
Set-Location $proto
pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$out" $Src | Out-Null
Push-Location $out
bibtex $job | Out-Null
Pop-Location
pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$out" $Src | Out-Null
pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$out" $Src | Out-Null
$log = Get-Content (Join-Path $out "$job.log") -Raw
$pages = (pdfinfo (Join-Path $out "$job.pdf") | Select-String "Pages").ToString()
"pages: $pages"
"errors: " + ([regex]::Matches($log, "(?m)^! ").Count)
"undefined: " + ([regex]::Matches($log, "undefined").Count)
"overfull: " + ([regex]::Matches($log, "Overfull").Count)
"underfull: " + ([regex]::Matches($log, "Underfull").Count)
"float too large: " + ([regex]::Matches($log, "Float too large").Count)
$blg = Join-Path $out "$job.blg"
if (Test-Path $blg) { "bibtex warnings: " + ([regex]::Matches((Get-Content $blg -Raw), "Warning").Count) }
