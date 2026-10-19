# REPRODUCE-PACKAGE COPY. Do not edit numerics here; edit the original and re-copy.
# Original: paper_jsv/reproduce_plan/rescued_scratch/build.ps1  (sha256 7f68dd6b9cb48fba)
# Changed: original lines 2-3 ($proto, $root absolute paths) -> derived from REPRO_ROOT ($env:REPRO_ROOT if set, else
#          (Resolve-Path "$PSScriptRoot\..\..").Path): $root = <REPRO_ROOT>\paper_jsv, $proto = <REPRO_ROOT>\paper_jsv\prototype;
#          original line 4 ($out = $env:TEMP\trim_build<Sub>) -> $env:REPRO_OUT (fallback "$PSScriptRoot\..\out")\manuscript_build<Sub>.
#          The pdflatex/bibtex sequence, TEXINPUTS/BIBINPUTS and the printed counters are unchanged. The build reads the manuscript, the
#          figures (pdf/png/tex/data) and the .bib files from the paper tree (read only) and writes only into the output folder above.
param([string]$Src = "manuscript_proto_voice.tex", [string]$Sub = "")
if ($env:REPRO_ROOT) { $reproRoot = $env:REPRO_ROOT } else { $reproRoot = (Resolve-Path "$PSScriptRoot\..\..").Path }
if ($env:REPRO_OUT) { $outBase = $env:REPRO_OUT } else { $outBase = [System.IO.Path]::GetFullPath("$PSScriptRoot\..\out") }
$root = Join-Path $reproRoot "paper_jsv"
$proto = Join-Path $root "prototype"
$out = Join-Path $outBase ("manuscript_build" + $Sub)
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
