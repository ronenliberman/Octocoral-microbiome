# =============================================================================
# 11_run_alignment_and_tree.ps1 — MAFFT + IQ-TREE on Windows (replaces the Mac step)
# =============================================================================
# Run after phase 1 of 03 and 04 (which write *_for_tree.fasta into phylogeny\bacteria\).
# Tools (portable, no install), unpacked in $HOME\tools\ (or $env:OCTO_TOOLS):
#   MAFFT v7.526  (mafft-7.526-win64-signed.zip, mafft.cbrc.jp)
#   IQ-TREE v3.1.4 (iqtree-3.1.4-Windows.zip, github.com/iqtree/iqtree3)
# Settings match the Methods:
#   MAFFT L-INS-i (--localpair --maxiterate 1000)
#   IQ-TREE ModelFinder (-m MFP), 1000 ultrafast bootstraps, 1000 SH-aLRT
# Work is done in a folder without spaces (MAFFT's Windows package dislikes
# them); results are copied back into phylogeny\bacteria\, where 03/04 read them.
#
# Usage (PowerShell):  & ".\11_run_alignment_and_tree.ps1"
# =============================================================================

# "Continue", not "Stop": Windows PowerShell 5.1 turns any stderr line from a
# native program (MAFFT prints "Active code page: 65001") into a fatal error.
# Failures are caught by the explicit output checks below instead.
$ErrorActionPreference = "Continue"

# Tools folder: set $env:OCTO_TOOLS to override the default ($HOME\tools).
$tools = if ($env:OCTO_TOOLS) { $env:OCTO_TOOLS } else { Join-Path $HOME "tools" }
$mafft  = Join-Path $tools "mafft-win\mafft.bat"
$iqtree = Join-Path $tools "iqtree-3.1.4-Windows\bin\iqtree3.exe"
$macDir = (Resolve-Path (Join-Path $PSScriptRoot "..\..\phylogeny\bacteria")).Path
$work   = Join-Path $tools "work"
New-Item -ItemType Directory -Force $work | Out-Null

foreach ($prefix in @("Endo_r3", "Myco_r3")) {
    $in = Join-Path $macDir "${prefix}_for_tree.fasta"
    if (-not (Test-Path $in)) { Write-Host "SKIP $prefix : $in not found (run phase 1 of 03/04 first)"; continue }

    Write-Host "`n===== $prefix ====="
    Copy-Item $in (Join-Path $work "${prefix}_for_tree.fasta") -Force
    Push-Location $work
    try {
        # --- MAFFT L-INS-i ---
        Write-Host "MAFFT L-INS-i ..."
        & $mafft --localpair --maxiterate 1000 --quiet --out "${prefix}_mafft_aligned.fasta" "${prefix}_for_tree.fasta" 2>&1 |
            Where-Object { $_ -notmatch "Active code page" } | Out-Host
        if (-not (Test-Path "${prefix}_mafft_aligned.fasta") -or (Get-Item "${prefix}_mafft_aligned.fasta").Length -eq 0) {
            throw "MAFFT produced no alignment for $prefix"
        }
        $nIn  = (Select-String "${prefix}_for_tree.fasta" -Pattern '^>').Count
        $nOut = (Select-String "${prefix}_mafft_aligned.fasta" -Pattern '^>').Count
        Write-Host "  aligned $nOut / $nIn sequences"
        if ($nIn -ne $nOut) { throw "Sequence count changed during alignment" }

        # --- IQ-TREE ---
        Write-Host "IQ-TREE (ModelFinder, 1000 UFBoot, 1000 SH-aLRT) ..."
        Get-ChildItem -Filter "${prefix}_iqtree.*" | Remove-Item -Force
        & $iqtree -s "${prefix}_mafft_aligned.fasta" -m MFP -B 1000 -alrt 1000 -T AUTO `
                  --seed 123 --prefix "${prefix}_iqtree" --quiet
        if (-not (Test-Path "${prefix}_iqtree.treefile")) { throw "IQ-TREE produced no tree for $prefix" }
        Select-String "${prefix}_iqtree.iqtree" -Pattern "^Input data:|parsimony-informative|^Best-fit model according to BIC" |
            ForEach-Object { Write-Host "  $($_.Line.Trim())" }

        # --- copy back ---
        foreach ($f in @("${prefix}_mafft_aligned.fasta", "${prefix}_iqtree.treefile",
                         "${prefix}_iqtree.iqtree", "${prefix}_iqtree.log", "${prefix}_iqtree.contree")) {
            if (Test-Path $f) { Copy-Item $f $macDir -Force }
        }
        Write-Host "  results copied to phylogeny\bacteria\"
    }
    finally { Pop-Location }
}
Write-Host "`nDone. Now re-run 03_Endo_tree_revised.R and 04_Myco_tree_revised.R (phase 2)."
