# 01d (local) — type-strain BLAST of the 26 ASVs with significant ANCOM-BC2 results (Appendix 1B)
# Used instead of 01d_ANCOMBC2_typestrain_BLAST.R because NCBI web BLAST downloads of the multi-query
# hit table were incomplete (7 Oct 2026: 3/26, then 3/13 and 2/13 queries). Same database:
# NCBI 16S_ribosomal_RNA (type strains), local copy downloaded 1 Oct 2026; BLAST+ 2.17.0.
# megablast, max 10 targets; queries without a megablast hit are re-run with blastn.
# Query:  output\Editor_round3\ANCOMBC2_r3_ASVs.fasta (written by 01d_ANCOMBC2_typestrain_BLAST.R).
# Output: blast\ANCOMBC2_r3_16S_hits.csv (12-column hit table, no header) + SEARCH_INFO.
# Paths are relative to the repository; the BLAST+ install and the 16S database are
# looked for in <tools>\audit_refs.
$ErrorActionPreference = "Stop"
# Tools folder: set $env:OCTO_TOOLS to override the default ($HOME\tools).
$tools = if ($env:OCTO_TOOLS) { $env:OCTO_TOOLS } else { Join-Path $HOME "tools" }
$tools = Join-Path $tools "audit_refs"
$repo  = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path
$bin   = Join-Path (Get-ChildItem $tools -Directory | Where-Object Name -like "ncbi-blast*" | Select-Object -First 1).FullName "bin"
$db    = Join-Path $tools "db16S\16S_ribosomal_RNA"
$o     = Join-Path $repo "blast"
$q     = Join-Path $repo "output\Editor_round3\ANCOMBC2_r3_ASVs.fasta"
$fmt   = "6 qseqid saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore"
$out1  = Join-Path $env:TEMP "ancom_mega.tsv"; $out2 = Join-Path $env:TEMP "ancom_blastn.tsv"

& "$bin\blastn.exe" -task megablast -db $db -query $q -outfmt $fmt -max_target_seqs 10 -out $out1
if ($LASTEXITCODE -ne 0) { throw "megablast failed" }
$ids  = Get-Content $q | Where-Object { $_ -like ">*" } | ForEach-Object { $_.Substring(1).Split(" ")[0] }
$rows = @(Get-Content $out1)
$hit  = $rows | ForEach-Object { $_.Split("`t")[0] } | Sort-Object -Unique
$miss = $ids | Where-Object { $_ -notin $hit }
"megablast: $($hit.Count)/$($ids.Count) queries with hits"
if ($miss) {
  $lines = Get-Content $q; $sub = @()
  for ($i = 0; $i -lt $lines.Count; $i++) { if ($lines[$i] -like ">*" -and $lines[$i].Substring(1).Split(" ")[0] -in $miss) { $sub += $lines[$i]; $sub += $lines[$i + 1] } }
  $qs = Join-Path $env:TEMP "ancom_missing.fasta"; $sub | Set-Content $qs -Encoding ascii
  & "$bin\blastn.exe" -task blastn -db $db -query $qs -outfmt $fmt -max_target_seqs 10 -out $out2
  if ($LASTEXITCODE -ne 0) { throw "blastn failed" }
  $rows += @(Get-Content $out2)
  "blastn re-run for: $($miss -join ' ')"
}
$hit = $rows | ForEach-Object { $_.Split("`t")[0] } | Sort-Object -Unique
"total: $($hit.Count)/$($ids.Count) queries with hits; rows: $($rows.Count)"
$rows | ForEach-Object { $_ -replace "`t", "," } | Set-Content (Join-Path $o "ANCOMBC2_r3_16S_hits.csv") -Encoding ascii
@("Searched: $(Get-Date -Format 'yyyy-MM-dd HH:mm')",
  "Local BLAST+ ($bin), database NCBI 16S_ribosomal_RNA downloaded 1 Oct 2026",
  "megablast, max_target_seqs 10; blastn for queries without a megablast hit: $($miss -join ' ')",
  "Queries: $($ids.Count) | with hits: $($hit.Count)") | Set-Content (Join-Path $o "ANCOMBC2_r3_16S_hits_SEARCH_INFO.txt") -Encoding utf8
