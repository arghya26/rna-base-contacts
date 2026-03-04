# examples/run_pipeline.ps1
#
# End-to-end example run of the rna-base-contacts pipeline on Windows.
#
# Usage (in PowerShell):
#   .\examples\run_pipeline.ps1 -Pdb 1msy.pdb -Pairs GU,GC,AU -NTop 10
#
# Defaults: Pairs=GU, Cut=3.5, NTop=10

param(
    [string]$Pdb    = "1msy.pdb",
    [string]$Pairs  = "GU",
    [float] $Cut    = 3.5,
    [int]   $NTop   = 10
)

$PdbId   = [System.IO.Path]::GetFileNameWithoutExtension($Pdb)
$Results = "results\$PdbId"

Write-Host "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"
Write-Host "  rna-base-contacts pipeline"
Write-Host "  PDB   : $Pdb"
Write-Host "  Pairs : $Pairs"
Write-Host "  Cut   : $Cut A"
Write-Host "  TopN  : $NTop"
Write-Host "  Out   : $Results\"
Write-Host "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━"

New-Item -ItemType Directory -Force -Path $Results | Out-Null

# ── Step 1: find contacts, write mini-PDBs ────────────────────────────────
Write-Host "`n[1/4] Finding H...O contacts ..."
python scripts\find_rna_base_contacts_NHO.py `
    --pdbin      $Pdb `
    --pairs      $Pairs `
    --cut        $Cut `
    --summary    "$Results\NO_summary.txt" `
    --summary_HO "$Results\HO_summary.txt"

$miniPdbs = Get-ChildItem -Path $Results -Filter "*.pdb" -ErrorAction SilentlyContinue
Write-Host "    Mini-PDBs written: $($miniPdbs.Count)"

if ($miniPdbs.Count -eq 0) {
    Write-Host "    No pairs found — check --cut or --pairs."
    exit 0
}

# ── Step 2: score mini-PDBs ───────────────────────────────────────────────
Write-Host "`n[2/4] Scoring mini-PDBs ..."
$pdbList = ($miniPdbs | ForEach-Object { $_.FullName }) -join " "
Invoke-Expression "python scripts\compute_rna_base_miniPDBs.py --pairs $Pairs --cut $Cut $pdbList > `"$Results\scored.txt`""
Write-Host "    Scored output: $Results\scored.txt"

# ── Step 3: rank and select top N ─────────────────────────────────────────
Write-Host "`n[3/4] Ranking top $NTop candidates ..."
python scripts\rank_rna_base_contact_miniPDBs.py `
    "$Results\scored.txt" $NTop --pairs $Pairs `
    | Out-File -Encoding utf8 "$Results\ranked_top$NTop.txt"
Write-Host "    Ranked output: $Results\ranked_top$NTop.txt"

# ── Step 4: plot H···O statistics ─────────────────────────────────────────
Write-Host "`n[4/4] Plotting H...O contact statistics ..."
python scripts\plot_gu_HO.py "$Results\HO_summary.txt"
Write-Host "    Figures written to current directory."

Write-Host "`nDone. Results are in $Results\"
