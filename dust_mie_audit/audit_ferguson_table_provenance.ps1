# Read-only comparison. Streams one tar member to memory; no extraction or source edits.
# Run: powershell -NoProfile -File .\audit_ferguson_table_provenance.ps1
param(
    [string]$Archive = 'C:\Users\hila\Dropbox\RMHD\IDORT\AGB\f05.gs98.tar.gz',
    [string]$Member = 'g98.7.02.tron',
    [string]$RossTable = 'C:\Users\hila\Dropbox\RMHD\IDORT\RSGg_m10lgl4.6_res6ib03n3s_r6s_new4_kap\RSGg_m10lgl4.6_res6ib03n3s_r6s_new4_kap\gs98_z0.02_x0.7_ma_ross.data'
)
$ErrorActionPreference = 'Stop'
$raw = @(tar -xOzf $Archive $Member)
if ($LASTEXITCODE -ne 0) { throw 'Unable to read archive member' }
$target = Get-Content -LiteralPath $RossTable
$sourceRows = @($raw | Where-Object { $_ -match '^\s*[234]\.\d+\s' } |
    ForEach-Object { ,([double[]]($_.Trim() -split '\s+')) })
$targetRows = @($target[7..186] | ForEach-Object { ,([double[]]($_.Trim() -split '\s+')) })
$sourceR = [double[]](($raw[3].Trim() -split '\s+')[2..20])
# The archive's header starts with the two tokens 'log T'; parse its 19 logR coordinates.
$targetR = [double[]]($target[5].Trim() -split '\s+')
$comparison = @()
foreach ($source in $sourceRows) {
    $matches = @($targetRows | Where-Object { [Math]::Abs($_[0] - $source[0]) -lt 1e-8 })
    if ($matches.Count -ne 1) { continue }
    for ($j = 0; $j -lt $sourceR.Count; $j++) {
        $column = [Array]::IndexOf($targetR, $sourceR[$j])
        if ($column -lt 0) { continue }
        $difference = $matches[0][$column + 1] - $source[$j + 1]
        $comparison += [pscustomobject]@{ LogT=$source[0]; LogR=$sourceR[$j];
            Difference=$difference; AbsDifference=[Math]::Abs($difference) }
    }
}
$low = @($comparison | Where-Object { $_.LogT -le 3.5 })
[pscustomobject]@{
    Archive=$Archive
    ArchiveMember=$Member
    RossTable=$RossTable
    ArchiveHeader=$raw[0].ToString()
    RossHeader=$target[0].ToString()
    ArchiveTemperatureRows=$sourceRows.Count
    SharedTemperatureValues=@($comparison.LogT | Sort-Object -Unique).Count
    SharedCells=$comparison.Count
    ExactMatches=@($comparison | Where-Object { $_.Difference -eq 0 }).Count
    MaximumAbsoluteDifference=($comparison.AbsDifference | Measure-Object -Maximum).Maximum
    LowTemperatureLimitLogT=3.5
    LowTemperatureSharedCells=$low.Count
    LowTemperatureExactMatches=@($low | Where-Object { $_.Difference -eq 0 }).Count
    ArchiveSHA256=(Get-FileHash -LiteralPath $Archive -Algorithm SHA256).Hash
    RossTableSHA256=(Get-FileHash -LiteralPath $RossTable -Algorithm SHA256).Hash
    Scope='Only these accessible local files were compared. No current cluster-copy hashes or Planck generation provenance were verified.'
} | ConvertTo-Json -Depth 4
