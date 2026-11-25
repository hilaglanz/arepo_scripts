$ErrorActionPreference = 'Stop'
$targetRoot = (Resolve-Path -LiteralPath 'C:\Users\hila\GitRepository\arepo').Path
$stageRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot 'staged')).Path
$manifest = Get-Content -Raw -LiteralPath (Join-Path $PSScriptRoot 'original_hashes.json') | ConvertFrom-Json
$relativePaths = @('src/mesa_opacity.c', 'documentation/bowen_opacity.md', 'tests/check_bowen_opacity.py')

# Check every target before any write; do not overwrite concurrent user edits.
foreach ($relative in $relativePaths) {
    $target = [IO.Path]::GetFullPath((Join-Path $targetRoot $relative))
    $source = [IO.Path]::GetFullPath((Join-Path $stageRoot $relative))
    if (-not $target.StartsWith($targetRoot + '\', [StringComparison]::OrdinalIgnoreCase) -or
        -not $source.StartsWith($stageRoot + '\', [StringComparison]::OrdinalIgnoreCase)) {
        throw "Path outside the reviewed roots: $relative"
    }
    $expected = $manifest.PSObject.Properties[$relative].Value
    if (-not $expected -or (Get-FileHash -Algorithm SHA256 -LiteralPath $target).Hash -ne $expected) {
        throw "File changed since review; refusing to overwrite: $target"
    }
    if (-not (Test-Path -LiteralPath $source -PathType Leaf)) {
        throw "Missing reviewed source: $source"
    }
}

foreach ($relative in $relativePaths) {
    $source = Join-Path $stageRoot $relative
    $target = Join-Path $targetRoot $relative
    Copy-Item -LiteralPath $source -Destination $target -Force
    if ((Get-FileHash -Algorithm SHA256 -LiteralPath $source).Hash -ne
        (Get-FileHash -Algorithm SHA256 -LiteralPath $target).Hash) {
        throw "Verification failed after writing: $target"
    }
    Write-Output "Updated $relative"
}
