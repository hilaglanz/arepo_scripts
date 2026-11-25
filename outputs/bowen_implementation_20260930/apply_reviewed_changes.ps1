$ErrorActionPreference = 'Stop'
$targetRoot = (Resolve-Path -LiteralPath 'C:\Users\hila\GitRepository\arepo').Path
$stageRoot = (Resolve-Path -LiteralPath (Join-Path $PSScriptRoot 'staged')).Path
$manifest = Get-Content -Raw -LiteralPath (Join-Path $PSScriptRoot 'original_hashes.json') | ConvertFrom-Json
$newFiles = @('documentation\bowen_opacity.md', 'tests\check_bowen_opacity.py')

# Check every original before changing any file, so concurrent user edits win.
foreach ($entry in $manifest) {
    $target = [IO.Path]::GetFullPath((Join-Path $targetRoot $entry.path))
    if (-not $target.StartsWith($targetRoot + '\', [StringComparison]::OrdinalIgnoreCase)) {
        throw "Target is outside the AREPO repository: $target"
    }
    if ((Get-FileHash -Algorithm SHA256 -LiteralPath $target).Hash -ne $entry.sha256) {
        throw "Source changed since review; refusing to overwrite: $target"
    }
}
foreach ($relative in $newFiles) {
    if (Test-Path -LiteralPath (Join-Path $targetRoot $relative)) {
        throw "New-file destination already exists; refusing to overwrite: $relative"
    }
}
$relativePaths = @($manifest | ForEach-Object { $_.path }) + $newFiles
foreach ($relative in $relativePaths) {
    $source = [IO.Path]::GetFullPath((Join-Path $stageRoot $relative))
    $target = [IO.Path]::GetFullPath((Join-Path $targetRoot $relative))
    if (-not $source.StartsWith($stageRoot + '\', [StringComparison]::OrdinalIgnoreCase) -or
        -not $target.StartsWith($targetRoot + '\', [StringComparison]::OrdinalIgnoreCase)) {
        throw "Path outside the reviewed roots: $relative"
    }
    if (-not (Test-Path -LiteralPath $source -PathType Leaf)) {
        throw "Missing reviewed source: $relative"
    }
}
foreach ($relative in $relativePaths) {
    $target = Join-Path $targetRoot $relative
    New-Item -ItemType Directory -Path (Split-Path -Parent $target) -Force | Out-Null
    Copy-Item -LiteralPath (Join-Path $stageRoot $relative) -Destination $target -Force
    Write-Output "Updated $relative"
}
