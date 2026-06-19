$ErrorActionPreference = "Stop"

$root = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$stage = Join-Path $root "build\wheel-win"
$buildRoot = Join-Path $root "build"
$outDir = Join-Path $root "dist\windows"
$python = if ($env:PYTHON) { $env:PYTHON } else { "python" }

if ((Test-Path $stage) -and -not ((Resolve-Path $stage).Path.StartsWith((Resolve-Path $buildRoot).Path))) {
    throw "Refusing to remove unexpected path: $stage"
}

if (Test-Path $stage) {
    Remove-Item -LiteralPath $stage -Recurse -Force
}

New-Item -ItemType Directory -Path $stage | Out-Null
New-Item -ItemType Directory -Path $outDir -Force | Out-Null
Copy-Item -LiteralPath `
    (Join-Path $root "bfee_docking"), `
    (Join-Path $root "pyproject.toml"), `
    (Join-Path $root "setup.py"), `
    (Join-Path $root "MANIFEST.in"), `
    (Join-Path $root "README.md"), `
    (Join-Path $root "LICENSE") `
    -Destination $stage -Recurse

Remove-Item -LiteralPath (Join-Path $stage "bfee_docking\third_party\vina\vina") -Force -ErrorAction SilentlyContinue
Remove-Item -LiteralPath (Join-Path $stage "bfee_docking\third_party\smina\smina") -Force -ErrorAction SilentlyContinue
Remove-Item -LiteralPath (Join-Path $stage "bfee_docking\third_party\qvina\qvina2") -Force -ErrorAction SilentlyContinue
Remove-Item -LiteralPath (Join-Path $stage "bfee_docking\third_party\qvina\qvina2_split") -Force -ErrorAction SilentlyContinue
Remove-Item -LiteralPath (Join-Path $stage "bfee_docking\third_party\qvina\qvinaw") -Force -ErrorAction SilentlyContinue
Remove-Item -LiteralPath (Join-Path $stage "bfee_docking\third_party\qvina\qvinaw_split") -Force -ErrorAction SilentlyContinue
Remove-Item -LiteralPath (Join-Path $stage "bfee_docking\third_party\qvina\vina") -Force -ErrorAction SilentlyContinue
Remove-Item -LiteralPath (Join-Path $stage "bfee_docking\third_party\qvina\vina_split") -Force -ErrorAction SilentlyContinue
Get-ChildItem -LiteralPath (Join-Path $stage "bfee_docking\third_party\qvina") -Filter "*.so*" -File -ErrorAction SilentlyContinue |
    Remove-Item -Force

$env:BFEE_DOCKING_PLATFORM_TAG = "win_amd64"
& $python -m build --wheel --no-isolation --outdir $outDir $stage
Remove-Item Env:\BFEE_DOCKING_PLATFORM_TAG -ErrorAction SilentlyContinue
