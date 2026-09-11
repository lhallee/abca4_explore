$ErrorActionPreference = "Stop"
$projectRoot = (Resolve-Path (Join-Path $PSScriptRoot "..")).Path
$env:UV_CACHE_DIR = Join-Path $projectRoot ".uv-cache"

Push-Location $projectRoot
try {
    uv sync --python 3.12 --all-groups
    uv run python scripts/run_analysis.py all --resume
}
finally {
    Pop-Location
}
