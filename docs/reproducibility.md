# Reproducibility and verification

## Environment

The project targets Python 3.12 and is locked by `uv.lock`. VCF support is not
required; the workflow uses TSV and Parquet to avoid platform-specific `pysam`
build failures.

```powershell
$env:UV_CACHE_DIR = "$PWD\.uv-cache"
uv sync --python 3.12 --all-groups --locked
```

## Tests

Run the offline suite:

```powershell
uv run pytest -m "not live"
uv run ruff check .
uv run scripts/verify_no_secret_leaks.py
```

The offline tests cover secret redaction, the gene-interval parser regression,
GRCh38 allele selection, negative-strand handling, clinical classification,
deduplication, splice-distance tiers, checkpoint resume, metrics, Atlas URLs, and
300 dpi plot output.

Live smoke tests are opt-in because they use external API quotas:

```powershell
uv run pytest -m live
```

## Generated and tracked artifacts

Raw ClinVar responses, normalized tables, final score tables, reports, and plots
are tracked for inspection. Local API checkpoints, downloaded annotation caches,
and the externally hosted AlphaGenome Atlas methods PDF are ignored. The complete
workflow recreates those ignored files.

`results/provenance_manifest.json` records SHA-256 hashes for analysis artifacts
that belong in the repository. A report-only rerun refreshes metrics, plots,
acceptance checks, and the manifest without repeating API calls:

```powershell
uv run scripts/run_analysis.py report
uv run scripts/make_showcase_plots.py
```

## Credential handling

`.secrets.env` and `.env` are ignored. The loader maps `NCBI` to
`NCBI_API_KEY` and `ALPHA_GENOME` to `ALPHAGENOME_API_KEY`, while preserving
already-set process environment variables. Credential values are never included
in logs, tables, plots, or the provenance manifest.
