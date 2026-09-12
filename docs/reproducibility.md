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

## Atlas-PPI extension

Run from the companion `synth` checkout. The app name is isolated from
`synth-atlas-dev` and `synth-atlas-prod`.

```powershell
uv run modal run scripts/core/abca4_variant_proteome_screen.py --abca4-root "C:\Users\lhall\Desktop\Research\abca4_alpha_genome" --threshold 0.9 --resume
```

Use `--smoke` to screen WT and one early missense product against 128 cached
human proteins. `--resume` reuses a completed remote run keyed by the score-file
schema, model revision, selected threshold, sweep grid, and ordered sequence
hashes. The full run validates that a freshly embedded WT exactly reproduces
both cached human ABCA4 projection vectors before accepting any scores.

The original H100 worker performs its prespecified screen before export. The 152.9 MB
float32 score matrix is downloaded as a stream with an incremental SHA-256
check, then opened locally only in memory-mapped mode for shape and dtype
validation. This avoids a second in-memory copy and avoids repeated local scans.
The matrix is intentionally ignored by Git and remains a local analysis artifact.

The final sensitivity grid is generated without re-embedding or rescoring:

```powershell
uv run scripts/gnomad_abca4_controls.py --refresh
uv run modal run scripts/core/abca4_variant_proteome_screen.py::resweep --abca4-root "C:\Users\lhall\Desktop\Research\abca4_alpha_genome" --threshold-min 0.33 --threshold-max 0.99 --threshold-step 0.01
```

The first command retrieves one compact gnomAD v4 gene response and caches it.
The second command scans the saved matrix on a CPU-only Modal worker across 68
thresholds, including the exact 0.3380771279335022 operating point. The laptop
loads only small indexes, annotations, and result tables.

```python
import numpy as np

scores = np.load("results/atlas_ppi/atlas_scores.float32.npy", mmap_mode="r")
```

Use `query_reference` and `human_reference` to resolve matrix rows and columns.
The small `threshold_sweep` table can be inspected without opening the matrix.

Focused offline checks in the `synth` checkout cover HGVS normalization,
canonical-transcript preference, reference-residue validation, truncation,
deduplication, float32 score persistence, memory mapping, strict thresholding,
threshold selection, FDR correction, enrichment, sequence-level aggregation,
population-control labeling, secret argument rejection, and 300 dpi PNG output.
