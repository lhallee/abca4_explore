---
name: abca4-alphagenome
description: Run or maintain the Windows Python 3.12 ABCA4 ClinVar and AlphaGenome AVI workflow in this project, including secret-safe credentials, allele-level GRCh38 parsing, resumable scoring, cohort metrics, sentinel checks, and 300 dpi plots.
---

# ABCA4 AlphaGenome workflow

Use this skill for the ABCA4 project pipeline in this repository.

## Required practices

- Run with Python 3.12 through `uv`; do not add `pysam` to the required path.
- Load `.secrets.env` only through `abca4_avi.config.load_credentials`.
- Give canonical process variables precedence over aliases `NCBI` and `ALPHA_GENOME`.
- Never print or persist credential values.
- Treat ClinVar gene intervals as invalid variant coordinates. Use only allele-level GRCh38 VCF attributes, with dbSNP chromosome placement as a fallback.
- Query each unique SNV once and preserve JSONL checkpoints.
- Always pair a reported AVI score with its AlphaGenome Atlas URL.
- Frame AlphaGenome outputs as research-only predicted molecular effects, not diagnoses.
- Save plots as PNG at 300 dpi.

## Commands

Run the full workflow:

```powershell
.\scripts\run_all.ps1
```

Run a quick live smoke check:

```powershell
$env:UV_CACHE_DIR = "$PWD\.uv-cache"
uv run python scripts/run_analysis.py all --resume --max-records 200 --max-variants 5 --bootstrap-replicates 20
```

Run offline verification:

```powershell
uv run pytest -m "not live"
uv run ruff check .
```

Use `results/provenance_manifest.json` and `results/tables/acceptance_checks.json` as the completion gates.
