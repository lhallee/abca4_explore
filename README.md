# ABCA4 variant impact with AlphaGenome

This repository tests whether AlphaGenome Variant Impact (AVI) scores agree with
clinical classifications and published functional measurements for `ABCA4`, the
gene associated with Stargardt disease 1. The workflow retrieves ClinVar records,
resolves allele-level GRCh38 coordinates, scores eligible SNVs, and evaluates
coding and splicing variants separately.

## Main result

The analysis retrieved 4,670 ClinVar records and scored 4,024 unique, precise
GRCh38 SNVs. In the primary cohort of 2,521 variants with at least one ClinVar
review star, AVI separated pathogenic or likely pathogenic from benign or likely
benign variants with ROC-AUC 0.990 and average precision 0.991.

This is strong retrospective separation, but it is not a clean prospective test.
The base AlphaGenome model learned experimental functional-genomics tracks, not
ABCA4 ClinVar labels. The downstream AVI ensemble was trained on 37 million
gnomAD variants using allele frequency as a proxy label, and chromosome 1 was in
its training partition. The public score API does not expose the sampled training
coordinates, so overlap with ABCA4 positions cannot be excluded here. See the
[interpretation report](results/impressiveness_assessment.md) for the evidence and
limits of the claim.

## Figures

![ABCA4 AVI overview](plots/showcase/abca4_avi_overview.png)

**Figure 1. Clinical classification performance.** Cumulative AVI score
distributions, ROC and precision-recall curves, and cohort-specific performance
for precise GRCh38 SNVs. The primary cohort contains 1,092 pathogenic or likely
pathogenic and 1,429 benign or likely benign variants. Its ROC-AUC is 0.990 and
average precision is 0.991. The evaluation is retrospective and cannot exclude
positions that overlap the AVI training set.

![AVI signal sources](plots/showcase/abca4_avi_signal_sources.png)

**Figure 2. Signal-source analysis.** ROC-AUC from the complete AVI score and
sums of signed AVI feature attributions. AlphaMissense and conservation account
for most coding-missense separation. Splicing attribution is strongest for
splice-region and deep-intronic variants, with ROC-AUC 0.991 and 0.953,
respectively. Attribution-only rows diagnose the trained AVI model; they are not
independently retrained ablations.

![Functional assay correlations](plots/showcase/abca4_avi_functional_validation.png)

**Figure 3. Published functional-assay comparisons.** Study-specific Spearman
correlations between predicted impact and retained molecular function, with
fixed-seed 95% bootstrap intervals. Negative correlations are expected because
greater predicted impact should correspond to lower expression, ATPase activity,
or correctly spliced mRNA. The larger Garces and Sangermano panels show moderate,
mechanism-consistent associations. The Aslaksen estimates are uncertain because
only eight variants were matched.

## Reproduce the analysis

Requirements:

- Python 3.12
- [`uv`](https://docs.astral.sh/uv/)
- an AlphaGenome API key accepted under the
  [AlphaGenome Terms of Use](https://deepmind.google/science/alphagenome/terms)
- an optional NCBI API key for higher E-utilities rate limits

Create `.secrets.env` in the repository root:

```dotenv
NCBI=your_ncbi_key
ALPHA_GENOME=your_alphagenome_key
```

The file is ignored by Git. Credentials are loaded only at runtime, canonical
environment variables take precedence, and values are never printed or copied to
outputs.

Run the complete resumable workflow on Windows:

```powershell
.\scripts\run_all.ps1
```

Or run the equivalent commands:

```powershell
$env:UV_CACHE_DIR = "$PWD\.uv-cache"
uv sync --python 3.12 --all-groups --locked
uv run scripts/run_analysis.py all --resume
uv run scripts/make_showcase_plots.py
uv run scripts/verify_no_secret_leaks.py
```

Use `--max-variants 5` for an opt-in smoke run. The full run checkpoints every
successful AlphaGenome query and each ClinVar batch, so interrupted runs resume
without duplicate requests.

## Command-line entry points

```powershell
uv run scripts/clinvar_api.py search --query "ABCA4[gene]" --output data/raw/clinvar_ids.json
uv run scripts/alphagenome_atlas_avi.py annotate -i data/processed/variants.tsv -o results/tables/variants_scored.tsv --include_features
uv run scripts/alphagenome_atlas_links.py table --input results/tables/variants_scored.tsv --output results/tables/variants_with_links.tsv
uv run scripts/run_analysis.py all --resume
```

## Repository contents

- `src/abca4_avi/`: retrieval, normalization, scoring, statistics, and plotting
- `scripts/`: compatible entry points and the one-command PowerShell runner
- `tests/`: offline regression tests and opt-in live API smoke tests
- `data/raw/`: checkpointed ClinVar and literature responses
- `data/processed/`: normalized ClinVar records and unique scoring variants
- `data/reference/`: curated functional-assay inputs and source catalog
- `results/`: score tables, metrics, exclusions, provenance, and reports
- `plots/`: 300 dpi PNG figures
- `docs/`: methods, data dictionary, and reproducibility notes

Generated caches and the externally hosted AlphaGenome Atlas paper are not
tracked. A clone can regenerate them through the locked workflow.

## Documentation

- [Methods and cohort definitions](docs/methods.md)
- [Output data dictionary](docs/data_dictionary.md)
- [Reproducibility and verification](docs/reproducibility.md)
- [Concise results](results/summary.md)
- [Interpretation and training-overlap caveat](results/impressiveness_assessment.md)

## Sources and terms

The workflow uses [NCBI ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/),
[AlphaGenome Atlas](https://deepmind.google/science/alphagenome/), GENCODE v46,
and the published assays cataloged in
[`data/reference/functional_studies.tsv`](data/reference/functional_studies.tsv).
AlphaGenome API outputs and Atlas information are restricted to non-commercial
use unless the applicable terms provide otherwise, and must not be used to train
other machine-learning models. Software in this repository does not grant rights
to third-party data or model outputs. LOVD is not scraped.
