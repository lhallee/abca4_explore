# Output data dictionary

## Main score table

`results/tables/variants_scored.tsv` and its Parquet equivalent contain one row
per unique, precise GRCh38 SNV.

| Field group | Columns | Meaning |
|---|---|---|
| Variant identity | `variant`, `chromosome`, `position`, `ref`, `alt` | One-based GRCh38 allele identity. |
| ClinVar identity | `clinvar_variation_id`, `vcv_accession`, `rsid` | Source record identifiers. |
| Clinical evidence | `clinical_significance_raw`, `clinical_class`, `review_status`, `review_stars` | Original aggregate assertion and normalized analysis class. |
| Disease context | `conditions`, `condition_identifiers`, `is_stgd1` | Submitted conditions and STGD1 subset flag. |
| Molecular annotation | `hgvs_c`, `hgvs_p`, `molecular_consequences`, `region_class`, `nearest_exon_distance` | Transcript-relative notation and region class. |
| AVI summary | `avi_raw`, `avi_tail_probability`, `avi_phred`, `top_percentile`, `top_modality` | AlphaGenome Variant Impact outputs. |
| AVI attribution | 18 columns beginning with `fi_` | Signed feature attribution weights from the AVI ensemble. |
| Navigation | `atlas_url` | AlphaGenome Atlas deep link for the scored allele. |
| Audit | `clinvar_variation_ids`, `sentinel_name` | Deduplicated source records and sentinel annotation. |

Higher `avi_phred` indicates a more extreme predicted molecular impact. Phred 20
and 30 correspond to the top 1% and 0.1% of genome-wide AVI scores,
respectively. Attribution weights describe contributions within the trained AVI
ensemble and are not independently retrained ablations.

## Supporting tables

- `clinvar_exclusions`: unscored records and explicit exclusion reasons
- `sentinel_variants`: the four required ABCA4 sentinel variants
- `selected_track_attributions`: detailed tracks for sentinels and top candidates
- `functional_assays_scored`: published assay measurements joined to AVI scores
- `showcase_signal_diagnostics`: region-specific attribution discrimination
- `showcase_functional_correlations`: study-specific rank correlations
- `metrics.json`: cohort metrics, thresholds, and bootstrap confidence intervals
- `acceptance_checks.json`: machine-readable completion checks

All TSV tables have Parquet equivalents where practical. The provenance manifest
records file hashes, source URLs, package versions, retrieval time, cohort
definitions, and record counts.
