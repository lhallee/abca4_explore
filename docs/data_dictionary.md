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

## Atlas-PPI extension tables

Tables under `results/atlas_ppi/` have matching TSV and Parquet forms where
practical.

| Table | Grain | Key contents |
|---|---|---|
| `variant_protein_mapping` | Source SNV | Canonical protein annotation, mapping status, observability, sequence hash, and exclusion reason. |
| `unique_sequences` | Unique truncated product | Sequence identifier, SHA-256, length, and WT flag. Raw sequences are not included. |
| `query_reference` | Score-matrix row | Zero-based query index, sequence identifier, and sequence SHA-256. |
| `human_reference` | Score-matrix column | Zero-based human index, Atlas identifier, accession, gene symbol, and sequence SHA-256. |
| `atlas_scores.float32.npy` | Unique product by human protein | Complete merged Atlas probability matrix in NumPy NPY format. Load with `numpy.load(..., mmap_mode="r")`. |
| `thresholded_edges` | ABCA4 product and human partner | Atlas score above 0.9, WT score, score change, human accession, and gene symbol. |
| `threshold_sweep` | Atlas score threshold | WT partner count, network-disruption summaries, clinical ROC-AUC, and AVI Spearman associations. |
| `sequence_network_summary_joined` | Unique product | Partner count, gains, losses, Jaccard distance, score deltas, AVI, ClinVar, and protein consequence summaries. |
| `variant_disruption_ranking` | Unique product | Global enrichment evidence, effect size, and changed-partner ranking. |
| `ora_results` | Product and pathway | Partner overlap, genes, p-value, and within-screen BH q-value. |
| `differential_enrichment` | Product and pathway | Gained and lost genes, corrected log2 odds ratio, p-value, within-screen q-value, and global per-library q-value. |
| `top_partner_changes` | Product and partner | Largest absolute Atlas score changes relative to WT. |
| `training_overlap_sequences` | Unique product | Exact-hash presence and matching identifiers in the accessible training sequence universe. |

`results/tables/gnomad_population_control_candidates` contains gnomAD v4
missense SNVs, ancestry-aware maximum frequency, ABCA4 VCEP frequency tier,
ClinVar class, and existing Atlas-PPI sequence mapping. The label scope is
population-compatible control candidate, not clinical benign. The compressed
source response is `data/raw/gnomad_abca4_v4.json.gz`.

`results/atlas_ppi/provenance.json` records the immutable model revision, cache
compatibility check, source and reference hashes, row counts, audit limits, and
the raw score matrix format, shape, indexes, and hash. Query embeddings are not
persisted.
