# ABCA4 AlphaGenome AVI results

AlphaGenome outputs are research-only molecular predictions, not clinical diagnoses.

Retrieved ClinVar records: 4,670.
Unique scored GRCh38 SNVs: 4,024.
Primary cohort: 2,521 variants (1,092 P/LP, 1,429 B/LB).
ROC-AUC: 0.990.
Average precision: 0.991.
Trapezoidal PR-AUC: 0.991.

## Clinical cohort performance

| Cohort | n | P/LP | B/LB | ROC-AUC | Average precision |
|---|---:|---:|---:|---:|---:|
| primary at least one star | 2,521 | 1,092 | 1,429 | 0.990 | 0.991 |
| all review statuses | 2,521 | 1,092 | 1,429 | 0.990 | 0.991 |
| stgd1 at least one star | 442 | 411 | 31 | 0.967 | 0.997 |
| primary coding missense | 690 | 659 | 31 | 0.966 | 0.998 |
| primary canonical splice | 119 | 119 | 0 | NA | NA |
| primary splice region | 499 | 30 | 469 | 0.953 | 0.878 |
| primary deep intronic | 236 | 21 | 215 | 0.870 | 0.603 |

## Sentinel variants

| Variant | GRCh38 allele | AVI Phred | Top modality | Atlas |
|---|---|---:|---|---|
| c.4539+2001G>A | chr1:94027444:C>T | 9.802 | MERGED_SPLICING | [open](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:94027444:C%3ET&m=variant&lItems=avi,section:RNA_SEQ,section:SPLICE_JUNCTIONS,section:SPLICE_SITE_USAGE,section:SPLICE_SITES) |
| c.4539+2028C>T | chr1:94027417:G>A | 5.304 | MERGED_SPLICING | [open](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:94027417:G%3EA&m=variant&lItems=avi,section:RNA_SEQ,section:SPLICE_JUNCTIONS,section:SPLICE_SITE_USAGE,section:SPLICE_SITES) |
| c.5461-10T>C | chr1:94011395:A>G | 21.550 | CACTUS_241_WAY | [open](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:94011395:A%3EG&m=variant&lItems=avi,section:RNA_SEQ,section:SPLICE_JUNCTIONS,section:SPLICE_SITE_USAGE,section:SPLICE_SITES) |
| c.769-784C>T | chr1:94084225:G>A | 8.694 | MERGED_SPLICING | [open](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:94084225:G%3EA&m=variant&lItems=avi,section:RNA_SEQ,section:SPLICE_JUNCTIONS,section:SPLICE_SITE_USAGE,section:SPLICE_SITES) |

## Published functional assay coverage

Assay scales are analyzed within each study and are not pooled.

| Study | Assay rows | Variants represented | Rows with AVI |
|---|---:|---:|---:|
| aslaksen_2024 | 40 | 10 | 32 |
| garces_2021 | 114 | 38 | 105 |
| sangermano_2018 | 47 | 47 | 32 |
| wang_2025 | 7 | 7 | 2 |

## Acceptance checks

- both_primary_classes_present: True
- all_18_attribution_columns_present: True
- all_four_sentinels_scored: True
- every_score_has_atlas_url: True
- no_failed_scores: True

See `metrics.json` for fixed-seed 95% bootstrap intervals and threshold metrics.
