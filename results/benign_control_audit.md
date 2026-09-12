# ABCA4 benign-control audit

## Outcome

The current ClinVar retrieval already exhausts the 4,670 records matching
`ABCA4[gene]`. It contains 1,429 benign or likely benign precise SNVs, but these
do not become 1,429 independent Atlas-PPI controls. After deterministic protein
construction and truncation, 710 map to the WT input, 692 are not protein
scorable, and only 27 distinct non-WT B/LB protein products remain in the
primary Atlas-PPI cohort.

The imbalance is therefore a protein-sequence problem, not a shortage of
benign genomic records. Counting synonymous or WT-equivalent variants as
independent observations would be pseudoreplication.

## Population-frequency candidates

A compact [gnomAD v4](https://gnomad.broadinstitute.org/) gene query returned
3,262 ABCA4 missense SNVs. The audit applied the
[ClinGen ABCA4 VCEP specifications](https://cspec.genome.network/cspec/ui/svi/doc/GN164?version=1.0):

- group-maximum allele frequency greater than 0.0163 for a strong-frequency
  proxy;
- group-maximum allele frequency greater than 0.00163 for a
  supporting-frequency proxy;
- gnomAD v4 groupmax ancestry exclusions; and
- removal of the seven variants for which the VCEP disallows BS1 because of
  hypomorphic or reduced-penetrance behavior.

This produced 54 population-compatible missense candidates. Fifty-two already
map to sequences in the completed Atlas-PPI screen. The source labels are 13
B/LB, 36 conflicting, two VUS, one P/LP, and two absent from the ClinVar
snapshot. No ClinVar label was changed.

The two unscreened candidates are:

| Variant | Protein | Groupmax AF | Group | Tier |
|---|---|---:|---|---|
| chr1:94045779:A>T | p.Asp961Glu | 0.059990 | AFR | Strong-frequency proxy |
| chr1:94045879:T>C | p.Gln928Arg | 0.014246 | EAS | Supporting-frequency proxy |

These two sequences were not embedded because the existing 52 candidates are
enough to evaluate whether the broader control definition materially changes
the result.

## Sensitivity result

At the unique-sequence level:

| Negative cohort | Negative sequences | Best exploratory ROC-AUC | Threshold | ROC-AUC at 0.338077128 | ROC-AUC at 0.9 |
|---|---:|---:|---:|---:|---:|
| ClinVar B/LB | 27 | 0.735 | 0.44 | 0.704 | 0.665 |
| B/LB plus strong-frequency proxy | 31 | 0.702 | 0.88 | 0.660 | 0.650 |
| B/LB plus supporting-frequency proxy | 67 | 0.676 | 0.62 | 0.632 | 0.626 |

![Control-cohort sensitivity](../plots/atlas_ppi/control_cohort_sensitivity.png)

**Figure.** Atlas-PPI clinical discrimination across the 0.33 to 0.99 threshold
grid for the original ClinVar negative set and two nested population-frequency
proxy sets. The vertical line is the published Atlas-PPI operating threshold.

Adding more weakly labeled controls narrows the class imbalance but does not
improve discrimination. This is a useful robustness result: the original
ClinVar-only estimate is sensitive to negative-cohort composition, and the
larger proxy cohorts should not be described as stronger clinical validation.

## Functional negatives

The Aslaksen et al. assay panel contains five normal or mild missense products:
p.Tyr106Phe, p.Gly172Ser, p.Pro940Arg, p.Lys1164Arg, and p.Val2050Leu. They are
valuable as a separate molecular-function endpoint. They are not pooled into
the clinical benign cohort because the ABCA4 VCEP states that currently
available functional assays cannot dependably rule out pathogenicity.

## Recommended interpretation

Retain the 27-sequence ClinVar cohort as the primary clinical analysis. Report
the 31- and 67-sequence population-proxy cohorts as sensitivity analyses. The
best route to genuinely stronger validation is more expert-curated B/LB
missense products or blinded functional WT-like variants, not synonymous
duplicates or frequency-only relabeling.

The complete candidate table is
[`gnomad_population_control_candidates.tsv`](tables/gnomad_population_control_candidates.tsv),
with a machine-readable count summary in
[`gnomad_population_control_candidates_summary.json`](tables/gnomad_population_control_candidates_summary.json).
