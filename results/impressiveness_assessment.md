# How impressive is the ABCA4 AVI result?

The result is impressive as a retrospective ABCA4 prioritization result, especially for splice-region and deep-intronic variants. The headline 0.990 ROC-AUC is not yet an independent estimate of prospective performance.

## What AlphaGenome and AVI trained on

The base AlphaGenome model learned to predict experimental genome tracks from human and mouse DNA sequence. It was not trained on ABCA4 ClinVar classifications.

AVI is a separate 18-feature ensemble. According to the [Atlas methods paper](https://storage.googleapis.com/deepmind-media/DeepMind.com/Blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/alphagenome-atlas.pdf), it was trained on 37 million gnomAD v4.1 variants. Variants with filtering allele frequency below 0.001 were proxy-impactful; more common variants were proxy-benign. Chromosome 1 was one of the AVI training chromosomes.

Therefore AVI did not simply memorize ClinVar labels, but some ABCA4 positions may have appeared in its gnomAD-derived training data. The exact overlap cannot be measured from the public score API because the sampled 37-million-variant coordinate list is not returned.

The official Atlas ClinVar benchmark removed every evaluation position overlapping AVI training or validation. This ABCA4 analysis cannot apply that filter, so it should be called retrospective separation, not held-out generalization. Rarity is also partly circular because population frequency contributes to clinical classification.

## What remains genuinely persuasive

- The primary cohort contains 2,521 variants and achieves ROC-AUC 0.990 and average precision 0.991. The bootstrap intervals are narrow.
- Among 583 records last evaluated after the paper's June 15, 2025 ClinVar snapshot, ROC-AUC remains 0.982 and average precision remains 0.990. Last-evaluated date is not a first-submission date, so this is a sensitivity analysis, not a strict temporal holdout.
- The splicing attribution alone ranks deep-intronic P/LP versus B/LB variants at ROC-AUC 0.953 and average precision 0.850. In the splice-region subset it reaches ROC-AUC 0.991. This is the most biologically distinctive result because AlphaMissense and conservation carry little deep-intronic signal.
- Published ABCA4 experiments provide a smaller, more independent check. In 32 Sangermano variants, the splicing attribution correlates with correctly spliced mRNA at Spearman rho -0.429. In 35 Garces missense variants, full AVI correlations with retained protein expression and ATPase measures range from -0.377 to -0.470.

## What would make the claim publication-grade

Obtain the AVI train/validation coordinate list and exclude every overlapping ABCA4 position, then repeat the analysis in rare-only benign and pathogenic strata. A truly prospective ClinVar release or new blinded functional assays would test generalization without the rarity proxy.

## Sentinel variants

| Variant | AVI Phred | Leading attribution | Atlas |
|---|---:|---|---|
| c.4539+2001G>A | 9.802 | MERGED_SPLICING | [open](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:94027444:C%3ET&m=variant&lItems=avi,section:RNA_SEQ,section:SPLICE_JUNCTIONS,section:SPLICE_SITE_USAGE,section:SPLICE_SITES) |
| c.4539+2028C>T | 5.304 | MERGED_SPLICING | [open](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:94027417:G%3EA&m=variant&lItems=avi,section:RNA_SEQ,section:SPLICE_JUNCTIONS,section:SPLICE_SITE_USAGE,section:SPLICE_SITES) |
| c.5461-10T>C | 21.550 | CACTUS_241_WAY | [open](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:94011395:A%3EG&m=variant&lItems=avi,section:RNA_SEQ,section:SPLICE_JUNCTIONS,section:SPLICE_SITE_USAGE,section:SPLICE_SITES) |
| c.769-784C>T | 8.694 | MERGED_SPLICING | [open](https://deepmind.google.com/science/alphagenome/atlas?q=chr1:94084225:G%3EA&m=variant&lItems=avi,section:RNA_SEQ,section:SPLICE_JUNCTIONS,section:SPLICE_SITE_USAGE,section:SPLICE_SITES) |

AlphaGenome and AVI outputs are research-only molecular predictions, not clinical diagnoses.
