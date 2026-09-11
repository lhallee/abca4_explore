# Methods

## Scope

The pipeline evaluates precise, biallelic GRCh38 SNVs associated with `ABCA4`.
Every retrieved ClinVar record is retained for descriptive counts. Variants that
cannot be scored are recorded with an explicit exclusion reason.

## ClinVar retrieval and normalization

The ClinVar query is `ABCA4[gene]`. E-utilities responses are requested in
batches of 200 and checkpointed after each batch. The client limits requests to
8 per second with an NCBI API key and 2.5 per second without one, and retries
transient HTTP 429 and 5xx responses with bounded exponential backoff.

Coordinate parsing selects allele-level GRCh38 `SequenceLocation` elements with
VCF attributes. Gene-level intervals are rejected. dbSNP GRCh38 placement is used
only when ClinVar lacks a precise allele. Indels, CNVs, haplotypes, missing
alleles, conflicting records, and unresolved coordinates remain in the
descriptive tables but are excluded from AVI scoring.

The GENCODE v46 MANE Select model for the negative-strand `ABCA4` transcript is
cached locally. Variants are classified as coding missense, canonical splice
positions 1 to 2, splice-region positions 3 to 20, deep intronic positions over
20 bases from an exon boundary, or other.

## AlphaGenome scoring

Each unique SNV is queried once at a default limit of five requests per second.
Successful responses are checkpointed immediately. The final score table stores
the AVI raw score, tail probability, Phred score, top percentile, top modality,
all 18 signed attribution weights, and an AlphaGenome Atlas URL. Detailed track
attributions are requested only for the four required sentinel variants and the
20 highest-scoring candidates.

The required sentinel variants are `c.4539+2001G>A`, `c.5461-10T>C`,
`c.4539+2028C>T`, and `c.769-784C>T`.

## Clinical cohorts

The primary cohort contains precise SNVs classified as pathogenic or likely
pathogenic versus benign or likely benign and reviewed at one ClinVar star or
higher. VUS and conflicting classifications are retained for descriptive plots
but excluded from binary metrics.

Sensitivity analyses use all review statuses and an STGD1 subset identified by
MedGen `C1855465`, MONDO `0009549`, OMIM `248200`, Orphanet `827`, or matching
synonyms. Coding missense, canonical splice, splice-region, and deep-intronic
strata are evaluated separately. ROC-AUC is omitted when a stratum contains fewer
than 10 positive or 10 negative variants.

Reported metrics are ROC-AUC, average precision, trapezoidal PR-AUC,
point-biserial association, rank-biserial effect size, and performance at Phred
20 and 30. Confidence intervals use 2,000 fixed-seed stratified bootstrap samples.

## Published functional assays

Functional measurements are analyzed within each study because assay scales are
not comparable across publications. Protein expression and ATPase measurements
are compared with full AVI and AlphaMissense attribution. Correctly spliced mRNA
is compared with full AVI and `MERGED_SPLICING` attribution.

The source studies are:

- Garces et al. 2021, [DOI 10.3390/ijms22010185](https://doi.org/10.3390/ijms22010185)
- Sangermano et al. 2018, [DOI 10.1101/gr.226621.117](https://doi.org/10.1101/gr.226621.117)
- Aslaksen et al. 2024, [DOI 10.1167/iovs.65.10.2](https://doi.org/10.1167/iovs.65.10.2)
- Wang et al. 2025, [DOI 10.1167/iovs.66.1.65](https://doi.org/10.1167/iovs.66.1.65)

## Interpretation

AlphaGenome and AVI produce molecular-impact predictions. ClinVar classifications
are categorical clinical evidence, not quantitative functional measurements.
The current analysis does not establish disease causality, penetrance, prognosis,
or clinical utility.

AVI used gnomAD frequency-derived proxy labels and included chromosome 1 in its
training partition. The public API does not expose the sampled training
coordinates, so this repository cannot reproduce the official position-overlap
filter. The reported clinical discrimination is therefore retrospective.
