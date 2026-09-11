from __future__ import annotations

import pandas as pd

from abca4_avi.gene_model import Exon, GeneModel, extract_gene_model
from abca4_avi.variants import classify_variant_region


def _model() -> GeneModel:
    return GeneModel(
        gene_name="ABCA4",
        gene_id="ENSG1",
        transcript_id="ENST1",
        chromosome="chr1",
        strand="-",
        start=100,
        end=500,
        exons=(
            Exon(number=1, start=400, end=450),
            Exon(number=2, start=200, end=250),
        ),
    )


def test_minus_strand_exons_are_in_transcript_order() -> None:
    gtf = pd.DataFrame(
        [
            {
                "gene_name": "ABCA4",
                "gene_id": "ENSG1",
                "transcript_id": "ENST1",
                "Feature": "transcript",
                "Chromosome": "chr1",
                "Strand": "-",
                "Start": 100,
                "End": 500,
                "tag": "MANE_Select",
            },
            {
                "gene_name": "ABCA4",
                "gene_id": "ENSG1",
                "transcript_id": "ENST1",
                "Feature": "exon",
                "Chromosome": "chr1",
                "Strand": "-",
                "Start": 200,
                "End": 250,
                "tag": "MANE_Select",
                "exon_number": 2,
            },
            {
                "gene_name": "ABCA4",
                "gene_id": "ENSG1",
                "transcript_id": "ENST1",
                "Feature": "exon",
                "Chromosome": "chr1",
                "Strand": "-",
                "Start": 400,
                "End": 450,
                "tag": "MANE_Select",
                "exon_number": 1,
            },
        ]
    )

    model = extract_gene_model(gtf)

    assert model.strand == "-"
    assert [exon.start for exon in model.exons] == [400, 200]


def test_splice_distance_tiers() -> None:
    model = _model()

    assert (
        classify_variant_region(pd.Series({"position": 399}), model)
        == "canonical_splice"
    )
    assert (
        classify_variant_region(pd.Series({"position": 390}), model) == "splice_region"
    )
    assert (
        classify_variant_region(pd.Series({"position": 300}), model) == "deep_intronic"
    )
    assert (
        classify_variant_region(pd.Series({"position": 420}), model) == "other_exonic"
    )
    assert (
        classify_variant_region(
            pd.Series({"position": 420, "molecular_consequences": "missense variant"}),
            model,
        )
        == "coding_missense"
    )
    assert classify_variant_region(pd.Series({"position": pd.NA}), model) == "other"
