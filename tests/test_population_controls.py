from __future__ import annotations

from abca4_avi.population_controls import (
    ABCA4_BA1_AF,
    ABCA4_BS1_STRONG_AF,
    ABCA4_BS1_SUPPORTING_AF,
    frequency_evidence_tier,
    normalize_gnomad_missense,
)


def test_frequency_tiers_use_strict_vcep_boundaries() -> None:
    assert (
        frequency_evidence_tier(ABCA4_BS1_SUPPORTING_AF, excluded_by_vcep=False)
        == "below_supporting_threshold"
    )
    assert (
        frequency_evidence_tier(0.00164, excluded_by_vcep=False)
        == "bs1_supporting_candidate"
    )
    assert (
        frequency_evidence_tier(ABCA4_BS1_STRONG_AF, excluded_by_vcep=False)
        == "bs1_supporting_candidate"
    )
    assert (
        frequency_evidence_tier(0.0164, excluded_by_vcep=False)
        == "bs1_strong_candidate"
    )
    assert (
        frequency_evidence_tier(ABCA4_BA1_AF, excluded_by_vcep=False)
        == "bs1_strong_candidate"
    )
    assert (
        frequency_evidence_tier(0.164, excluded_by_vcep=False)
        == "ba1_candidate"
    )


def test_vcep_exclusion_overrides_frequency() -> None:
    assert (
        frequency_evidence_tier(0.5, excluded_by_vcep=True)
        == "excluded_vcep_reduced_penetrance_or_hypomorphic"
    )


def test_normalization_selects_missense_snvs_and_computes_population_max() -> None:
    response = {
        "data": {
            "gene": {
                "variants": [
                    {
                        "variant_id": "1-10-A-G",
                        "chrom": "1",
                        "pos": 10,
                        "ref": "A",
                        "alt": "G",
                        "rsids": ["rs1"],
                        "transcript_consequence": {
                            "consequence_terms": ["missense_variant"],
                            "hgvsc": "ENST1:c.100A>G",
                            "hgvsp": "ENSP1:p.Lys34Arg",
                            "is_mane_select": True,
                            "transcript_id": "ENST1",
                        },
                        "joint": {
                            "ac": 3,
                            "an": 1_000,
                            "homozygote_count": 0,
                            "populations": [
                                {"id": "afr", "ac": 1, "an": 1_000},
                                {"id": "nfe", "ac": 2, "an": 1_000},
                                {"id": "fin", "ac": 500, "an": 1_000},
                            ],
                        },
                    },
                    {
                        "variant_id": "1-11-A-T",
                        "chrom": "1",
                        "pos": 11,
                        "ref": "A",
                        "alt": "T",
                        "transcript_consequence": {
                            "consequence_terms": ["synonymous_variant"]
                        },
                    },
                ]
            }
        }
    }

    result = normalize_gnomad_missense(response)

    assert result["variant"].tolist() == ["chr1:10:A>G"]
    assert result.loc[0, "max_population_af"] == 0.002
    assert result.loc[0, "max_population"] == "nfe"
    assert result.loc[0, "population_evidence_tier"] == "bs1_supporting_candidate"


def test_known_hypomorphic_variant_is_flagged() -> None:
    response = {
        "data": {
            "gene": {
                "variants": [
                    {
                        "variant_id": "1-10-A-G",
                        "chrom": "1",
                        "pos": 10,
                        "ref": "A",
                        "alt": "G",
                        "transcript_consequence": {
                            "consequence_terms": ["missense_variant"],
                            "hgvsc": "ENST1:c.5603A>T",
                        },
                        "joint": {
                            "populations": [{"id": "nfe", "ac": 500, "an": 1_000}]
                        },
                    }
                ]
            }
        }
    }

    result = normalize_gnomad_missense(response)

    assert bool(result.loc[0, "vcep_frequency_exclusion"])
    assert result.loc[0, "population_evidence_tier"].startswith("excluded_vcep")
