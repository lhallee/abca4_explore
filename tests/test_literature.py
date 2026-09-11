from __future__ import annotations

from abca4_avi.literature import SANGERMANO_CORRECT_MRNA, WANG_DEEP_INTRONIC, _mean
from abca4_avi.pipeline import _coding_key, _protein_key


def test_published_numeric_mean_parsing() -> None:
    assert _mean("34,3#") == 34.3
    assert _mean("−0.02") == -0.02
    assert _mean("102 ± 15") == 102.0


def test_hgvs_keys_match_one_and_three_letter_protein_notation() -> None:
    assert _coding_key("NM_000350.3:c.317A>T") == "317A>T"
    assert _protein_key("p.R653C") == "R653C"
    assert _protein_key("NP_000341.2:p.(Arg653Cys)") == "R653C"


def test_published_splicing_subsets_have_expected_sizes() -> None:
    assert len(SANGERMANO_CORRECT_MRNA) == 47
    assert len(WANG_DEEP_INTRONIC) == 7
