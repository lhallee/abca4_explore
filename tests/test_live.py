from __future__ import annotations

from pathlib import Path

import pytest

from abca4_avi.atlas import score_variant
from abca4_avi.clinvar import ClinVarClient, parse_vcv_xml
from abca4_avi.config import load_credentials


@pytest.mark.live
def test_one_live_clinvar_record() -> None:
    credentials = load_credentials(Path.cwd())
    client = ClinVarClient(api_key=credentials.ncbi_api_key)
    search = client.search("ABCA4[gene]", retmax=1)
    payload = client.fetch_vcv_batch(search["variant_ids"])

    assert len(parse_vcv_xml(payload)) == 1


@pytest.mark.live
def test_one_live_alphagenome_snv() -> None:
    from alphagenome.atlas import atlas

    credentials = load_credentials(Path.cwd())
    if not credentials.alphagenome_api_key:
        pytest.skip("AlphaGenome key unavailable")
    result = score_variant(
        atlas.create(credentials.alphagenome_api_key), "chr1:94023381:C>A"
    )

    assert result["avi_phred"] >= 0
    assert result["atlas_url"]
