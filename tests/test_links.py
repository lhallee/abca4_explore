from __future__ import annotations

import urllib.parse

from abca4_avi.links import build_variant_url


def test_atlas_url_contains_variant_and_splicing_modalities() -> None:
    url = build_variant_url("1:94023381:c>a")
    parsed = urllib.parse.urlparse(url)
    query = urllib.parse.parse_qs(parsed.query)

    assert query["q"] == ["chr1:94023381:C>A"]
    assert "section:RNA_SEQ" in query["lItems"][0]
    assert "section:SPLICE_JUNCTIONS" in query["lItems"][0]
    assert "section:SPLICE_SITE_USAGE" in query["lItems"][0]
    assert "section:SPLICE_SITES" in query["lItems"][0]
