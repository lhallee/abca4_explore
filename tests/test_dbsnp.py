from __future__ import annotations

import json

from pathlib import Path

from abca4_avi import dbsnp
from abca4_avi.dbsnp import fill_missing_coordinates, parse_refsnp_grch38


def test_dbsnp_selects_grch38_chromosome_snv() -> None:
    payload = {
        "primary_snapshot_data": {
            "placements_with_allele": [
                {
                    "placement_annot": {
                        "is_chromosome": True,
                        "seq_id_traits_by_assembly": [{"assembly_name": "GRCh38.p14"}],
                    },
                    "alleles": [
                        {
                            "allele": {
                                "spdi": {
                                    "seq_id": "NC_000001.11",
                                    "position": 94023380,
                                    "deleted_sequence": "C",
                                    "inserted_sequence": "C",
                                }
                            }
                        },
                        {
                            "allele": {
                                "spdi": {
                                    "seq_id": "NC_000001.11",
                                    "position": 94023380,
                                    "deleted_sequence": "C",
                                    "inserted_sequence": "A",
                                }
                            }
                        },
                    ],
                }
            ]
        }
    }

    result = parse_refsnp_grch38(payload)

    assert result is not None
    assert result["variant"] == "chr1:94023381:C>A"


def test_dbsnp_resolution_cache_prevents_duplicate_calls(
    tmp_path: Path, monkeypatch
) -> None:
    calls: list[str] = []

    def fake_resolve(rsid: str) -> dict[str, object]:
        calls.append(rsid)
        return {
            "chromosome": "1",
            "position": 10,
            "ref": "A",
            "alt": "G",
            "variant": "chr1:10:A>G",
        }

    monkeypatch.setattr(dbsnp, "resolve_rsid", fake_resolve)
    cache_path = tmp_path / "dbsnp.json"
    records = [{"variant": None, "rsid": "rs1"}]

    first = fill_missing_coordinates(records, cache_path=cache_path)
    second = fill_missing_coordinates(records, cache_path=cache_path)

    assert first == second
    assert calls == ["rs1"]
    assert json.loads(cache_path.read_text(encoding="utf-8"))["rs1"] is not None
