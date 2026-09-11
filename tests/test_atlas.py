from __future__ import annotations

import json

from pathlib import Path

import abca4_avi.atlas as atlas_module


def test_checkpoint_resume_skips_completed_variant(tmp_path: Path, monkeypatch) -> None:
    checkpoint = tmp_path / "scores.jsonl"
    checkpoint.write_text(
        json.dumps({"variant": "chr1:100:A>G", "avi_phred": 12.0}) + "\n",
        encoding="utf-8",
    )
    calls: list[str] = []
    monkeypatch.setattr(atlas_module.atlas, "create", lambda _: object())

    def fake_score(client: object, variant: str) -> dict[str, object]:
        del client
        calls.append(variant)
        return {"variant": variant, "avi_phred": 15.0}

    monkeypatch.setattr(atlas_module, "score_variant", fake_score)

    result = atlas_module.score_variants(
        ["chr1:100:A>G", "chr1:101:C>T"],
        api_key="hidden",
        checkpoint_path=checkpoint,
        failure_path=tmp_path / "failures.jsonl",
        requests_per_second=1_000,
    )

    assert calls == ["chr1:101:C>T"]
    assert set(result["variant"]) == {"chr1:100:A>G", "chr1:101:C>T"}


def test_fully_complete_checkpoint_does_not_create_client(
    tmp_path: Path, monkeypatch
) -> None:
    checkpoint = tmp_path / "scores.jsonl"
    checkpoint.write_text(
        json.dumps({"variant": "chr1:100:A>G", "avi_phred": 12.0}) + "\n",
        encoding="utf-8",
    )

    def unexpected_client(api_key: str) -> object:
        del api_key
        raise AssertionError("client should not be created")

    monkeypatch.setattr(atlas_module.atlas, "create", unexpected_client)

    result = atlas_module.score_variants(
        ["chr1:100:A>G"],
        api_key="hidden",
        checkpoint_path=checkpoint,
        failure_path=tmp_path / "failures.jsonl",
    )

    assert len(result) == 1
