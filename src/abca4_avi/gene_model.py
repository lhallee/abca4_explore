"""AlphaGenome GENCODE v46 gene-model caching."""

from __future__ import annotations

import json
import os
import urllib.request

from dataclasses import dataclass
from pathlib import Path
from collections.abc import Sequence

import pandas as pd


GTF_URL = (
    "https://storage.googleapis.com/alphagenome/reference/gencode/"
    "hg38/gencode.v46.annotation.gtf.gz.feather"
)


@dataclass(frozen=True, kw_only=True)
class Exon:
    """One 1-based closed exon interval."""

    number: int
    start: int
    end: int


@dataclass(frozen=True, kw_only=True)
class GeneModel:
    """A MANE Select transcript model."""

    gene_name: str
    gene_id: str
    transcript_id: str
    chromosome: str
    strand: str
    start: int
    end: int
    exons: tuple[Exon, ...]


def _download_file(url: str, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = path.with_suffix(path.suffix + ".part")
    urllib.request.urlretrieve(url, temporary_path)
    os.replace(temporary_path, path)


def _serialize_model(model: GeneModel, path: Path) -> None:
    payload = {
        "gene_name": model.gene_name,
        "gene_id": model.gene_id,
        "transcript_id": model.transcript_id,
        "chromosome": model.chromosome,
        "strand": model.strand,
        "start": model.start,
        "end": model.end,
        "exons": [
            {"number": exon.number, "start": exon.start, "end": exon.end}
            for exon in model.exons
        ],
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _deserialize_model(path: Path) -> GeneModel:
    payload = json.loads(path.read_text(encoding="utf-8"))
    exons = tuple(Exon(**exon) for exon in payload.pop("exons"))
    return GeneModel(**payload, exons=exons)


def extract_gene_model(gtf: pd.DataFrame, gene_name: str = "ABCA4") -> GeneModel:
    """Extract the MANE Select transcript for one gene."""
    gene_rows = gtf[gtf["gene_name"].fillna("").str.upper() == gene_name.upper()]
    if "tag" in gene_rows.columns:
        mane_rows = gene_rows[gene_rows["tag"].fillna("").str.contains("MANE_Select")]
        if not mane_rows.empty:
            gene_rows = mane_rows
    transcript_rows = gene_rows[gene_rows["Feature"] == "transcript"]
    if transcript_rows.empty:
        raise ValueError(f"No transcript found for {gene_name}")
    transcript = transcript_rows.iloc[0]
    transcript_id = str(transcript["transcript_id"])
    transcript_data = gene_rows[gene_rows["transcript_id"] == transcript_id]
    exon_rows = transcript_data[transcript_data["Feature"] == "exon"].copy()
    if exon_rows.empty:
        raise ValueError(f"No exons found for {transcript_id}")
    strand = str(transcript["Strand"])
    exon_rows = exon_rows.sort_values("Start", ascending=strand == "+")
    exons = tuple(
        Exon(
            number=int(row.get("exon_number", index + 1)),
            start=int(row["Start"]),
            end=int(row["End"]),
        )
        for index, (_, row) in enumerate(exon_rows.iterrows())
    )
    return GeneModel(
        gene_name=str(transcript["gene_name"]),
        gene_id=str(transcript["gene_id"]),
        transcript_id=transcript_id,
        chromosome=str(transcript["Chromosome"]),
        strand=strand,
        start=int(transcript["Start"]),
        end=int(transcript["End"]),
        exons=exons,
    )


def load_abca4_gene_model(cache_dir: Path) -> GeneModel:
    """Load a cached gene model or derive it from AlphaGenome GENCODE v46."""
    model_path = cache_dir / "abca4_gencode_v46_mane.json"
    if model_path.exists():
        return _deserialize_model(model_path)

    feather_path = cache_dir / "gencode.v46.annotation.gtf.gz.feather"
    if not feather_path.exists():
        _download_file(GTF_URL, feather_path)
    gtf = pd.read_feather(feather_path)
    model = extract_gene_model(gtf)
    _serialize_model(model, model_path)
    return model


def nearest_exon_distance(position: int, exons: Sequence[Exon]) -> int:
    """Return distance from an intronic position to the nearest exon edge."""
    distances: list[int] = []
    for exon in exons:
        if exon.start <= position <= exon.end:
            return 0
        if position < exon.start:
            distances.append(exon.start - position)
        else:
            distances.append(position - exon.end)
    if not distances:
        raise ValueError("gene model has no exons")
    return min(distances)
