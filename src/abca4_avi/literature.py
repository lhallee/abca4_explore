"""Reproducible extraction of open-access ABCA4 functional assay tables."""

from __future__ import annotations

import re
import urllib.request
import xml.etree.ElementTree as ET

from pathlib import Path

import pandas as pd


EUROPE_PMC_XML = "https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
STUDIES = {
    "garces_2021": {
        "pmcid": "PMC7796138",
        "doi": "https://doi.org/10.3390/ijms22010185",
        "expected_variants": 38,
    },
    "sangermano_2018": {
        "pmcid": "PMC5749174",
        "doi": "https://doi.org/10.1101/gr.226621.117",
        "expected_variants": 47,
    },
    "aslaksen_2024": {
        "pmcid": "PMC11305421",
        "doi": "https://doi.org/10.1167/iovs.65.10.2",
        "expected_variants": 10,
    },
    "wang_2025": {
        "pmcid": "PMC11781324",
        "doi": "https://doi.org/10.1167/iovs.66.1.65",
        "expected_variants": 7,
    },
}
SANGERMANO_CORRECT_MRNA = (
    ("c.160+5G>C", 34.3, "Moderate"),
    ("c.161G>T", 0.0, "Severe"),
    ("c.302+4A>C", 0.0, "Severe"),
    ("c.303-3C>G", 0.0, "Severe"),
    ("c.768G>T", 0.0, "Severe"),
    ("c.859-9T>C", 75.7, "Mild"),
    ("c.1100-6T>A", 0.0, "Severe"),
    ("c.1937+13T>G", 14.0, "Severe"),
    ("c.2382+5G>C", 47.9, "Moderate"),
    ("c.2588G>C", 60.0, "Mild"),
    ("c.2919-10T>C", 61.1, "Moderate"),
    ("c.2919-6C>A", 79.6, "Mild"),
    ("c.3050+5G>A", 0.0, "Severe"),
    ("c.3522+5del", 53.0, "Moderate"),
    ("c.3607G>A", 10.9, "Severe"),
    ("c.3607+3A>T", 0.0, "Severe"),
    ("c.3812A>G", 0.0, "Severe"),
    ("c.3813G>C", 0.0, "Severe"),
    ("c.3862+3A>G", 53.4, "Moderate"),
    ("c.4128G>A", 0.0, "Severe"),
    ("c.4253+4C>T", 7.8, "Severe"),
    ("c.4253+5G>A", 0.0, "Severe"),
    ("c.4253+5G>T", 5.4, "Severe"),
    ("c.4538A>G", 0.0, "Severe"),
    ("c.4538A>C", 4.3, "Severe"),
    ("c.4634G>A", 100.0, "Benign"),
    ("c.4667G>C", 0.0, "Severe"),
    ("c.4773G>C", 0.0, "Severe"),
    ("c.4773+3A>G", 24.6, "Severe"),
    ("c.4773+5G>A", 29.1, "Severe"),
    ("c.5196+3_5196+6del", 0.0, "Severe"),
    ("c.5196+3_5196+8del", 100.0, "Benign"),
    ("c.5312+3A>T", 0.0, "Severe"),
    ("c.5313-3C>G", 0.0, "Severe"),
    ("c.5460+5G>A", 0.0, "Severe"),
    ("c.5461-10T>C", 0.0, "Severe"),
    ("c.5461-8T>G", 0.0, "Severe"),
    ("c.5584G>C", 0.0, "Severe"),
    ("c.5584+5G>A", 0.0, "Severe"),
    ("c.5584+6T>C", 0.0, "Severe"),
    ("c.5585-10T>C", 100.0, "Benign"),
    ("c.5714+5G>A", 39.8, "Moderate"),
    ("c.5836-3C>A", 0.0, "Severe"),
    ("c.5898+5del", 4.5, "Severe"),
    ("c.6478A>G", 55.2, "Moderate"),
    ("c.6479+4A>G", 0.0, "Severe"),
    ("c.6729+5_6729+19del", 0.0, "Severe"),
)
WANG_DEEP_INTRONIC = (
    "c.161-395G>A",
    "c.2919-826T>A",
    "c.4539+1100A>G",
    "c.4634+741A>G",
    "c.5461-1321A>G",
    "c.1937+188A>G",
    "c.3329-551C>G",
)


def _fetch_xml(pmcid: str, cache_dir: Path) -> ET.Element:
    path = cache_dir / f"{pmcid}.xml"
    if not path.exists():
        request = urllib.request.Request(
            EUROPE_PMC_XML.format(pmcid=pmcid),
            headers={"User-Agent": "abca4-alpha-genome/0.1"},
        )
        with urllib.request.urlopen(request, timeout=120) as response:
            payload = response.read()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(payload)
    return ET.fromstring(path.read_bytes())


def _table_rows(root: ET.Element, table_id: str) -> list[list[str]]:
    table = next(
        item
        for item in root.findall(".//table-wrap")
        if item.attrib.get("id") == table_id
    )
    return [
        [" ".join("".join(cell.itertext()).split()) for cell in list(row)]
        for row in table.findall(".//tr")
    ]


def _mean(value: str) -> float:
    normalized = value.replace("−", "-").replace(",", ".")
    match = re.search(r"[-+]?\d+(?:\.\d+)?", normalized)
    if not match:
        raise ValueError(f"No numeric assay value in {value!r}")
    return float(match.group())


def _protein_row(
    *,
    study_key: str,
    doi: str,
    protein_change: str,
    coding_change: str | None,
    effect_name: str,
    effect_value: float,
    source_location: str,
) -> dict[str, object]:
    return {
        "study_key": study_key,
        "doi": doi,
        "assay_group": f"protein_{effect_name}",
        "variant": None,
        "hgvs_c": coding_change,
        "hgvs_p": protein_change,
        "effect_name": effect_name,
        "effect_value": effect_value,
        "effect_unit": "reported relative percent",
        "direction_interpretation": "higher value indicates more retained function",
        "source_location": source_location,
    }


def _extract_garces(root: ET.Element) -> list[dict[str, object]]:
    study = STUDIES["garces_2021"]
    output: list[dict[str, object]] = []
    for table_id in ("ijms-22-00185-t001", "ijms-22-00185-t002"):
        rows = _table_rows(root, table_id)
        for row in rows[1:]:
            if not row or row[0] == "WT":
                continue
            protein_change = row[0] if row[0].startswith("p.") else f"p.{row[0]}"
            measurements = {
                "expression": row[2],
                "basal_atpase": row[3],
                "ret_pe_induced_atpase": row[4],
            }
            for effect_name, value in measurements.items():
                output.append(
                    _protein_row(
                        study_key="garces_2021",
                        doi=str(study["doi"]),
                        protein_change=protein_change,
                        coding_change=None,
                        effect_name=effect_name,
                        effect_value=_mean(value),
                        source_location=table_id,
                    )
                )
    return output


def _extract_aslaksen(root: ET.Element) -> list[dict[str, object]]:
    study = STUDIES["aslaksen_2024"]
    rows = _table_rows(root, "tbl3")
    output: list[dict[str, object]] = []
    for row in rows[1:]:
        if not row or row[0] == "WT":
            continue
        match = re.match(r"(c\.\S+)\s+(p\.\(\S+\))", row[0])
        if not match:
            continue
        measurements = {
            "expression": row[1],
            "basal_atpase": row[3],
            "ret_pe_induced_atpase": row[4],
            "f_index": row[5],
        }
        for effect_name, value in measurements.items():
            output.append(
                _protein_row(
                    study_key="aslaksen_2024",
                    doi=str(study["doi"]),
                    protein_change=match.group(2),
                    coding_change=match.group(1),
                    effect_name=effect_name,
                    effect_value=_mean(value),
                    source_location="tbl3",
                )
            )
    return output


def _extract_sangermano() -> list[dict[str, object]]:
    study = STUDIES["sangermano_2018"]
    return [
        {
            "study_key": "sangermano_2018",
            "doi": study["doi"],
            "assay_group": "splicing_correct_mrna",
            "variant": None,
            "hgvs_c": coding_change,
            "hgvs_p": None,
            "effect_name": "correctly_spliced_mrna",
            "effect_value": correct_mrna,
            "effect_unit": "percent",
            "direction_interpretation": (
                "higher value indicates more correctly spliced transcript"
            ),
            "source_location": "Table 1",
            "reported_class": reported_class,
        }
        for coding_change, correct_mrna, reported_class in SANGERMANO_CORRECT_MRNA
    ]


def _extract_wang() -> list[dict[str, object]]:
    study = STUDIES["wang_2025"]
    tested_here = {
        "c.161-395G>A",
        "c.4634+741A>G",
        "c.1937+188A>G",
        "c.3329-551C>G",
    }
    return [
        {
            "study_key": "wang_2025",
            "doi": study["doi"],
            "assay_group": "splicing_minigene",
            "variant": None,
            "hgvs_c": coding_change,
            "hgvs_p": None,
            "effect_name": "aberrant_splicing",
            "effect_value": None,
            "effect_unit": "qualitative",
            "direction_interpretation": (
                "study reports a splice effect; exact band fractions are "
                "not tabulated in machine-readable text"
            ),
            "source_location": (
                "Results and Figure 2" if coding_change in tested_here else "Figure 1B"
            ),
            "reported_class": (
                "minigene tested in this study; aberrant bands exceed 50%"
                if coding_change in tested_here
                else "previously reported splice effect"
            ),
        }
        for coding_change in WANG_DEEP_INTRONIC
    ]


def prepare_functional_assays(reference_dir: Path, raw_dir: Path) -> pd.DataFrame:
    """Extract comparable quantitative tables and record study coverage."""
    roots = {
        key: _fetch_xml(str(study["pmcid"]), raw_dir) for key, study in STUDIES.items()
    }
    rows = [
        *_extract_garces(roots["garces_2021"]),
        *_extract_sangermano(),
        *_extract_aslaksen(roots["aslaksen_2024"]),
        *_extract_wang(),
    ]
    frame = pd.DataFrame(rows)
    reference_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(reference_dir / "functional_assays.tsv", sep="\t", index=False)
    coverage = pd.DataFrame(
        [
            {
                "study_key": key,
                "doi": study["doi"],
                "expected_variant_count": study["expected_variants"],
                "assay_rows_extracted": int((frame["study_key"] == key).sum()),
                "unique_variants_extracted": int(
                    frame.loc[frame["study_key"] == key, ["hgvs_c", "hgvs_p"]]
                    .drop_duplicates()
                    .shape[0]
                ),
                "status": (
                    "quantitative_table_extracted"
                    if key in {"garces_2021", "sangermano_2018", "aslaksen_2024"}
                    else (
                        "all variants represented; exact quantitative band "
                        "fractions require figure digitization"
                    )
                ),
            }
            for key, study in STUDIES.items()
        ]
    )
    coverage.to_csv(
        reference_dir / "functional_study_coverage.tsv", sep="\t", index=False
    )
    return frame
