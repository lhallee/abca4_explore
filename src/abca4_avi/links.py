"""AlphaGenome Atlas URL construction."""

from __future__ import annotations

import urllib.parse


ATLAS_URL = "https://deepmind.google.com/science/alphagenome/atlas"
DEFAULT_MODALITIES = (
    "RNA_SEQ",
    "SPLICE_JUNCTIONS",
    "SPLICE_SITE_USAGE",
    "SPLICE_SITES",
)


def normalize_variant(variant: str) -> str:
    """Validate and normalize a 1-based chr:pos:ref>alt variant string."""
    chromosome, remainder = variant.strip().split(":", maxsplit=1)
    position_text, alleles = remainder.split(":", maxsplit=1)
    reference, alternate = alleles.split(">", maxsplit=1)
    position = int(position_text)
    if position < 1 or not reference or not alternate:
        raise ValueError(f"Invalid variant: {variant}")
    if not chromosome.lower().startswith("chr"):
        chromosome = f"chr{chromosome}"
    return f"{chromosome}:{position}:{reference.upper()}>{alternate.upper()}"


def build_variant_url(
    variant: str,
    *,
    biosample: str | None = None,
    modalities: tuple[str, ...] = DEFAULT_MODALITIES,
) -> str:
    """Build an Atlas variant link with RNA and splicing sections."""
    normalized = normalize_variant(variant)
    layout = ["avi", *(f"section:{name}" for name in modalities)]
    query: dict[str, str] = {
        "q": normalized,
        "m": "variant",
        "lItems": ",".join(layout),
    }
    if biosample:
        query["f"] = f"BIOSAMPLE_NAME:{biosample}"
    return f"{ATLAS_URL}?" + urllib.parse.urlencode(
        query,
        quote_via=urllib.parse.quote,
        safe=":,-",
    )
