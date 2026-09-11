from __future__ import annotations

import xml.etree.ElementTree as ET

from pathlib import Path

from abca4_avi.clinvar import (
    extract_allele_location,
    parse_vcv_xml,
    save_summary_batches,
)


ALLELE_XML = b"""<?xml version="1.0"?>
<VariationArchive VariationID="4887205" Accession="VCV004790316"
  VariationType="single nucleotide variant" RecordStatus="current">
  <IncludedRecord>
    <GeneList><Gene><Location>
      <SequenceLocation Assembly="GRCh38" Chr="1" start="93992834"
        stop="94121148" />
    </Location></Gene></GeneList>
  </IncludedRecord>
  <ClassifiedRecord>
    <SimpleAllele>
      <Name>NP_000341.2:p.(Arg1514Gln)</Name>
      <ProteinChange>R1514Q</ProteinChange>
      <Location>
        <SequenceLocation Assembly="GRCh37" Chr="1" positionVCF="94556337"
          referenceAlleleVCF="C" alternateAlleleVCF="A" />
        <SequenceLocation Assembly="GRCh38" Chr="1" positionVCF="94023381"
          start="94023381" stop="94023381" referenceAlleleVCF="C"
          alternateAlleleVCF="A" Accession="NC_000001.11" forDisplay="true" />
      </Location>
      <HGVSlist><HGVS><Expression>NM_000350.3:c.4539+2001G&gt;A</Expression></HGVS></HGVSlist>
      <MolecularConsequence Type="intron variant" ID="SO:0001627" />
    </SimpleAllele>
    <Classifications><GermlineClassification>
      <Description>Pathogenic</Description>
      <ReviewStatus>criteria provided, single submitter</ReviewStatus>
    </GermlineClassification></Classifications>
  </ClassifiedRecord>
</VariationArchive>
"""


def test_gene_interval_is_never_selected_as_variant() -> None:
    location = extract_allele_location(ET.fromstring(ALLELE_XML))

    assert location is not None
    assert location.position == 94_023_381
    assert location.position != 93_992_834
    assert location.variant == "chr1:94023381:C>A"


def test_parser_preserves_grch38_allele_and_review_stars() -> None:
    record = parse_vcv_xml(ALLELE_XML)[0]

    assert record["variant"] == "chr1:94023381:C>A"
    assert record["coordinate_source"] == "ClinVar allele VCF"
    assert record["review_stars"] == 1
    assert record["hgvs_c"] == "NM_000350.3:c.4539+2001G>A"
    assert record["hgvs_p"] == "p.R1514Q"
    assert record["molecular_consequences"] == "intron variant"


def test_summary_batches_resume_independently(tmp_path: Path) -> None:
    class FakeClient:
        def __init__(self) -> None:
            self.calls: list[list[str]] = []

        def summary(self, variant_ids: list[str]) -> list[dict[str, object]]:
            self.calls.append(variant_ids)
            return [{"variant_id": value} for value in variant_ids]

    client = FakeClient()
    identifiers = ["1", "2", "3"]

    paths = save_summary_batches(client, identifiers, tmp_path, batch_size=2)
    resumed = save_summary_batches(client, identifiers, tmp_path, batch_size=2)

    assert paths == resumed
    assert client.calls == [["1", "2"], ["3"]]
