import shutil
from pathlib import Path

import pytest

from osteosarc import SNAPSHOT_SOURCES, TIMELINE_SOURCES, Cache, Dataset

DATA = Path(__file__).parent / "data"
FILES = {"bams": "bams.json", "bucket": "bucket_listing.json", "variant_index": "variants.html",
         "vafs": "vafs.tsv", "vaf_columns": "vafs-columns.tsv", "vaccine_overlap": "vaccine_overlap.json",
         "source_variants": "source-variants.json", "bam_metadata": "bam-metadata.tsv", "data_page": "data.html",
         "events": "events.json", "events_sheet": "timeline.csv", "mrd": "mrd.json",
         "specimens": "samples-consolidated.tsv", "timepoint_summary": "samples.json",
         "fastqs": "fastqs-consolidated.tsv", "flow": "flow-manifest.json", "imaging": "dicom-studies.json",
         "pathology": "pathology-slides.json", "labs": "lab_results.tsv", "cytometry": "cytometry.tsv"}


@pytest.fixture
def dataset(tmp_path):
    cache = Cache(tmp_path / "cache", offline=True)
    for key, name in FILES.items():
        cache.import_file(DATA / name, {**SNAPSHOT_SOURCES, **TIMELINE_SOURCES}[key])
    return Dataset.sync("fixture", cache=cache)


@pytest.fixture
def bam(tmp_path):
    pysam = pytest.importorskip("pysam")
    if shutil.which("samtools") is None:
        pytest.skip("samtools required for indexed extraction integration tests")
    path = tmp_path / "input.bam"
    header = dict(HD={"VN": "1.6", "SO": "coordinate"},
                  SQ=[dict(SN="chr1", LN=248956422), dict(SN="chr2", LN=242193529)],
                  RG=[dict(ID="rg1", SM="sample"), dict(ID="rg2", SM="another-library")])
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        # Exact repeat records and secondary/supplementary/duplicate-marked reads
        # are deliberately distinct from duplicate emission across query intervals.
        for name, pos, flag, mapq, rg, barcode in [
            ("repeated", 100, 0, 60, "rg1", "A"), ("repeated", 100, 0, 60, "rg1", "A"),
            ("repeated", 100, 0, 60, "rg2", "B"), ("secondary", 110, 256, 20, "rg1", "A"),
            ("duplicate", 115, 1024, 10, "rg1", "B"), ("supplementary", 120, 2048, 60, "rg1", "A"),
            ("outside", 1000, 0, 60, "rg1", "A"),
        ]:
            read = pysam.AlignedSegment(out.header)
            read.query_name = name
            read.query_sequence = "ACGT" * 10
            read.query_qualities = pysam.qualitystring_to_array("I" * 40)
            read.reference_id = 0
            read.reference_start = pos
            read.flag = flag
            read.mapping_quality = mapq
            read.cigarstring = "40M"
            read.set_tag("RG", rg)
            read.set_tag("CB", barcode)
            read.set_tag("UB", "original-umi")
            read.set_tag("xf", 1.23456789, value_type="f")
            out.write(read)
    pysam.index(str(path))
    return path
