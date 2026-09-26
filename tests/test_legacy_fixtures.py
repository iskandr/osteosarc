import json
import zipfile

import pysam
import pytest

from osteosarc import cohort_bundle
from osteosarc import legacy_fixtures as legacy
from osteosarc.cohort_bundle import record_digest, select_records


def test_historical_segment_policies_are_explicit(bam):
    windows = [("chr1", 100, 120), ("chr1", 120, 150)]
    with pysam.AlignmentFile(bam) as inp:
        kept, n, _ = legacy.select_window_segments(inp, windows)
        once, n_once, _ = legacy.select_window_segments(inp, windows, duplicate_policy="identical-SAM-lines-once")
    assert len(kept) == len(once) + 1
    assert n == n_once


def test_cohort_float_identity_and_exact_multiplicity(bam):
    from pathlib import Path

    from osteosarc import ReadSubset
    with pysam.AlignmentFile(bam) as inp:
        records = list(inp)
    cohort = dict(path="all.bam", format="bam", records=[record_digest(r) for r in records])
    selected = select_records(ReadSubset(bam, Path(str(bam) + ".bai"), {}), cohort)
    assert [record_digest(r) for r in selected] == cohort["records"]


def test_cohort_size_budget_includes_manifest_before_extraction(tmp_path, monkeypatch):
    archive, destination = tmp_path / "bundle.zip", tmp_path / "out"
    destination.mkdir()
    monkeypatch.setattr(cohort_bundle, "DEFAULT_SIZE_BUDGET", 1024)
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as out:
        out.writestr("bundle.json", json.dumps(dict(schema_version=1, files={}, padding="x" * 1200)))
    with pytest.raises(ValueError, match="size budget"):
        cohort_bundle.extract_bundle(archive, destination)
    assert not list(destination.iterdir())


def test_fusion_metadata_mismatch_does_not_silently_drop_records(tmp_path):
    cohort = dict(path="fusion.json.gz", format="fusion", record_metadata=[{}])
    with pytest.raises(ValueError, match="metadata count differs"):
        cohort_bundle.write_cohort(tmp_path, cohort, [None, None], tmp_path / "missing-recipe")
    assert not (tmp_path / "fusion.json.gz").exists()


def test_the_helpers_isovar_imports(bam, tmp_path):
    import hashlib
    with pysam.AlignmentFile(bam) as inp:
        header, records = inp.header.to_dict(), [r.to_string() for r in inp]
    assert legacy.sam_digest(records[0]) == hashlib.sha256(records[0].encode("ascii")).hexdigest()
    path = tmp_path / "x.json.gz"
    legacy.write_json(path, dict(a=1))
    assert legacy.read_json(path) == dict(a=1)
    minimal = legacy.minimal_header(header, records)
    assert [row["SN"] for row in minimal["SQ"]] == ["chr1"] and {row["ID"] for row in minimal["RG"]} == {"rg1", "rg2"}
    regions = legacy.sam_regions(["chr1:101-200"], "GRCh38")  # one-based, inclusive
    assert (regions[0].contig, regions[0].start, regions[0].end, regions[0].assembly) == ("chr1", 100, 200, "GRCh38")
    with pysam.AlignmentFile(bam) as inp:
        read = next(iter(inp))
    assert legacy.segment_key(read) == ("rg1", "repeated", 0)
