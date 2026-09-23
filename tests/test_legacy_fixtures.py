from collections import Counter

import pysam
import pytest

from osteosarc import legacy_fixtures as legacy
from osteosarc.cohort_bundle import record_digest, retrieval_regions, select_records


def test_historical_sam_selection_and_duplicate_references(bam):
    with pysam.AlignmentFile(bam) as inp:
        records = list(inp)
    required = Counter(legacy.sam_digest(r.to_string()) for r in records)
    selected = legacy.select_records(records, required)
    assert len(selected) == len(required)
    with pytest.raises(ValueError, match="missing"):
        legacy.select_records(records[1:], required)
    assert legacy.required_counts(dict(fixtures={"a": dict(source="s", records=["r", "r"]),
           "b": dict(source="s", records=["r"] * 10, record_references=True)}), "s") == {"r": 2}


def test_historical_segment_and_qname_policies_are_explicit(bam, tmp_path):
    windows = [("chr1", 100, 120), ("chr1", 120, 150)]
    with pysam.AlignmentFile(bam) as inp:
        kept, n, _ = legacy.select_window_segments(inp, windows)
        once, n_once, _ = legacy.select_window_segments(inp, windows, duplicate_policy="identical-SAM-lines-once")
    assert len(kept) == len(once) + 1
    assert n == n_once
    assert legacy.select_assigned_names({"alt": {"a": "2", "b": "1"}}, ["rare"], limit=1) == {"rare", "b"}
    rows = legacy.write_selected_names(bam, tmp_path / "selected.bam", ("chr1", 100, 150), {"repeated"})
    assert len(rows) == 3


def test_cohort_float_identity_and_exact_multiplicity(bam):
    from pathlib import Path

    from osteosarc import ReadSubset
    with pysam.AlignmentFile(bam) as inp:
        records = list(inp)
    assert retrieval_regions(records, "GRCh38")
    cohort = dict(path="all.bam", format="bam", records=[record_digest(r) for r in records])
    selected = select_records(ReadSubset(bam, Path(str(bam) + ".bai"), {}), cohort)
    assert [record_digest(r) for r in selected] == cohort["records"]
