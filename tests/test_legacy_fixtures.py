import gzip
import json
import zipfile
from collections import Counter

import pysam
import pytest

from osteosarc import cohort_bundle
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


def test_cohort_size_budget_includes_manifest_before_extraction(tmp_path, monkeypatch):
    archive, destination = tmp_path / "bundle.zip", tmp_path / "out"
    destination.mkdir()
    monkeypatch.setattr(cohort_bundle, "DEFAULT_SIZE_BUDGET", 1024)
    with zipfile.ZipFile(archive, "w", compression=zipfile.ZIP_DEFLATED) as out:
        out.writestr("bundle.json", json.dumps(dict(schema_version=1, files={}, padding="x" * 1200)))
    with pytest.raises(ValueError, match="size budget"):
        cohort_bundle.extract_bundle(archive, destination)
    assert not list(destination.iterdir())


@pytest.mark.parametrize("symlink", [False, True])
def test_cohort_existing_output_refused_before_acquisition(tmp_path, symlink):
    output = tmp_path / "existing.zip"
    if symlink:
        output.symlink_to(tmp_path / "missing-target")
    else:
        output.write_text("preserved")
    with pytest.raises(FileExistsError):
        cohort_bundle.generate_cohort_bundle(tmp_path / "missing-recipe", None, output)
    assert output.is_symlink() if symlink else output.read_text() == "preserved"


def test_cohort_publication_never_overwrites_concurrent_output(tmp_path, monkeypatch):
    recipe = tmp_path / "recipe"
    recipe.mkdir()
    (recipe / "support").mkdir()
    (recipe / "support/data.txt").write_text("original")
    (recipe / "selection.json.gz").write_bytes(gzip.compress(json.dumps(dict(snapshot_id="snapshot", cohorts=[])).encode()))
    (recipe / "catalog.json").write_text(json.dumps(dict(snapshot=dict(id="snapshot"), corrections=False, variants={})))
    # Exercise publication independently of the historical manifest format.
    monkeypatch.setattr(cohort_bundle, "update_manifests", lambda root: None)
    output = tmp_path / "output.zip"
    link = cohort_bundle.os.link

    def concurrent_writer(staged, target):
        target.write_text("another writer")
        return link(staged, target)

    monkeypatch.setattr(cohort_bundle.os, "link", concurrent_writer)
    with pytest.raises(FileExistsError):
        cohort_bundle.generate_cohort_bundle(recipe, None, output)
    assert output.read_text() == "another writer"


def test_fusion_metadata_mismatch_does_not_silently_drop_records(tmp_path):
    cohort = dict(path="fusion.json.gz", format="fusion", record_metadata=[{}])
    with pytest.raises(ValueError, match="metadata count differs"):
        cohort_bundle.write_cohort(tmp_path, cohort, [None, None], tmp_path / "missing-recipe")
    assert not (tmp_path / "fusion.json.gz").exists()
