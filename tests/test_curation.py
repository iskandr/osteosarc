"""Central corrections: optional, re-checked against sources, and never partial."""

import json
import re
from pathlib import Path

import pytest

from osteosarc import Dataset
from osteosarc.cli import main
from osteosarc.curation import (
    CORRECTIONS,
    Change,
    Correction,
    CurationWarning,
    glob,
    normalize_tissue,
)

SMC5 = "SMC5-chr9-70298024"


def reopen(dataset, corrections):
    return Dataset.open("fixture", cache=dataset.cache, corrections=corrections)


def status(data, correction_id):
    return next(r for r in data.corrections if r["id"] == correction_id)["status"]


def test_edit_and_flag_are_applied_and_traceable(dataset):
    rename = Correction("rename", "fixture rename", (
        Change("variant_index", {"id": SMC5}, expect={"gene": "SMC5"}, set={"gene": "SMC5-renamed"}),))
    flag = Correction("flag", "fixture flag", (Change("vafs", {"variant_id": SMC5}),))
    data = reopen(dataset, [rename, flag])
    variant = data.variants()[SMC5]
    assert variant.gene == "SMC5-renamed"
    assert variant.annotations["corrections"] == ("rename", "flag")
    rows = data.vafs.select(variant_id=SMC5)
    assert rows.rows and all(r["corrections"] == "flag" for r in rows)
    assert data.vafs.columns[-1] == "corrections"
    assert {status(data, "rename"), status(data, "flag")} == {"applied"}
    assert data.variants().source["corrections"] == ("rename", "flag")
    # The published sources are untouched and remain available.
    raw = reopen(dataset, False)
    assert raw.variants()[SMC5].gene == "SMC5"
    assert "corrections" not in raw.vafs.columns
    assert status(raw, "viewer-label-BG009368") == "disabled"


def test_upstream_fix_is_recognized(dataset):
    already = Correction("already", "source already agrees", (
        Change("variant_index", {"id": SMC5}, expect={"gene": "OLD"}, set={"gene": "SMC5"}),))
    data = reopen(dataset, [already])
    assert status(data, "already") == "fixed_upstream"
    assert data.variants()[SMC5].annotations["corrections"] == ()


def test_changed_or_missing_source_is_stale_and_atomic(dataset):
    changed = Correction("changed", "expectation no longer holds", (
        Change("variant_index", {"id": SMC5}, set={"gene": "EDITED"}),
        Change("vafs", {"variant_id": SMC5}, expect={"ref": "T"}, set={"alt": "C"})))
    missing = Correction("missing", "record removed upstream", (
        Change("variant_index", {"id": "GONE-chr1-1"}, set={"gene": "X"}),))
    data = reopen(dataset, [changed, missing])
    with pytest.warns(CurationWarning) as caught:
        variants = data.variants()
    assert sorted(str(w.message).split("'")[1] for w in caught) == ["changed", "missing"]
    # Neither change is made when one of them is stale.
    assert variants[SMC5].gene == "SMC5"
    assert variants[SMC5].allele == ("chr9", 70298024, "G", "A")
    report = {r["id"]: r for r in data.corrections}
    assert report["changed"]["status"] == report["missing"]["status"] == "stale"
    assert report["changed"]["changes"][1]["differing"] == ["ref"]


def test_asset_flags_use_bucket_globs(dataset):
    flag = Correction("bam-flag", "every BAM", (Change("bucket", {"key": glob("*.bam")}),))
    data = reopen(dataset, [flag])
    flagged = [a for a in data.assets if "bam-flag" in a.metadata.get("corrections", ())]
    assert flagged and all(a.key.endswith(".bam") for a in flagged)


def test_duplicate_ids_are_rejected(dataset):
    twice = Correction("same", "", (Change("vafs", {"variant_id": SMC5}),))
    with pytest.raises(ValueError, match="unique"):
        reopen(dataset, [twice, twice])


def test_unrecognized_labels_are_reported(dataset):
    assert not reopen(dataset, ()).unrecognized
    flag = Correction("new-assay", "simulate an upstream assay label",
                      (Change("bam_metadata", {"assay": "RNA"}, set={"assay": "Spatial"}),))
    unknown = reopen(dataset, [flag]).unrecognized
    assert [(r["source"], r["field"], r["value"]) for r in unknown] == [("bam_metadata", "assay", "Spatial")]


def test_cli_reports_and_strict_mode(dataset, capsys):
    root = str(dataset.cache.root)
    # The fixture excerpt lacks most corrected records, so they are stale here,
    # and a drift check fails even when this run does not apply corrections.
    with pytest.warns(CurationWarning):
        assert main(["--cache", root, "curation", "fixture", "--strict"]) == 1
    assert "stale corrections" in capsys.readouterr().err
    assert main(["--cache", root, "--no-corrections", "curation", "fixture", "--strict"]) == 1
    report = json.loads(capsys.readouterr().out)
    assert {r["status"] for r in report["corrections"]} == {"disabled"}
    assert main(["--cache", root, "--no-corrections", "curation", "fixture"]) == 0


def test_normal_and_blood_are_one_tissue():
    assert normalize_tissue("Normal") == normalize_tissue("Blood Normal") == "blood"
    assert normalize_tissue("Tumor + Normal") is None
    assert normalize_tissue("Bone marrow") == "bone marrow"  # kept and reported, not guessed


def test_registry_is_complete_and_documented():
    from osteosarc.catalog import TIMELINE_SOURCES
    sources = {"vafs", "bam_metadata", "source_variants", "vaccine_overlap", "variant_index", "bams",
               "bucket", *TIMELINE_SOURCES}
    docs = (Path(__file__).parents[1] / "docs" / "curation.md").read_text()
    assert len({c.id for c in CORRECTIONS}) == len(CORRECTIONS)
    for correction in CORRECTIONS:
        assert correction.summary and correction.evidence and correction.verified
        assert f"`{correction.id}`" in docs, correction.id
        for change in correction.changes:
            assert change.source in sources
            if change.source in ("source_variants", "vafs") and "ref" in change.set:
                assert re.fullmatch("[ACGT]+", change.set["ref"]) and re.fullmatch("[ACGT]+", change.set["alt"])


def test_cached_records_are_not_exposed_for_mutation(dataset):
    data = reopen(dataset, True)
    record = data.variants()[SMC5].annotations["source_record"]
    for row in data.annotations:
        row["detection"].clear()
    for row in data.vaccines:
        row["vaccines"] = {}
    assert data.annotations.select(id=SMC5).rows[0]["detection"]
    assert record["detection"] and data.pipeline_names
    assert any(row["vaccines"] for row in data.vaccines)


def test_source_revision_pins_every_repository_source():
    from osteosarc.catalog import SOURCE_REPO, TIMELINE_SOURCES
    from osteosarc.cli import pinned_sources
    revision = "0" * 40
    pinned = pinned_sources(revision)
    assert {"source_variants", "bam_metadata", "specimens", "events_sheet", "imaging"} <= set(pinned)
    assert all(f"/-/raw/{revision}/" in url for url in pinned.values())
    assert all(url.startswith(SOURCE_REPO) is False for url in pinned.values())
    assert set(pinned) == {k for k, v in {**TIMELINE_SOURCES}.items() if v.startswith(SOURCE_REPO)} | {
        "source_variants", "bam_metadata"}
    with pytest.raises(ValueError):
        pinned_sources("main")


def test_malformed_bucket_rows_are_schema_errors(tmp_path):
    from conftest import DATA, FILES

    from osteosarc import SNAPSHOT_SOURCES, TIMELINE_SOURCES, Cache, SchemaError
    cache = Cache(tmp_path / "cache", offline=True)
    urls = {**SNAPSHOT_SOURCES, **TIMELINE_SOURCES}
    for key, name in FILES.items():
        path = DATA / name
        if key == "bucket":
            listing = json.loads(path.read_text())
            listing["files"].append(["truncated-row"])
            path = tmp_path / name
            path.write_text(json.dumps(listing))
        cache.import_file(path, urls[key])
    with pytest.raises(SchemaError, match="Malformed bucket object"):
        Dataset.sync("broken", cache=cache)
