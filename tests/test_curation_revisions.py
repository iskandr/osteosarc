"""Reviewed upstream edits must not disable historical or future drift checks."""

import copy
import json
import warnings
from pathlib import Path

import pytest

from osteosarc.curation import CORRECTIONS, Change, Correction, Curation, CurationWarning

DATA = json.loads((Path(__file__).parent / "data" / "curation_snapshots.json").read_text())
RELOCATED = {
    "CABLES1": ("CABLES1-chr18-23135500", "CABLES1-chr18-23135764"),
    "CCDC40": ("CCDC40-chr17-80058951", "CCDC40-chr17-80090148"),
    "DCHS2": ("DCHS2-chr4-154322488", "DCHS2-chr4-154323273"),
    "GAPVD1": ("GAPVD1-chr9-125299105", "GAPVD1-chr9-125301980"),
    "GOLGA6L2": ("GOLGA6L2-chr15-23441121", "GOLGA6L2-chr15-23440197"),
}
IDS = {"allele-" + old for old, _ in RELOCATED.values()} | {
    "tempus-grch37-counts", "transcript-DCHS2", "ush2a-transposed-duplicate"}
RULES = {c.id: c for c in CORRECTIONS if c.id in IDS}
COUNT_FIELDS = ("ref_reads", "alt_reads", "other_reads", "total_reads", "vaf")


def sources(version):
    return copy.deepcopy(DATA[version]["sources"])


def evaluate(raw, rule_ids=IDS, enabled=True):
    return Curation([RULES[key] for key in sorted(rule_ids)], raw.get, enabled=enabled)


def by_id(rows):
    return {r["id"]: r for r in rows}


@pytest.mark.parametrize("version", ["historical", "current"])
def test_reviewed_snapshots_have_no_stale_rules(version):
    raw = sources(version)
    before = copy.deepcopy(raw)
    curated = evaluate(raw)
    with warnings.catch_warnings():
        warnings.simplefilter("error", CurationWarning)
        report = {r["id"]: r["status"] for r in curated.report()}
        for source in raw:
            curated.records(source)
    assert set(report) == IDS
    assert report == {key: "fixed_upstream" if version == "current" and key == "tempus-grch37-counts"
                     else "applied" for key in IDS}
    assert raw == before


def test_relocated_ids_keep_verified_alleles_and_do_not_reuse_placeholder_counts():
    old, new = evaluate(sources("historical")), evaluate(sources("current"))
    old_records = by_id(old.records("source_variants")[0])
    new_records, marks = new.records("source_variants")
    current = by_id(new_records)
    for gene, (old_id, new_id) in RELOCATED.items():
        assert old_id not in current
        for field in ("chr", "pos", "ref", "alt", "variant_type", "genomic_location"):
            assert current[new_id][field] == old_records[old_id][field]
        index = next(i for i, r in enumerate(new_records) if r["id"] == new_id)
        assert "allele-" + old_id in marks[index]
        rows = [r for r in new.records("vafs")[0] if r["variant_id"] == new_id]
        assert len(rows) == 2
        assert all(row[field] == "" for row in rows for field in COUNT_FIELDS)
        assert all((row["ref"], row["alt"]) == (current[new_id]["ref"], current[new_id]["alt"])
                   for row in rows)
    assert (current[RELOCATED["CABLES1"][1]]["ref"], current[RELOCATED["CABLES1"][1]]["alt"]) == ("T", "TGGCGGC")
    assert current[RELOCATED["DCHS2"][1]]["refseq_id"] == "NM_001142552.1"
    # Published zeroes are not evidence against the corrected literal alleles.
    assert all(row["total_reads"] == "0" for row in DATA["current"]["sources"]["vafs"]
               if row["gene"] in RELOCATED)


def test_historical_tempus_counts_are_cleared_and_other_counts_are_preserved():
    raw = sources("historical")
    corrected = evaluate(raw, ["tempus-grch37-counts"]).records("vafs")[0]
    for before, after in zip(raw["vafs"], corrected):
        if "5GQLV9WSXQ" in before["bam_file"]:
            assert all(after[field] == "" for field in COUNT_FIELDS)
        else:
            assert after == before
    raw = sources("current")
    assert evaluate(raw, ["tempus-grch37-counts"]).records("vafs")[0] == raw["vafs"]


@pytest.mark.parametrize("drift", ["empty_vafs", "empty_metadata", "restored_viewer", "renamed_bam"])
def test_tempus_removal_requires_both_absence_and_surviving_evidence(drift):
    raw = sources("current")
    if drift == "empty_vafs":
        raw["vafs"] = []
    elif drift == "empty_metadata":
        raw["bam_metadata"] = []
    elif drift == "restored_viewer":
        raw["bams"].extend(r for r in sources("historical")["bams"] if "5GQLV9WSXQ" in r["url"])
    else:
        raw["vafs"].append(dict(raw["vafs"][0], bam_file="renamed-TL-24-5GQLV9WSXQ.bam"))
    correction = evaluate(raw, ["tempus-grch37-counts"])
    with pytest.warns(CurationWarning):
        assert correction.report()[0]["status"] == "stale"
    assert correction.records("vafs")[0] == raw["vafs"]


@pytest.mark.parametrize("gene", RELOCATED)
@pytest.mark.parametrize("source, key, field, value", [
    ("source_variants", "id", "ref", "N"),
    ("vafs", "variant_id", "pos", "1"),
    ("variant_index", "id", "location", "chr1:1"),
])
def test_unexpected_relocated_records_are_stale_and_atomic(gene, source, key, field, value):
    raw = sources("current")
    old_id, new_id = RELOCATED[gene]
    next(r for r in raw[source] if r[key] == new_id)[field] = value
    curated = evaluate(raw, ["allele-" + old_id])
    with pytest.warns(CurationWarning):
        assert curated.report()[0]["status"] == "stale"
    for name in raw:
        assert curated.records(name)[0] == raw[name]


@pytest.mark.parametrize("drift", ["both_ids", "mixed_sources", "unknown_id"])
def test_id_changes_require_a_complete_known_layout(drift):
    raw = sources("current")
    old = sources("historical")
    old_id, new_id = RELOCATED["CABLES1"]
    if drift == "both_ids":
        for name, key in [("source_variants", "id"), ("variant_index", "id"), ("vafs", "variant_id")]:
            raw[name].extend(r for r in old[name] if r[key] == old_id)
    elif drift == "mixed_sources":
        raw["source_variants"] = old["source_variants"]
    else:
        for name, key in [("source_variants", "id"), ("variant_index", "id"), ("vafs", "variant_id")]:
            for row in raw[name]:
                if row[key] == new_id:
                    row[key] = "CABLES1-chr18-unknown"
    curated = evaluate(raw, ["allele-" + old_id])
    with pytest.warns(CurationWarning):
        assert curated.report()[0]["status"] == "stale"
    assert curated.records("source_variants")[0] == raw["source_variants"]


def test_ush2a_history_and_current_merge_remain_distinct():
    old_id, new_id = "USH2A-chr1-215560752", "USH2A-chr1-215650752"
    historical = evaluate(sources("historical"), ["ush2a-transposed-duplicate"])
    old_records = by_id(historical.records("source_variants")[0])
    assert old_records[old_id]["alt"] == "not_reported"
    assert old_records[old_id]["allele_resolution"]["status"] == "source_identity_unconfirmed"
    assert old_records[new_id]["ref"] == "C"
    raw = sources("current")
    current = evaluate(raw, ["ush2a-transposed-duplicate"])
    new_records = by_id(current.records("source_variants")[0])
    assert old_id not in new_records
    assert (new_records[new_id]["pos"], new_records[new_id]["ref"], new_records[new_id]["alt"]) == (215650752, "C", "A")
    assert new_records[new_id]["genomic_location"] == "chr1:215650752"
    assert new_records[new_id]["upstream_merge"]["retired_id"] == old_id
    assert new_records[new_id]["upstream_merge"]["original_source_identity"] == "unconfirmed"
    assert current.records("vafs")[0] == raw["vafs"]


def test_changed_ush2a_allele_is_not_accepted_as_the_reviewed_merge():
    raw = sources("current")
    next(r for r in raw["source_variants"] if r["gene"] == "USH2A")["alt"] = "T"
    curated = evaluate(raw, ["ush2a-transposed-duplicate"])
    with pytest.warns(CurationWarning):
        assert curated.report()[0]["status"] == "stale"
    assert curated.records("source_variants")[0] == raw["source_variants"]


def test_ambiguous_versions_do_not_choose_an_edit():
    correction = Correction("ambiguous", "two layouts fit", (
        Change("rows", {"id": "one"}, set={"value": 1}),), alternatives=((
            Change("rows", {"id": "one"}, set={"value": 2}),),))
    raw = [{"id": "one", "value": 0}]
    curated = Curation([correction], lambda _: raw)
    with pytest.warns(CurationWarning, match="ambiguous_versions"):
        assert curated.report()[0]["status"] == "stale"
    assert curated.records("rows")[0] == raw


def test_corrected_field_does_not_hide_changed_identity_guard():
    correction = Correction("guard", "guard still matters", (
        Change("rows", {"id": "one"}, expect={"gene": "GENE", "value": 0}, set={"value": 1}),))
    curated = Curation([correction], lambda _: [{"id": "one", "gene": "OTHER", "value": 1}])
    with pytest.warns(CurationWarning):
        assert curated.report()[0]["status"] == "stale"


@pytest.mark.parametrize("guard", [[], [{"id": "witness", "gene": "OTHER"}]])
def test_corrected_values_do_not_hide_a_missing_or_changed_witness(guard):
    correction = Correction("guard", "check the whole layout", (
        Change("rows", {"id": "one"}, expect={"value": 0}, set={"value": 1}),
        Change("rows", {"id": "witness"}, expect={"gene": "GENE"})))
    curated = Curation([correction], lambda _: [{"id": "one", "value": 1}, *guard])
    with pytest.warns(CurationWarning):
        assert curated.report()[0]["status"] == "stale"


def test_cleared_historical_tempus_counts_are_fixed_upstream():
    raw = sources("historical")
    for row in raw["vafs"]:
        if "5GQLV9WSXQ" in row["bam_file"]:
            row.update(dict.fromkeys(COUNT_FIELDS, ""))
    curated = evaluate(raw, ["tempus-grch37-counts"])
    assert curated.report()[0]["status"] == "fixed_upstream"
    assert curated.records("vafs")[0] == raw["vafs"]


@pytest.mark.parametrize("version", ["historical", "current"])
def test_disabling_corrections_preserves_each_published_layout(version):
    raw = sources(version)
    curated = evaluate(raw, enabled=False)
    assert {r["status"] for r in curated.report()} == {"disabled"}
    for source in raw:
        assert curated.records(source)[0] == raw[source]


def test_absence_assertions_cannot_edit_or_expect_values():
    with pytest.raises(ValueError, match="absence"):
        Change("rows", {}, absent=True, set={"gene": "X"})
    with pytest.raises(ValueError, match="absence"):
        Change("rows", {}, absent=True, expect={"gene": "X"})
