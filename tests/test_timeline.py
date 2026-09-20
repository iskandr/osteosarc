"""Timeline events, specimens, and the plain-text explorer, on public excerpts."""

import io
import json

import pytest

from osteosarc import Dataset, SchemaError
from osteosarc.cli import main
from osteosarc.explore import Explorer, specimen_view, specimens_view
from osteosarc.timeline import normalize_date


def test_events_keep_sources_precision_and_stated_timepoints(dataset):
    timeline = dataset.timeline
    assert {e.source for e in timeline} == {"events", "mrd", "flow", "imaging", "pathology", "specimens",
                                            "labs", "cytometry"}
    points = timeline.select(lane="Time points")
    assert [(e.date, e.timepoint) for e in points] == [
        ("2022-12-16", "T0"), ("2024-06-06", "T1"), ("2025-01-28", "T2"), ("2025-04-17", "T3")]
    # "T4/5" is a vertebral level, not a timepoint.
    biopsy = timeline.select(lane="Biopsy", contains="T4/5")[0]
    assert biopsy.timepoint is None
    proton = timeline.select(contains="Proton therapy")[0]
    assert (proton.date, proton.end, proton.open_end) == ("2024-06-12", "2024-08-02", False)
    assert timeline.select(kind="measurement", source="mrd")
    assert {e.details["value_kind"] for e in timeline.select(source="mrd")} >= {"numeric", "not_detected"}
    ongoing = [e for e in timeline if e.open_end]
    assert ongoing and all(e.end == "2026-09-16" for e in ongoing)
    assert list(timeline) == sorted(timeline, key=lambda e: (e.first_day, e.lane, e.label, e.id))


def test_vaccine_doses_carry_sheet_dose_numbers(dataset):
    doses = dataset.timeline.select(lane="Cancer vaccines")
    assert doses and any(e.value and "#" in e.value for e in doses)


def test_windows_and_neighbourhoods(dataset):
    timeline = dataset.timeline
    june = timeline.select(since="2024-06", until="2024-06")
    assert june and all(e.last_day.isoformat() >= "2024-06-01" and e.first_day.isoformat() <= "2024-06-30"
                        for e in june)
    assert "2024-06-12..2024-08-02" in timeline.around("2024-06-11", 3).listing()
    assert normalize_date("06/11/2024") == "2024-06-11" and normalize_date("20240611") == "2024-06-11"


def test_render_draws_lanes_axis_and_legend(dataset):
    chart = dataset.timeline.render(width=100)
    lines = chart.splitlines()
    assert "2023" in lines[0] and "|" in lines[1]
    assert any(line.startswith("Time points") and "T1" in line for line in lines)
    assert any(line.startswith("Treatments: Radiation") and "=" in line for line in lines)
    assert any(line.startswith("MRD: ") for line in lines)
    assert "one column =" in chart and max(len(line) for line in lines if not line.startswith("20")) <= 100
    assert dataset.timeline.select(lane="nothing").render() == "(no events)"


def test_specimens_are_corrected_linked_and_compared(dataset):
    specimens = {r["sample_id"]: r for r in dataset.specimens}
    t2 = specimens["T2_tumor"]
    assert (t2["date"], t2["site"], t2["corrections"]) == ("2025-01-28", "UCLA", ("specimen-T2-date-site",))
    assert not t2["disagreements"]
    t1 = specimens["T1_tumor"]
    assert "rna-seq/reprocessed/BG009368/BG009368.Aligned.sortedByCoord.out.md.bam" in t1["assets"]
    assert t1["fastq_folders"] and t1["tissue"] == "tumor"
    assert specimens["blood_2025-06-26"]["corrections"] == ("pbmc-capture-dates",)
    raw = {r["sample_id"]: r for r in Dataset.open("fixture", cache=dataset.cache, corrections=False).specimens}
    assert {(d["source"], d["field"]) for d in raw["T2_tumor"]["disagreements"]} == {
        ("timepoint_summary", "date"), ("timepoint_summary", "site"), ("events", "date")}


def test_measurements_keep_raw_values(dataset):
    rows = dataset.measurements
    assert {r["source"] for r in rows} == {"mrd", "labs", "cytometry"}
    assert all(isinstance(r["value"], str) for r in rows if r["source"] != "mrd")
    assert {r["kind"] for r in rows.select(source="mrd")} >= {"numeric", "not_detected"}


def test_timeline_corrections_flag_events(dataset):
    tempus = dataset.timeline.select(contains="Tempus")
    assert len(tempus) == 3
    statuses = {r["id"]: r["status"] for r in dataset.corrections}
    for name in ("tempus-timepoint", "events-duplicate-rows", "apheresis-date", "specimen-T3-site",
                 "reyagel-end-date"):
        assert statuses[name] == "applied"
    reyagel = dataset.timeline.select(contains="ReyaGel")[0]
    assert (reyagel.end, reyagel.corrections) == (None, ("reyagel-end-date",))


def test_views_and_interactive_explorer(dataset):
    assert "T2_tumor" in specimens_view(dataset, width=200)
    detail = specimen_view(dataset, "T1_tumor")
    assert "BG009368" in detail and "corrected by specimen-T1-site" in detail
    with pytest.raises(KeyError, match="known"):
        specimen_view(dataset, "T9_tumor")
    output = io.StringIO()
    script = io.StringIO("lanes\nzoom 2024-05 2024-09\nonly MRD\nevents Proton\nreset\n"
                         "on 2024-06-11 2\nspecimen T2_tumor\nvariants SMC5\ncorrections tempus-timepoint\n"
                         "specimen nope\nbogus\nquit\n")
    Explorer(dataset, width=100, stdin=script, stdout=output).cmdloop()
    text = output.getvalue()
    assert "Time points" in text and "2024-06" in text and "MRD: " in text
    assert "Proton therapy" in text  # events search all lanes despite 'only MRD'
    assert "SMC5-chr9-70298024" in text and "tempus-timepoint [applied]" in text
    assert "error: " in text and "Unknown syntax: bogus" in text


def test_cli_timeline_specimens_and_on(dataset, capsys):
    root = str(dataset.cache.root)
    assert main(["--cache", root, "timeline", "fixture", "--since", "2024", "--until", "2025", "--width", "90"]) == 0
    assert "Time points" in capsys.readouterr().out
    assert main(["--cache", root, "timeline", "fixture", "--lane", "MRD", "--json"]) == 0
    assert all(e["category"] == "MRD" for e in json.loads(capsys.readouterr().out))
    assert main(["--cache", root, "on", "fixture", "2025-01-28", "--days", "0"]) == 0
    assert "T2" in capsys.readouterr().out
    assert main(["--cache", root, "specimens", "fixture", "T3_tumor"]) == 0
    assert "MSKCC" in capsys.readouterr().out


def test_snapshots_without_timeline_sources_still_open(dataset):
    manifest = dict(dataset.manifest)
    from osteosarc.cache import stable_id
    from osteosarc.catalog import TIMELINE_SOURCES
    manifest["sources"] = {k: v for k, v in manifest["sources"].items() if k not in TIMELINE_SOURCES}
    manifest["id"] = stable_id(manifest["sources"])
    older = Dataset(dataset.cache, manifest)
    assert len(older.variants()) == len(dataset.variants())
    with pytest.raises(SchemaError, match="predates the timeline"):
        older.timeline
    statuses = {r["id"]: r["status"] for r in older.corrections}
    assert statuses["specimen-T2-date-site"] == "unavailable"


def test_timeline_events_carry_their_corrections(dataset):
    tempus = dataset.timeline.select(contains="Tempus")
    assert all(e.corrections == ("tempus-timepoint",) for e in tempus)
    specimen = dataset.timeline.select(lane="Specimens", contains="blood_2025-06-26")[0]
    assert specimen.corrections == ("pbmc-capture-dates",)
    assert "(corrections: tempus-timepoint)" in tempus.listing()


def test_mrd_measurements_are_raw_strings_with_units(dataset):
    rows = dataset.measurements.select(source="mrd")
    assert all(isinstance(r["value"], str) and "{" not in r["value"] and r["value"] != "None" for r in rows)
    below = [r for r in rows if r["kind"] == "below_loq"]
    assert below and all(float(r["value"]) > 0 for r in below)
    assert {r["unit"] for r in rows if r["measurement"] == "Signatera"} == {"MTM/mL"}
    assert all(r["value"] == "not_detected" for r in rows if r["kind"] == "not_detected")


def test_invalid_dates_are_clear_errors_everywhere(dataset, capsys):
    for bad in ("2024/06", "June", "2024-13"):
        with pytest.raises(ValueError):
            dataset.timeline.select(since=bad)
        with pytest.raises(ValueError):
            dataset.timeline.around(bad)
    assert main(["--cache", str(dataset.cache.root), "timeline", "fixture", "--since", "2024/06"]) == 1
    assert "Expected a date" in capsys.readouterr().err
    output = io.StringIO()
    session = io.StringIO("zoom 2024-05 2024-09\nzoom June\nzoom 2024-13\nevents Proton\nquit\n")
    explorer = Explorer(dataset, width=100, stdin=session, stdout=output)
    explorer.cmdloop()
    text = output.getvalue()
    assert text.count("error: ") == 2 and "Proton therapy" in text
    assert (explorer.since, explorer.until) == ("2024-05", "2024-09")  # bad zooms left it unchanged


def test_unreadable_source_dates_are_reported_not_fatal(dataset):
    from osteosarc import CORRECTIONS, Change, Correction
    drift = Correction("drift", "simulate upstream date formats", (
        Change("flow", {"draw_date": "2025-06-24"}, set={"draw_date": "06/24/2025"}),
        Change("flow", {"draw_date": "2025-07-22"}, set={"draw_date": "soon"}),
        Change("events", {"title": "Trabectedin"}, set={"date": "2024-02-30"})))
    data = Dataset.open("fixture", cache=dataset.cache, corrections=[*CORRECTIONS, drift])
    draws = data.timeline.select(source="flow")
    assert "2025-06-24" in [e.date for e in draws]
    undated = data.timeline.source["undated"]
    assert {(u["source"], u["record"].get("draw_date") or u["record"].get("date")) for u in undated} == {
        ("flow", "soon"), ("events", "2024-02-30")}
    assert len(data.specimens) == len(dataset.specimens)


def test_explorer_options_convert_booleans_and_integers(dataset):
    from osteosarc.explore import options
    assert options(["include_conflicts=false", "limit=5", "assay=rna-seq"]) == dict(
        include_conflicts=False, limit=5, assay="rna-seq")
    with pytest.raises(ValueError, match="true or false"):
        options(["include_conflicts=maybe"])
    output = io.StringIO()
    Explorer(dataset, width=120, stdin=io.StringIO("variants limit=1\nassets include_inferred=no kind=bam\nquit\n"),
             stdout=output).cmdloop()
    assert "... " in output.getvalue() and "error" not in output.getvalue()


def test_table_limits_are_validated():
    from osteosarc.explore import table
    rows = [dict(a=str(i)) for i in range(4)]
    assert table(rows, ("a",), limit=0).endswith("... 4 more")
    assert "... 2 more" in table(rows, ("a",), limit=2)
    with pytest.raises(ValueError):
        table(rows, ("a",), limit=-2)


def test_sample_overview_keeps_full_assay_labels_and_filters(dataset, capsys):
    text = dataset.describe_samples(timepoint="T1", tissue="tumor", width=80)
    assert "T1_tumor" in text and "T0_tumor" not in text and "T1_blood" not in text
    assert "sequencing" in text and "FASTQ_folders" in text
    specimen = next(r for r in dataset.specimens if r["sample_id"] == "T1_tumor")
    for assay in specimen["assays"]:
        assert assay in text
    assert dataset.describe_samples(timepoint="missing") == "(no matching samples)"
    assert main(["--cache", str(dataset.cache.root), "samples", "fixture", "--timepoint", "T1", "--tissue", "tumor"]) == 0
    assert "T1_tumor" in capsys.readouterr().out
    assert main(["--cache", str(dataset.cache.root), "samples", "fixture", "--json"]) == 0
    assert json.loads(capsys.readouterr().out) == list(dataset.samples)


def test_sample_asset_selection_includes_only_linked_files(dataset):
    from osteosarc import Asset, Assets
    from osteosarc.cache import stable_id

    sample = next(r for r in dataset.specimens if r["sample_id"] == "T1_tumor")
    folder = dataset._object_key(sample["fastq_folders"][0]).rstrip("/")
    files = []
    for key in (folder + "/reads.fastq.gz", folder + "-other/reads.fastq.gz"):
        url = "https://example.test/" + key
        files.append(Asset(stable_id(url), key, url, "reads", "fastq"))
    dataset.assets = Assets([*dataset.assets, *files])
    selected = dataset.assets_for_sample("T1_tumor")
    assert files[0] in selected and files[1] not in selected
    assert selected[files[0].key] is selected[files[0].url] is selected[files[0].id]
    with pytest.raises(KeyError):
        selected[files[1].key]
    assert any(a.key in sample["assets"] for a in selected)
    assert all(a.kind == "alignment" for a in dataset.assets_for_sample("T1_tumor", kind="alignment"))
    with pytest.raises(KeyError, match="Unknown sample"):
        dataset.assets_for_sample("typo")
