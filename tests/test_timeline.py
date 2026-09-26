"""Timeline events, samples, and their plain-text views, on public excerpts."""

import json

import pytest

from osteosarc import Dataset, SchemaError
from osteosarc.cli import main
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


def test_the_chart_names_treatments_and_hides_frequent_records(dataset):
    chart = dataset.timeline.render(width=100)
    lines = chart.splitlines()
    assert "2023" in lines[0] and "J" in lines[1]  # years over month initials
    assert any(line.startswith("Time points") and "T1" in line for line in lines)
    # Each treatment is a named row under its group, with ranges drawn as =.
    assert "Radiation" in lines and any(line.startswith("  Proton therapy") and "=" in line for line in lines)
    assert "MRD" in lines and any(line.startswith("  Signatera") for line in lines)
    assert not any(c.isdigit() for line in lines[2:] for c in line.split("  ")[-1] if line.startswith("  "))
    legend = " ".join(chart.split())  # the legend wraps to the chart's width
    assert "Each column is a month." in legend and "Not shown:" in legend and "Lab draws" in legend
    assert max(len(line) for line in lines) <= 100
    everything = dataset.timeline.render(width=100, everything=True)
    assert "Records" in everything.splitlines() and "Not shown" not in everything
    # A short window gets several columns per month, labeled by name.
    zoomed = dataset.timeline.render(width=100, since="2024-05", until="2024-09")
    assert "May" in zoomed.splitlines()[1] and "Each month is" in zoomed
    assert dataset.timeline.select(lane="nothing").render() == "(no events)"
    # Asking only for hidden lanes still charts them.
    assert "Lab draws" in dataset.timeline.select(lane="Lab").render(width=100)


def test_samples_are_corrected_linked_and_compared(dataset):
    samples = dataset.samples
    t2 = samples["T2_tumor"]
    assert (t2.date, t2.site, t2.corrections) == ("2025-01-28", "UCLA", ("specimen-T2-date-site",))
    assert not t2.disagreements
    t1 = samples["T1_tumor"]
    assert "rna-seq/reprocessed/BG009368/BG009368.Aligned.sortedByCoord.out.md.bam" in t1.bams
    assert t1.fastq_folders and t1.tissue == "tumor" and "rna-seq" in t1.assays
    assert samples["blood_2025-06-26"].corrections == ("pbmc-capture-dates",)
    raw = Dataset.open("fixture", cache=dataset.cache, corrections=False).samples
    assert {(d["source"], d["field"]) for d in raw["T2_tumor"].disagreements} == {
        ("timepoint_summary", "date"), ("timepoint_summary", "site"), ("events", "date")}
    with pytest.raises(KeyError, match="No sample 'T9_tumor'; samples: T1_tumor"):
        samples["T9_tumor"]


def test_sequencing_also_comes_from_the_fastq_table(dataset):
    from osteosarc import CORRECTIONS, Change, Correction
    # Like the site's 2026 blood draws: nothing in the registry, but FASTQ folders that say what they hold.
    blank = Correction("blank", "simulate a registry row without assays", (
        Change("specimens", {"sample_id": "T2_tumor"}, set={"assays_run": "", "vendors_involved": ""}),))
    data = Dataset.open("fixture", cache=dataset.cache, corrections=[*CORRECTIONS, blank])
    t2 = data.samples["T2_tumor"]
    assert t2.details["registry"]["assays_run"] == ""
    assert t2.sequencing and set(t2.assays) <= set(dataset.samples["T2_tumor"].assays)
    assert t2.providers


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


def test_one_sample_shows_its_files_and_how_to_get_them(dataset):
    from osteosarc.views import sample_view
    t1 = dataset.samples["T1_tumor"]
    shell, python = sample_view(t1, dataset), sample_view(t1, dataset, python=True)
    assert shell.startswith("T1_tumor: ") and "Corrected by specimen-T1-site" in shell
    assert "BAMs, aligned reads (1)" in shell and t1.bams[0] in shell
    assert "FASTQ folders, raw reads" in shell and t1.fastq_folders[0] + "/" in shell
    assert f"osteosarc download {t1.bams[0]} --to ." in shell and "aws s3 cp --recursive --no-sign-request" in shell
    assert f'data.download("{t1.bams[0]}", to=".")' in python and "osteosarc download" not in python
    assert repr(t1) == python
    empty = dataset.samples["blood_2025-06-26"]
    assert "BAMs: none" in sample_view(empty, dataset) and "(no files)" in sample_view(empty, dataset)


def test_cli_timeline_samples_and_on(dataset, capsys):
    root = str(dataset.cache.root)
    assert main(["--cache", root, "timeline", "--snapshot", "fixture", "--since", "2024", "--until", "2025", "--width", "90"]) == 0
    assert "Time points" in capsys.readouterr().out
    assert main(["--cache", root, "timeline", "--snapshot", "fixture", "--lane", "MRD", "--json"]) == 0
    assert all(e["category"] == "MRD" for e in json.loads(capsys.readouterr().out))
    assert main(["--cache", root, "timeline", "--snapshot", "fixture", "--days", "3"]) == 1
    assert "--around DATE" in capsys.readouterr().err
    assert main(["--cache", root, "timeline", "--snapshot", "fixture", "--around", "2025-01-28", "--days", "0"]) == 0
    assert "T2" in capsys.readouterr().out
    assert main(["--cache", root, "samples", "--snapshot", "fixture", "T3_tumor"]) == 0
    assert "MSKCC" in capsys.readouterr().out
    assert main(["--cache", root, "timeline", "--snapshot", "fixture", "--list", "--contains", "Proton"]) == 0
    assert "Proton therapy" in capsys.readouterr().out
    assert main(["--cache", root, "timeline", "--snapshot", "fixture", "--all", "--width", "90"]) == 0
    assert "Lab draws" in capsys.readouterr().out


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
    drawn = dataset.timeline.select(lane="Samples", contains="blood_2025-06-26")[0]
    assert drawn.corrections == ("pbmc-capture-dates",)
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
    assert main(["--cache", str(dataset.cache.root), "timeline", "--snapshot", "fixture", "--since", "2024/06"]) == 1
    assert "Expected a date" in capsys.readouterr().err


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
    assert len(data.samples) == len(dataset.samples)


def test_table_limits_are_validated():
    from osteosarc.views import table
    rows = [dict(a=str(i)) for i in range(4)]
    assert table(rows, ("a",), limit=0).endswith("... 4 more")
    assert "... 2 more" in table(rows, ("a",), limit=2)
    with pytest.raises(ValueError):
        table(rows, ("a",), limit=-2)


def test_sample_overview_shows_filter_names_and_filters_by_assay(dataset, capsys):
    from osteosarc.curation import sequencing_pairs
    from osteosarc.views import sequencing_text
    assert sequencing_text(sequencing_pairs(["PacBio", "RNA", "WES", "WGS", "scRNA", "scRNA_ONT"])) == \
        "rna-seq, wes, wgs, scrna-seq (ont, pacbio)"
    assert sequencing_text(sequencing_pairs(["CITE", "New label", "scRNA_TCRgd", "bulk RNA"])) == \
        "rna-seq, scrna-seq, cite-seq, New label"
    text = repr(dataset.samples.select(timepoint="T1", tissue="tumor"))
    assert "T1_tumor" in text and "T0_tumor" not in text and "T1_blood" not in text
    assert "sequencing" in text and "FASTQ folders" in text
    t1 = dataset.samples["T1_tumor"]
    assert sequencing_text(t1.sequencing) in repr(dataset.samples)
    assert repr(dataset.samples.select(timepoint="missing")) == "(no matching samples)"
    with_rna = [s.id for s in dataset.samples if "rna-seq" in s.assays]
    assert with_rna and [s.id for s in dataset.samples.select(assay="rna-seq")] == with_rna
    assert not dataset.samples.select(assay="cite-seq", platform="pacbio")
    with pytest.raises(ValueError, match="registry label"):
        dataset.samples.select(assay="scRNA_ONT")
    root = str(dataset.cache.root)
    assert main(["--cache", root, "samples", "--snapshot", "fixture", "--timepoint", "T1", "--tissue", "tumor"]) == 0
    assert "T1_tumor" in capsys.readouterr().out
    assert main(["--cache", root, "samples", "--snapshot", "fixture", "--assay", "rna-seq"]) == 0
    shown = capsys.readouterr().out
    assert all(s in shown for s in with_rna) and "osteosarc samples T1_tumor" in shown
    assert main(["--cache", root, "samples", "--snapshot", "fixture", "--json"]) == 0
    assert [r["id"] for r in json.loads(capsys.readouterr().out)] == [s.id for s in dataset.samples]
    assert main(["--cache", root, "samples", "--snapshot", "fixture", "--files"]) == 0
    assert "== T1_tumor: " in (shown := capsys.readouterr().out) and t1.bams[0] in shown
    assert main(["--cache", root, "samples", "--snapshot", "fixture", "T1_tumor", "--tissue", "blood"]) == 1
    assert "leave them out to see one" in capsys.readouterr().err


def test_a_samples_files_are_its_bams_and_its_fastq_folders(dataset):
    from osteosarc import File, Files
    from osteosarc.cache import stable_id

    t1 = dataset.samples["T1_tumor"]
    folder = t1.fastq_folders[0]
    files = []
    for key in (folder + "/reads.fastq.gz", folder + "-other/reads.fastq.gz"):
        url = "https://example.test/" + key
        files.append(File(stable_id(url), key, url, "reads", "fastq"))
    dataset._tag_samples(files)
    dataset.files = Files([*dataset.files, *files])
    selected = t1.files
    assert files[0] in selected and files[1] not in selected
    assert files[0].samples == ("T1_tumor",) and files[1].samples == ()
    assert selected[files[0].key] is selected[files[0].url] is selected[files[0].id]
    with pytest.raises(KeyError):
        selected[files[1].key]
    assert set(t1.bams) <= {f.key for f in selected}
    assert all(f.kind == "alignment" for f in dataset.files.select(sample="T1_tumor", kind="alignment"))
    assert not dataset.files.select(sample="typo")
