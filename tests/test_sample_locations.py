"""Current published sample locations, including filename prefixes and row aliases."""

import html
import json
from collections import Counter
from pathlib import Path

from osteosarc import File, Files, Table
from osteosarc.catalog import build_files, fastq_prefix, matches_fastq, parse_data_paths
from osteosarc.curation import ASSAYS
from osteosarc.models import SampleClaim
from osteosarc.views import _fastq_folders, get_data_hints

REVIEW = json.loads((Path(__file__).parent / "data/sample-locations-2026-10-07.json").read_text())


def test_current_locations_link_every_fastq_without_mixing_tumor_and_normal(dataset, monkeypatch):
    sources = [dict(row=row, corrections=(), bams=(), missing=(), fastqs=tuple(
        r for r in REVIEW["fastqs"] if r["sample_id"] == row["sample_id"])) for row in REVIEW["specimens"]]
    monkeypatch.setattr(dataset, "_sample_sources", sources)
    files = build_files(dict(files=REVIEW["bucket"]), dict(categories=[]), Table([]), Table([]),
                        (("vendor/cegat/P116686_2_S000048", SampleClaim("data_page", assay="wes")),),
                        fastqs=REVIEW["fastqs"])
    dataset._tag_samples(files)
    monkeypatch.setattr(dataset, "files", files)
    assert len(dataset.samples) == 31
    for row in REVIEW["fastqs"]:
        sample = dataset.samples[row["sample_id"]]
        matched = [f for f in sample.files if matches_fastq(f.key, row["s3_folder"])]
        assert len(matched) == int(row["file_count"]) == 2
        assert sum(f.size for f in matched) == int(row["total_size_bytes"])
        assert all(f.samples == (sample.id,) for f in matched)
        assert all(f.resolved("assay") == ASSAYS[row["assay"]][0] for f in matched)
        assert all(row in f.metadata["fastq_rows"] for f in matched)
        shown = next(r for r in _fastq_folders(sample, sample.files)
                     if r["published_location"] == row["s3_folder"])
        assert shown["files"] == 2 and not shown["directory"]
        assert shown["folder"] == fastq_prefix(row["s3_folder"])
        assert row["s3_folder"] in sample.fastq_folders  # original alias stays visible
        commands = get_data_hints(dataset, [], [shown])
        assert "--exclude '*'" in commands
        for f in matched:
            assert "--include " + Path(f.key).name in commands
    assert not any(f.samples for f in files if f.kind != "reads")
    samples = dataset.samples
    assert Counter(s.tissue for s in samples) == {"blood": 24, "tumor": 5, "organoid": 2}
    assert [(samples[s].date, samples[s].site) for s in
            ("T0_tumor", "T1_tumor", "T2_tumor", "T3_tumor")] == [
                ("2022-12-16", "UCSF"), ("2024-06-06", "UCLA"),
                ("2025-01-28", "UCLA"), ("2025-04-17", "MSKCC")]
    assert samples["T2_organoid"].assays == ("wgs",)
    assert "another person's" in samples["T2_organoid"].notes
    assert "Composite of separate draws" in samples["T0_blood"].notes
    assert samples["blood_2024-06-11"].date == "2024-06-11"
    assert samples["blood_2024-09-21"].date == "2024-09-21"


def test_prefix_boundaries_sidecars_and_exact_files():
    prefix = "vendor/cegat/P116686_1"
    assert matches_fastq(prefix + ".1.fastq.gz", prefix)
    assert matches_fastq(prefix + "_R1.fq.gz", prefix)
    assert matches_fastq(prefix + "/lane/read.fastq.gz", prefix)
    assert matches_fastq(prefix + ".fastq.gz", prefix + ".fastq.gz")
    assert not matches_fastq(prefix + "0.1.fastq.gz", prefix)
    assert not matches_fastq(prefix + ".1.fastq.gz.seqkit-stats.tsv", prefix)
    assert not matches_fastq(prefix + ".bam", prefix)
    assert not matches_fastq("other/read.fastq.gz", prefix)
    assert not matches_fastq("other/read.fastq.gz", "")
    # Only the reviewed Tempus labels have aliases; no generic /T or /N guessing.
    assert fastq_prefix("another/FastQ/T") == "another/FastQ/T"


def test_fastq_library_metadata_preserves_conflicting_exact_file_claims():
    key = "library/read.fastq.gz"
    files = build_files(dict(files=[[key, 10, 0]]), dict(categories=[]), Table([]), Table([]),
                        ((key, SampleClaim("data_page", assay="wes")),), fastqs=[
                            dict(s3_folder="library", assay="RNA", display_name="RNA library")])
    file = next(f for f in files if f.key == key)
    assert file.conflicts["assay"] == ("rna-seq", "wes")
    assert file.resolved("assay") is None


def test_directory_counts_exclude_indexes_and_keep_nested_fastqs():
    from osteosarc.models import Sample

    sample = Sample("s", fastq_folders=("reads/run",))
    files = Files([File(str(i), key, key, kind, fmt, size=size) for i, (key, kind, fmt, size) in enumerate([
        ("reads/run/lane/R1.fastq.gz", "reads", "fastq", 10),
        ("reads/run/lane/R2.fastq.gz", "reads", "fastq", 20),
        ("reads/run/report.tsv", "table", "tsv", 100),
        ("reads/run2/R1.fastq.gz", "reads", "fastq", 200)])])
    row, = _fastq_folders(sample, files)
    assert row["files"] == 2 and row["size"] == "30 B"
    assert row["directory"] and row["folder"] == "reads/run/"


def test_data_page_notes_do_not_create_provider_or_assay_conflicts():
    claims = []
    for entry in REVIEW["data_rows"]:
        row = entry["row"]
        page = ("<h2>" + html.escape(entry["context"]) + "</h2><table><thead><tr>"
                + "".join("<th>" + html.escape(k) + "</th>" for k in row)
                + "</tr></thead><tbody><tr>"
                + "".join("<td>" + html.escape(v) + "</td>" for v in row.values())
                + "</tr></tbody></table>")
        claims.extend(parse_data_paths(page))
    assert claims
    for path, claim in claims:
        if "normal_bg24" in path:
            assert claim.provider == "BostonGene" and claim.assay == "wes"
        elif "5GQLV9WSXQ" in path:
            assert claim.provider == "Tempus" and claim.assay == "panel"
        else:
            assert "oncoanalyser" in path and claim.assay is None
