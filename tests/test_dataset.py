import json
from collections import Counter
from pathlib import Path

import pytest

from osteosarc import Dataset, IntegrityError
from osteosarc.cache import stable_id
from osteosarc.cli import main


def test_real_source_identity_and_unresolved_variants(dataset):
    variants = dataset.variants()
    assert len(variants) == 5
    assert variants["DYNC1H1-chr14-101980529"].allele == ("chr14", 101980529, "G", "A")
    assert variants["COL3A1-Splice"].status == "ready"
    assert len(dataset.variants(status="ready")) == 5
    deletion = variants["GTF3C5-chr9-133057893"].region()
    # The original anchored REF has 13 bases: anchor plus twelve deleted bases.
    assert (deletion.start, deletion.end) == (133057892, 133057905)
    assert len(variants.regions()) == 5
    raw = Dataset.open("fixture", cache=dataset.cache, corrections=False).variants()
    assert raw["COL3A1-Splice"].status == "missing_literal_allele"
    with pytest.raises(ValueError):
        raw.regions()


def test_vaccine_peptides_annotations_and_missing_states(dataset):
    smc5 = dataset.variants(gene="SMC5")[0]
    assert {"mRNA", "JLF V3"} <= set(smc5.vaccines)
    assert dataset.variants(vaccine="mRNA", gene="SMC5")[0] == smc5
    assert not dataset.variants(vaccine="mRNA", vaccine_source="source_variants", gene="SMC5")
    assert "vaccine_membership_conflict" in smc5.annotations
    row = next(r for r in dataset.vaccines if r["gene"] == "SMC5")
    assert row["elispot_status"] == "not_tested"
    assert row["elispot_response"] is None
    assert len(dataset.vaccine_peptides("mRNA")) > 0
    assert smc5.annotations["source_record"]["consequence"] == "missense_variant"


def test_sample_conflicts_do_not_silently_become_replicates(dataset):
    raw = Dataset.open("fixture", cache=dataset.cache, corrections=False)
    bams = raw.files.select(kind="alignment")
    assert len(bams) == 3
    source = bams.select(contains="BG009368")[0]
    assert set(source.values("timepoint")) == {"T0", "T1"}
    assert source.resolved("timepoint") is None
    assert not bams.select(contains="BG009368", timepoint="T0")
    assert len(bams.select(contains="BG009368", timepoint="T1", include_conflicts=True)) == 1
    assert source.metadata["catalog_genome_assertion"] == "hg38"
    assert "assembly" not in source.metadata
    # The viewer label's stale date is surfaced alongside its timepoint.
    assert source.conflicts["date"] == ("2022-12", "2024-06")
    # The central correction for that stale label resolves it, traceably.
    corrected = dataset.file(source.key)
    assert corrected.resolved("timepoint") == "T1" and not corrected.conflicts
    assert corrected.metadata["corrections"] == ("viewer-label-BG009368",)


def test_catalog_normalization_does_not_invent_conflicts():
    from osteosarc.catalog import BUCKET, build_files
    from osteosarc.models import SampleClaim
    from osteosarc.parsing import Table
    cite, rna = "kamil/blood/Nov2025_CITE/possorted_genome_bam.bam", "vendor/cegat/P2/P3.bam"
    listing = dict(files=[[cite, 1, 0], [rna, 1, 0], ["vendor/cegat/P2/P2.1.fastq.gz", 1, 0]])
    bams = dict(baseUrl=BUCKET, genome="hg38", categories=[dict(name="Blood scRNA", bams=[
        dict(name="Hudson Lab Blood CITE 2025-11 (Pool 1)", url=BUCKET + cite, tissue="blood")])])
    metadata = Table([dict(s3_path=BUCKET + cite, display_name="CITE pool", assay="CITE",
                           timepoint="", sample_date="2025-11-06", tissue="Blood", provider="Hudson Lab")])
    vafs = Table([], columns=("bam_file",))
    # A directory row for WES FASTQs must not override the RNA BAM's own row.
    paths = (("vendor/cegat/P2", SampleClaim("data_page", "WES FASTQ", assay="wes")),
             (rna, SampleClaim("data_page", "RNA BAM", assay="rna-seq")))
    files = {f.key: f for f in build_files(listing, bams, metadata, vafs, paths)}
    assert files[cite].values("assay") == ("cite-seq",)
    assert files[cite].values("provider") == ("Hudson Lab",)
    assert files[cite].values("date") == ("2025-11-06",)  # month precision agrees
    assert not files[cite].conflicts
    assert files[rna].values("assay") == ("rna-seq",)
    assert files["vendor/cegat/P2/P2.1.fastq.gz"].values("assay") == ("wes",)


def test_reopen_is_offline_and_uses_pinned_receipts(dataset):
    reopened = Dataset.open("fixture", cache=dataset.cache)
    assert reopened.cache.offline
    assert reopened.id == dataset.id
    assert reopened.variants().to_records() == dataset.variants().to_records()
    path = reopened.download("vafs")
    assert len(reopened.parse("vafs")) > 0
    assert all(isinstance(r["alt_reads"], str) for r in reopened.parse("vafs"))
    assert Counter(v.status for v in reopened.variants()) == Counter(v.status for v in dataset.variants())
    path.write_text("corrupt")
    with pytest.raises(IntegrityError):
        Dataset.open("fixture", cache=dataset.cache)


def test_snapshot_and_cli_use_same_selection(dataset, capsys):
    root = str(dataset.cache.root)
    assert main(["--cache", root, "variants", "--snapshot", "fixture", "--gene", "SMC5", "--json"]) == 0
    assert json.loads(capsys.readouterr().out)[0]["id"] == dataset.variants(gene="SMC5")[0].id
    # Text by default: a table, and one variant in detail with its published counts.
    assert main(["--cache", root, "variants", "--snapshot", "fixture", "--gene", "SMC5"]) == 0
    text = capsys.readouterr().out
    assert text.startswith("1 variant") and dataset.variants(gene="SMC5")[0].id in text
    assert main(["--cache", root, "variants", "--snapshot", "fixture", dataset.variants(gene="SMC5")[0].id]) == 0
    assert "Read counts published by the site" in (text := capsys.readouterr().out) and "Vaccines:" in text
    assert main(["--cache", root, "variants", "--snapshot", "fixture", "NOPE-chr1-1"]) == 1
    assert "No variant 'NOPE-chr1-1'" in capsys.readouterr().err
    assert main(["--cache", root, "files", "--snapshot", "fixture", "--kind", "alignment", "--limit", "1", "--json"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["total"] == 3
    assert len(result["files"]) == 1
    with pytest.raises(FileExistsError):
        Dataset.sync("fixture", cache=dataset.cache, refresh=True)


def test_cli_version_and_sample_files_match_the_api(dataset, capsys):
    import osteosarc
    with pytest.raises(SystemExit) as exited:
        main(["--version"])
    assert exited.value.code == 0
    assert capsys.readouterr().out.strip() == f"osteosarc {osteosarc.__version__}"
    root = str(dataset.cache.root)
    assert main(["--cache", root, "files", "--snapshot", "fixture", "--sample", "T1_tumor", "--kind", "alignment",
                 "--json"]) == 0
    result = json.loads(capsys.readouterr().out)
    expected = dataset.samples["T1_tumor"].files.select(kind="alignment")
    assert result["total"] == len(expected) == 1
    assert [f["key"] for f in result["files"]] == [f.key for f in expected]
    assert main(["--cache", root, "files", "--snapshot", "fixture", "--sample", "T1_tumor", "--assay", "wgs",
                 "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["total"] == 0
    # Without --json, a readable table with the complete key.
    assert main(["--cache", root, "files", "--snapshot", "fixture", "--sample", "T1_tumor", "--kind", "alignment"]) == 0
    text = capsys.readouterr().out
    assert text.startswith("1 file (") and expected[0].key in text and "rna-seq" in text
    assert main(["--cache", root, "files", "--snapshot", "fixture", "--kind", "alignment", "--limit", "1"]) == 0
    assert "... 2 more" in (text := capsys.readouterr().out) and "--limit N" in text
    # With no filters, an overview by kind and folder.
    assert main(["--cache", root, "files", "--snapshot", "fixture"]) == 0
    assert "Largest top-level folders" in (text := capsys.readouterr().out) and "alignment" in text
    # Registry labels and unknown values are errors, not empty selections.
    assert main(["--cache", root, "files", "--snapshot", "fixture", "--assay", "scRNA_ONT"]) == 1
    assert "select assay 'scrna-seq' with platform 'ont'" in capsys.readouterr().err
    assert main(["--cache", root, "files", "--snapshot", "fixture", "--platform", "nanopore"]) == 1
    assert "choose from: illumina, ont, pacbio" in capsys.readouterr().err
    assert main(["--cache", root, "files", "--snapshot", "fixture", "--sample", "T9_tumor"]) == 1
    assert "No sample 'T9_tumor'" in capsys.readouterr().err


def test_cli_reads_accepts_catalogue_variants(dataset, capsys, monkeypatch, tmp_path):
    from osteosarc import ReadSubset
    calls = []

    def extract_reads(self, asset, regions=None, **kwargs):
        calls.append(dict(asset=asset, regions=regions, **kwargs))
        return ReadSubset(tmp_path / "reads.bam", tmp_path / "reads.bam.bai", {"records": 0})
    monkeypatch.setattr(Dataset, "extract_reads", extract_reads)
    root, source = str(dataset.cache.root), dataset.files.select(format="bam")[0].key
    ids = ["DYNC1H1-chr14-101980529", "SMC5-chr9-70298024"]
    command = ["--cache", root, "reads", "--snapshot", "fixture", source, "--variant", ids[0], "--variant", ids[1]]
    assert main(command + ["--padding", "100"]) == 0
    assert capsys.readouterr().out.strip() == str(tmp_path / "reads.bam")
    calls.pop()
    assert main(command + ["--padding", "100", "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["path"] == str(tmp_path / "reads.bam")
    call = calls.pop()
    assert [v.id for v in call["variants"]] == sorted(ids)
    assert call["regions"] is None and call["padding"] == 100
    # The same selection through Python yields the same padded one-based regions.
    assert call["variants"].regions(padding=100) == dataset.variants(ids=ids).regions(padding=100)

    expected = [("chr14", 101980528, 101980530, "GRCh38"), ("chr14", 101980599, 101980600, "GRCh38")]
    # Regions may precede or follow options, as with argparse on every supported Python.
    for arguments in (["chr14:101980529-101980530", "chr14:101980600-101980600", "--assembly", "GRCh38"],
                      ["--assembly", "GRCh38", "chr14:101980529-101980530", "chr14:101980600-101980600"],
                      ["chr14:101980529-101980530", "--min-mapq", "0", "--assembly", "GRCh38",
                       "chr14:101980600-101980600"]):
        assert main(["--cache", root, "reads", "--snapshot", "fixture", source, *arguments]) == 0
        capsys.readouterr()
        assert [(r.contig, r.start, r.end, r.assembly) for r in calls.pop()["regions"]] == expected

    for arguments, message in [
        (["--variant", "NOT-A-VARIANT"], "Unknown variant ID"),
        (["chr14:1-2", "--variant", ids[0]], "either regions (with --assembly) or --variant"),
        (["--variant", ids[0], "chr14:1-2"], "either regions (with --assembly) or --variant"),
        (["--variant", ids[0], "--assembly", "GRCh38"], "either regions (with --assembly) or --variant"),
        (["chr14:1-2"], "--assembly is required"),
    ]:
        assert main(["--cache", root, "reads", "--snapshot", "fixture", source, *arguments]) == 1
        assert message in capsys.readouterr().err
    assert not calls
    # Dataset.extract_reads enforces the remaining rules, with the same messages as in Python.
    monkeypatch.undo()
    for arguments, message in [
        (["chr14:1-2", "--assembly", "GRCh38", "--padding", "5"], "padding requires variants"),
        ([], "nonempty sequence of regions or ready variants"),
    ]:
        assert main(["--cache", root, "reads", "--snapshot", "fixture", source, *arguments]) == 1
        assert message in capsys.readouterr().err
    with pytest.raises(SystemExit):
        main(["--cache", root, "reads", "--snapshot", "fixture", source, "--bogus"])
    with pytest.raises(SystemExit):
        main(["--cache", root, "samples", "--snapshot", "fixture", "T1_tumor", "extra"])


def test_acquired_data_remains_pinned_after_url_refresh(dataset, tmp_path):
    from osteosarc import File
    url = "https://example.test/results.tsv"
    path = tmp_path / "results.tsv"
    path.write_text("id\tcount\na\t0\n")
    asset = File(stable_id(url), "results.tsv", url, "table", "tsv")
    old = dataset.cache.import_file(path, url)
    pinned = dataset.download(asset)
    path.write_text("id\tcount\na\t9\n")
    new = dataset.cache.import_file(path, url)
    assert old.sha256 != new.sha256
    assert dataset.download(asset) == pinned
    assert dataset.parse(asset).rows[0]["count"] == "0"
    with pytest.raises(ValueError, match="new snapshot"):
        dataset.download(asset, refresh=True)


def test_unsupported_parser_refuses_before_downloading_alignment(dataset, monkeypatch):
    asset = dataset.files.select(format="bam")[0]
    monkeypatch.setattr(dataset, "download", lambda *a: pytest.fail("downloaded unsupported format"))
    with pytest.raises(ValueError, match="No built-in parser"):
        dataset.parse(asset)


def test_first_download_must_match_the_inventory_time():
    from osteosarc import File
    from osteosarc.cache import Receipt
    from osteosarc.dataset import _check_inventory_time
    asset = File("id", "a.tsv", "https://example.test/a.tsv", "table", "tsv", modified=1784586023)
    same = Receipt("https://example.test/a.tsv", "0" * 64, 1, "a.tsv", "now",
                   last_modified="Mon, 20 Jul 2026 22:20:23 GMT")
    _check_inventory_time(asset, same)
    newer = Receipt(**{**same.to_dict(), "last_modified": "Tue, 21 Jul 2026 09:00:00 GMT"})
    with pytest.raises(IntegrityError, match="new snapshot"):
        _check_inventory_time(asset, newer)


def test_extraction_binds_the_listed_index_to_the_snapshot(dataset, monkeypatch):
    import osteosarc.reads
    from osteosarc import Region
    source = dataset.file("rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam")
    assert source.index_urls
    # A source from an earlier session must also work on a freshly reopened
    # dataset, before its catalog and curation state have been populated.
    dataset = Dataset.open("fixture", cache=dataset.cache, offline=False)
    downloads, calls = [], []
    monkeypatch.setattr(osteosarc.reads, "require_samtools", lambda **k: None)
    monkeypatch.setattr(dataset, "download", lambda asset: downloads.append(asset) or "/pinned/index.bai")
    monkeypatch.setattr(osteosarc.reads, "extract_reads", lambda *a, **k: calls.append(k))
    dataset.extract_reads(source, [Region("chr1", 100, 140, "GRCh38")])
    assert downloads == [source.index_urls[0]] and calls[0]["index"] == "/pinned/index.bai"


def test_mirrored_site_tables_stay_pinned(tmp_path):
    from conftest import DATA, FILES

    from osteosarc import Cache
    from osteosarc.catalog import SNAPSHOT_SOURCES, TIMELINE_SOURCES
    mirror = "https://mirror.example.test/variant_vafs_long.tsv"
    cache = Cache(tmp_path / "cache", offline=True)
    urls = {**SNAPSHOT_SOURCES, **TIMELINE_SOURCES, "vafs": mirror}
    for key, name in FILES.items():
        cache.import_file(DATA / name, urls[key])
    data = Dataset.sync("mirrored", cache=cache, sources={"vafs": mirror})
    assert data.file("vafs").url == mirror
    assert len(data.parse("vafs")) == len(data.vafs) > 0   # offline: the pinned copy, not a download


def test_choose_snapshot_selects_exact_names_download_dates_or_ids():
    from osteosarc.dataset import choose_snapshot, downloaded_at
    rows = [dict(name="2026-10-02", downloaded="2026-10-02T09:00:00+00:00", id="c" * 64),
            dict(name="2026-09-24.2", downloaded="2026-09-24T18:00:00+00:00", id="b" * 64),
            dict(name="2026-09-24", downloaded="2026-09-24T08:00:00+00:00", id="a" * 64),
            dict(name="baseline", downloaded="2026-09-18T16:07:17+00:00", id="a" * 64)]
    assert choose_snapshot(rows) == "2026-10-02"
    assert choose_snapshot(rows, "2026-09-24") == "2026-09-24"  # Names are exact
    assert choose_snapshot(rows, date="2026-09-24") == "2026-09-24.2"  # Newest downloaded that day
    assert choose_snapshot(rows, date="2026-09") == "2026-09-24.2"
    assert choose_snapshot(rows, date="2026-09-18") == "baseline"  # Dates are download dates
    assert choose_snapshot(rows, date="2026") == "2026-10-02"
    assert choose_snapshot(rows, "bbbbbb") == "2026-09-24.2"
    assert choose_snapshot(rows, "AAAAAA") == "2026-09-24"  # Same content: the newest name
    for name in ("2026-09", "zzzzzz", "bbb", "missing"):
        with pytest.raises(FileNotFoundError, match="No snapshot named"):
            choose_snapshot(rows, name)
    with pytest.raises(FileNotFoundError, match="No snapshot downloaded in 2025"):
        choose_snapshot(rows, date="2025")
    with pytest.raises(ValueError, match="YYYY-MM"):
        choose_snapshot(rows, date="2026-9")
    with pytest.raises(ValueError, match="not both"):
        choose_snapshot(rows, "baseline", date="2026")
    with pytest.raises(FileNotFoundError, match="Dataset.sync"):
        choose_snapshot([])
    ambiguous = [dict(rows[0], id="abcdef" + "0" * 58), dict(rows[1], id="abcdef" + "1" * 58)]
    with pytest.raises(FileNotFoundError, match="ambiguous"):
        choose_snapshot(ambiguous, "abcdef")
    assert downloaded_at(dict(sources={}, created_at="2026-09-24T00:00:00+00:00")) == "2026-09-24T00:00:00+00:00"


@pytest.fixture
def source_cache(tmp_path):
    """An offline cache holding every fixture source, with no snapshot yet."""
    from conftest import DATA, FILES

    from osteosarc import Cache
    from osteosarc.catalog import SNAPSHOT_SOURCES, TIMELINE_SOURCES
    cache = Cache(tmp_path / "cache", offline=True)
    for key, name in FILES.items():
        cache.import_file(DATA / name, {**SNAPSHOT_SOURCES, **TIMELINE_SOURCES}[key])
    return cache


def test_dated_sync_reuses_today_and_open_uses_the_newest(source_cache, tmp_path, monkeypatch):
    from osteosarc import Cache, OfflineError
    with pytest.raises(FileNotFoundError, match="No saved snapshots"):
        Dataset.open(cache=source_cache)
    # Offline, sync builds a snapshot from cached sources, dated by their download.
    first = Dataset.sync(cache=source_cache)
    day = first.downloaded[:10]
    assert first.name == day
    assert Dataset.sync(cache=source_cache).name == day  # Same day: reopened
    with pytest.raises(OfflineError):
        Dataset.sync(cache=source_cache, refresh=True)
    with pytest.raises(ValueError, match="Name a snapshot"):
        Dataset.sync(cache=source_cache, sources={"bams": "https://example.test/bams.json"})

    # Online, refresh=True downloads again under the next name for that day.
    fetch = Cache.fetch
    monkeypatch.setattr(Cache, "fetch", lambda self, url, refresh=False, **kw: fetch(self, url, **kw))
    second = Dataset.sync(cache=Cache(source_cache.root), refresh=True)
    assert second.name == f"{day}.2" and second.id == first.id
    named = Dataset.sync("named", cache=source_cache)

    rows = list(Dataset.snapshots(cache=source_cache))
    assert [r["name"] for r in rows] == ["named", f"{day}.2", day]  # Same download: newest creation first
    assert Dataset.open(cache=source_cache).name == "named"
    assert Dataset.open(day, cache=source_cache).name == day
    assert Dataset.open(date=day, cache=source_cache).name == "named"
    assert Dataset.open(date=day[:7], cache=source_cache).name == "named"
    assert Dataset.open(first.id[:8], cache=source_cache).name == "named"
    assert Dataset.open(named.name, cache=source_cache).name == "named"
    with pytest.raises(FileNotFoundError, match="No snapshot named 'missing'"):
        Dataset.open("missing", cache=source_cache)
    with pytest.raises(FileNotFoundError, match="No snapshot named '2026-01'"):
        Dataset.open("2026-01", cache=source_cache)  # A missing name is never read as a date
    assert Dataset.open(cache=source_cache).cache.offline
    assert list(Dataset.snapshots(cache=Cache(tmp_path / "empty"))) == []

    (source_cache.workspace / "snapshots" / "broken.json").write_text("{")
    with pytest.warns(UserWarning, match="broken.json"):
        assert len(Dataset.snapshots(cache=source_cache)) == 3
    with pytest.warns(UserWarning):
        assert Dataset.open(cache=source_cache).name == "named"


def test_cli_uses_the_newest_snapshot_unless_one_is_chosen(source_cache, capsys, monkeypatch):
    root = str(source_cache.root)
    cli = ["--offline", "--cache", root]
    assert main([*cli, "snapshots"]) == 0
    assert "No saved snapshots" in capsys.readouterr().out
    assert main([*cli, "sync", "--json"]) == 0
    day = json.loads(capsys.readouterr().out)["snapshot"]
    assert main([*cli, "sync", "--source-revision", "0" * 40]) == 1
    assert "Name a snapshot" in capsys.readouterr().err
    assert main([*cli, "sync", "older"]) == 0
    assert capsys.readouterr().out.startswith("Saved snapshot older (")
    assert main([*cli, "snapshots"]) == 0
    listing = capsys.readouterr().out
    assert day in listing and "older" in listing and "Commands use older" in listing
    assert main([*cli, "snapshots", "--json"]) == 0
    assert [r["name"] for r in json.loads(capsys.readouterr().out)] == ["older", day]

    expected = [v.id for v in Dataset.open(cache=source_cache).variants(gene="SMC5")]
    for snapshot in ([], ["--snapshot", day], ["--snapshot", day[:7]], ["--snapshot", "older"]):
        assert main([*cli, "variants", *snapshot, "--gene", "SMC5", "--json"]) == 0
        assert [v["id"] for v in json.loads(capsys.readouterr().out)] == expected
    assert main([*cli, "variants", "--snapshot", "missing"]) == 1
    assert "No snapshot named 'missing'" in capsys.readouterr().err
    assert main([*cli, "variants", "--snapshot", "1999"]) == 1
    assert "No snapshot downloaded in 1999" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        main([*cli, "snapshots", "older"])  # A snapshot is chosen only with --snapshot
    capsys.readouterr()

    assert main([*cli, "samples", "T2_tumor", "--snapshot", day]) == 0
    assert capsys.readouterr().out.startswith("T2_tumor")
    assert main([*cli, "timeline", "--around", "2025-01-28", "--days", "1"]) == 0
    capsys.readouterr()

    from osteosarc import ReadSubset
    calls = []

    def extract_reads(self, asset, regions=None, **kwargs):
        calls.append((self.name, getattr(asset, "key", asset)))
        return ReadSubset(source_cache.root / "reads.bam", source_cache.root / "reads.bam.bai", {})
    monkeypatch.setattr(Dataset, "extract_reads", extract_reads)
    key = Dataset.open(cache=source_cache).files.select(format="bam")[0].key
    for arguments, snapshot in (([key, "--variant", "SMC5-chr9-70298024"], "older"),
                                # A date-shaped value is a download date: the newest that day.
                                ([key, "--variant", "SMC5-chr9-70298024", "--snapshot", day], "older")):
        assert main([*cli, "reads", *arguments]) == 0
        capsys.readouterr()
        assert calls.pop() == (snapshot, key)
    assert main([*cli, "reads", "--snapshot", day, key, "--variant", "MISSING"]) == 1
    assert f"variants --set all --snapshot {day}" in capsys.readouterr().err


def test_file_lists_never_shorten_keys(dataset):
    from osteosarc.views import files_view
    expected = dataset.samples["T1_tumor"].files.select(kind="alignment")
    shown = files_view(expected, dataset, width=40)
    assert shown.startswith(f"{len(expected)} file")
    assert expected[0].key in shown  # Never truncated, even when narrower than the key
    assert "--limit N" not in shown
    assert "--limit N" in files_view(dataset.files.select(kind="alignment"), dataset, limit=1)


def test_file_filters_reject_labels_and_values_no_source_uses(dataset):
    from osteosarc import File, Files
    from osteosarc.models import SampleClaim
    with pytest.raises(ValueError, match="registry label; select assay 'rna-seq'"):
        dataset.files.select(assay="RNA")
    with pytest.raises(ValueError, match="Unknown tissue 'normal'"):
        dataset.files.select(tissue="normal")
    assert len(dataset.files.select(assay="cite-seq")) == 0  # Known but absent: simply empty
    # A value outside the vocabulary that the data does use still selects it.
    marrow = Files([File("x", "x.bam", "https://example.test/x.bam", "alignment", "bam",
                           claims=(SampleClaim("bams", tissue="marrow"),))])
    assert len(marrow.select(tissue="marrow")) == 1


def test_a_first_look_is_readable_in_python(dataset):
    from osteosarc.display import Text
    shown = repr(dataset)
    assert shown.startswith("Osteosarc snapshot fixture (")
    for name in ("data.samples", "data.files", "data.variants()", "data.timeline", "data.download(key)",
                 "data.downloads()"):
        assert name in shown
    summary = dataset.summary()
    assert isinstance(summary, Text) and repr(summary) == str(summary)  # Shown without quotes
    assert "variants:" in summary and "samples:" in summary and "data.samples" in summary

    variants = repr(dataset.variants())
    assert variants.startswith("5 variants, 5 with a ready allele") and "chr14:101980529 G>A" in variants
    files = repr(dataset.files)
    assert files.startswith(f"{len(dataset.files):,} files") and dataset.files[0].key in files
    samples = repr(dataset.samples)
    assert samples.startswith(f"{len(dataset.samples)} samples") and "T1_tumor" in samples
    assert 'data.samples["T1_tumor"]' in samples
    one = repr(dataset.samples["T1_tumor"])
    assert one.startswith("T1_tumor: ") and "BAMs, aligned reads" in one and "data.download(" in one
    timeline = repr(dataset.timeline)
    assert timeline.startswith(f"Timeline: {len(dataset.timeline)} events") and ".render()" in timeline
    assert repr(dataset.timeline.render()) == str(dataset.timeline.render())
    assert repr(dataset.timeline.select(lane="no such lane")) == "Timeline: no events"
    # Single records stay short: the long source records aren't printed.
    assert "annotations=" not in repr(dataset.variants()[0])
    assert "details=" not in repr(dataset.timeline[0]) and "metadata=" not in repr(dataset.files[0])

    # Headers survive a snapshot without a download time; alleles show whenever there's one.
    from types import SimpleNamespace

    from osteosarc import Variant, Variants
    from osteosarc.views import snapshot_line
    stub = SimpleNamespace(name="old", id="a" * 64, downloaded=None)
    assert snapshot_line(stub) == "snapshot old (aaaaaaaaaaaa), download time unknown"
    odd = Variants([Variant("x", "G", "GRCh38", (("chr1", 5, "A" * 20, "dup"),), "non_literal_allele")])
    assert "chr1:5 AAAAAAAAAA..>dup" in repr(odd) and "0 with a ready allele" in repr(odd)


def test_a_first_look_is_guided_on_the_command_line(source_cache, capsys, monkeypatch):
    from osteosarc import Cache
    root = str(source_cache.root)
    monkeypatch.setenv("OSTEOSARC_CACHE", root)
    assert main([]) == 0
    out = capsys.readouterr().out
    assert "No data yet: start with `osteosarc sync`" in out and "Browse:" in out and "samples [SAMPLE]" in out
    with pytest.raises(SystemExit):
        main(["--help"])
    assert "Get data:" in capsys.readouterr().out
    import re

    from osteosarc import NoSnapshotsError
    with pytest.raises(NoSnapshotsError, match=re.escape(f"No saved snapshots in {root}")) as raised:
        Dataset.open(cache=Cache(root))
    assert str(raised.value.root) == root
    # Without a snapshot, commands say where they looked and what to run.
    assert main(["--offline", "--cache", root, "samples"]) == 1
    err = capsys.readouterr().err
    assert f"No snapshot yet in {root}" in err and "osteosarc sync" in err
    # In a terminal, the REPL offers to download one first (served here from the test files).
    import osteosarc.cli as cli
    fetch = Cache.fetch
    monkeypatch.setattr(Cache, "fetch", lambda self, url, refresh=False, **kw: fetch(self, url, **kw))
    monkeypatch.setattr("sys.stdin.isatty", lambda: True)
    opened = []
    monkeypatch.setattr(cli, "repl", opened.append)

    def closed(prompt):  # Ctrl-D at the prompt means no, without a traceback
        raise EOFError
    monkeypatch.setattr("builtins.input", closed)
    assert main(["--cache", root, "repl"]) == 1
    assert "osteosarc sync" in capsys.readouterr().err and not opened
    monkeypatch.setattr("builtins.input", lambda prompt: "")
    assert main(["--cache", root, "repl"]) == 0
    assert "Saved snapshot" in capsys.readouterr().out
    assert len(Dataset.snapshots(cache=source_cache)) == 1 and opened[0].name == Dataset.open(cache=source_cache).name
    assert main([]) == 0 and f"Using snapshot {opened[0].name}" in capsys.readouterr().out


def test_the_repl_starts_python_with_the_data_loaded(dataset, monkeypatch, capsys):
    import builtins
    import code

    import osteosarc.cli as cli
    real_import, seen = builtins.__import__, {}

    def without_ipython(name, *args, **kwargs):
        if name == "IPython":
            raise ImportError(name)
        return real_import(name, *args, **kwargs)
    monkeypatch.setattr(builtins, "__import__", without_ipython)
    monkeypatch.setattr(code, "interact", lambda banner, local, exitmsg: seen.update(banner=banner, local=local))
    cli.repl(dataset)
    assert seen["local"]["data"] is dataset and "data.samples" in seen["banner"]
    assert "help(data)" in seen["banner"]


def test_downloads_are_listed_found_and_placed_under_their_own_names(dataset, tmp_path, capsys):
    from osteosarc import File, Files
    base = "https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/"
    key = "calls/tumor.vcf.gz"
    source = tmp_path / "tumor.vcf.gz"
    source.write_bytes(b"not really a vcf")
    (tmp_path / "tumor.vcf.gz.tbi").write_bytes(b"index")
    vcf = File(stable_id(base + key), key, base + key, "variants", "vcf", size=16,
               index_urls=(base + key + ".tbi",))
    index = File(stable_id(base + key + ".tbi"), key + ".tbi", base + key + ".tbi", "index", "tbi", size=5)
    dataset.files = Files([*dataset.files, vcf, index])
    assert dataset.local_path(vcf) is None and not dataset.downloads()
    dataset.cache.import_file(source, vcf.url)
    dataset.cache.import_file(tmp_path / "tumor.vcf.gz.tbi", index.url)
    cached = dataset.download(key)
    assert dataset.local_path(vcf) == cached and dataset.local_path(index) is not None
    # A local copy whose size doesn't match this snapshot's listing isn't this file.
    assert dataset.local_path(File(vcf.id, key, vcf.url, "variants", "vcf", size=99)) is None
    # to= puts the file and its index in a folder under their own names, without a second copy.
    placed = dataset.download(key, to=tmp_path / "out")
    assert placed == tmp_path / "out" / "tumor.vcf.gz" and placed.read_bytes() == source.read_bytes()
    assert (tmp_path / "out" / "tumor.vcf.gz.tbi").read_bytes() == b"index"
    assert dataset.download(key, to=tmp_path / "out") == placed  # again: nothing changes
    rows = {r["key"]: r for r in dataset.downloads()}
    assert rows[key]["kind"] == "file" and rows[key]["path"] == str(cached)
    assert "vafs" not in rows and not any(k.startswith("https://osteosarc.com") for k in rows)
    # Extracted reads are listed with their source file and regions.
    derived = dataset.cache.workspace / "derived" / "example"
    derived.mkdir(parents=True)
    (derived / "reads.bam").write_bytes(b"bam")
    (derived / "receipt.json").write_text(json.dumps(dict(request=dict(source=vcf.url, regions=[
        dict(contig="chr1", start=9, end=20), dict(contig="chr2", start=0, end=5)]))))
    reads = next(r for r in dataset.downloads() if r["kind"] == "reads")
    assert (reads["key"], reads["regions"]) == (key, "2 regions from chr1:10-20")
    from osteosarc.views import downloads_view, files_view
    shown = downloads_view(list(dataset.downloads()), dataset.cache.root)
    assert "Downloaded files (2," in shown and "Extracted reads (1)" in shown and str(cached) in shown
    assert "yes" in files_view(Files([vcf]), dataset)
    # A fresh open doesn't list these files, but still shows the bucket's downloads, by URL.
    assert main(["--cache", str(dataset.cache.root), "downloads", "--snapshot", "fixture", "--json"]) == 0
    assert {r["key"] for r in json.loads(capsys.readouterr().out)} == {vcf.url, index.url}


def test_text_commands_for_vaccines_tables_and_sync(dataset, capsys):
    root = str(dataset.cache.root)
    assert main(["--cache", root, "vaccines", "--snapshot", "fixture"]) == 0
    text = capsys.readouterr().out
    assert "vaccine targets across" in text and "ELISPOT" in text and "SMC5" in text
    assert main(["--cache", root, "vaccines", "--snapshot", "fixture", "--json"]) == 0
    assert json.loads(capsys.readouterr().out)[0]["gene"]
    assert main(["--cache", root, "files", "--snapshot", "fixture", "--downloaded", "--prefix", "site/"]) == 0
    assert "site/vafs" in capsys.readouterr().out  # the snapshot's own tables are on this computer
    assert main(["--cache", root, "--offline", "sync", "fixture"]) == 0
    assert capsys.readouterr().out.startswith("Saved snapshot fixture (")


def test_review_fixes_for_samples_links_and_offline_downloads(dataset, tmp_path, monkeypatch, capsys):
    import os
    import stat

    from osteosarc import CORRECTIONS, Change, Correction, OfflineError
    from osteosarc.cache import place
    # Samples with disagreements are hashable, and providers are named one way.
    raw = Dataset.open("fixture", cache=dataset.cache, corrections=False)
    assert raw.samples["T2_tumor"].disagreements and len(set(raw.samples)) == len(raw.samples)
    boston = Correction("spelling", "simulate the registry's spelling", (
        Change("specimens", {"sample_id": "T1_tumor"}, set={"vendors_involved": "Boston Gene"}),))
    spelled = Dataset.open("fixture", cache=dataset.cache, corrections=[*CORRECTIONS, boston])
    assert "Boston Gene" not in spelled.samples["T1_tumor"].providers
    # A registry problem leaves files unlinked, with a warning, instead of breaking every file operation.
    broken = Correction("broken", "simulate a registry row without an ID", (
        Change("specimens", {"sample_id": "T3_tumor"}, set={"sample_id": ""}),))
    unlinked = Dataset.open("fixture", cache=dataset.cache, corrections=[*CORRECTIONS, broken])
    with pytest.warns(UserWarning, match="aren't linked to samples"):
        assert len(unlinked.files) == len(dataset.files)
    # A hard-linked copy is read-only, so editing it can't corrupt the cache.
    cached = tmp_path / "object.bam"
    cached.write_bytes(b"bytes")
    placed = place(cached, tmp_path / "out" / "copy.bam")
    assert placed.samefile(cached) and not os.stat(cached).st_mode & (stat.S_IWUSR | stat.S_IWGRP | stat.S_IWOTH)
    assert place(cached, placed) == placed
    # Offline, a download that needs the network says how to allow it.
    key = dataset.files.select(format="bam")[0].key
    with pytest.raises(OfflineError, match="opened offline; reopen it with Dataset.open"):
        dataset.download(key)
    # The REPL opens its snapshot able to download; --offline keeps it offline.
    import osteosarc.cli as cli
    opened = []
    monkeypatch.setattr(cli, "repl", opened.append)
    root = str(dataset.cache.root)
    assert main(["--cache", root, "repl", "--snapshot", "fixture"]) == 0
    assert main(["--cache", root, "--offline", "repl", "--snapshot", "fixture"]) == 0
    assert [d.cache.offline for d in opened] == [False, True]
    assert main(["--cache", root, "--offline", "download", "--snapshot", "fixture", key, "--refresh"]) == 1
    assert "Cannot refresh in offline mode" in capsys.readouterr().err


def test_sample_hints_promise_only_what_download_does(dataset):
    from osteosarc.views import get_data_hints
    rows = [dict(key="a.bam", size="1 GB", indexed=False, local="")]
    text = get_data_hints(dataset, rows, [])
    assert "osteosarc download a.bam --to ." in text and "its index" not in text
    assert "osteosarc reads" not in text  # an unindexed BAM can't be read by region
    folders = [dict(folder="x/y/")]
    assert "aws s3 cp --recursive --no-sign-request s3://sid-sijbrandij-osteosarc-dataset/x/y/ y/" in \
        get_data_hints(dataset, [], folders)


def test_reads_for_a_sample_land_in_a_folder_with_readable_names(dataset, bam, tmp_path, monkeypatch, capsys):
    import osteosarc.reads
    from osteosarc import ReadSubset
    subset = ReadSubset(bam, Path(str(bam) + ".bai"), {"records": 7})
    monkeypatch.setattr(osteosarc.reads, "extract_reads", lambda *a, **k: subset)
    monkeypatch.setattr(osteosarc.reads, "require_samtools", lambda **k: None)
    source = dataset.samples["T1_tumor"].files.select(kind="alignment")[0]
    monkeypatch.setattr(Dataset, "_download", lambda self, file, **k: bam)
    variant = dataset.variants(status="ready")[0]
    one = (v for v in dataset.variants(ids=variant.id))  # any iterable of variants will do
    placed = dataset.extract_reads(source, variants=one, to=tmp_path / "tests")
    stem = dataset.short_name(source)
    assert placed.path == tmp_path / "tests" / f"{stem}.{variant.id}.bam" and placed.path.read_bytes() == bam.read_bytes()
    assert placed.index_path.name == f"{stem}.{variant.id}.bam.bai" and placed.receipt == subset.receipt
    # From the command line, a sample ID means each of its indexed BAMs, filtered by assay.
    root = str(dataset.cache.root)
    assert main(["--cache", root, "reads", "--snapshot", "fixture", "T1_tumor", "--assay", "rna-seq",
                 "--variant", variant.id, "--to", str(tmp_path / "cli")]) == 0
    assert capsys.readouterr().out.strip() == str(tmp_path / "cli" / f"{stem}.{variant.id}.bam")
    assert main(["--cache", root, "reads", "--snapshot", "fixture", "T1_tumor", "--assay", "wgs",
                 "--variant", variant.id]) == 1
    assert "no indexed BAMs of that kind" in capsys.readouterr().err
    assert main(["--cache", root, "reads", "--snapshot", "fixture", source.key, "--assay", "rna-seq",
                 "--variant", variant.id]) == 1
    assert "give a sample ID" in capsys.readouterr().err
    assert main(["--cache", root, "reads", "--snapshot", "fixture", "T1_tumor", "--index", "x.bai",
                 "--variant", variant.id]) == 1
    assert "--index belongs to one BAM" in capsys.readouterr().err
    # A sample's JSON is always a list, even with one BAM.
    assert main(["--cache", root, "reads", "--snapshot", "fixture", "T1_tumor", "--variant", variant.id,
                 "--json"]) == 0
    assert isinstance(json.loads(capsys.readouterr().out), list)
    # Different reads never silently replace a file already in the folder.
    (tmp_path / "tests" / "other.bam").write_bytes(b"not these reads")
    with pytest.raises(FileExistsError, match="different contents"):
        dataset.extract_reads(source, variants=dataset.variants(ids=variant.id), to=tmp_path / "tests",
                              name="other")
    # A sample whose BAMs are all skipped (here, the only one) fails rather than printing nothing.
    from osteosarc import CoordinateError

    def wrong_build(*args, **kwargs):
        raise CoordinateError("the BAM is GRCh37")
    monkeypatch.setattr(Dataset, "extract_reads", wrong_build)
    assert main(["--cache", root, "reads", "--snapshot", "fixture", "T1_tumor", "--variant", variant.id]) == 1
    assert "GRCh37" in capsys.readouterr().err


def test_short_names_tell_same_named_files_apart(dataset):
    from osteosarc import File, Files
    keys = ["kamil/blood/output/Pool_1/outs/possorted_genome_bam.bam",
            "kamil/blood/output/Pool_2/outs/possorted_genome_bam.bam",
            "rna-seq/tempus/TL/RNA/x_sorted.bam", "vendor/tempus/TL/RNA/x_sorted.bam", "solo/unique.bam"]
    dataset.files = Files([File(k, k, "https://example.test/" + k, "alignment", "bam") for k in keys])
    assert [dataset.short_name(k) for k in keys] == [
        "Pool_1.possorted_genome_bam", "Pool_2.possorted_genome_bam", "rna-seq.x_sorted", "vendor.x_sorted", "unique"]


def test_a_sample_whose_bams_are_all_skipped_is_an_error(dataset, monkeypatch, capsys):
    import osteosarc.cli as cli
    from osteosarc import CoordinateError
    bams = list(dataset.files.select(kind="alignment"))[:2]
    monkeypatch.setattr(cli, "read_sources", lambda data, args: bams)

    def wrong_build(*args, **kwargs):
        raise CoordinateError("the BAM is GRCh37")
    monkeypatch.setattr(Dataset, "extract_reads", wrong_build)
    variant = dataset.variants(status="ready")[0].id
    assert main(["--cache", str(dataset.cache.root), "reads", "--snapshot", "fixture", "T1_tumor",
                 "--variant", variant]) == 1
    err = capsys.readouterr().err
    assert err.count("skipping") == 2 and "every BAM of T1_tumor was skipped" in err
