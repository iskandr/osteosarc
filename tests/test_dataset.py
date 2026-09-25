import json
from collections import Counter

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
    assets = raw.assets.select(kind="alignment")
    assert len(assets) == 3
    source = assets.select(contains="BG009368")[0]
    assert set(source.values("timepoint")) == {"T0", "T1"}
    assert source.resolved("timepoint") is None
    assert not assets.select(contains="BG009368", timepoint="T0")
    assert len(assets.select(contains="BG009368", timepoint="T1", include_conflicts=True)) == 1
    assert source.metadata["catalog_genome_assertion"] == "hg38"
    assert "assembly" not in source.metadata
    assert len(dataset.samples) > len(assets)
    # The viewer label's stale date is surfaced alongside its timepoint.
    assert source.conflicts["date"] == ("2022-12", "2024-06")
    # The central correction for that stale label resolves it, traceably.
    corrected = dataset.asset(source.key)
    assert corrected.resolved("timepoint") == "T1" and not corrected.conflicts
    assert corrected.metadata["corrections"] == ("viewer-label-BG009368",)


def test_catalog_normalization_does_not_invent_conflicts():
    from osteosarc import SampleClaim
    from osteosarc.catalog import BUCKET, build_assets
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
    assets = {a.key: a for a in build_assets(listing, bams, metadata, vafs, paths)}
    assert assets[cite].values("assay") == ("cite-seq",)
    assert assets[cite].values("provider") == ("Hudson Lab",)
    assert assets[cite].values("date") == ("2025-11-06",)  # month precision agrees
    assert not assets[cite].conflicts
    assert assets[rna].values("assay") == ("rna-seq",)
    assert assets["vendor/cegat/P2/P2.1.fastq.gz"].values("assay") == ("wes",)


def test_reopen_is_offline_and_uses_pinned_receipts(dataset):
    reopened = Dataset.open("fixture", cache=dataset.cache)
    assert reopened.cache.offline
    assert reopened.id == dataset.id
    assert reopened.variants().to_records() == dataset.variants().to_records()
    path = reopened.download("vafs")
    assert len(reopened.table("vafs")) > 0
    assert all(isinstance(r["alt_reads"], str) for r in reopened.table("vafs"))
    assert Counter(v.status for v in reopened.variants()) == Counter(v.status for v in dataset.variants())
    path.write_text("corrupt")
    with pytest.raises(IntegrityError):
        Dataset.open("fixture", cache=dataset.cache)


def test_snapshot_and_cli_use_same_selection(dataset, capsys):
    assert main(["--cache", str(dataset.cache.root), "variants", "--snapshot", "fixture", "--gene", "SMC5"]) == 0
    assert json.loads(capsys.readouterr().out)[0]["id"] == dataset.variants(gene="SMC5")[0].id
    assert main(["--cache", str(dataset.cache.root), "assets", "--snapshot", "fixture", "--kind", "alignment", "--limit", "1", "--json"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert result["total"] == 3
    assert len(result["assets"]) == 1
    with pytest.raises(FileExistsError):
        Dataset.sync("fixture", cache=dataset.cache, refresh=True)


def test_cli_version_and_sample_assets_match_the_api(dataset, capsys):
    import osteosarc
    with pytest.raises(SystemExit) as exited:
        main(["--version"])
    assert exited.value.code == 0
    assert capsys.readouterr().out.strip() == f"osteosarc {osteosarc.__version__}"
    root = str(dataset.cache.root)
    assert main(["--cache", root, "assets", "--snapshot", "fixture", "--sample", "T1_tumor", "--kind", "alignment",
                 "--json"]) == 0
    result = json.loads(capsys.readouterr().out)
    expected = dataset.assets_for_sample("T1_tumor", kind="alignment")
    assert result["total"] == len(expected) == 1
    assert [a["key"] for a in result["assets"]] == [a.key for a in expected]
    assert main(["--cache", root, "assets", "--snapshot", "fixture", "--sample", "T1_tumor", "--assay", "wgs",
                 "--json"]) == 0
    assert json.loads(capsys.readouterr().out)["total"] == 0
    # Without --json, a readable table with the complete key.
    assert main(["--cache", root, "assets", "--snapshot", "fixture", "--sample", "T1_tumor", "--kind", "alignment"]) == 0
    text = capsys.readouterr().out
    assert text.startswith("1 files") and expected[0].key in text and "rna-seq" in text
    assert main(["--cache", root, "assets", "--snapshot", "fixture", "--kind", "alignment", "--limit", "1"]) == 0
    assert "... 2 more" in (text := capsys.readouterr().out) and "--limit N" in text
    # Registry labels and unknown values are errors, not empty selections.
    assert main(["--cache", root, "assets", "--snapshot", "fixture", "--assay", "scRNA_ONT"]) == 1
    assert "select assay 'scrna-seq' with platform 'ont'" in capsys.readouterr().err
    assert main(["--cache", root, "assets", "--snapshot", "fixture", "--platform", "nanopore"]) == 1
    assert "choose from: illumina, ont, pacbio" in capsys.readouterr().err
    assert main(["--cache", root, "assets", "--snapshot", "fixture", "--sample", "T9_tumor"]) == 1
    assert "Unknown sample" in capsys.readouterr().err


def test_cli_reads_accepts_catalogue_variants(dataset, capsys, monkeypatch, tmp_path):
    from osteosarc import ReadSubset
    calls = []

    def extract_reads(self, asset, regions=None, **kwargs):
        calls.append(dict(asset=asset, regions=regions, **kwargs))
        return ReadSubset(tmp_path / "reads.bam", tmp_path / "reads.bam.bai", {"records": 0})
    monkeypatch.setattr(Dataset, "extract_reads", extract_reads)
    root, source = str(dataset.cache.root), dataset.assets.select(format="bam")[0].key
    ids = ["DYNC1H1-chr14-101980529", "SMC5-chr9-70298024"]
    command = ["--cache", root, "reads", "--snapshot", "fixture", source, "--variant", ids[0], "--variant", ids[1]]
    assert main(command + ["--padding", "100"]) == 0
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
        main(["--cache", root, "variants", "--snapshot", "fixture", "extra"])


def test_acquired_data_remains_pinned_after_url_refresh(dataset, tmp_path):
    from osteosarc import Asset
    url = "https://example.test/results.tsv"
    path = tmp_path / "results.tsv"
    path.write_text("id\tcount\na\t0\n")
    asset = Asset(stable_id(url), "results.tsv", url, "table", "tsv")
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
    asset = dataset.assets.select(format="bam")[0]
    monkeypatch.setattr(dataset, "download", lambda *a: pytest.fail("downloaded unsupported format"))
    with pytest.raises(ValueError, match="No built-in parser"):
        dataset.parse(asset)


def test_vcf_annotations_and_indexed_subsetting(dataset, tmp_path):
    import pysam

    from osteosarc import Asset, Assets
    source = tmp_path / "calls.vcf"
    source.write_text('##fileformat=VCFv4.2\n'
                      '##contig=<ID=chr1,length=1000>\n'
                      '##INFO=<ID=GENE,Number=1,Type=String,Description="gene">\n'
                      '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ttumor\n'
                      'chr1\t11\tfirst\tA\tC,G\t.\tPASS\tGENE=ONE\tGT\t1/2\n'
                      'chr1\t50\tsecond\tT\tG\t.\tPASS\tGENE=TWO\tGT\t0/1\n')
    compressed = tmp_path / "calls.vcf.gz"
    pysam.tabix_compress(str(source), str(compressed))
    pysam.tabix_index(str(compressed), preset="vcf")
    url = "https://example.test/calls.vcf.gz"
    index_url = url + ".tbi"
    asset = Asset(stable_id(url), "calls.vcf.gz", url, "variants", "vcf", index_urls=(index_url,))
    index = Asset(stable_id(index_url), "calls.vcf.gz.tbi", index_url, "index", "tbi")
    dataset.assets = Assets([asset, index])
    dataset.cache.import_file(compressed, url)
    dataset.cache.import_file(str(compressed) + ".tbi", index_url)
    with dataset.open_variants(asset) as calls:
        rows = list(calls.fetch("chr1", 10, 11))
    assert len(rows) == 1
    assert rows[0].alts == ("C", "G")
    assert rows[0].info["GENE"] == "ONE"
    assert rows[0].samples["tumor"]["GT"] == (1, 2)


def test_first_download_must_match_the_inventory_time():
    from osteosarc import Asset, Receipt
    from osteosarc.dataset import _check_inventory_time
    asset = Asset("id", "a.tsv", "https://example.test/a.tsv", "table", "tsv", modified=1784586023)
    same = Receipt("https://example.test/a.tsv", "0" * 64, 1, "a.tsv", "now",
                   last_modified="Mon, 20 Jul 2026 22:20:23 GMT")
    _check_inventory_time(asset, same)
    newer = Receipt(**{**same.to_dict(), "last_modified": "Tue, 21 Jul 2026 09:00:00 GMT"})
    with pytest.raises(IntegrityError, match="new snapshot"):
        _check_inventory_time(asset, newer)


def test_extraction_binds_the_listed_index_to_the_snapshot(dataset, monkeypatch):
    import osteosarc.reads
    from osteosarc import Region
    source = dataset.asset("rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam")
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

    from osteosarc import SNAPSHOT_SOURCES, TIMELINE_SOURCES, Cache
    mirror = "https://mirror.example.test/variant_vafs_long.tsv"
    cache = Cache(tmp_path / "cache", offline=True)
    urls = {**SNAPSHOT_SOURCES, **TIMELINE_SOURCES, "vafs": mirror}
    for key, name in FILES.items():
        cache.import_file(DATA / name, urls[key])
    data = Dataset.sync("mirrored", cache=cache, sources={"vafs": mirror})
    assert data.asset("vafs").url == mirror
    assert len(data.table("vafs")) == len(data.vafs) > 0   # offline: the pinned copy, not a download


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
    with pytest.raises(FileNotFoundError, match="osteosarc sync"):
        choose_snapshot([])
    ambiguous = [dict(rows[0], id="abcdef" + "0" * 58), dict(rows[1], id="abcdef" + "1" * 58)]
    with pytest.raises(FileNotFoundError, match="ambiguous"):
        choose_snapshot(ambiguous, "abcdef")
    assert downloaded_at(dict(sources={}, created_at="2026-09-24T00:00:00+00:00")) == "2026-09-24T00:00:00+00:00"


@pytest.fixture
def source_cache(tmp_path):
    """An offline cache holding every fixture source, with no snapshot yet."""
    from conftest import DATA, FILES

    from osteosarc import SNAPSHOT_SOURCES, TIMELINE_SOURCES, Cache
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
    assert main([*cli, "sync"]) == 0
    day = json.loads(capsys.readouterr().out)["snapshot"]
    assert main([*cli, "sync", "--source-revision", "0" * 40]) == 1
    assert "Name a snapshot" in capsys.readouterr().err
    assert main([*cli, "sync", "older"]) == 0
    capsys.readouterr()
    assert main([*cli, "snapshots"]) == 0
    listing = capsys.readouterr().out
    assert day in listing and "older" in listing and "Commands use older" in listing
    assert main([*cli, "snapshots", "--json"]) == 0
    assert [r["name"] for r in json.loads(capsys.readouterr().out)] == ["older", day]

    expected = [v.id for v in Dataset.open(cache=source_cache).variants(gene="SMC5")]
    for snapshot in ([], ["--snapshot", day], ["--snapshot", day[:7]], ["--snapshot", "older"]):
        assert main([*cli, "variants", *snapshot, "--gene", "SMC5"]) == 0
        assert [v["id"] for v in json.loads(capsys.readouterr().out)] == expected
    assert main([*cli, "variants", "--snapshot", "missing"]) == 1
    assert "No snapshot named 'missing'" in capsys.readouterr().err
    assert main([*cli, "variants", "--snapshot", "1999"]) == 1
    assert "No snapshot downloaded in 1999" in capsys.readouterr().err
    with pytest.raises(SystemExit):
        main([*cli, "variants", "older"])  # A snapshot is chosen only with --snapshot
    capsys.readouterr()

    assert main([*cli, "specimens", "T2_tumor", "--snapshot", day]) == 0
    assert capsys.readouterr().out.startswith("T2_tumor")
    assert main([*cli, "on", "2025-01-28", "--days", "1"]) == 0
    capsys.readouterr()

    from osteosarc import ReadSubset
    calls = []

    def extract_reads(self, asset, regions=None, **kwargs):
        calls.append((self.name, asset))
        return ReadSubset(source_cache.root / "reads.bam", source_cache.root / "reads.bam.bai", {})
    monkeypatch.setattr(Dataset, "extract_reads", extract_reads)
    key = Dataset.open(cache=source_cache).assets.select(format="bam")[0].key
    for arguments, snapshot in (([key, "--variant", "SMC5-chr9-70298024"], "older"),
                                # A date-shaped value is a download date: the newest that day.
                                ([key, "--variant", "SMC5-chr9-70298024", "--snapshot", day], "older")):
        assert main([*cli, "reads", *arguments]) == 0
        capsys.readouterr()
        assert calls.pop() == (snapshot, key)
    assert main([*cli, "reads", "--snapshot", day, key, "--variant", "MISSING"]) == 1
    assert f"variants --set all --snapshot {day}" in capsys.readouterr().err


def test_explorer_lists_a_samples_files_with_complete_keys(dataset):
    from osteosarc.explore import assets_view
    expected = dataset.assets_for_sample("T1_tumor", kind="alignment")
    shown = assets_view(dataset, sample="T1_tumor", kind="alignment", width=40)
    assert shown.startswith(f"{len(expected)} files")
    assert expected[0].key in shown  # Never truncated, even when narrower than the key
    assert "limit=N" not in shown
    assert "limit=N" in assets_view(dataset, kind="alignment", limit=1)


def test_asset_filters_reject_labels_and_values_no_source_uses(dataset):
    from osteosarc import Asset, Assets, SampleClaim
    with pytest.raises(ValueError, match="registry label; select assay 'rna-seq'"):
        dataset.assets.select(assay="RNA")
    with pytest.raises(ValueError, match="Unknown tissue 'normal'"):
        dataset.assets.select(tissue="normal")
    assert len(dataset.assets.select(assay="cite-seq")) == 0  # Known but absent: simply empty
    # A value outside the vocabulary that the data does use still selects it.
    marrow = Assets([Asset("x", "x.bam", "https://example.test/x.bam", "alignment", "bam",
                           claims=(SampleClaim("bams", tissue="marrow"),))])
    assert len(marrow.select(tissue="marrow")) == 1


def test_a_first_look_is_readable_in_python(dataset, monkeypatch):
    from osteosarc.display import Text
    assert repr(dataset).startswith("Osteosarc snapshot fixture (") and "data.summary()" in repr(dataset)
    summary = dataset.summary()
    assert isinstance(summary, Text) and repr(summary) == str(summary)  # Shown without quotes
    assert "variants:" in summary and "data.explore()" in summary
    assert repr(dataset.describe_samples()) == str(dataset.describe_samples())

    variants = repr(dataset.variants())
    assert variants.startswith("5 variants, 5 ready") and "chr14:101980529 G>A" in variants
    assets = repr(dataset.assets)
    assert assets.startswith(f"{len(dataset.assets):,} files") and dataset.assets[0].key in assets
    table = repr(dataset.specimens)
    assert table.startswith(f"Table: {len(dataset.specimens)} rows") and "sample_id" in table
    timeline = repr(dataset.timeline)
    assert timeline.startswith(f"Timeline: {len(dataset.timeline)} events") and ".render()" in timeline
    assert repr(dataset.timeline.render()) == str(dataset.timeline.render())
    assert repr(dataset.timeline.select(lane="no such lane")) == "Timeline: no events"
    # Single records stay short: the long source records aren't printed.
    assert "annotations=" not in repr(dataset.variants()[0])
    assert "details=" not in repr(dataset.timeline[0]) and "metadata=" not in repr(dataset.assets[0])

    import osteosarc.explore as explore
    opened = []
    monkeypatch.setattr(explore.Explorer, "cmdloop", lambda self: opened.append(self.data))
    dataset.explore()
    assert opened == [dataset]


def test_a_first_look_is_guided_on_the_command_line(source_cache, capsys, monkeypatch):
    from osteosarc import Cache
    assert main([]) == 0
    assert "osteosarc sync" in capsys.readouterr().out
    root = str(source_cache.root)
    with pytest.raises(FileNotFoundError, match=f"No saved snapshots in {root}"):
        Dataset.open(cache=Cache(root))
    # Without a snapshot, commands say where they looked and what to run.
    assert main(["--offline", "--cache", root, "samples"]) == 1
    err = capsys.readouterr().err
    assert f"No snapshot yet in {root}" in err and "osteosarc sync" in err
    # In a terminal, the explorer offers to download one first (served here from the test files).
    import osteosarc.explore as explore
    fetch = Cache.fetch
    monkeypatch.setattr(Cache, "fetch", lambda self, url, refresh=False, **kw: fetch(self, url, **kw))
    monkeypatch.setattr("sys.stdin.isatty", lambda: True)
    monkeypatch.setattr("builtins.input", lambda prompt: "")
    monkeypatch.setattr(explore.Explorer, "cmdloop", lambda self: None)
    assert main(["--cache", root, "explore"]) == 0
    assert len(Dataset.snapshots(cache=source_cache)) == 1
