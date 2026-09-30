from collections import Counter
from dataclasses import replace

import pytest

from osteosarc import Cache, CoordinateError, IntegrityError, ReadFilter, Region, extract_reads
from osteosarc.reads import inspect_alignment, resolve_regions


def records(path):
    import pysam
    with pysam.AlignmentFile(path) as bam:
        return [r.to_string() for r in bam]


def test_indexed_union_preserves_original_records_and_tags(bam, tmp_path):
    cache = Cache(tmp_path / "cache", offline=True)
    subset = extract_reads(bam, [Region("1", 100, 140, "GRCh38"), Region("chr1", 115, 160, "hg38")], cache=cache)
    expected = [r for r in records(bam) if not r.startswith("outside\t")]
    assert Counter(records(subset.path)) == Counter(expected)
    assert len(expected) == subset.receipt["records"] == 6
    assert subset.receipt["scope"] == "regional_records"
    with subset.open() as alignments:
        assert sum(1 for _ in alignments.fetch("chr1", 100, 160)) == 6
        for read in alignments.fetch("chr1", 100, 160):
            assert read.get_tag("UB") == "original-umi"


def test_filters_have_effect_and_are_part_of_cache_identity(bam, tmp_path):
    cache = Cache(tmp_path / "cache")
    regions = [Region("chr1", 100, 160, "GRCh38")]
    unfiltered = extract_reads(bam, regions, cache=cache)
    filtered = extract_reads(bam, regions, cache=cache, filters=ReadFilter(exclude_flags=256 | 1024 | 2048))
    assert len(records(unfiltered.path)) == 6
    assert len(records(filtered.path)) == 3
    assert filtered.path != unfiltered.path
    barcode = extract_reads(bam, regions, cache=cache, filters=ReadFilter(barcodes=("B",)))
    assert len(records(barcode.path)) == 2
    assert all("CB:Z:B" in r for r in records(barcode.path))


def test_empty_results_are_valid_but_empty_requests_raise(bam, tmp_path):
    cache = Cache(tmp_path / "cache")
    with pytest.raises(CoordinateError):
        extract_reads(bam, [], cache=cache)
    subset = extract_reads(bam, [Region("chr2", 10, 20, "GRCh38")], cache=cache)
    assert subset.receipt["records"] == 0
    assert records(subset.path) == []


def test_cached_reads_verified_and_reused_offline(bam, tmp_path, monkeypatch):
    cache = Cache(tmp_path / "cache")
    regions = [Region("chr1", 100, 160, "GRCh38")]
    subset = extract_reads(bam, regions, cache=cache)
    import osteosarc.reads
    monkeypatch.setattr(osteosarc.reads, "_run", lambda *a: pytest.fail("cache reuse executed a command"))
    assert extract_reads(bam, regions, cache=Cache(cache.root, offline=True)).path == subset.path
    subset.path.write_bytes(b"modified")
    with pytest.raises(IntegrityError):
        extract_reads(bam, regions, cache=cache)


def test_too_many_records_is_remembered_too(bam, tmp_path, monkeypatch):
    import osteosarc.reads
    from osteosarc.errors import RecordLimitError
    cache = Cache(tmp_path / "cache")
    regions = [Region("chr1", 100, 160, "GRCh38")]
    with pytest.raises(RecordLimitError):
        extract_reads(bam, regions, cache=cache, max_records=5)  # it holds 6
    # Asked again, even offline, the answer comes from the cache without reading.
    monkeypatch.setattr(osteosarc.reads, "_run_bounded", lambda *a: pytest.fail("read the records again"))
    with pytest.raises(RecordLimitError):
        extract_reads(bam, regions, cache=Cache(cache.root, offline=True), max_records=5)
    monkeypatch.undo()
    assert extract_reads(bam, regions, cache=cache, max_records=6).receipt["records"] == 6


def test_cached_header_reuse_and_corruption(bam, tmp_path, monkeypatch):
    cache = Cache(tmp_path / "cache", offline=True)
    info = inspect_alignment(bam, cache=cache)
    assert info.assembly == "GRCh38"
    import osteosarc.reads
    monkeypatch.setattr(osteosarc.reads, "_run", lambda *a: pytest.fail("cached header executed a command"))
    assert inspect_alignment(bam, cache=cache).header == info.header
    info.path.write_text("corrupted")
    with pytest.raises(IntegrityError):
        inspect_alignment(bam, cache=cache)


@pytest.mark.parametrize("limit", [None, 6])
def test_regional_counts_preserve_multiplicity_without_storing_bams(bam, tmp_path, limit):
    from osteosarc import count_reads
    cache = Cache(tmp_path / "cache", offline=True)
    regions = [Region("1", 100, 140, "GRCh38"), Region("chr1", 115, 160, "GRCh38")]
    count = count_reads(bam, regions, cache=cache, max_records=limit)
    assert count.records == 6  # overlapping regions count once; identical records remain twice
    assert count.receipt["scope"] == "regional_record_count"
    assert count.receipt["request"]["operation"] == "count_reads"
    assert "count.json" in count.receipt["files"]
    assert not list((cache.workspace / "derived").rglob("*.bam*"))
    assert count_reads(bam, regions, cache=cache, filters=ReadFilter(exclude_flags=256 | 1024 | 2048)).records == 3
    assert count_reads(bam, [Region("chr2", 10, 20, "GRCh38")], cache=cache).records == 0
    with pytest.raises(CoordinateError):
        count_reads(bam, [], cache=cache)


def test_count_cache_and_overflow_are_reusable_and_verified(bam, tmp_path, monkeypatch):
    import json

    import osteosarc.reads as reads
    from osteosarc import RecordLimitError, count_reads
    cache = Cache(tmp_path / "cache")
    regions = [Region("chr1", 100, 160, "GRCh38")]
    count = count_reads(bam, regions, cache=cache, max_records=6)
    with pytest.raises(RecordLimitError):
        count_reads(bam, regions, cache=cache, max_records=5)
    assert not list((cache.workspace / "derived").rglob("*.bam*"))
    monkeypatch.setattr(reads, "_run_bounded", lambda *a: pytest.fail("repeated the count query"))
    offline = Cache(cache.root, offline=True)
    assert count_reads(bam, regions, cache=offline, max_records=6) == count
    with pytest.raises(RecordLimitError):
        count_reads(bam, regions, cache=offline, max_records=5)
    receipt_path, = (cache.workspace / "derived").glob("*/receipt.json")
    receipt = json.loads(receipt_path.read_text())
    receipt["records"] += 1
    receipt_path.write_text(json.dumps(receipt))
    with pytest.raises(IntegrityError, match="count differs"):
        count_reads(bam, regions, cache=offline, max_records=6)


def test_counts_reuse_older_extractions_and_overflow_results(bam, tmp_path, monkeypatch):
    import osteosarc.reads as reads
    from osteosarc import RecordLimitError, count_reads
    cache = Cache(tmp_path / "cache")
    regions = [Region("chr1", 100, 160, "GRCh38")]
    subset = extract_reads(bam, regions, cache=cache, max_records=6)
    with pytest.raises(RecordLimitError):
        extract_reads(bam, regions, cache=cache, max_records=5)
    monkeypatch.setattr(reads, "_run_bounded", lambda *a: pytest.fail("repeated the old probe"))
    monkeypatch.setattr(reads, "_run", lambda *a: pytest.fail("started a process for cached reads"))
    count = count_reads(bam, regions, cache=Cache(cache.root, offline=True), max_records=6)
    assert count.records == 6 and count.receipt == subset.receipt
    with pytest.raises(RecordLimitError):
        count_reads(bam, regions, cache=cache, max_records=5)
    assert len(list((cache.workspace / "derived").glob("*/receipt.json"))) == 1


def test_dataset_counts_share_extraction_identity(bam, dataset):
    from osteosarc import File
    file = File("local", "source.bam", str(bam), "alignment", "bam")
    regions = [Region("chr1", 100, 160, "GRCh38")]
    subset = dataset.extract_reads(file, regions, max_records=6)
    count = dataset.count_reads(file, regions, max_records=6)
    assert count.records == 6 and count.receipt == subset.receipt


@pytest.mark.parametrize("records", [5, 7, -1, True, "6", None])
def test_counts_reject_incorrect_legacy_record_metadata(bam, tmp_path, records):
    import json

    from osteosarc import count_reads
    cache = Cache(tmp_path / "cache")
    regions = [Region("chr1", 100, 160, "GRCh38")]
    subset = extract_reads(bam, regions, cache=cache, max_records=6)
    receipt = dict(subset.receipt, records=records)
    subset.receipt_path.write_text(json.dumps(receipt))
    with pytest.raises(IntegrityError, match="BAM count differs"):
        count_reads(bam, regions, cache=Cache(cache.root, offline=True), max_records=6)


@pytest.mark.parametrize("unpinned", ["reads.bam", "reads.bam.bai"])
def test_counts_require_legacy_bam_and_index_integrity_pins(bam, tmp_path, unpinned):
    import json

    from osteosarc import count_reads
    cache = Cache(tmp_path / "cache")
    regions = [Region("chr1", 100, 160, "GRCh38")]
    subset = extract_reads(bam, regions, cache=cache, max_records=6)
    receipt = subset.receipt
    del receipt["files"][unpinned]
    subset.receipt_path.write_text(json.dumps(receipt))
    with pytest.raises(IntegrityError, match="missing from its receipt"):
        count_reads(bam, regions, cache=Cache(cache.root, offline=True), max_records=6)


@pytest.mark.parametrize("legacy", [False, True])
def test_counts_reuse_compact_receipts_and_overflow_offline(bam, tmp_path, monkeypatch, legacy):
    import osteosarc.reads as reads
    from osteosarc import RecordLimitError, count_reads
    from osteosarc.read_receipts import compact_read_cache
    cache = Cache(tmp_path / "cache")
    regions = [Region("chr1", 100, 160, "GRCh38")]
    filters = ReadFilter(query_names=("repeated", "supplementary"))
    acquire = extract_reads if legacy else count_reads
    acquire(bam, regions, cache=cache, filters=filters, max_records=4)
    with pytest.raises(RecordLimitError):
        acquire(bam, regions, cache=cache, filters=filters, max_records=3)
    compact_read_cache(cache)
    monkeypatch.setattr(reads, "_run", lambda *a: pytest.fail("cached query started SAMtools"))
    monkeypatch.setattr(reads, "_run_bounded", lambda *a: pytest.fail("cached query acquired records"))
    offline = Cache(cache.root, offline=True)
    count = count_reads(bam, regions, cache=offline, filters=filters, max_records=4)
    assert count.records == 4
    with pytest.raises(RecordLimitError):
        count_reads(bam, regions, cache=offline, filters=filters, max_records=3)
    asset, = (cache.workspace / "query-names").glob("*.txt")
    asset.write_text("changed\n")
    with pytest.raises(IntegrityError, match="query-name asset"):
        count_reads(bam, regions, cache=offline, filters=filters, max_records=4)


@pytest.mark.parametrize("limit", [5, 6])
def test_count_queries_do_not_cache_results_from_changed_remote_sources(bam, tmp_path, monkeypatch, limit):
    import osteosarc.reads as reads
    from osteosarc import RecordLimitError, count_reads
    url = "https://example.test/alignment.bam"
    identity = {"etag": '"original"', "content-length": str(bam.stat().st_size), "last-modified": None}
    run, bounded = reads._run, reads._run_bounded
    monkeypatch.setattr(reads, "_remote_identity", lambda *a, **kw: dict(identity))
    monkeypatch.setattr(reads, "_run", lambda command, timeout: run(
        [str(bam) if c == url else c for c in command], timeout))

    def changed(command, output, maximum, timeout):
        try:
            return bounded([str(bam) if c == url else c for c in command], output, maximum, timeout)
        finally:
            identity["etag"] = '"replacement"'

    monkeypatch.setattr(reads, "_run_bounded", changed)
    cache = Cache(tmp_path / "cache")
    with pytest.raises(RecordLimitError if limit == 5 else IntegrityError):
        count_reads(url, [Region("chr1", 100, 160, "GRCh38")], cache=cache,
                    index=str(bam) + ".bai", snapshot_id="snapshot", max_records=limit)
    assert not list((cache.workspace / "derived").glob("*/receipt.json"))
    assert not list((cache.workspace / "derived").glob("*.over-limit.json"))


def test_remote_extraction_refuses_changed_header_source(bam, tmp_path, monkeypatch):
    import subprocess

    import osteosarc.reads as reads
    real_run = reads._run
    url = "https://example.test/alignment.bam"
    identity = {"etag": '"original"', "content-length": str(bam.stat().st_size), "last-modified": None}
    monkeypatch.setattr(reads, "_remote_identity", lambda *a, **k: dict(identity))

    def remote(command, timeout):  # the local BAM stands in for the remote one
        return real_run([str(bam) if c == url else c for c in command], timeout)

    monkeypatch.setattr(reads, "_run", remote)
    cache = Cache(tmp_path / "cache")
    info = inspect_alignment(url, cache=cache, snapshot_id="snapshot")
    assert info.assembly == "GRCh38"
    offline = Cache(cache.root, offline=True)
    assert inspect_alignment(url, cache=offline, snapshot_id="snapshot").path == info.path
    identity["etag"] = '"replacement"'
    with pytest.raises(IntegrityError, match="since header inspection"):
        extract_reads(url, [Region("chr1", 100, 160, "GRCh38")], cache=cache,
                      index=str(bam) + ".bai", snapshot_id="snapshot")
    assert not list((cache.workspace / "derived").glob("*/receipt.json"))  # nothing kept

    # An object replaced while samtools reads it is what's reported when the read fails.
    identity["etag"] = '"original"'

    def replaced_while_read(command, timeout):
        if command[:4] == ["samtools", "view", "--no-PG", "-b"]:
            identity["etag"] = '"replacement"'
            raise subprocess.CalledProcessError(1, command, stderr=b"[E::bgzf_read] Read block operation failed")
        return remote(command, timeout)
    monkeypatch.setattr(reads, "_run", replaced_while_read)
    with pytest.raises(IntegrityError, match="changed during extraction"):
        extract_reads(url, [Region("chr1", 100, 170, "GRCh38")], cache=cache,
                      index=str(bam) + ".bai", snapshot_id="snapshot")

    # Too many records is remembered only while the object is still the one read.
    from osteosarc.errors import RecordLimitError
    monkeypatch.setattr(reads, "_run", remote)

    def too_many(*args, changed):
        if changed:
            identity["etag"] = '"replacement"'
        raise RecordLimitError("Acquisition exceeds record limit")
    for changed in (True, False):
        identity["etag"] = '"original"'
        monkeypatch.setattr(reads, "_run_bounded", lambda *a, changed=changed: too_many(*a, changed=changed))
        with pytest.raises(RecordLimitError):
            extract_reads(url, [Region("chr1", 100, 180, "GRCh38")], cache=cache,
                          index=str(bam) + ".bai", snapshot_id="snapshot", max_records=1)
        assert len(list((cache.workspace / "derived").glob("*.over-limit.json"))) == (0 if changed else 1)


def test_samtools_version_ignores_non_utf8_distribution_build_flags(monkeypatch):
    from subprocess import CompletedProcess

    import osteosarc.reads as reads
    monkeypatch.setattr(reads, "_run", lambda *a: CompletedProcess([], 0, b"samtools 1.19.2\nflags: \xab\n"))
    assert reads._samtools_version() == "samtools 1.19.2"


def test_assembly_mismatch_and_missing_index_fail_before_query(bam, tmp_path):
    with pytest.raises(CoordinateError, match="header establishes GRCh38"):
        extract_reads(bam, [Region("chr1", 1, 2, "GRCh37")], cache=tmp_path / "cache")
    (bam.parent / (bam.name + ".bai")).unlink()
    with pytest.raises(ValueError, match="index"):
        extract_reads(bam, [Region("chr1", 1, 2, "GRCh38")], cache=tmp_path / "cache")


def test_mitochondrial_lengths_and_ambiguous_contigs():
    header = {"SQ": [dict(SN="1", LN=249250621), dict(SN="2", LN=243199373),
                     dict(SN="MT", LN=16569)]}
    region = Region("chrM", 1, 2, "GRCh37")
    with pytest.raises(CoordinateError, match="reference_length"):
        resolve_regions([region], header)
    with pytest.raises(CoordinateError, match="length mismatch"):
        resolve_regions([replace(region, reference_length=16571)], header)
    assert resolve_regions([replace(region, reference_length=16569)], header)[0].contig == "MT"
    header["SQ"].append(dict(SN="M", LN=16569))
    with pytest.raises(CoordinateError, match="ambiguous"):
        resolve_regions([replace(region, reference_length=16569)], header)


def test_cram_requires_explicit_reference_and_decodes_regional_records(bam, tmp_path):
    import pysam
    # The small reference covers the queried decoy only; unused primary contigs
    # establish assembly without creating a half-gigabyte reference fixture.
    reference = tmp_path / "reference.fa"
    reference.write_text(">decoy\n" + "ACGT" * 500 + "\n")
    pysam.faidx(str(reference))
    with pysam.AlignmentFile(bam) as handle:
        header = handle.header.to_dict()
    header["SQ"].append(dict(SN="decoy", LN=2000))
    cram = tmp_path / "input.cram"
    with pysam.AlignmentFile(cram, "wc", header=header, reference_filename=str(reference)) as out:
        read = pysam.AlignedSegment(out.header)
        read.query_name = "cram-read"
        read.query_sequence = "ACGT" * 10
        read.query_qualities = pysam.qualitystring_to_array("I" * 40)
        read.reference_id = 2
        read.reference_start = 100
        read.mapping_quality = 60
        read.cigarstring = "40M"
        read.set_tag("CB", "A")
        out.write(read)
    pysam.index(str(cram))
    region = Region("decoy", 100, 140, "GRCh38")
    with pytest.raises(CoordinateError, match="reference FASTA"):
        extract_reads(cram, [region], cache=tmp_path / "cache")
    subset = extract_reads(cram, [region], cache=tmp_path / "cache", reference=reference)
    with subset.open() as handle:
        read = next(handle)
    assert read.query_sequence == "ACGT" * 10
    assert read.get_tag("CB") == "A"


def test_a_single_barcode_string_is_not_split_into_characters():
    from osteosarc import ReadFilter
    assert ReadFilter(barcodes="AAACCTGAGAAACCAT").barcodes == ("AAACCTGAGAAACCAT",)
    assert ReadFilter(barcodes=["A1", "B2"]).barcodes == ("A1", "B2")


@pytest.mark.parametrize("help_text,options,missing", [
    (b"--no-PG -M -X", {"fetch_pairs": True}, "--fetch-pairs"),
    (b"--no-PG -M", {}, "-X"),
    (b"--no-PG -M -X", {"filters": ReadFilter(barcodes=("A",))}, "-D"),
    (b"--no-PG -M -X", {"filters": ReadFilter(query_names=("q",))}, "-N"),
    (b"--no-PG -M -X --fetch-pairs", {"fetch_pairs": True, "unplaced_mates": False}, "-N"),
])
def test_missing_samtools_options_fail_before_any_acquisition(dataset, monkeypatch, help_text, options, missing):
    from subprocess import CompletedProcess

    import osteosarc.reads as reads
    from osteosarc import OsteosarcError

    dataset.cache = Cache(dataset.cache.root)
    source = dataset.files.select(format="bam")[0]
    assert source.index_urls
    def run(command, timeout):
        assert command == ["samtools", "view", "--help"]
        return CompletedProcess(command, 0, help_text, b"")
    monkeypatch.setattr(reads, "_run", run)
    monkeypatch.setattr("osteosarc.cache.http_identity", lambda *a: pytest.fail("contacted source before capability check"))
    with pytest.raises(OsteosarcError, match=missing):
        dataset.extract_reads(source, [Region("chr1", 100, 140, "GRCh38")], **options)


def test_variant_selection_produces_the_same_complete_regional_records(dataset, bam, monkeypatch):
    import osteosarc.reads as reads
    from osteosarc import File, Variant, Variants
    from osteosarc.cache import stable_id

    # A local file URL exercises the real Dataset path without remote data.
    source = File(stable_id(str(bam)), "toy.bam", str(bam), "alignment", "bam")
    selected = Variants([Variant("toy", "GENE", "GRCh38", (("chr1", 106, "A", "C"),), "ready")])
    subset = dataset.extract_reads(source, variants=selected, padding=40)
    expected = extract_reads(source, selected.regions(padding=40), cache=dataset.cache, snapshot_id=dataset.id)
    assert subset.path == expected.path
    assert Counter(records(subset.path)) == Counter(r for r in records(bam) if not r.startswith("outside\t"))
    monkeypatch.setattr(reads, "_run", lambda *a: pytest.fail("cached request needed SAMtools"))
    assert dataset.extract_reads(source, variants=selected, padding=40).path == subset.path
    with pytest.raises(ValueError, match="either"):
        dataset.extract_reads(source, selected.regions(), variants=selected)
    with pytest.raises(CoordinateError, match="nonempty"):
        dataset.extract_reads(source, variants=[])


@pytest.mark.parametrize("names", [("bad name",), ("bad\tname",), ("bad\nname",), ("",), ("@bad",), (123,), ("x" * 255,)])
def test_invalid_query_name_filters(names):
    with pytest.raises(ValueError, match="query name"):
        ReadFilter(query_names=names)


def test_query_name_filter_is_canonical_and_scopes_cache(bam, tmp_path):
    regions = [Region("chr1", 100, 160, "GRCh38")]
    all_reads = extract_reads(bam, regions, cache=tmp_path)
    name = records(all_reads.path)[0].split("\t")[0]
    chosen = extract_reads(bam, regions, cache=tmp_path, filters=ReadFilter(query_names=name))
    assert Counter(records(chosen.path)) == Counter(r for r in records(all_reads.path) if r.split("\t")[0] == name)
    assert chosen.path != all_reads.path
    assert extract_reads(bam, regions, cache=tmp_path, filters=ReadFilter(query_names=[name, name])).path == chosen.path
    assert ReadFilter(query_names=["b", "a", "a"]).query_names == ("a", "b")
    absent = extract_reads(bam, regions, cache=tmp_path, filters=ReadFilter(query_names="absent"))
    assert absent.receipt["records"] == 0


def test_mates_with_no_position_can_be_left_out(bam, tmp_path):
    import pysam
    paired = tmp_path / "paired.bam"
    with pysam.AlignmentFile(bam) as template, pysam.AlignmentFile(paired, "wb", template=template) as output:
        # (name, flag, contig, start, mate contig, mate start); -1 is no position.
        for name, flag, contig, start, mate_contig, mate_start in [
            ("placed", 99, 0, 100, 0, 1000), ("placed", 147, 0, 1000, 0, 100),
            ("other-contig", 65, 0, 130, 1, 500), ("other-contig", 129, 1, 500, 0, 130),
            ("beside", 73, 0, 120, 0, 120), ("beside", 133, 0, 120, 0, 120),  # unmapped, at its mate's place
            ("straddling", 97, 0, 140, 0, 70), ("straddling", 145, 0, 70, 0, 140),  # starts before, runs in
            # Both mates here, and a supplementary where another read's mate is: not a mate.
            ("inside", 99, 0, 105, 0, 125), ("inside", 147, 0, 125, 0, 105), ("inside", 2145, 0, 1000, 0, 105),
            ("single", 0, 0, 150, -1, -1), ("single", 2048, 0, 1000, -1, -1),
            ("unknown", 65, 0, 155, 1, -1),  # its mate's contig, but no position
            # samtools asks for these names too (its mate has no position, or starts at a
            # region's start), so their records at another read's mate position come along.
            ("nowhere", 2121, 0, 1000, -1, -1), ("unknown", 2113, 0, 1000, 1, -1),
            ("at-start", 97, 0, 105, 0, 100), ("at-start", 145, 0, 100, 0, 105), ("at-start", 2145, 0, 1000, 0, 105),
            ("nowhere", 73, 0, 110, -1, -1), ("nowhere", 133, -1, -1, -1, -1),  # unmapped, no position
            ("bystander", 0, 0, 1000, -1, -1),  # at a mate's position, but not asked for
        ]:
            read = pysam.AlignedSegment(output.header)
            read.query_name, read.flag = name, flag
            read.query_sequence = "ACGT" * 10
            read.query_qualities = pysam.qualitystring_to_array("I" * 40)
            read.reference_id, read.reference_start = contig, start
            read.next_reference_id, read.next_reference_start = mate_contig, mate_start
            if not flag & 4:
                read.mapping_quality, read.cigarstring = 60, "40M"
            read.set_tag("RG", "rg1")  # records are kept byte for byte: tags stay in order
            read.set_tag("NH", 1)
            output.write(read)
    pysam.sort("-o", str(tmp_path / "sorted.bam"), str(paired))
    paired = tmp_path / "sorted.bam"
    pysam.index(str(paired))
    cache = Cache(tmp_path / "cache")
    regions = [Region("chr1", 100, 160, "GRCh38")]
    every = extract_reads(paired, regions, cache=cache, fetch_pairs=True)
    placed = extract_reads(paired, regions, cache=cache, fetch_pairs=True, unplaced_mates=False)
    unplaced = [r for r in records(every.path) if r.split("\t")[2] == "*"]
    assert len(unplaced) == 1 and unplaced[0].startswith("nowhere\t")
    assert Counter(records(placed.path)) == Counter(r for r in records(every.path) if r not in unplaced)
    assert {r.split("\t")[0] for r in records(placed.path)} == {"placed", "other-contig", "beside", "nowhere",
                                                                "straddling", "inside", "single", "unknown", "at-start"}
    assert not any(r.startswith(("inside\t2145", "single\t2048")) for r in records(placed.path))
    assert all(any(r.startswith(f"{name}\t{flag}\t") for r in records(placed.path))
               for name, flag in [("nowhere", 2121), ("unknown", 2113), ("at-start", 2145)])
    assert placed.receipt["scope"] == "regional_records_and_placed_mates"
    assert "mates.bed" in placed.receipt["files"]
    # Asking for every mate is the request it always was; leaving some out is another.
    assert "unplaced_mates" not in every.receipt["request"]
    assert placed.receipt["request"]["unplaced_mates"] is False and placed.path != every.path
    assert "--fetch-pairs" not in placed.receipt["command"] and placed.receipt["mates_command"]
    # With every mate already in the regions (or no reads), one read of the BAM does.
    around_every_mate = [Region("chr1", 60, 1100, "GRCh38"), Region("chr2", 490, 540, "GRCh38")]
    for regions, count in [(around_every_mate, 21), ([Region("chr2", 10, 20, "GRCh38")], 0)]:
        subset = extract_reads(paired, regions, cache=cache, fetch_pairs=True, unplaced_mates=False)
        assert subset.receipt["records"] == count and "mates_command" not in subset.receipt


def test_a_region_with_no_assembly_is_in_the_bams_own_build(bam, tmp_path):
    subset = extract_reads(bam, [Region("chr1", 100, 160, None)], cache=tmp_path)
    assert subset.receipt["records"] == 6 and subset.receipt["resolved_regions"][0]["assembly"] == "GRCh38"
    with pytest.raises(CoordinateError):
        Region("chr1", 100, 160, "")
