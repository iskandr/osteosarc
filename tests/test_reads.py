from collections import Counter
from dataclasses import replace

import pytest

from osteosarc import (
    Cache,
    CoordinateError,
    IntegrityError,
    ReadFilter,
    Region,
    extract_reads,
    inspect_alignment,
    resolve_regions,
    subset_templates,
)


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


def test_remote_extraction_refuses_changed_header_source(bam, tmp_path, monkeypatch):
    import osteosarc.reads as reads
    real_run = reads._run
    url = "https://example.test/alignment.bam"
    identity = {"etag": '"original"', "content-length": str(bam.stat().st_size), "last-modified": None}
    monkeypatch.setattr(reads, "_remote_identity", lambda *a: dict(identity))

    def remote_header(command, timeout):
        assert command[:4] == ["samtools", "view", "--no-PG", "-H"] or command in (["samtools", "--version"], ["samtools", "view", "--help"])
        return real_run([str(bam) if c == url else c for c in command], timeout)

    monkeypatch.setattr(reads, "_run", remote_header)
    cache = Cache(tmp_path / "cache")
    info = inspect_alignment(url, cache=cache, snapshot_id="snapshot")
    assert info.assembly == "GRCh38"
    offline = Cache(cache.root, offline=True)
    assert inspect_alignment(url, cache=offline, snapshot_id="snapshot").path == info.path
    identity["etag"] = '"replacement"'
    with pytest.raises(IntegrityError, match="since header inspection"):
        extract_reads(url, [Region("chr1", 100, 160, "GRCh38")], cache=cache,
                      index=str(bam) + ".bai", snapshot_id="snapshot")


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


def test_template_sampling_keeps_full_record_multiplicity_and_read_groups(bam, tmp_path):
    cache = Cache(tmp_path / "cache")
    fixture = subset_templates(bam, 2, cache=cache, seed="example")
    selected = {tuple(t) for t in fixture.receipt["selected_templates"]}
    import pysam
    with pysam.AlignmentFile(bam) as handle:
        expected = [r.to_string() for r in handle if (r.get_tag("RG"), r.query_name) in selected]
    assert Counter(records(fixture.path)) == Counter(expected)
    assert fixture.receipt["templates"] == 2
    assert fixture.receipt["suitable_for_vaf"] is False
    assert subset_templates(bam, 2, cache=cache, seed="example").path == fixture.path


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
])
def test_missing_samtools_options_fail_before_any_acquisition(dataset, monkeypatch, help_text, options, missing):
    from subprocess import CompletedProcess

    import osteosarc.reads as reads
    from osteosarc import OsteosarcError

    dataset.cache = Cache(dataset.cache.root)
    source = dataset.assets.select(format="bam")[0]
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
    from osteosarc import Asset, Variant, Variants
    from osteosarc.cache import stable_id

    # A local file URL exercises the real Dataset path without remote data.
    source = Asset(stable_id(str(bam)), "toy.bam", str(bam), "alignment", "bam")
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
