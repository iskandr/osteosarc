from collections import Counter

import pysam
import pytest

from osteosarc import Cache, IntegrityError, RecoveryPolicy, Region, extract_reads, record_multiset


@pytest.fixture
def split_bam(tmp_path):
    path = tmp_path / "split.bam"
    header = dict(HD=dict(VN="1.6", SO="coordinate"),
                  SQ=[dict(SN="chr1", LN=248956422), dict(SN="chr2", LN=242193529)],
                  RG=[dict(ID="a", SM="sample", LB="library"), dict(ID="b", SM="different")])
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        def add(name, contig, pos, flag=0, rg="a", cigar="10M", sa=None, seq=True, mate=None):
            r = pysam.AlignedSegment(out.header)
            r.query_name, r.flag, r.reference_id, r.reference_start = name, flag, contig, pos
            r.mapping_quality, r.cigarstring = 60, cigar
            if seq:
                r.query_sequence = "ACGTACGTAA"
                r.query_qualities = [30] * 10
            r.set_tag("RG", rg)
            r.set_tag("CB", "shared-cell")
            r.set_tag("UB", "shared-umi")
            r.set_tag("NM", 0)
            if sa:
                r.set_tag("SA", sa)
            if mate:
                r.next_reference_id, r.next_reference_start = mate
            out.write(r)
        add("split", 0, 100, 65, sa="chr2,501,+,10H10M,60,0;", mate=(1, 800))
        add("split", 0, 100, 65, sa="chr2,501,+,10H10M,60,0;", mate=(1, 800))
        add("competitor", 0, 110)
        add("missing", 0, 120, sa="chr2,901,+,10M,60,0;")
        add("conflict", 0, 125, sa="chr2,701,+,5M5S,60,0;")
        add("split", 1, 500, 2113, cigar="10H10M", sa="chr1,101,+,10M,60,0;", mate=(1, 800))
        add("split", 1, 500, 2113, rg="b", cigar="10H10M")
        add("conflict", 1, 700, 2048)
        add("split", 1, 800, 129, seq=False, mate=(0, 100))
    pysam.index(str(path))
    return path


def test_bounded_recovery_keeps_original_records_duplicates_and_context(split_bam, tmp_path):
    regions = [Region("chr1", 100, 140, "GRCh38"), Region("chr1", 110, 150, "GRCh38")]
    cache = Cache(tmp_path / "cache", offline=True)
    regional = extract_reads(split_bam, regions, cache=cache)
    recovered = extract_reads(split_bam, regions, cache=cache, recovery=RecoveryPolicy())
    assert regional.receipt["records"] == 5
    assert recovered.receipt["records"] == 7
    assert record_multiset(recovered.path) <= record_multiset(split_bam)
    with recovered.open() as bam:
        reads = list(bam)
    assert Counter(r.query_name for r in reads)["split"] == 4
    assert all(r.get_tag("RG") == "a" for r in reads)
    assert any(r.cigarstring == "10H10M" for r in reads)
    assert any(r.query_sequence is None for r in reads)
    assert recovered.receipt["complete_template"] is False
    assert recovered.receipt["repeated_or_cyclic_leads"]
    assert len(recovered.receipt["unresolved"]) == 2
    assert any(o["reason"] == "missing SEQ" for o in recovered.receipt["observations"])
    assert len(recovered.receipt["visited_intervals"]) > 1


def test_cap_is_distinct_from_empty_and_cache_reuse_needs_no_commands(split_bam, tmp_path, monkeypatch):
    cache = Cache(tmp_path / "cache", offline=True)
    regions = [Region("chr1", 100, 150, "GRCh38")]
    policy = RecoveryPolicy(max_intervals=1)
    result = extract_reads(split_bam, regions, cache=cache, recovery=policy)
    assert result.receipt["status"] == "truncated"
    assert result.receipt["limits"] == ["max_intervals"]
    assert result.receipt["records"] == 5
    import osteosarc.reads
    monkeypatch.setattr(osteosarc.reads, "_run", lambda *a: pytest.fail("cached recovery ran a command"))
    assert extract_reads(split_bam, regions, cache=cache, recovery=policy).receipt == result.receipt


def test_record_limit_fails_instead_of_publishing_partial_success(split_bam, tmp_path):
    with pytest.raises(IntegrityError, match="limit"):
        extract_reads(split_bam, [Region("chr1", 100, 150, "GRCh38")], cache=tmp_path / "cache",
                      recovery=RecoveryPolicy(max_records=2))


def test_changed_source_across_acquisitions_fails(split_bam, tmp_path, monkeypatch):
    import osteosarc.recovery as recovery
    original = recovery.extract_reads
    calls = []

    def changed(*args, **kwargs):
        result = original(*args, **kwargs)
        calls.append(result)
        if len(calls) > 1:
            result.receipt["request"]["source_sha256"] = "changed-source"
        return result

    monkeypatch.setattr(recovery, "extract_reads", changed)
    with pytest.raises(IntegrityError, match="changed across"):
        extract_reads(split_bam, [Region("chr1", 100, 150, "GRCh38")], cache=tmp_path / "cache",
                      recovery=RecoveryPolicy())


def test_interruption_reuses_only_verified_acquisitions(split_bam, tmp_path, monkeypatch):
    import osteosarc.recovery as recovery
    original = recovery.extract_reads
    calls = []

    def interrupted(*args, **kwargs):
        calls.append(args)
        if len(calls) == 2:
            raise OSError("interrupted transfer")
        return original(*args, **kwargs)

    monkeypatch.setattr(recovery, "extract_reads", interrupted)
    regions = [Region("chr1", 100, 150, "GRCh38")]
    with pytest.raises(OSError, match="interrupted"):
        extract_reads(split_bam, regions, cache=tmp_path / "cache", recovery=RecoveryPolicy())
    monkeypatch.setattr(recovery, "extract_reads", original)
    result = extract_reads(split_bam, regions, cache=tmp_path / "cache", recovery=RecoveryPolicy())
    assert result.receipt["records"] == 7
    assert record_multiset(result.path) <= record_multiset(split_bam)


def test_missing_index_never_falls_back_to_scan(split_bam, tmp_path):
    from pathlib import Path
    Path(str(split_bam) + ".bai").unlink()
    with pytest.raises(ValueError, match="No known index"):
        extract_reads(split_bam, [Region("chr1", 100, 150, "GRCh38")], cache=tmp_path / "cache",
                      recovery=RecoveryPolicy())


def test_recipe_partner_retention_has_effect(split_bam, tmp_path):
    from osteosarc import select_fixtures
    subset = extract_reads(split_bam, [Region("chr1", 100, 150, "GRCh38")], cache=tmp_path / "cache",
                           recovery=RecoveryPolicy())
    recipe = dict(schema_version=1, id="split", targets={"sv": dict(
        kind="sv", assembly="GRCh38", reference={"assembly": "GRCh38"}, coordinates="zero-based-interbase",
        breakends=[dict(contig="chr1", position=110, orientation="+"), dict(contig="chr2", position=500, orientation="+")])},
        sources={"rna": dict(identity={"id": "original"}, assembly="GRCh38", sample=None, library=None, product=None)},
        members={"split": dict(source="rna", target="sv", policy=dict(kind="regional", version=1),
                               regions=[dict(contig="chr1", start=100, end=110, assembly="GRCh38")])})
    seed = select_fixtures(recipe, {"rna": subset}).members["split"]
    recipe["members"]["split"]["retain_partners"] = True
    linked = select_fixtures(recipe, {"rna": subset}).members["split"]
    assert seed["record_count"] == 2
    assert linked["record_count"] == 4
    assert linked["acquisition_status"] == "bounded"
    assert any("SA-linked partner" in reasons for reasons in linked["reasons"].values())
