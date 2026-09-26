"""Allele-balanced selection and required-record matching for the shared test data."""

import pytest

from osteosarc import IntegrityError, validate_recipe
from osteosarc.alleles import allele_window
from osteosarc.records import FixtureRecord
from osteosarc.shared import match_required, sam_span, select_allele_balanced, template_order

pysam = pytest.importorskip("pysam")

REFERENCE = "GATTACAGCTAAAAAGCGTCAGTCAGTTGACCATGGCATG"
HEADER = pysam.AlignmentHeader.from_dict(dict(SQ=[dict(SN="chr1", LN=1000)]))


def record(name, base, *, quality=30, start=2):
    read = pysam.AlignedSegment(HEADER)
    read.query_name, read.flag, read.reference_id, read.reference_start = name, 0, 0, start
    sequence = REFERENCE[start:7] + base + REFERENCE[8:start + 10]
    read.cigarstring, read.query_sequence = f"{len(sequence)}M", sequence
    read.query_qualities = pysam.qualitystring_to_array(chr(33 + quality) * len(sequence))
    return FixtureRecord(read, name + base)


def test_selection_is_balanced_capped_and_ordered_by_hash():
    window = allele_window("chr1", 8, "G", "T", REFERENCE, 0)
    templates = {}
    for i in range(30):
        templates[(None, f"alt{i}")] = [record(f"alt{i}", "T", quality=10 + i)]
        templates[(None, f"ref{i}")] = [record(f"ref{i}", "G")]
    templates[(None, "other")] = [record("other", "A")]
    templates[(None, "far")] = [record("far", "G", start=20)]  # doesn't reach the window
    chosen, seen = select_allele_balanced(templates, window, caps=dict(alt=3, ref=2), low_quality_alt=2)
    assert seen == dict(alt=30, ref=30, other=1, uncallable=0)
    by_class = {}
    for template, (name, reason) in chosen.items():
        by_class.setdefault(name, []).append((template, reason))
    hashed = sorted((t for t in templates if t[1].startswith("alt")), key=template_order)
    assert [t for t, why in by_class["alt"] if "hash" in why] == hashed[:3]
    # The two lowest-quality alt templates the hash order didn't already pick.
    low = sorted((t for t in hashed[3:]), key=lambda t: int(t[1][3:]))[:2]
    assert sorted(t for t, why in by_class["alt"] if "quality" in why) == sorted(low)
    assert len(by_class["ref"]) == 2 and len(by_class["other"]) == 1 and (None, "far") not in chosen


def test_required_records_match_by_sam_text_with_multiplicity():
    a, b = record("a", "T"), record("b", "G")
    duplicate = FixtureRecord(a.read, a.digest)
    counts = match_required([a, duplicate, b], [a.read.to_string()] * 2 + [b.read.to_string()])
    assert counts == {a.digest: 2, b.digest: 1}
    with pytest.raises(IntegrityError, match="required records"):
        match_required([a, b], [a.read.to_string()] * 2)
    assert sam_span(a.read.to_string()) == ("chr1", 2, 12)
    assert sam_span("q\t4\t*\t0\t0\t*\t*\t0\t0\tA\tI") is None


def test_recipes_accept_library_fixture_targets():
    recipe = dict(schema_version=1, id="x", targets=dict(f=dict(kind="fixture", assembly="GRCh38",
                  reference=dict(source="library"), description="Isovar's six-locus fixture")),
                  sources=dict(s=dict(identity=dict(id="s"), assembly="GRCh38", sample=None, library=None,
                                      product=None)),
                  members=dict(m=dict(target="f", source="s", policy=dict(version=1, kind="exact", records={}))))
    validate_recipe(recipe)
    del recipe["targets"]["f"]["description"]
    with pytest.raises(Exception, match="description"):
        validate_recipe(recipe)


def test_fusion_templates_must_reach_every_breakend_window():
    from osteosarc.shared import select_breakend_templates
    header = pysam.AlignmentHeader.from_dict(dict(SQ=[dict(SN="chr1", LN=10000), dict(SN="chr2", LN=10000)]))

    def aligned(name, contig, start, flag=0):
        read = pysam.AlignedSegment(header)
        read.query_name, read.flag, read.reference_id, read.reference_start = name, flag, contig, start
        read.cigarstring, read.query_sequence = "10M", "ACGTACGTAC"
        return FixtureRecord(read, f"{name}{contig}{start}")
    templates = {
        (None, "split"): [aligned("split", 0, 100), aligned("split", 1, 500, flag=2048)],
        (None, "pair"): [aligned("pair", 0, 150, flag=65), aligned("pair", 1, 520, flag=129)],
        (None, "one-sided"): [aligned("one-sided", 0, 120)],
    }
    windows = [("chr1", 50, 250), ("chr2", 400, 600)]
    chosen, joined = select_breakend_templates(templates, windows, cap=5)
    assert joined == 2 and set(chosen) == {(None, "split"), (None, "pair")}
    capped, _ = select_breakend_templates(templates, windows, cap=1)
    assert list(capped) == sorted([(None, "split"), (None, "pair")], key=template_order)[:1]


def test_required_reads_can_be_named_instead_of_listed():
    from osteosarc.shared import match_named, required_spans
    a, b = record("a", "T"), record("b", "G")
    assert match_named([a, b, FixtureRecord(a.read, a.digest)], {"a"}) == {a.digest: 2}
    with pytest.raises(IntegrityError, match="required reads"):
        match_named([b], {"a"})
    assert required_spans(dict(names=["a"], regions=[["chr1", 0, 50]])) == [("chr1", 0, 50)]
    assert required_spans(dict(sam=[a.read.to_string()])) == [("chr1", 2, 12)]


def test_required_files_need_sam_or_named_reads(tmp_path):
    import json

    from osteosarc.errors import SchemaError
    from osteosarc.shared import load_required
    path = tmp_path / "required.json"
    path.write_text(json.dumps(dict(consumer="x", subsets={"x/named": dict(source="s", names=["r"])})))
    with pytest.raises(SchemaError, match="names with regions"):
        load_required([path])


def test_the_record_index_finds_what_a_scan_would():
    from osteosarc.shared import RecordIndex
    records = [record(f"r{i}", "T", start=i) for i in range(0, 30, 3)]
    index = RecordIndex(records)
    for start, end in ((0, 5), (7, 9), (25, 40), (100, 200)):
        scanned = [r for r in records if r.read.reference_start < end and r.read.reference_end > start]
        assert index.overlapping("chr1", start, end) == scanned
    assert set(index.touching("chr1", 7, 9)) == {r.template for r in records if r.read.reference_start < 9}
    # Required lines match through the index as they do against the whole list.
    lines = [records[2].read.to_string(), records[5].read.to_string()]
    assert match_required(index, lines) == match_required(records, lines)
