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

    def aligned(name, contig, start, flag=0, cigar="10M"):
        read = pysam.AlignedSegment(header)
        read.query_name, read.flag, read.reference_id, read.reference_start = name, flag, contig, start
        read.cigarstring = cigar
        read.query_sequence = "A" * read.infer_query_length()
        return FixtureRecord(read, f"{name}{contig}{start}")
    templates = {
        (None, "split"): [aligned("split", 0, 100), aligned("split", 1, 500, flag=2048)],
        (None, "pair"): [aligned("pair", 0, 150, flag=65), aligned("pair", 1, 520, flag=129)],
        (None, "one-sided"): [aligned("one-sided", 0, 120)],
    }
    breakends = [("chr1", 150), ("chr2", 500)]
    chosen, joined = select_breakend_templates(templates, breakends, pad=100, cap=5)
    assert joined == 2 and set(chosen) == {(None, "split"), (None, "pair")}
    capped, _ = select_breakend_templates(templates, breakends, pad=100, cap=1)
    assert list(capped) == sorted([(None, "split"), (None, "pair")], key=template_order)[:1]

    # Breakends on one contig: reads that merely run across them, or are spliced from
    # an exon near one to an exon near the other, don't join them; a junction does.
    close = [("chr1", 1100), ("chr1", 1400)]
    ordinary = {
        (None, "spans"): [aligned("spans", 0, 1150, cigar="200M")],
        (None, "spliced"): [aligned("spliced", 0, 1150, cigar="5M200N5M")],
        (None, "intron-over-both"): [aligned("intron-over-both", 0, 900, cigar="5M700N5M")],
        (None, "deletion-elsewhere"): [aligned("deletion-elsewhere", 0, 1150, cigar="5M200D5M")],
        (None, "proper"): [aligned("proper", 0, 1100, flag=67), aligned("proper", 0, 1400, flag=131)],
    }
    junctions = {
        (None, "spliced-junction"): [aligned("spliced-junction", 0, 1050, cigar="50M300N50M")],
        (None, "deletion"): [aligned("deletion", 0, 1052, cigar="50M296D50M")],
        (None, "discordant"): [aligned("discordant", 0, 1100, flag=65), aligned("discordant", 0, 1400, flag=129)],
        (None, "supplementary"): [aligned("supplementary", 0, 1100), aligned("supplementary", 0, 1400, flag=2048)],
    }
    chosen, joined = select_breakend_templates({**ordinary, **junctions}, close, pad=100, cap=10)
    assert set(chosen) == set(junctions) and joined == 4


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


def test_every_shipped_spec_has_a_published_bundle():
    from osteosarc.cache import stable_id
    from osteosarc.shared import BUNDLES, published, read_json
    specs = sorted(BUNDLES.glob("*.spec.json"))
    assert specs
    for path in specs:
        name = path.name.removesuffix(".spec.json")
        spec = read_json(path)
        assert spec["id"] == name
        release = published(name)
        assert release["spec_sha256"] == stable_id(spec), f"rebuild {name}: its spec changed since its release"
        assert release["url"].endswith(f"/releases/download/{name}/{name}.tar.gz")
        assert len(release["sha256"]) == len(release["manifest_sha256"]) == 64 and release["size_bytes"] > 0


def test_member_names_must_be_file_names_and_unique():
    from osteosarc.shared import _add_member
    members = {}
    _add_member(members, "T2_rna.MAP2-chr2-209694768-historical", {})
    with pytest.raises(IntegrityError, match="share the name"):
        _add_member(members, "T2_rna.MAP2-chr2-209694768-historical", {})
    with pytest.raises(IntegrityError, match="Unsafe"):
        _add_member(members, "T2_rna.MAP2-chr2-209694768:historical", {})


def test_fixtures_kept_in_json_are_read_in_each_layout(tmp_path):
    import json

    from osteosarc import SchemaError
    from osteosarc.shared import _local_sam
    line = "r1\t0\tchr1\t3\t60\t4M\t*\t0\t0\tACGT\tIIII"
    partner = "r1\t2048\tchr2\t9\t60\t4M\t*\t0\t0\tACGT\tIIII"
    (tmp_path / "f.json").write_text(json.dumps({
        "lines": [line, line], "single": line, "by_digest": {"a" * 64: line},
        # Isovar keeps each split read with its partner; field names don't matter.
        "pairs": [dict(sam=line, partner_sam=partner, read_id="r1"), dict(sam=line, partner_sam=None)],
        "lost_sam": [dict(sam=None, partner_sam=partner)],
        "records_by_digest": {"b" * 64: dict(sam=line, partner_sam=partner)},
        "other_field": [dict(read=line, id="r1")],
        "with_metadata": dict(version="1.0", source="https://example.test/x.bam", records=[line]),
        "joined": line + "\n" + partner + "\r\n", "a/b": [line], "numbers": [1, 2, None]}))
    keys = [("/lines", [line, line]), ("/single", [line]), ("/by_digest", [line]),
            ("/pairs", [line, partner, line]), ("/lost_sam", [partner]), ("/records_by_digest", [line, partner]),
            ("/other_field", [line]), ("/with_metadata", [line]), ("/joined", [line, partner]),
            ("/a~1b", [line]), ("/numbers", [])]
    for pointer, expected in keys + [("/pairs/0", [line, partner])]:
        assert _local_sam(dict(json="f.json", pointer=pointer), tmp_path) == expected, pointer
    whole = _local_sam(dict(json="f.json", pointer=""), tmp_path)
    assert whole == _local_sam(dict(json="f.json", pointer="/"), tmp_path) == _local_sam(dict(json="f.json"), tmp_path)
    assert whole == [x for _, expected in keys for x in expected]
    for pointer, where in [("/record", "/record"), ("/lines/2", "/lines/2"), ("/single/x", "/single/x"),
                           ("/a/b", "/a"), ("lines", None)]:
        with pytest.raises(SchemaError, match="nothing at " + where if where else "must start with"):
            _local_sam(dict(json="f.json", pointer=pointer), tmp_path)


def _dync1h1_bam(path):
    """A local BAM at DYNC1H1-chr14-101980529 (G>A): 3 alt and 2 ref reads, and reference context."""
    position = 101980529
    start = position - 1 - 20
    reference = {i: "ACGT"[(i * 7) % 4] for i in range(start - 200, start + 260)}
    reference[position - 1] = "G"
    header = dict(HD={"VN": "1.6", "SO": "coordinate"},
                  SQ=[dict(SN="chr1", LN=248956422), dict(SN="chr2", LN=242193529), dict(SN="chr14", LN=107043718)])
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for name, base in [("alt1", "A"), ("alt2", "A"), ("alt3", "A"), ("ref1", "G"), ("ref2", "G")]:
            read = pysam.AlignedSegment(out.header)
            sequence = "".join(base if i == position - 1 else reference[i] for i in range(start, start + 40))
            read.query_name, read.flag, read.reference_id, read.reference_start = name, 0, 2, start
            read.cigarstring, read.query_sequence = "40M", sequence
            read.query_qualities = pysam.qualitystring_to_array("I" * 40)
            read.mapping_quality = 60
            out.write(read)
    pysam.index(str(path))
    return lambda contig, lo, hi, assembly, **kwargs: "".join(reference.get(i, "N") for i in range(lo, hi))


def test_a_bundle_is_made_from_variants_and_files_and_read_back(dataset, tmp_path, monkeypatch):
    import os
    import shutil

    import osteosarc.shared as shared
    from osteosarc import bundle_file, extract_reads, inspect_alignment, list_bundle
    from osteosarc.records import record_multiset
    if shutil.which("samtools") is None:
        pytest.skip("samtools required")
    local = tmp_path / "rna.bam"
    monkeypatch.setattr(shared, "reference_sequence", _dync1h1_bam(local))
    # The snapshot's BAM, read from the local copy.
    monkeypatch.setattr(dataset, "inspect_alignment", lambda file, **kw: inspect_alignment(str(local), cache=dataset.cache))
    monkeypatch.setattr(dataset, "extract_reads", lambda file, regions, **kw: extract_reads(
        str(local), regions, cache=dataset.cache, fetch_pairs=kw.get("fetch_pairs", False)))
    key = "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam"
    folder = dataset.make_bundle(tmp_path / "dync1h1", variants=["DYNC1H1-chr14-101980529"], files=[key])
    member = "BG003082.Aligned.sortedByCoord.out.md.DYNC1H1-chr14-101980529"
    assert list(list_bundle(folder)) == [member]
    # A test reads the member as a local, read-only BAM, exported once and then reused offline.
    bam = bundle_file(folder, member, cache=tmp_path / "cache")
    assert sum(record_multiset(bam).values()) == 5 and not os.access(bam, os.W_OK)
    assert bundle_file(folder, member, cache=tmp_path / "cache", offline=True) == bam
    with pytest.raises(KeyError, match="did you mean"):
        bundle_file(folder, member[:-1], cache=tmp_path / "cache")
    with pytest.raises(FileExistsError, match="new folder"):
        dataset.make_bundle(folder, variants=["DYNC1H1-chr14-101980529"], files=[key])


def test_a_bundle_spec_says_what_it_needs(dataset):
    from osteosarc.shared import bundle_spec
    key = "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam"
    spec = bundle_spec(dataset, "mine", variants=dataset.variants(ids=["DYNC1H1-chr14-101980529"]),
                       files=dataset.file(key), svs="SV0461", caps=dict(alt=3))
    assert spec["targets"]["ids"] == ["DYNC1H1-chr14-101980529"]
    assert spec["targets"]["structural"] == [dict(name="SV0461", label="sv", **{"from": {"sv_candidates": "SV0461"}})]
    assert spec["sources"] == dict(all_targets=[dataset.file(key).url], structural=[dataset.file(key).url],
                                   observed=False)
    assert spec["selection"]["caps"]["alt"] == 3 and spec["selection"]["caps"]["ref"] == 10
    assert bundle_spec(dataset, "x", svs=["GABBR1-SLC29A1"], files=[key])["targets"]["structural"][0]["from"] == {
        "panel": "sv-regressions-v1", "id": "GABBR1-SLC29A1"}
    with pytest.raises(ValueError, match="did you mean DYNC1H1-chr14-101980529"):
        bundle_spec(dataset, "x", variants=["DYNC1H1-chr14-101980528"], files=[key])
    with pytest.raises(ValueError, match="No SV"):
        bundle_spec(dataset, "x", svs=["SV9999"], files=[key])
    with pytest.raises(ValueError, match="variants"):
        bundle_spec(dataset, "x", files=[key])
    with pytest.raises(ValueError, match="BAMs"):
        bundle_spec(dataset, "x", variants=["DYNC1H1-chr14-101980529"])


def test_the_cli_makes_a_bundle_from_a_sample_and_variants(dataset, tmp_path, monkeypatch, capsys):
    import osteosarc.cli as cli
    import osteosarc.shared as shared
    from osteosarc import File
    made = {}

    def make_bundle(data, to, **kwargs):
        made.update(kwargs, to=to)
        raise shared.OfflineError("stop here")
    monkeypatch.setattr(shared, "make_bundle", make_bundle)
    monkeypatch.setattr(cli, "open_snapshot", lambda args, cache, online: dataset)
    assert cli.main(["--offline", "test-data", "make", str(tmp_path / "b"), "rna-seq/reprocessed/BG003082/"
                     "BG003082.Aligned.sortedByCoord.out.md.bam", "--variant", "DYNC1H1-chr14-101980529",
                     "--sv", "SV0461"]) == 1
    assert "stop here" in capsys.readouterr().err
    assert made["variants"] == ["DYNC1H1-chr14-101980529"] and made["svs"] == ["SV0461"]
    assert [f.key for f in made["files"]] == ["rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam"]
    assert all(isinstance(f, File) for f in made["files"])
