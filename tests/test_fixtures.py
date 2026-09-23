import copy
import json

import pysam
import pytest

from osteosarc import (
    IntegrityError,
    load_panel,
    read_records,
    record_multiset,
    select_fixture_records,
    select_fixtures,
)
from osteosarc.cli import main


def recipe_for(bam):
    return dict(schema_version=1, id="synthetic-v1", targets={"snv": dict(
        kind="small_variant", assembly="GRCh38", reference={"id": "GRCh38"},
        coordinates="one-based", contig="chr1", position=101, ref="A", alt="C")},
        sources={"rna": dict(identity={"url": "https://example.test/source.bam"},
                             assembly="GRCh38", sample="T0", library="rna", product="original")},
        members={"alt": dict(source="rna", target="snv", policy=dict(
            kind="exact", version=1, records=dict(record_multiset(bam))))})


def test_api_dataset_cli_and_input_order_agree(bam, dataset, tmp_path, capsys):
    recipe = recipe_for(bam)
    recipe["members"]["background"] = dict(source="rna", target="snv", regions=[dict(
        contig="chr1", start=90, end=160, assembly="GRCh38")], policy=dict(kind="regional", version=1, cap=2, seed="pinned"))
    selected = select_fixtures(recipe, {"rna": bam})
    assert selected.manifest == dataset.select_fixtures(recipe, {"rna": bam}).manifest
    path = tmp_path / "recipe.json"
    path.write_text(json.dumps(recipe))
    assert main(["--cache", str(tmp_path / "cache"), "--offline", "fixtures", "select", str(path), "--source", f"rna={bam}"]) == 0
    assert json.loads(capsys.readouterr().out) == selected.manifest
    assert selected.members["alt"]["records"] == record_multiset(bam)
    records = list(read_records(bam))
    policy = dict(kind="stratified", version=1, seed="0", cap=1, assignments=[dict(
        selector=dict(rg=r.template[0], qname=r.template[1]), reason="reference control",
        producer=dict(name="historical-classifier", version="1.2"), required=False) for r in records])
    assert select_fixture_records(records, policy) == select_fixture_records(records[::-1], policy)


def test_scoped_witnesses_duplicates_controls_caps_and_missing(bam):
    recipe = recipe_for(bam)
    recipe["sources"]["other"] = copy.deepcopy(recipe["sources"]["rna"])
    policy = dict(kind="stratified", version=1, cap=0, assignments=[dict(
        selector=dict(rg="rg1", qname="repeated", segment=0), reason="rare alt witness",
        producer=dict(name="isovar", version="historical")), dict(
        selector=dict(rg="rg2", qname="repeated"), reason="reference control",
        producer=dict(name="isovar", version="historical"))])
    recipe["members"]["alt"]["policy"] = policy
    recipe["members"]["separate-source"] = dict(source="other", target="snv", policy=dict(kind="empty", version=1))
    selected = select_fixtures(recipe, dict(rna=bam, other=bam))
    assert selected.members["alt"]["record_count"] == 3
    assert sorted(selected.members["alt"]["records"].values()) == [1, 2]
    assert selected.members["separate-source"]["status"] == "empty"
    policy["assignments"][0]["selector"]["rg"] = "absent"
    with pytest.raises(IntegrityError, match="Missing pinned witness"):
        select_fixtures(recipe, dict(rna=bam, other=bam))


def test_required_exact_multiplicity_and_explicit_states(bam):
    recipe = recipe_for(bam)
    counts = recipe["members"]["alt"]["policy"]["records"]
    counts[next(iter(counts))] += 1
    with pytest.raises(IntegrityError, match="duplicate"):
        select_fixtures(recipe, {"rna": bam})
    recipe["targets"]["pending"] = dict(kind="unresolved", reason="unreviewed breakends")
    recipe["members"] = {"unknown": dict(target="pending", source="rna", policy=dict(kind="empty", version=1)),
                         "omit": dict(target="snv", source="rna", policy=dict(kind="omitted", version=1, reason="size budget"))}
    assert {v["status"] for v in select_fixtures(recipe, {}).members.values()} == {"unresolved", "omitted"}


def test_binary_identity_keeps_float_bits_and_tag_types(bam, tmp_path):
    records = list(read_records(bam))
    original = records[0].read
    path = tmp_path / "float.bam"
    with pysam.AlignmentFile(bam) as inp, pysam.AlignmentFile(path, "wb", header=inp.header) as out:
        # These two float32 values print identically at SAM's limited precision.
        original.set_tag("xf", 1.23456788, value_type="f")
        out.write(original)
        sam_before = original.to_string()
        original.set_tag("xf", 1.23456800, value_type="f")
        assert original.to_string() == sam_before
        out.write(original)
    assert len(record_multiset(path)) == 2
    assert sum(record_multiset(bam).values()) == 7


def test_context_is_retained_without_support_claim(bam):
    recipe = recipe_for(bam)
    member = recipe["members"]["alt"]
    member["policy"]["records"] = {}
    member["context_regions"] = [dict(contig="chr1", start=1000, end=1040, assembly="GRCh38")]
    result = select_fixtures(recipe, {"rna": bam}).members["alt"]
    assert result["record_count"] == 1
    assert list(result["reasons"].values()) == [["assembly context"]]


def test_panels_are_offline_copies():
    panel = load_panel("vaccine-loci-v1")
    assert any(key.startswith("NTF3-") for key in panel)
    panel.clear()
    assert load_panel("vaccine-loci-v1")
    assert load_panel("sv-regressions-v1")["SV0461"]["kind"] == "sv"


def test_streamed_records_and_regions_match_materialized_selection(bam):
    from osteosarc import Region
    records = list(read_records(bam))
    regions = [Region("chr1", 90, 160, "GRCh38")]
    context = [Region("chr1", 1000, 1040, "GRCh38")]
    policy = dict(kind="regional", version=1, cap=2)
    expected = select_fixture_records(records, policy, regions=regions, context_regions=context)
    actual = select_fixture_records(read_records(bam), policy, regions=iter(regions), context_regions=iter(context))
    assert actual == expected
    assert sum(actual[0].values()) > 1


@pytest.mark.parametrize("path,value", [
    (("schema_version",), True),
    (("targets", "snv"), None),
    (("sources", "rna", "identity"), "not-an-object"),
    (("members", "alt", "source"), ["rna"]),
    (("members", "alt", "policy"), None),
    (("members", "alt", "policy", "records"), None),
    (("members", "alt", "policy", "strata"), []),
    (("members", "alt", "policy", "assignments"), [None]),
    (("members", "alt", "regions"), None),
])
def test_malformed_recipes_fail_before_reading_sources(bam, path, value):
    from osteosarc import SchemaError
    recipe = recipe_for(bam)
    parent = recipe
    for key in path[:-1]:
        parent = parent[key]
    parent[path[-1]] = value
    with pytest.raises(SchemaError):
        select_fixtures(recipe, {})


def test_exact_selection_requires_an_explicit_record_mapping(bam):
    from osteosarc import SchemaError
    recipe = recipe_for(bam)
    del recipe["members"]["alt"]["policy"]["records"]
    with pytest.raises(SchemaError, match="explicit records"):
        select_fixtures(recipe, {"rna": bam})


def test_equivalent_acquired_archive_and_direct_input_agree(bam, tmp_path):
    from osteosarc import Region, digest, extract_reads
    recipe = recipe_for(bam)
    recipe["sources"]["rna"]["archive_sha256"] = digest(bam)
    direct = select_fixtures(recipe, {"rna": bam})
    subset = extract_reads(bam, [Region("chr1", 0, 2000, "GRCh38")], cache=tmp_path / "cache")
    acquired = select_fixtures(recipe, {"rna": subset})
    assert acquired.members == direct.members
    subset.receipt["request"]["source"] = "https://example.test/different-source.bam"
    with pytest.raises(IntegrityError, match="source identity"):
        select_fixtures(recipe, {"rna": subset})
