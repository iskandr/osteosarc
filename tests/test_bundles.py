import copy
import json
import shutil
from pathlib import Path

import pytest

from osteosarc import (
    IntegrityError,
    digest,
    export_bundle,
    generate_bundle,
    list_bundle,
    pack_bundle,
    record_multiset,
    select_fixtures,
    verify_bundle,
)
from osteosarc.cli import main


def bundle_recipe(bam):
    return dict(schema_version=1, id="portable-test", redistribution=dict(license="CC0-1.0"), targets={
        "snv": dict(kind="small_variant", assembly="GRCh38", reference=dict(id="GRCh38"),
                    coordinates="one-based", contig="chr1", position=101, ref="A", alt="C")},
        sources={"rna": dict(identity=dict(url="https://example.test/source.bam"), assembly="GRCh38",
                             sample="T0", library="rna", product="original", archive_sha256=digest(bam))},
        members={"duplicates": dict(source="rna", target="snv", policy=dict(kind="exact", version=1,
                                                                           records=dict(record_multiset(bam))))})


def test_generate_pack_export_verify_fresh_offline_directory(bam, tmp_path, monkeypatch):
    import osteosarc.reads
    recipe = bundle_recipe(bam)
    recipe["members"]["shared"] = copy.deepcopy(recipe["members"]["duplicates"])
    recipe["members"]["empty"] = dict(source="rna", target="snv", policy=dict(kind="empty", version=1))
    first = tmp_path / "bundle"
    manifest = generate_bundle(recipe, first, sources={"rna": bam})
    assert manifest["sources"]["rna"]["record_count"] == 7  # shared, not 14
    packed = tmp_path / "packed"
    assert pack_bundle(select_fixtures(recipe, {"rna": bam}), packed) == manifest
    assert digest(first / "manifest.json") == digest(packed / "manifest.json")
    bam.unlink()
    Path(str(bam) + ".bai").unlink()
    monkeypatch.setattr(osteosarc.reads, "_run", lambda *a: pytest.fail("offline export ran a command"))
    fresh = tmp_path / "fresh"
    shutil.copytree(first, fresh)
    exported = export_bundle(fresh, tmp_path / "exported")
    verify_bundle(tmp_path / "exported", sha256=digest(tmp_path / "exported/manifest.json"))
    assert list_bundle(fresh) == manifest["members"]
    assert len(exported["exports"]) == 3
    assert record_multiset(tmp_path / "exported/members/empty.bam") == {}
    assert record_multiset(tmp_path / "exported/members/duplicates.bam") == recipe["members"]["duplicates"]["policy"]["records"]
    sam = export_bundle(fresh, tmp_path / "sam", members=["duplicates"], format="sam.gz")
    assert sam["exports"]["duplicates"]["fidelity"] == "sam-text-v1"


def test_cli_and_dataset_produce_same_bundle(bam, dataset, tmp_path, capsys):
    recipe = bundle_recipe(bam)
    path = tmp_path / "recipe.json"
    path.write_text(json.dumps(recipe))
    expected = dataset.generate_bundle(recipe, tmp_path / "api", sources={"rna": bam})
    assert main(["--cache", str(tmp_path / "cache"), "--offline", "fixtures", "generate", str(path),
                 str(tmp_path / "cli"), "--source", f"rna={bam}"]) == 0
    assert json.loads(capsys.readouterr().out) == expected
    assert main(["--cache", str(tmp_path / "cache"), "--offline", "fixtures", "verify", str(tmp_path / "cli")]) == 0


@pytest.mark.parametrize("corruption", ["record", "index", "recipe", "nested", "duplicate", "traversal"])
def test_corruption_is_actionable(bam, tmp_path, corruption):
    recipe = bundle_recipe(bam)
    directory = tmp_path / "bundle"
    manifest = generate_bundle(recipe, directory, sources={"rna": bam})
    source = manifest["sources"]["rna"]
    if corruption in ("record", "index"):
        path = directory / source["bam" if corruption == "record" else "index"]
        path.write_bytes(path.read_bytes()[:-12])
    elif corruption == "nested":
        (directory / source["original_header"]).unlink()
    elif corruption == "recipe":
        manifest["recipe_sha256"] = "0" * 64
    elif corruption == "duplicate":
        key = next(k for k, n in source["records"].items() if n == 2)
        source["records"][key] = 1
    else:
        manifest["files"]["../escape"] = dict(sha256="0" * 64, size_bytes=0)
    (directory / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(IntegrityError, match="fixture|Fixture|recipe|multiset|Unsafe"):
        verify_bundle(directory)


def test_atomic_failure_size_budget_and_no_overwrite(bam, tmp_path):
    recipe = bundle_recipe(bam)
    target = tmp_path / "bundle"
    with pytest.raises(IntegrityError, match="size budget"):
        generate_bundle(recipe, target, sources={"rna": bam}, size_budget=100)
    assert not target.exists()
    generate_bundle(recipe, target, sources={"rna": bam})
    with pytest.raises(FileExistsError):
        generate_bundle(recipe, target, sources={"rna": bam})
    before = digest(target / "manifest.json")
    recipe["members"]["../bad"] = recipe["members"].pop("duplicates")
    other = tmp_path / "unsafe"
    generate_bundle(recipe, other, sources={"rna": bam})
    with pytest.raises(IntegrityError, match="Unsafe"):
        export_bundle(other, tmp_path / "export")
    assert not (tmp_path / "export").exists()
    assert digest(target / "manifest.json") == before


def test_compact_keeps_verified_assembly_and_full_header_provenance(bam, tmp_path):
    from osteosarc import Region, resolve_regions
    manifest = generate_bundle(bundle_recipe(bam), tmp_path / "bundle", sources={"rna": bam}, header_policy="compact")
    source = manifest["sources"]["rna"]
    assert resolve_regions([Region("chr1", 100, 150, "GRCh38")], source["exported_header"])
    assert (tmp_path / "bundle" / source["original_header"]).is_file()
    verify_bundle(tmp_path / "bundle")


def test_swapped_index_detected_even_if_file_hash_is_relisted(bam, tmp_path):
    import pysam
    directory = tmp_path / "bundle"
    manifest = generate_bundle(bundle_recipe(bam), directory, sources={"rna": bam})
    source = manifest["sources"]["rna"]
    empty = tmp_path / "empty.bam"
    with pysam.AlignmentFile(bam) as inp, pysam.AlignmentFile(empty, "wb", header=inp.header):
        pass
    pysam.index(str(empty))
    path = directory / source["index"]
    path.write_bytes(Path(str(empty) + ".bai").read_bytes())
    manifest["files"][source["index"]] = dict(sha256=digest(path), size_bytes=path.stat().st_size)
    manifest["total_size_bytes"] = sum(v["size_bytes"] for v in manifest["files"].values())
    (directory / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(IntegrityError, match="Index does not enumerate"):
        verify_bundle(directory)


def test_compact_retains_pg_ancestry_and_unresolved_metadata():
    import pysam

    from osteosarc import compact_header
    from osteosarc.records import FixtureRecord
    header = dict(SQ=[dict(SN="chr1", LN=248956422), dict(SN="chr2", LN=242193529)],
                  RG=[dict(ID="rg", SM="sample", LB="library", PG="aligned")],
                  PG=[dict(ID="base"), dict(ID="aligned", PP="base"), dict(ID="unrelated")],
                  CO=["unresolved source annotation"])
    h = pysam.AlignmentHeader.from_dict(header)
    r = pysam.AlignedSegment(h)
    r.query_name = "q"
    r.set_tag("RG", "rg")
    r.set_tag("PG", "aligned")
    compact = compact_header(header, [FixtureRecord(r, "digest")])
    assert [p["ID"] for p in compact["PG"]] == ["base", "aligned"]
    assert compact["SQ"] == header["SQ"]
    assert compact["CO"] == header["CO"]
    r.set_tag("PG", None)
    assert compact_header(header, [FixtureRecord(r, "digest")])["PG"] == header["PG"]
