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


@pytest.mark.parametrize("format", ["sam", "sam.gz"])
def test_sam_exports_are_coordinate_sorted_and_indexable(bam, tmp_path, format):
    import pysam
    source = tmp_path / "bundle"
    generate_bundle(bundle_recipe(bam), source, sources={"rna": bam})
    dest = tmp_path / "export"
    export_bundle(source, dest, format=format)
    converted = tmp_path / "converted.bam"
    with pysam.AlignmentFile(dest / ("members/duplicates." + format), "r") as sam:
        assert sam.header.to_dict()["HD"]["SO"] == "coordinate"
        with pysam.AlignmentFile(converted, "wb", template=sam) as out:
            positions = []
            for read in sam:
                positions.append((read.reference_id, read.reference_start))
                out.write(read)
    assert positions == sorted(positions)
    pysam.index(str(converted))


def test_reexport_replaces_export_set_without_stale_files(bam, tmp_path):
    source = tmp_path / "bundle"
    generate_bundle(bundle_recipe(bam), source, sources={"rna": bam})
    first, second, third = [tmp_path / name for name in ("first", "second", "third")]
    export_bundle(source, first)
    export_bundle(first, second, members=["duplicates", "duplicates"])
    assert record_multiset(first / "members/duplicates.bam") == record_multiset(second / "members/duplicates.bam")
    changed = export_bundle(second, third, format="sam.gz")
    assert sorted(p.name for p in (third / "members").iterdir()) == ["duplicates.sam.gz"]
    assert changed["parent"]["manifest_sha256"] == digest(second / "manifest.json")
    verify_bundle(third)
    with pytest.raises(FileExistsError):
        export_bundle(first, second)


@pytest.mark.parametrize("corruption", ["record", "format"])
def test_sam_export_semantics_checked_after_relisting_bytes(bam, tmp_path, corruption):
    source, dest = tmp_path / "bundle", tmp_path / "export"
    generate_bundle(bundle_recipe(bam), source, sources={"rna": bam})
    manifest = export_bundle(source, dest, format="sam")
    if corruption == "format":
        manifest["exports"]["duplicates"]["format"] = "unknown"
    else:
        name = "members/duplicates.sam"
        path = dest / name
        path.write_text(path.read_text().replace("repeated", "modified", 1))
        manifest["files"][name] = dict(sha256=digest(path), size_bytes=path.stat().st_size)
        manifest["total_size_bytes"] = sum(v["size_bytes"] for v in manifest["files"].values())
    (dest / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(IntegrityError, match="SAM export record|Unsupported export format"):
        verify_bundle(dest)


def test_exact_pins_remain_required_when_context_is_declared(bam, tmp_path):
    recipe = bundle_recipe(bam)
    recipe["members"]["duplicates"]["context_regions"] = [dict(contig="chr1", start=100, end=1100, assembly="GRCh38")]
    selection = select_fixtures(recipe, {"rna": bam})
    member = selection.members["duplicates"]
    key = next(iter(member["records"]))
    member["record_count"] -= member["records"].pop(key)
    member["reasons"].pop(key)
    with pytest.raises(IntegrityError, match="pinned exact recipe"):
        pack_bundle(selection, tmp_path / "invalid")
    assert not (tmp_path / "invalid").exists()


@pytest.mark.parametrize("field,value", [("record_count", 999), ("record_count", True),
                                        ("status", "empty"), ("status", "omitted"), ("status", "unknown")])
def test_member_counts_and_status_are_verified(bam, tmp_path, field, value):
    directory = tmp_path / "bundle"
    manifest = generate_bundle(bundle_recipe(bam), directory, sources={"rna": bam})
    manifest["members"]["duplicates"][field] = value
    (directory / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(IntegrityError, match="count|status"):
        list_bundle(directory)


def test_unavailable_exact_member_stays_explicit_without_fabricated_records(bam, tmp_path):
    recipe = bundle_recipe(bam)
    recipe["targets"]["snv"] = dict(kind="unresolved", reason="No resolved coordinates")
    directory = tmp_path / "bundle"
    manifest = generate_bundle(recipe, directory)
    assert manifest["members"]["duplicates"]["record_count"] == 0
    assert not manifest["sources"]
    assert not export_bundle(directory, tmp_path / "export")["exports"]


@pytest.mark.parametrize("url", ["https://example.test/source.bam", "https://example.test/alignment"])
def test_direct_acquisition_preserves_inventory_size_pin(bam, tmp_path, monkeypatch, url):
    import osteosarc.reads as reads
    recipe = bundle_recipe(bam)
    source = recipe["sources"]["rna"]
    source["identity"]["url"] = url
    source.pop("archive_sha256")
    source["identity"]["size"] = bam.stat().st_size + 100
    source["index"] = str(bam) + ".bai"
    source["regions"] = [dict(contig="chr1", start=100, end=1100, assembly="GRCh38")]
    monkeypatch.setattr(reads, "_remote_identity", lambda *a: {"content-length": str(bam.stat().st_size)})
    original = reads._run
    url = source["identity"]["url"]

    def local_transport(command, timeout):
        return original([str(bam) if arg == url else arg for arg in command], timeout)

    monkeypatch.setattr(reads, "_run", local_transport)
    with pytest.raises(IntegrityError, match="pinned inventory"):
        generate_bundle(recipe, tmp_path / "rejected", cache=tmp_path / "cache")
    assert not (tmp_path / "rejected").exists()
    source["identity"]["size"] = bam.stat().st_size
    generate_bundle(recipe, tmp_path / "accepted", cache=tmp_path / "cache")
    acquired = json.loads((tmp_path / "accepted/acquisition.json").read_text())["rna"]
    assert acquired["request"]["source_size"] == bam.stat().st_size


@pytest.mark.parametrize("lookup", ["key", "id"])
def test_dataset_acquisition_accepts_identity_without_url(bam, dataset, tmp_path, monkeypatch, lookup):
    import osteosarc.reads as reads
    dataset.cache.offline = False
    asset = next(a for a in dataset.assets if a.format == "bam")
    recipe = bundle_recipe(bam)
    source = recipe["sources"]["rna"]
    source.pop("archive_sha256")
    source["identity"] = {lookup: getattr(asset, lookup)}
    source["index"] = str(bam) + ".bai"
    source["regions"] = [dict(contig="chr1", start=100, end=1100, assembly="GRCh38")]
    monkeypatch.setattr(reads, "_remote_identity", lambda *a: {"content-length": str(asset.size) if asset.size is not None else None})
    original = reads._run
    monkeypatch.setattr(reads, "_run", lambda command, timeout: original(
        [str(bam) if arg == asset.url else arg for arg in command], timeout))
    manifest = dataset.generate_bundle(recipe, tmp_path / "bundle")
    assert manifest["members"]["duplicates"]["records"] == dict(record_multiset(bam))
    acquired = json.loads((tmp_path / "bundle/acquisition.json").read_text())["rna"]
    assert acquired["request"]["source"] == asset.url
    assert acquired["request"]["snapshot_id"] == dataset.id
