import copy
import json
import os
import shutil
from pathlib import Path

import pytest

from osteosarc import IntegrityError, export_bundle, generate_bundle, list_bundle, verify_bundle
from osteosarc.bundles import pack_bundle
from osteosarc.cache import digest
from osteosarc.cli import main
from osteosarc.fixtures import select_fixtures
from osteosarc.records import read_records, record_multiset


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
    assert list_bundle(fresh) == manifest["members"]
    assert exported == {name: tmp_path / "exported" / f"{name}.bam" for name in ("duplicates", "empty", "shared")}
    assert record_multiset(tmp_path / "exported/empty.bam") == {}
    assert record_multiset(tmp_path / "exported/duplicates.bam") == recipe["members"]["duplicates"]["policy"]["records"]
    assert (tmp_path / "exported/duplicates.bam.bai").is_file()
    sam = export_bundle(fresh, tmp_path / "sam", members=["duplicates"], format="sam.gz")
    assert sam == {"duplicates": tmp_path / "sam/duplicates.sam.gz"}


def test_cli_and_dataset_produce_same_bundle(bam, dataset, tmp_path, capsys):
    recipe = bundle_recipe(bam)
    path = tmp_path / "recipe.json"
    path.write_text(json.dumps(recipe))
    expected = generate_bundle(recipe, tmp_path / "api", sources={"rna": bam}, dataset=dataset)
    assert main(["--cache", str(tmp_path / "cache"), "--offline", "test-data", "generate", str(path),
                 str(tmp_path / "cli"), "--source", f"rna={bam}", "--json"]) == 0
    assert json.loads(capsys.readouterr().out) == expected
    assert main(["--cache", str(tmp_path / "cache"), "--offline", "test-data", "verify", str(tmp_path / "cli")]) == 0
    assert "1 member (1 with reads), 7 records from 1 BAM" in capsys.readouterr().out


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
    from osteosarc import Region
    from osteosarc.reads import resolve_regions
    manifest = generate_bundle(bundle_recipe(bam), tmp_path / "bundle", sources={"rna": bam}, header_policy="compact")
    source = manifest["sources"]["rna"]
    import pysam
    with pysam.AlignmentFile(tmp_path / "bundle" / source["bam"]) as inp:
        assert resolve_regions([Region("chr1", 100, 150, "GRCh38")], inp.header.to_dict())
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

    from osteosarc.bundles import compact_header
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
    with pysam.AlignmentFile(dest / ("duplicates." + format), "r") as sam:
        assert sam.header.to_dict()["HD"]["SO"] == "coordinate"
        with pysam.AlignmentFile(converted, "wb", template=sam) as out:
            positions = []
            for read in sam:
                positions.append((read.reference_id, read.reference_start))
                out.write(read)
    assert positions == sorted(positions)
    pysam.index(str(converted))


def test_export_into_an_existing_folder_keeps_identical_files_and_refuses_different_ones(bam, tmp_path):
    source = tmp_path / "bundle"
    generate_bundle(bundle_recipe(bam), source, sources={"rna": bam})
    to = tmp_path / "tests/data"
    to.mkdir(parents=True)
    (to / "other.txt").write_text("kept")
    first = export_bundle(source, to, members=["duplicates", "duplicates"])
    assert export_bundle(source, to) == first
    assert sorted(p.name for p in to.iterdir()) == ["duplicates.bam", "duplicates.bam.bai", "other.txt"]
    # The same records in differently compressed bytes are the same export; a lost index comes back.
    import pysam
    with pysam.AlignmentFile(str(to / "duplicates.bam")) as inp:
        header, reads = inp.header, list(inp)
    (to / "duplicates.bam.bai").unlink()
    with pysam.AlignmentFile(str(to / "duplicates.bam"), "wb0", header=header) as out:
        for read in reads:
            out.write(read)
    assert export_bundle(source, to) == first and (to / "duplicates.bam.bai").is_file()
    # A conflict anywhere writes nothing: not the SAM, and not a BAM beside a stale index.
    (to / "duplicates.sam").write_text("different")
    with pytest.raises(FileExistsError, match="different contents"):
        export_bundle(source, to, format="sam")
    assert (to / "duplicates.sam").read_text() == "different"
    other = tmp_path / "other"
    other.mkdir()
    (other / "duplicates.bam.bai").write_text("stale")
    with pytest.raises(FileExistsError, match="duplicates.bam.bai"):
        export_bundle(source, other)
    assert sorted(p.name for p in other.iterdir()) == ["duplicates.bam.bai"]
    assert sorted(p.name for p in to.iterdir()) == ["duplicates.bam", "duplicates.bam.bai", "duplicates.sam",
                                                    "other.txt"]


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
    assert export_bundle(directory, tmp_path / "export") == {}


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
    asset = next(a for a in dataset.files if a.format == "bam")
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
    manifest = generate_bundle(recipe, tmp_path / "bundle", dataset=dataset)
    assert manifest["members"]["duplicates"]["records"] == dict(record_multiset(bam))
    acquired = json.loads((tmp_path / "bundle/acquisition.json").read_text())["rna"]
    assert acquired["request"]["source"] == asset.url
    assert acquired["request"]["snapshot_id"] == dataset.id


def test_frozen_additional_sv_panel_generates_offline(tmp_path):
    import osteosarc
    recipe = json.loads((Path(osteosarc.__file__).parent / "data/additional_sv_recipe.json").read_text())
    root = Path(__file__).parent / "data/additional_svs"
    sources = {sid: root / (sid + ".bam") for sid in recipe["sources"]}
    manifest = generate_bundle(recipe, tmp_path / "additional", sources=sources, size_budget=4_000_000)
    assert len(manifest["members"]) == 24
    assert sum(s["record_count"] for s in manifest["sources"].values()) == 251
    assert {m["target"] for m in manifest["members"].values()} == {"SV0055", "SV0175", "SV0402", "SV0461", "SV0499"}
    assert manifest["total_size_bytes"] < 4_000_000
    exported = export_bundle(tmp_path / "additional", tmp_path / "exported", members=["SV0461/T1-PacBio"])
    assert exported == {"SV0461/T1-PacBio": tmp_path / "exported/SV0461/T1-PacBio.bam"}
    assert record_multiset(exported["SV0461/T1-PacBio"]) == manifest["members"]["SV0461/T1-PacBio"]["records"]


def test_legacy_full_header_bundles_remain_readable_and_exportable(bam, tmp_path):
    import pysam
    directory = tmp_path / "legacy"
    manifest = generate_bundle(bundle_recipe(bam), directory, sources={"rna": bam})
    source = manifest["sources"]["rna"]
    with pysam.AlignmentFile(directory / source["bam"]) as inp:
        source["exported_header"] = inp.header.to_dict()
    del source["header_sha256"]
    (directory / "manifest.json").write_text(json.dumps(manifest))
    assert verify_bundle(directory)["members"] == manifest["members"]
    exported = export_bundle(directory, tmp_path / "export", format="sam")
    assert set(exported) == set(manifest["members"])
    source["exported_header"]["HD"]["SO"] = "unsorted"
    (directory / "manifest.json").write_text(json.dumps(manifest))
    with pytest.raises(IntegrityError, match="Exported header differs"):
        verify_bundle(directory)


def test_a_published_bundle_is_fetched_verified_and_checked_offline(bam, tmp_path, monkeypatch):
    import gzip
    import json

    import pysam

    import osteosarc.shared as shared
    from osteosarc import Cache
    bundle = tmp_path / "bundle"
    generate_bundle(bundle_recipe(bam), bundle, sources={"rna": bam})
    release = shared.pack_release(bundle, tmp_path / "panel.tar.gz")
    assert shared.pack_release(bundle, tmp_path / "again.tar.gz") == release  # reproducible
    published_dir = tmp_path / "bundles"
    published_dir.mkdir()
    url = "https://example.test/panel.tar.gz"
    (published_dir / "tiny-v1.release.json").write_text(json.dumps(dict(release, url=url)))
    monkeypatch.setattr(shared, "BUNDLES", published_dir)
    cache = Cache(tmp_path / "cache", offline=True)
    cache.import_file(tmp_path / "panel.tar.gz", url)
    fetched = shared.fetch_bundle("tiny-v1", cache=cache)
    assert fetched == shared.fetch_bundle("tiny-v1", cache=cache)  # reused
    with pytest.raises(KeyError, match="published: tiny-v1"):
        shared.published("missing")
    # A library's copy passes when its records are the member's; its file format doesn't matter.
    with pysam.AlignmentFile(str(bam)) as handle:
        lines = [read.to_string() for read in handle]
        text = str(handle.header) + "".join(line + "\n" for line in lines)
    (tmp_path / "copy.sam.gz").write_bytes(gzip.compress(text.encode()))
    (tmp_path / "evidence.json").write_text(json.dumps(dict(paths=[dict(records=lines[:-1])])))
    fixtures = {"duplicates": str(bam), "empty": "copy.sam.gz"}
    problems = shared.check_fixtures(fetched, {"duplicates": "copy.sam.gz"}, root=tmp_path)
    assert problems == {}
    problems = shared.check_fixtures(fetched, dict(fixtures, duplicates=dict(json="evidence.json",
                                     pointer="/paths/0/records")), root=tmp_path)
    assert problems["duplicates"] == dict(missing=1, extra=0) and "not a member" in problems["empty"]["error"]
    # The cached copy is read-only, and a change to it is noticed.
    manifest = fetched / "manifest.json"
    assert not os.access(manifest, os.W_OK)
    manifest.chmod(0o644)
    manifest.write_text(manifest.read_text() + " ")
    with pytest.raises(IntegrityError, match="changed since it was downloaded"):
        shared.fetch_bundle("tiny-v1", cache=cache)


def test_check_compares_members_without_reads_too(bam, tmp_path):
    import osteosarc.shared as shared
    recipe = bundle_recipe(bam)
    recipe["targets"]["pending"] = dict(kind="unresolved", reason="unreviewed breakends")
    recipe["sources"]["other"] = copy.deepcopy(recipe["sources"]["rna"])
    recipe["members"]["pending"] = dict(target="pending", source="other", policy=dict(kind="empty", version=1))
    bundle = tmp_path / "bundle"
    manifest = generate_bundle(recipe, bundle, sources={"rna": bam})
    assert "other" not in manifest["sources"]
    (tmp_path / "empty.sam").write_text("@HD\tVN:1.6\n")
    assert shared.check_fixtures(bundle, {"pending": "empty.sam"}, root=tmp_path) == {}
    assert shared.check_fixtures(bundle, {"pending": str(bam)}, root=tmp_path) == {"pending": dict(missing=0, extra=7)}
    # One unreadable fixture is reported without stopping the others.
    (tmp_path / "lines.json").write_text(json.dumps({"records": []}))
    problems = shared.check_fixtures(bundle, {"pending": dict(json="lines.json", pointer="/record"),
                                              "duplicates": str(bam)}, root=tmp_path)
    assert list(problems) == ["pending"] and "nothing at /record" in problems["pending"]["error"]


def test_the_cli_lists_and_checks_a_published_bundle(bam, tmp_path, monkeypatch, capsys):
    import json

    import osteosarc.shared as shared
    from osteosarc import Cache
    from osteosarc.cli import main
    bundle = tmp_path / "bundle"
    generate_bundle(bundle_recipe(bam), bundle, sources={"rna": bam})
    release = shared.pack_release(bundle, tmp_path / "panel.tar.gz")
    published_dir = tmp_path / "bundles"
    published_dir.mkdir()
    url = "https://example.test/panel.tar.gz"
    (published_dir / "tiny-v1.release.json").write_text(json.dumps(dict(release, url=url)))
    monkeypatch.setattr(shared, "BUNDLES", published_dir)
    root = tmp_path / "cache"
    Cache(root, offline=True).import_file(tmp_path / "panel.tar.gz", url)
    (tmp_path / "fixtures.json").write_text(json.dumps({"duplicates": str(bam)}))
    assert main(["--cache", str(root), "--offline", "test-data", "check", "tiny-v1", str(tmp_path / "fixtures.json")]) == 0
    assert "matches the bundle" in capsys.readouterr().out
    (tmp_path / "wrong.json").write_text(json.dumps({"shared": str(bam)}))
    assert main(["--cache", str(root), "--offline", "test-data", "check", "tiny-v1", str(tmp_path / "wrong.json")]) == 1
    assert "not a member" in capsys.readouterr().out
    assert main(["--cache", str(root), "--offline", "test-data", "list", "tiny-v1"]) == 0
    assert capsys.readouterr().out.split("\n")[2].split() == ["selected", "7", "duplicates"]
    assert main(["--cache", str(root), "--offline", "test-data", "list", "tiny-v1", "--json"]) == 0
    assert "duplicates" in json.loads(capsys.readouterr().out)
    out = tmp_path / "tests/data"
    assert main(["--cache", str(root), "--offline", "test-data", "export", "tiny-v1", str(out),
                 "--member", "duplicates"]) == 0
    assert capsys.readouterr().out.strip() == str(out / "duplicates.bam")
    assert record_multiset(out / "duplicates.bam") == record_multiset(bam)
    assert main(["--cache", str(root), "--offline", "test-data", "list", "no-such-bundle"]) == 1
    assert "published: tiny-v1" in capsys.readouterr().err


def test_a_bundles_library_fixtures_carry_forward_by_checksum(bam, tmp_path):
    from osteosarc import SchemaError
    from osteosarc.shared import bundle_fixtures, match_records
    recipe = bundle_recipe(bam)
    fixture = dict(kind="fixture", assembly="GRCh38", reference=dict(source="library", consumer="isovar"),
                   consumer="isovar", description="Isovar fixture x.sam")
    recipe["targets"]["fixture:isovar/x.sam"] = fixture
    recipe["targets"]["fixture:isovar/planned.sam"] = dict(fixture, regions=[["chr1", 0, 5000]])
    recipe["targets"]["fixture:isovar/none.sam"] = dict(fixture)
    exact = recipe["members"]["duplicates"]["policy"]
    for name in ("x.sam", "planned.sam"):
        recipe["members"]["isovar/" + name] = dict(target="fixture:isovar/" + name, source="rna",
                                                   policy=copy.deepcopy(exact))
    recipe["members"]["isovar/none.sam"] = dict(target="fixture:isovar/none.sam", source="rna",
                                                policy=dict(kind="exact", version=1, records={}))
    generate_bundle(recipe, tmp_path / "bundle", sources={"rna": bam})
    carried = bundle_fixtures(tmp_path / "bundle")
    assert list(carried) == ["isovar/none.sam", "isovar/planned.sam", "isovar/x.sam"]  # not selected members
    subset = carried["isovar/x.sam"]
    assert subset["records"] == exact["records"]  # the exact records, repeats kept, by checksum
    assert (subset["consumer"], subset["source"], subset["description"]) == (
        "isovar", "https://example.test/source.bam", "Isovar fixture x.sam")
    # Planned from the target's recorded regions, or else from the records' own spans.
    assert carried["isovar/planned.sam"]["regions"] == [["chr1", 0, 5000]]
    assert all(contig == "chr1" for contig, _, _ in subset["regions"]) and subset["regions"]
    assert carried["isovar/none.sam"] == dict(subset, records={}, regions=[])
    records = list(read_records(bam))
    assert match_records(records, subset["records"]) == subset["records"]
    with pytest.raises(IntegrityError, match="pinned records"):
        match_records(records[:1], subset["records"])
    recipe["targets"]["fixture:isovar/x.sam"].pop("consumer")
    generate_bundle(recipe, tmp_path / "anonymous", sources={"rna": bam})
    with pytest.raises(SchemaError, match="Can't carry isovar/x.sam"):
        bundle_fixtures(tmp_path / "anonymous")


def test_fresh_lists_replace_a_librarys_carried_fixtures(tmp_path):
    import gzip

    from osteosarc import SchemaError
    from osteosarc.shared import merge_required
    carried = {"isovar/a": dict(consumer="isovar", source="s", records={}, regions=[]),
               "varcode/b": dict(consumer="varcode", source="s", records={}, regions=[])}

    def required(name, consumer, subsets):
        path = tmp_path / name
        path.write_bytes(gzip.compress(json.dumps(dict(consumer=consumer, subsets=subsets)).encode()))
        return path
    line = "r\t0\tchr1\t1\t60\t4M\t*\t0\t0\tACGT\tIIII"
    fresh = required("isovar.json.gz", "isovar", {"isovar/c": dict(source="s", sam=[line])})
    assert set(merge_required(carried, [fresh])) == {"isovar/c", "varcode/b"}
    nothing = required("varcode.json.gz", "varcode", {})  # a library that needs no reads any more
    assert set(merge_required(carried, [nothing])) == {"isovar/a"}
    clash = required("topiary.json.gz", "topiary", {"varcode/b": dict(source="s", sam=[line])})
    with pytest.raises(SchemaError, match="share names"):
        merge_required(carried, [clash])
