"""Portable, atomic, size-bounded fixture bundles and offline indexed export."""

from __future__ import annotations

import copy
import gzip
import json
import os
import shutil
import tempfile
from collections import Counter
from contextlib import contextmanager
from pathlib import Path

from .cache import Cache, digest, file_lock, stable_id, write_json
from .errors import IntegrityError
from .fixtures import select_fixtures, validate_recipe
from .models import Region
from .reads import extract_reads
from .records import RECORD_ENCODING, read_records, record_multiset

DEFAULT_SIZE_BUDGET = 64 * 1024 * 1024


def safe_path(root, name):
    """Resolve a portable relative filename, refusing traversal and symlinks."""
    if (not isinstance(name, str) or not name or any(c in name for c in "\\:\x00")
            or any(p in ("", ".", "..") for p in name.split("/"))):
        raise IntegrityError(f"Unsafe relative fixture path: {name!r}")
    root = Path(root)
    path = root
    for part in name.split("/"):
        path = path / part
        if path.is_symlink():
            raise IntegrityError(f"Symlink in fixture path: {name}")
    return path


def verify_digest(path, expected, *, label=None):
    """Read a small pinned artifact and fail on changed bytes (also under -O)."""
    import hashlib
    data = Path(path).read_bytes()
    if hashlib.sha256(data).hexdigest() != expected:
        raise IntegrityError(f"Fixture checksum mismatch: {label or path}")
    return data


def verify_gzip_digests(path, compressed, uncompressed, *, label=None):
    """Check both a stored gzip artifact and its decompressed historical content."""
    import hashlib
    data = gzip.decompress(verify_digest(path, compressed, label=label))
    if hashlib.sha256(data).hexdigest() != uncompressed:
        raise IntegrityError(f"Uncompressed fixture checksum mismatch: {label or path}")
    return data


def verify_manifest_files(root, files, *, label=None):
    """Verify safe nested {relative filename: sha256} or {filename: metadata}."""
    for name, expected in files.items():
        path = safe_path(root, name)
        checksum = expected if isinstance(expected, str) else expected["sha256"]
        if not path.is_file():
            raise IntegrityError(f"Missing fixture file: {label or ''}/{name}")
        if isinstance(expected, dict) and path.stat().st_size != expected["size_bytes"]:
            raise IntegrityError(f"Fixture size mismatch: {name}")
        if digest(path) != checksum:
            raise IntegrityError(f"Fixture checksum mismatch: {label or ''}/{name}")


@contextmanager
def _publication(destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with file_lock(destination.parent / ("." + destination.name + ".lock")):
        if destination.exists() or destination.is_symlink():
            raise FileExistsError(destination)
        with tempfile.TemporaryDirectory(dir=destination.parent, prefix=".bundle-") as temporary:
            work = Path(temporary)
            yield work
            os.rename(work, destination)


def compact_header(header, records):
    """Keep assembly-defining SQs, comments, RG identity and proven PG ancestry.

    If any selected record has no explicit PG, producer lineage is unresolved:
    preserve all PG entries rather than guessing which aligner produced it.
    Missing source RG/PG declarations stay missing. Original header is archived.
    """
    result = copy.deepcopy(header)
    groups = {r.read.get_tag("RG") for r in records if r.read.has_tag("RG")}
    result["RG"] = [r for r in result.get("RG", []) if r["ID"] in groups]
    if records and all(r.read.has_tag("PG") for r in records):
        programs = {r.read.get_tag("PG") for r in records}
        programs.update(row["PG"] for row in result["RG"] if "PG" in row)
        by_id = {row["ID"]: row for row in result.get("PG", [])}
        pending = list(programs)
        while pending:
            previous = by_id.get(pending.pop(), {}).get("PP")
            if previous and previous not in programs:
                programs.add(previous)
                pending.append(previous)
        result["PG"] = [row for row in result.get("PG", []) if row["ID"] in programs]
    return result


def _write_bam(path, header, records, counts):
    import pysam
    header = copy.deepcopy(header)
    header.setdefault("HD", {})["SO"] = "coordinate"
    header["HD"].pop("SS", None)
    header["HD"].pop("GO", None)
    by_id = {r.digest: r for r in records}
    if missing := Counter(counts) - Counter(r.digest for r in records):
        raise IntegrityError(f"Missing export record occurrences: {dict(missing)}")
    ordered = sorted(counts, key=lambda k: (by_id[k].read.reference_id if by_id[k].read.reference_id >= 0 else len(header["SQ"]),
                                           by_id[k].read.reference_start, k))
    path.parent.mkdir(parents=True, exist_ok=True)
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for key in ordered:
            for _ in range(counts[key]):
                out.write(by_id[key].read)
    pysam.index(str(path))
    if record_multiset(path) != Counter(counts):
        raise IntegrityError("BAM export changed original record bytes/multiplicity")
    with pysam.AlignmentFile(path) as bam:
        return bam.header.to_dict()


def _inventory(work):
    return {p.relative_to(work).as_posix(): dict(sha256=digest(p), size_bytes=p.stat().st_size)
            for p in sorted(work.rglob("*")) if p.is_file() and p != work / "manifest.json"}


def _finish(work, manifest, budget):
    manifest["files"] = _inventory(work)
    manifest["size_budget_bytes"] = budget
    manifest["total_size_bytes"] = sum(f["size_bytes"] for f in manifest["files"].values())
    write_json(work / "manifest.json", manifest)
    if manifest["total_size_bytes"] + (work / "manifest.json").stat().st_size > budget:
        raise IntegrityError(f"Fixture bundle exceeds {budget} byte size budget")
    verify_bundle(work)
    return manifest


def pack_bundle(selection, destination, *, header_policy="full", size_budget=DEFAULT_SIZE_BUDGET, parent=None):
    """Store shared records once per source; atomically publish a verified bundle.

    The manifest embeds the recipe digest, selected multiset, original header,
    acquisition receipts, scope, reasons, source license and toolchain. No source
    BAM or sibling repository is needed to verify, list or export it afterward.
    """
    import pysam
    if header_policy not in ("full", "compact"):
        raise ValueError("header_policy must be full or compact")
    from . import __version__
    with _publication(destination) as work:
        write_json(work / "recipe.json", selection.recipe)
        write_json(work / "acquisition.json", selection.receipts)
        manifest = dict(schema_version=1, kind="osteosarc-fixture-bundle", recipe_sha256=stable_id(selection.recipe),
                        record_encoding=RECORD_ENCODING, members=selection.members, sources={}, exports={},
                        parent=parent, header_policy=header_policy, suitable_for_abundance=False,
                        toolchain=dict(osteosarc=__version__, pysam=pysam.__version__, samtools=pysam.__samtools_version__),
                        redistribution=selection.recipe.get("redistribution", {"license": "unresolved"}))
        for sid, header in sorted(selection.headers.items()):
            counts = Counter()
            for member in selection.members.values():
                if member["source"] == sid:
                    counts |= Counter(member["records"])
            stem = "sources/" + stable_id(sid)[:24]
            path = safe_path(work, stem + ".bam")
            header_path = safe_path(work, stem + ".header.json")
            records = selection.records[sid]
            chosen = [r for r in records if r.digest in counts]
            exported = compact_header(header, chosen) if header_policy == "compact" else header
            exported = _write_bam(path, exported, records, counts)
            write_json(header_path, header)
            manifest["sources"][sid] = dict(identity=selection.recipe["sources"][sid], bam=stem + ".bam",
                index=stem + ".bam.bai", original_header=stem + ".header.json", exported_header=exported,
                records=dict(sorted(counts.items())), record_multiset_sha256=stable_id(dict(counts)),
                record_count=sum(counts.values()))
        return _finish(work, manifest, size_budget)


def _verify_index(path, index):
    import pysam
    try:
        with pysam.AlignmentFile(path, index_filename=str(index)) as bam:
            observed = Counter(r.to_string() for r in bam if r.reference_id >= 0)
            indexed = Counter(r.to_string() for contig in bam.references for r in bam.fetch(contig))
            if observed != indexed:
                raise IntegrityError(f"Index does not enumerate the stored mapped records: {index.name}")
    except (ValueError, OSError) as error:
        raise IntegrityError(f"Invalid fixture BAM/index {path.name}: {error}") from error


def verify_bundle(directory, *, sha256=None):
    """Verify bytes, recipe/membership, original record multiplicity and indexes offline.

    Pin ``sha256`` to the manifest hash when loading an externally supplied data
    version. Internal checks detect corruption; an unpinned manifest is not an
    authenticity signature. This function needs neither SAMtools nor network.
    """
    import pysam
    root = Path(directory)
    if sha256 and digest(root / "manifest.json") != sha256:
        raise IntegrityError("Pinned bundle manifest checksum mismatch")
    manifest = json.loads((root / "manifest.json").read_text())
    if manifest.get("schema_version") != 1 or manifest.get("kind") != "osteosarc-fixture-bundle":
        raise IntegrityError("Unsupported fixture bundle schema")
    files = manifest["files"]
    verify_manifest_files(root, files)
    actual = {p.relative_to(root).as_posix() for p in root.rglob("*") if p.is_file() or p.is_symlink()}
    if actual != set(files) | {"manifest.json"}:
        raise IntegrityError("Bundle contains unlisted or missing files")
    recipe = validate_recipe(json.loads(safe_path(root, "recipe.json").read_text()))
    if stable_id(recipe) != manifest["recipe_sha256"]:
        raise IntegrityError("Bundle recipe digest mismatch")
    if set(manifest["members"]) != set(recipe["members"]):
        raise IntegrityError("Bundle member set differs from recipe")
    if manifest.get("record_encoding") != RECORD_ENCODING:
        raise IntegrityError("Unsupported record identity encoding")
    total = sum(p.stat().st_size for p in root.rglob("*") if p.is_file())
    if total > manifest["size_budget_bytes"] or sum(f["size_bytes"] for f in files.values()) != manifest["total_size_bytes"]:
        raise IntegrityError("Bundle size budget or total differs")
    needed_sources = {m["source"] for m in manifest["members"].values() if m["status"] not in ("unresolved", "omitted")}
    if set(manifest["sources"]) != needed_sources:
        raise IntegrityError("Bundle source set differs from selected members")
    for sid, source in manifest["sources"].items():
        if source["identity"] != recipe["sources"][sid]:
            raise IntegrityError(f"Bundle source identity differs from recipe: {sid}")
        counts = Counter(source["records"])
        if sum(counts.values()) != source["record_count"] or stable_id(dict(counts)) != source["record_multiset_sha256"]:
            raise IntegrityError(f"Source multiset/count digest mismatch: {sid}")
        path, index = safe_path(root, source["bam"]), safe_path(root, source["index"])
        if record_multiset(path) != counts:
            raise IntegrityError(f"Source record multiset/multiplicity differs: {sid}")
        with pysam.AlignmentFile(path) as bam:
            if bam.header.to_dict() != source["exported_header"]:
                raise IntegrityError(f"Exported header differs: {sid}")
        original = json.loads(safe_path(root, source["original_header"]).read_text())
        expected_header = (compact_header(original, list(read_records(path))) if manifest["header_policy"] == "compact" else original)
        expected_header = copy.deepcopy(expected_header)
        expected_header.setdefault("HD", {})["SO"] = "coordinate"
        for key in ("SS", "GO"):
            expected_header["HD"].pop(key, None)
        # pysam omits empty header sections.
        expected_header = {k: v for k, v in expected_header.items() if v}
        if expected_header != source["exported_header"]:
            raise IntegrityError(f"Source header semantics differ: {sid}")
        union = Counter()
        for member in manifest["members"].values():
            if member["source"] == sid:
                union |= Counter(member["records"])
        if union != counts:
            raise IntegrityError(f"Shared source union differs from members: {sid}")
        _verify_index(path, index)
    for name, member in manifest["members"].items():
        declared = recipe["members"][name]
        if member["source"] != declared["source"] or member["target"] != declared["target"]:
            raise IntegrityError(f"Member target/source differs: {name}")
        if set(member["reasons"]) != set(member["records"]):
            raise IntegrityError(f"Member reasons differ from records: {name}")
        if declared["policy"]["kind"] == "exact" and declared["policy"].get("encoding", RECORD_ENCODING) == RECORD_ENCODING and not declared.get("context_regions") and not declared.get("retain_partners"):
            expected = Counter(declared["policy"].get("records", {}))
            if declared["policy"].get("duplicate_policy") == "identical-record-once":
                expected = Counter(dict.fromkeys(expected, 1))
            if Counter(member["records"]) != expected:
                raise IntegrityError(f"Member differs from pinned exact recipe: {name}")
    for name, exported in manifest.get("exports", {}).items():
        path = safe_path(root, exported["path"])
        if exported["format"] == "bam":
            if record_multiset(path) != Counter(manifest["members"][name]["records"]):
                raise IntegrityError(f"Export record multiset differs: {name}")
            _verify_index(path, safe_path(root, exported["index"]))
    return manifest


def list_bundle(directory):
    """Return member status, scope, record counts and reasons after verification."""
    return verify_bundle(directory)["members"]


def export_bundle(directory, destination, *, members=None, format="bam", size_budget=DEFAULT_SIZE_BUDGET):
    """Publish a self-contained copy with named indexed BAM or legacy SAM exports.

    SAM export preserves SAM field values but cannot promise binary float/tag
    fidelity. BAM export verifies the lossless record multiset. Parent lineage
    and the original portable source pool travel with every export.
    """
    import pysam
    if format not in ("bam", "sam", "sam.gz"):
        raise ValueError("Export format must be bam, sam or sam.gz")
    root = Path(directory)
    original = verify_bundle(root)
    members = sorted(original["members"] if members is None else members)
    with _publication(destination) as work:
        for filename in original["files"]:
            target = safe_path(work, filename)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(safe_path(root, filename), target)
        manifest = copy.deepcopy(original)
        manifest["parent"] = dict(manifest_sha256=digest(root / "manifest.json"), parent=original.get("parent"))
        for name in members:
            member = manifest["members"][name]
            if member["status"] in ("unresolved", "omitted"):
                continue
            source = manifest["sources"][member["source"]]
            records = list(read_records(safe_path(root, source["bam"])))
            relative = f"members/{name}.{format}"
            target = safe_path(work, relative)
            if target.exists():
                raise FileExistsError(target)
            if format == "bam":
                _write_bam(target, source["exported_header"], records, member["records"])
                manifest["exports"][name] = dict(format=format, path=relative, index=relative + ".bai", fidelity=RECORD_ENCODING)
            else:
                by_id = {r.digest: r for r in records}
                text = str(pysam.AlignmentHeader.from_dict(source["exported_header"]))
                text += "".join((by_id[key].read.to_string() + "\n") * n for key, n in sorted(member["records"].items()))
                target.parent.mkdir(parents=True, exist_ok=True)
                target.write_bytes(gzip.compress(text.encode(), mtime=0) if format == "sam.gz" else text.encode())
                manifest["exports"][name] = dict(format=format, path=relative, fidelity="sam-text-v1")
        return _finish(work, manifest, size_budget)


def generate_bundle(recipe, destination, *, sources=None, cache=None, dataset=None, **pack_options):
    """Acquire declared pinned archives or indexed windows, select and pack.

    Explicit ``sources`` take precedence and are checked against archive hashes.
    A source ``archive`` requires url/sha256/size_bytes; it is a bounded historical
    input, not permission to download the original full BAM. Otherwise source
    ``identity.url`` and ``index`` plus declared regions drive indexed acquisition.
    """
    recipe = validate_recipe(recipe)
    cache = cache if isinstance(cache, Cache) else dataset.cache if dataset is not None else Cache(cache)
    inputs = dict(sources or {})
    archive_receipts = {}
    for sid, source in recipe["sources"].items():
        active = [m for m in recipe["members"].values() if m["source"] == sid
                  and m["policy"]["kind"] != "omitted" and recipe["targets"][m["target"]]["kind"] != "unresolved"]
        if not active or sid in inputs:
            continue
        if "archive" in source:
            archive = source["archive"]
            receipt = cache.fetch(archive["url"], sha256=archive["sha256"], size=archive["size_bytes"],
                                  max_bytes=pack_options.get("size_budget", DEFAULT_SIZE_BUDGET))
            inputs[sid] = cache.path(receipt)
            archive_receipts[sid] = receipt.to_dict()
        else:
            regions = [Region(**r) for r in source.get("regions", [])]
            if not regions:
                regions = [Region(**r) for m in active for r in m.get("regions", []) + m.get("context_regions", [])]
            options = dict(source.get("acquisition", {}))
            if "filters" in options:
                from .reads import ReadFilter
                options["filters"] = ReadFilter(**options["filters"])
            if dataset is None:
                inputs[sid] = extract_reads(source["identity"]["url"], regions, cache=cache,
                                            index=source.get("index"), snapshot_id=source.get("snapshot_id"), **options)
            else:
                asset = dataset.asset(source["identity"].get("key", source["identity"]["url"]))
                for key in ("id", "key", "url", "size", "modified"):
                    if key in source["identity"] and getattr(asset, key) != source["identity"][key]:
                        raise IntegrityError(f"Snapshot source identity changed: {sid}/{key}")
                if source.get("snapshot_id", dataset.id) != dataset.id:
                    raise IntegrityError(f"Snapshot identity differs: {sid}")
                inputs[sid] = dataset.extract_reads(asset, regions, **options)
    selection = select_fixtures(recipe, inputs)
    for sid, receipt in archive_receipts.items():
        selection.receipts[sid]["archive_acquisition"] = receipt
    return pack_bundle(selection, destination, **pack_options)
