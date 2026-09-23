"""Historical exact-record cohort transport format, migrated from Vaxrank.

Original implementation: OpenVax Vaxrank contributors (Apache-2.0). Frozen
cohort assignments and reference/prediction metadata are recipe inputs; this
module never invokes reconstruction, prediction or ranking.
"""
import gzip
import hashlib
import json
import os
import shutil
import struct
import tempfile
import zipfile
from collections import Counter, defaultdict
from dataclasses import asdict
from pathlib import Path, PurePosixPath

import pysam

from . import Asset, Dataset, Region, SampleClaim, digest, extract_reads
from . import __version__ as osteosarc_version
from .bundles import DEFAULT_SIZE_BUDGET, safe_path


def record_digest(read, *, text_only=False):
    """Include native float bits which SAM text formatting rounds away."""
    identity = [read.to_string()]
    if not text_only:
        for tag, value, kind in read.get_tags(with_value_type=True):
            if kind in ("f", "d"):
                identity.append([tag, kind, struct.pack("<" + kind, value).hex()])
            elif kind == "B" and value.typecode in ("f", "d"):
                # Explicit little-endian encoding is portable between hosts.
                identity.append([tag, kind, value.typecode,
                                 struct.pack("<" + value.typecode * len(value), *value).hex()])
    return hashlib.sha256(json.dumps(identity, separators=(",", ":")).encode()).hexdigest()


def retrieval_regions(records, assembly, *, pin_lengths=True):
    """Minimal point cover of alignment spans, including explicit mate records.

    Select the last base of the earliest-ending uncovered interval, per contig.
    Interval stabbing greedily minimizes the number of retrieval points; it
    does not change the required read list or infer biological coverage.
    """
    intervals = defaultdict(list)
    lengths = {}
    for read in records:
        if read.is_unmapped or read.reference_end is None:
            raise ValueError("Indexed retrieval requires mapped selected records")
        intervals[read.reference_name].append((read.reference_start, read.reference_end))
        lengths[read.reference_name] = read.header.get_reference_length(read.reference_name) if pin_lengths else None
    regions = []
    for contig, spans in sorted(intervals.items()):
        anchor = -1
        for start, end in sorted(spans, key=lambda span: span[1]):
            if start > anchor:
                anchor = end - 1
                regions.append([contig, anchor, anchor + 1, assembly, lengths[contig]])
    return regions


def select_records(subset, cohort):
    """Select exact records with multiplicity, preserving the reviewed order."""
    wanted = Counter(cohort["records"])
    found = defaultdict(list)
    with subset.open() as source:
        for read in source:
            key = record_digest(read, text_only=cohort["format"] != "bam")
            if key in wanted:
                found[key].append(read)
    observed = Counter({key: len(values) for key, values in found.items()})
    if observed != wanted:
        raise ValueError("%s: selected record multiplicity changed (missing=%s, extra=%s)" % (
            cohort["path"], sum((wanted - observed).values()), sum((observed - wanted).values())))
    return [found[key].pop(0) for key in cohort["records"]]


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def write_cohort(root, cohort, records, recipe):
    path = safe_path(root, cohort["path"])
    path.parent.mkdir(parents=True, exist_ok=True)
    if cohort["format"] == "fusion":
        if len(records) != len(cohort["record_metadata"]):
            raise ValueError("Fusion record metadata count differs: " + cohort["path"])
        data = json.loads((recipe / "fusion" / path.name.removesuffix(".gz")).read_text())
        data["original_records"] = [dict(sam=r.to_string(), **metadata)
                                    for r, metadata in zip(records, cohort["record_metadata"])]
        raw = (json.dumps(data, indent=2, sort_keys=True) + "\n").encode()
        path.write_bytes(gzip.compress(raw, mtime=0))
    else:
        header = pysam.AlignmentHeader.from_text(gzip.decompress(
            (recipe / "headers" / (cohort["header"] + ".sam.gz")).read_bytes()).decode())
        # IDs are numeric in BAM: require an identical reference dictionary
        # before copying records into the reviewed fixture header.
        if records and records[0].header.references != header.references:
            raise ValueError("Source reference dictionary changed: " + cohort["path"])
        if records and records[0].header.lengths != header.lengths:
            raise ValueError("Source reference lengths changed: " + cohort["path"])
        if records:
            selected_groups = {r.get_tag("RG") for r in records if r.has_tag("RG")}
            original_groups = {g["ID"]: g for g in records[0].header.to_dict().get("RG", [])
                               if g["ID"] in selected_groups}
            fixture_groups = {g["ID"]: g for g in header.to_dict().get("RG", [])
                              if g["ID"] in selected_groups}
            if original_groups != fixture_groups:
                raise ValueError("Source read-group definitions changed: " + cohort["path"])
        if cohort["format"] == "bam":
            with pysam.AlignmentFile(path, "wb", header=header) as output:
                for read in records:
                    output.write(read)
            pysam.index(str(path))
        else:
            sam = str(header) + "".join(r.to_string() + "\n" for r in records)
            path.write_bytes(gzip.compress(sam.encode(), mtime=0))


def update_manifests(root):
    """Update transport checksums, retaining independent biological expectations."""
    base = root / "osteosarc"
    manifest = json.loads((base / "manifest.json").read_text())
    for data in manifest["datasets"].values():
        data["sam_sha256"] = hashlib.sha256(gzip.decompress((base / data["file"]).read_bytes())).hexdigest()
    write_json(base / "manifest.json", manifest)
    shared = base / "shared-v1"
    manifest = json.loads((shared / "manifest.json").read_text())
    manifest["assets"] = [dict(filename=p.name, sha256=digest(p), size_bytes=p.stat().st_size,
                               bundle_path="osteosarc/shared-v1/" + p.name)
                          for p in sorted(shared.iterdir()) if p.suffix in (".bam", ".bai")]
    for case in manifest["cases"]:
        with pysam.AlignmentFile(shared / case["bam"]) as bam:
            case["selected_record_count"] = sum(1 for _ in bam)
        case["historical_selection"] = {key: case.pop(key) for key in (
            "selection", "source_receipt_sha256", "source_region_sha256", "selection_row_sha256")}
        case["selection"] = "One reviewed template for indexed retrieval; NTF3 retains compound AG>GT evidence"
    manifest["historical_import"] = manifest.pop("upstream")
    manifest["acquisition"] = "osteosarc; current extraction receipts are in ../../provenance.json"
    manifest["data_version"] = "minimal-vaccine-rna-v2"
    manifest["scope"] = "49 indexed-retrieval cases; one explicit template per case; not coverage or VAF"
    write_json(shared / "manifest.json", manifest)
    selection = base / "selection_validation"
    manifest = json.loads((selection / "isovar/manifest.json").read_text())
    manifest["files"] = {name: digest(selection / "isovar" / name) for name in manifest["files"]}
    write_json(selection / "isovar/manifest.json", manifest)
    prediction = json.loads((selection / "predictions_manifest.json").read_text())
    # Keep the original prediction-generation input receipt. This separate
    # digest identifies the repackaged, record-equivalent test inputs.
    prediction["bundled_input_manifest_sha256"] = digest(selection / "isovar/manifest.json")
    write_json(selection / "predictions_manifest.json", prediction)


def generate_cohort_bundle(recipe, cache, output):
    recipe, output = Path(recipe), Path(output)
    if output.exists() or output.is_symlink():
        raise FileExistsError(output)
    plan = json.loads(gzip.decompress((recipe / "selection.json.gz").read_bytes()))
    catalog = json.loads((recipe / "catalog.json").read_text())
    if catalog["snapshot"]["id"] != plan["snapshot_id"] or catalog["corrections"] is not False:
        raise ValueError("Wrong osteosarc catalog or correction policy")
    native = catalog["variants"]
    groups = defaultdict(list)
    for cohort in plan["cohorts"]:
        groups[cohort["source"]].append(cohort)
    provenance = dict(schema_version=1, snapshot_id=plan["snapshot_id"], corrections=False,
                      osteosarc_version=osteosarc_version, recipe_files={p.relative_to(recipe).as_posix(): digest(p)
                                    for p in sorted(recipe.rglob("*")) if p.is_file()},
                      variants=native, sources={}, cohorts=[])
    with tempfile.TemporaryDirectory(prefix="sid-bundle-build-") as temporary:
        root = Path(temporary)
        shutil.copytree(recipe / "support", root, dirs_exist_ok=True)
        for number, (url, cohorts) in enumerate(groups.items(), 1):
            source = catalog["assets"][url]
            fields = dict(source["asset"])
            fields["claims"] = tuple(SampleClaim(**c) for c in fields["claims"])
            fields["index_urls"] = tuple(fields["index_urls"])
            asset = Asset(**fields)
            pinned_index = source["index_receipt"]
            index = cache.path(cache.fetch(pinned_index["url"], sha256=pinned_index["sha256"],
                                           size=pinned_index["size"]))
            regions = {tuple(r) for c in cohorts for r in c["regions"]}
            print("[%d/%d] %s: %d required records" % (
                number, len(groups), asset.key, sum(len(c["records"]) for c in cohorts)), flush=True)
            subset = extract_reads(asset, [Region(*r) for r in sorted(regions)], cache=cache,
                                   index=str(index), snapshot_id=plan["snapshot_id"], fetch_pairs=False)
            provenance["sources"][url] = dict(asset=asdict(asset), index=pinned_index, extraction=subset.receipt)
            for cohort in cohorts:
                records = select_records(subset, cohort)
                write_cohort(root, cohort, records, recipe)
                provenance["cohorts"].append({k: v for k, v in cohort.items()
                                             if k not in ("header", "record_metadata")})
        update_manifests(root)
        write_json(root / "provenance.json", provenance)
        members = sorted(p for p in root.rglob("*") if p.is_file())
        manifest = {p.relative_to(root).as_posix(): dict(sha256=digest(p), size=p.stat().st_size)
                    for p in members}
        write_json(root / "bundle.json", dict(schema_version=1, files=manifest))
        # Fixed member timestamps/order/mode; volatile acquisition receipts are
        # preserved honestly, so a fresh online acquisition may differ in bytes.
        output.parent.mkdir(parents=True, exist_ok=True)
        if output.exists():
            raise FileExistsError(output)
        with tempfile.TemporaryDirectory(prefix=".sid-package-", dir=output.parent) as package_directory:
            staged = Path(package_directory) / output.name
            with zipfile.ZipFile(staged, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
                for path in sorted(root.rglob("*")):
                    if path.is_file():
                        info = zipfile.ZipInfo(path.relative_to(root).as_posix(), (1980, 1, 1, 0, 0, 0))
                        info.compress_type = zipfile.ZIP_DEFLATED
                        info.external_attr = 0o100644 << 16
                        archive.writestr(info, path.read_bytes())
            if staged.stat().st_size > plan.get("size_budget_bytes", 64 * 1024 * 1024):
                raise ValueError("Cohort bundle exceeds size budget")
            # Same-filesystem hard linking publishes atomically and exclusively,
            # including when another writer creates output during acquisition.
            os.link(staged, output)
        print("Wrote %s (%d bytes)" % (output, output.stat().st_size), flush=True)



def _extract_bundle(archive, destination):
    """Verify an explicit allowlist before materializing into an empty directory."""
    destination = Path(destination)
    if destination.is_symlink() or not destination.is_dir() or any(destination.iterdir()):
        raise ValueError("Bundle destination must be an empty, real directory")
    with zipfile.ZipFile(archive) as source:
        if sum(info.file_size for info in source.infolist()) > DEFAULT_SIZE_BUDGET:
            raise ValueError(f"Bundle exceeds {DEFAULT_SIZE_BUDGET} byte size budget")
        names = source.namelist()
        if len(names) != len(set(names)):
            raise ValueError("Duplicate bundle member")
        manifest = json.loads(source.read("bundle.json"))
        if manifest.get("schema_version") != 1:
            raise ValueError("Unsupported Sid bundle schema")
        expected = manifest["files"]
        if sum(v["size"] for v in expected.values()) > 64 * 1024 * 1024:
            raise ValueError("Bundle exceeds 64 MiB size budget")
        if set(names) != set(expected) | {"bundle.json"}:
            raise ValueError("Missing or unexpected bundle member")
        for name in names:
            path = PurePosixPath(name)
            if (path.is_absolute() or ".." in path.parts or not path.parts
                    or "\\" in name or ":" in name or str(path) != name):
                raise ValueError("Unsafe bundle path: " + name)
            info = source.getinfo(name)
            if info.is_dir() or (info.external_attr >> 16) & 0o170000 == 0o120000:
                raise ValueError("Bundle contains a non-file member")
            if name != "bundle.json" and info.file_size != expected[name]["size"]:
                raise ValueError("Bundle size mismatch: " + name)
        for name in names:
            path = destination / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(source.read(name))
            if name != "bundle.json" and digest(path) != expected[name]["sha256"]:
                raise ValueError("Bundle checksum mismatch: " + name)
    return destination



def pin_catalog(recipe, cache):
    plan = json.loads(gzip.decompress((recipe / "selection.json.gz").read_bytes()))
    dataset = Dataset.open(plan["snapshot_name"], cache=cache, offline=cache.offline, corrections=False)
    if dataset.id != plan["snapshot_id"]:
        raise ValueError("Snapshot differs from the reviewed selection recipe")
    variants = dataset.variants("all")
    catalog = dict(snapshot=dataset.manifest, corrections=False, variants={}, assets={})
    for variant_id in sorted({v for c in plan["cohorts"] for v in c["variants"]}):
        catalog["variants"][variant_id] = {k: v for k, v in asdict(variants[variant_id]).items()
                                           if k != "annotations"}
    for url in sorted({c["source"] for c in plan["cohorts"]}):
        asset = dataset.asset(url)
        index = dataset.asset(asset.index_urls[0])
        path = dataset.download(index)
        receipt = cache.fetch(index.url, sha256=digest(path), size=path.stat().st_size)
        catalog["assets"][url] = dict(asset=asdict(asset), index=asdict(index), index_receipt=receipt.to_dict())
    (recipe / "catalog.json").write_text(json.dumps(catalog, indent=2, sort_keys=True) + "\n")



def extract_bundle(archive, destination):
    """Verify into private staging, then atomically replace an empty directory."""
    destination = Path(destination)
    if destination.is_symlink() or not destination.is_dir() or any(destination.iterdir()):
        raise ValueError("Bundle destination must be an empty, real directory")
    with tempfile.TemporaryDirectory(dir=destination.parent, prefix=".cohort-") as temporary:
        staging = Path(temporary) / "verified"
        staging.mkdir()
        _extract_bundle(archive, staging)
        if destination.is_symlink() or any(destination.iterdir()):
            raise ValueError("Bundle destination changed during extraction")
        os.replace(staging, destination)
    return destination
