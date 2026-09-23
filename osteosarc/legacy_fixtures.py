"""Historical exact-SAM fixture format, migrated from openvax/isovar.

These adapters preserve reviewed SAM membership/order and header policy. They
promise SAM-text fidelity, not BAM auxiliary-bit equality. New binary bundles
use osteosarc.bundles. No Isovar import or sibling checkout is required.
Original implementation: Isovar contributors, Apache-2.0.
"""

import gzip
import json
import os
import tempfile
from collections import Counter
from hashlib import sha256
from pathlib import Path

from .bundles import safe_path


def read_json(path):
    with gzip.open(path, "rt") if str(path).endswith(".gz") else open(path) as handle:
        return json.load(handle)


def write_json(path, value):
    """Atomically publish deterministic historical gzip JSON without overwrite."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=path.parent, prefix=".legacy-") as temporary:
        staged = Path(temporary) / "data.gz"
        with staged.open("wb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
                handle.write(json.dumps(value, sort_keys=True, separators=(",", ":")).encode())
        if staged.stat().st_size > 64 * 1024 * 1024:
            raise ValueError("Historical fixture exceeds 64 MiB size budget")
        os.link(staged, path)  # atomic publication, fails if a destination exists


def sam_digest(line):
    return sha256(line.encode("ascii")).hexdigest()


def sam_regions(regions, assembly, reference_lengths=None):
    """Convert explicit 1-based inclusive SAM intervals to osteosarc Regions."""
    from osteosarc import Region
    result = []
    for region in regions:
        contig, span = region.rsplit(":", 1)
        start, end = map(int, span.split("-"))
        result.append(Region(contig, start - 1, end, assembly,
                             reference_length=(reference_lengths or {}).get(contig)))
    if not result:
        raise ValueError("An explicit nonempty list of locus intervals is required")
    return result


def asset_identity(asset):
    return {key: getattr(asset, key) for key in ("id", "key", "url", "size", "modified")}


def required_counts(recipe, source_id):
    """Repeated fixture use is shared; true source duplicates remain required."""
    counts = Counter()
    for fixture in recipe["fixtures"].values():
        if fixture["source"] == source_id:
            needed = Counter(fixture["records"])
            if fixture.get("record_references", False):
                needed = Counter(dict.fromkeys(needed, 1))
            counts |= needed
    return counts


def select_records(records, required):
    """Fail if even one required original record or duplicate is absent."""
    observed = Counter()
    selected = {}
    for record in records:
        line = record.to_string()
        checksum = sam_digest(line)
        if checksum in required:
            observed[checksum] += 1
            selected[checksum] = line
    missing = required - observed
    if missing:
        raise ValueError("Source is missing %d required SAM records: %s" % (
            sum(missing.values()), dict(missing)))
    return selected


def minimal_header(header, records):
    """Retain decoding identity, without shipping thousands of irrelevant PGs."""
    references, groups, programs = set(), set(), set()
    for line in records:
        fields = line.split("\t")
        references.update(v for v in (fields[2], fields[6]) if v not in ("*", "="))
        groups.update(v[5:] for v in fields[11:] if v.startswith("RG:Z:"))
        programs.update(v[5:] for v in fields[11:] if v.startswith("PG:Z:"))
        for tag in fields[11:]:
            if tag.startswith("SA:Z:"):
                references.update(row.split(",", 1)[0] for row in tag[5:].split(";") if row)
    sq = [row for row in header.get("SQ", []) if row["SN"] in references]
    if {row["SN"] for row in sq} != references:
        raise ValueError("Source header is missing a required reference sequence")
    # Some source records themselves have RG tags absent from their header.
    # Preserve that source condition instead of inventing read-group metadata.
    result = {"HD": {"VN": "1.6", "SO": "unknown"}, "SQ": sq}
    rg = [row for row in header.get("RG", []) if row["ID"] in groups]
    if rg:
        result["RG"] = rg
    programs.update(row["PG"] for row in rg if "PG" in row)
    by_id = {row["ID"]: row for row in header.get("PG", [])}
    pending = list(programs)
    while pending:
        previous = by_id.get(pending.pop(), {}).get("PP")
        if previous and previous not in programs:
            programs.add(previous)
            pending.append(previous)
    pg = [row for row in header.get("PG", []) if row["ID"] in programs]
    if pg:
        result["PG"] = pg
    return result


def generate(recipe, output, dataset, *, source_ids=None):
    """Regenerate only the exact test selections from indexed osteosarc data.

    Separate source files make an interrupted acquisition resumable. Every
    reused file is checked against the recipe and source snapshot identity.
    """
    from osteosarc import Region
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    chosen = sorted(recipe["sources"] if source_ids is None else source_ids)
    for identity in chosen:
        source = recipe["sources"][identity]
        asset = dataset.asset(source["asset"]["key"])
        if asset_identity(asset) != source["asset"]:
            raise ValueError("Snapshot asset changed: " + asset.key)
        required = required_counts(recipe, identity)
        path = safe_path(output, identity + ".json.gz")
        if path.exists():
            cached = read_json(path)
            validate_source(cached, required)
            if cached["asset"] != source["asset"] or cached["snapshot_id"] != dataset.id:
                raise ValueError("Acquisition snapshot/asset mismatch: " + identity)
            continue
        regions = [Region(**region) for region in source["regions"]]
        if not regions:
            raise ValueError("No bounded query regions: " + identity)
        subset = dataset.extract_reads(asset, regions)
        with subset.open() as handle:
            records = select_records(handle, required)
            header = minimal_header(handle.header.to_dict(), records.values())
        value = dict(asset=source["asset"], snapshot_id=dataset.id,
                     header=header, records=records,
                     duplicate_counts={k: n for k, n in required.items() if n > 1},
                     acquisition=subset.receipt)
        write_json(path, value)
        print(identity, len(records), "selected records", flush=True)
    return chosen


def validate_source(source, required):
    records = source["records"]
    if set(records) != set(required) or source["duplicate_counts"] != {
            k: n for k, n in required.items() if n > 1}:
        raise ValueError("Selected records differ from the test recipe")
    if any(sam_digest(line) != checksum for checksum, line in records.items()):
        raise ValueError("Original SAM checksum mismatch")


def pack(recipe, acquired, output):
    sources = {}
    for identity in sorted(recipe["sources"]):
        source = read_json(safe_path(acquired, identity + ".json.gz"))
        validate_source(source, required_counts(recipe, identity))
        if source["asset"] != recipe["sources"][identity]["asset"]:
            raise ValueError("Source asset differs from recipe")
        sources[identity] = source
    bundle = dict(schema_version=1, sources=sources, recipe_sha256=recipe_digest(recipe))
    write_json(output, bundle)
    return bundle


def verify(recipe, bundle):
    if bundle["recipe_sha256"] != recipe_digest(recipe):
        raise ValueError("Bundle recipe checksum mismatch")
    if set(bundle["sources"]) != set(recipe["sources"]):
        raise ValueError("Bundle source set differs from recipe")
    for identity, source in bundle["sources"].items():
        validate_source(source, required_counts(recipe, identity))
        if source["asset"] != recipe["sources"][identity]["asset"]:
            raise ValueError("Bundle source identity differs from recipe")
    return dict(fixtures=len(recipe["fixtures"]), sources=len(bundle["sources"]),
                unique_records=sum(len(s["records"]) for s in bundle["sources"].values()))


def recipe_digest(recipe):
    return sha256(json.dumps(recipe, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def export_fixture(name, output, *, recipe, bundle):
    """Export original SAM records in fixture order, including repeated records."""
    import pysam
    verify(recipe, bundle)
    fixture = recipe["fixtures"][name]
    source = bundle["sources"][fixture["source"]]
    with Path(output).open("x") as handle:
        handle.write(str(pysam.AlignmentHeader.from_dict(source["header"])))
        for checksum in fixture["records"]:
            handle.write(source["records"][checksum] + "\n")



def select_window_segments(bam, windows, *, available_context=False, duplicate_policy="preserve"):
    """Select original segments touching every window, with explicit legacy scope.

    This is regional path candidacy, never a claim of junction/ORF support.
    ``available_context`` also retains their records outside the seed windows
    in the supplied local archive. This explicitly scans only that local input.
    """
    from collections import defaultdict
    windows = list(windows)
    if not windows:
        raise ValueError("At least one explicit window is required")
    touches, entries = defaultdict(set), []
    # One pass avoids retrieval duplication but preserves true source repeats.
    bam.reset()
    for r in bam:
        if r.is_unmapped:
            continue
        key = segment_key(r)
        sides = {i for i, (contig, start, end) in enumerate(windows)
                 if r.reference_name == contig and r.reference_start < end
                 and start < r.reference_end}
        if sides or available_context:
            entries.append((r.to_string(), key))
        touches[key].update(sides)
    selected = {key for key, sides in touches.items() if sides == set(range(len(windows)))}
    lines = [line for line, key in entries if key in selected]
    if duplicate_policy == "identical-SAM-lines-once":
        lines = list(set(lines))
    elif duplicate_policy != "preserve":
        raise ValueError("Unknown historical duplicate policy")
    return sorted(lines), len(selected), len(entries)


def select_assigned_names(groups, preferred_names=(), *, limit=16, preferred_limit=32):
    """Execute the pinned historical QName-only stratum ordering (one source).

    ``groups`` contains producer-supplied name -> record-score mappings. This
    compatibility policy deliberately preserves old membership across RGs; new
    recipes should use source/RG/segment assignments. No support is inferred.
    """
    selected = set(preferred_names[:preferred_limit])
    for values in groups.values():
        selected.update(sorted(values, key=lambda name: (values[name], name))[:limit])
    return selected


def write_selected_names(source_bam, output, region, selected_names, *, complete=False):
    """Write the historical source-scoped QName selection, preserving occurrences."""
    import pysam
    if Path(output).exists():
        raise FileExistsError(output)
    records = []
    with pysam.AlignmentFile(source_bam) as source, pysam.AlignmentFile(output, "wb", header=source.header) as target:
        for ordinal, read in enumerate(source.fetch(*region)):
            if complete or read.query_name in selected_names:
                target.write(read)
                records.append(dict(regional_ordinal=ordinal, sam_sha256=sam_digest(read.to_string()), name=read.query_name))
    pysam.index(str(output))
    return records


def segment_key(read):
    """Historical source-local RG/QNAME/paired-segment identity, without inference."""
    return (read.get_tag("RG") if read.has_tag("RG") else "", read.query_name,
            read.flag & 0xc0 if read.is_paired else 0)
