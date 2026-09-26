"""Helpers from Isovar's historical exact-SAM fixtures, kept because Isovar imports them.

read_json, write_json, sam_digest, sam_regions, minimal_header, segment_key and
select_window_segments; new test data uses osteosarc.bundles and osteosarc.shared.
Original implementation: Isovar contributors, Apache-2.0.
"""

import gzip
import json
import os
import tempfile
from hashlib import sha256
from pathlib import Path


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


def segment_key(read):
    """Historical source-local RG/QNAME/paired-segment identity, without inference."""
    return (read.get_tag("RG") if read.has_tag("RG") else "", read.query_name,
            read.flag & 0xc0 if read.is_paired else 0)
