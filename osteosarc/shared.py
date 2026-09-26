"""Test data shared by several libraries: one bundle, selected once from a spec.

A spec names the targets (every catalog variant of one snapshot, plus
extra alleles libraries test), sources that should cover every target, and
selection caps. Each library also supplies the records its current fixtures
hold (required records). build_shared_recipe() then writes a frozen fixture
recipe in which every member lists its exact records (bam-record-v1) and why
each was kept. generate_bundle() on that recipe rebuilds the same records from
the public BAMs, and fails if any of them changed upstream.

At each small-variant target, in each source that covers it:

1. Take the reads overlapping the target's allele window (osteosarc.alleles),
   with their mates.
2. Classify each template as alt, ref, other or uncallable.
3. Keep up to caps[class] templates of each class, in order of SHA-256 of the
   read group and read name, then the lowest-quality alt templates not
   already kept. A kept template keeps every record the extraction holds.

A source covers a target when the spec lists the source for every target, or
when a library's required records from that source overlap the target.

For a fusion or SV, a source keeps up to a cap of templates, in hash order,
that have an alignment within a window of every breakend: split reads,
discordant pairs and chimeric long reads. A source covers it when the spec
lists the source for structural targets, when a library's required records
lie within the window, or when osteosarc's SV catalogue saw the junction there.
Required records are matched by SAM text and kept exactly, as they are.
"""

from __future__ import annotations

import functools
import gzip
import hashlib
import json
import os
from collections import Counter, defaultdict
from pathlib import Path

from .alleles import CLASSES, allele_window, read_allele, template_allele
from .errors import CoordinateError, IntegrityError, OfflineError, OsteosarcError, SchemaError
from .models import Region
from .reads import resolve_regions
from .records import read_records
from .reference import reference_sequence

DEFAULT_CAPS = {"alt": 20, "ref": 10, "other": 5, "uncallable": 2}
#: Seconds one source's extraction may take; hundreds of regions with mates take minutes.
EXTRACTION_TIMEOUT = 3600
CONTEXT_FLANKS = (80, 600, 5000)


def read_json(path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as handle:
        return json.load(handle)


def template_order(template):
    """The hash order of a template: SHA-256 of its read group and read name."""
    rg, qname = template
    return hashlib.sha256(f"{rg or ''}\t{qname}".encode()).hexdigest()


def _window_quality(read, window):
    """The lowest base quality the read has across the window and its anchors."""
    qualities = read.query_qualities
    if qualities is None:
        return 255
    values = [qualities[q] for q, r in read.get_aligned_pairs(matches_only=True)
              if window.start - 1 <= r <= window.end]
    return min(values, default=255)


def _overlaps(read, window):
    end = read.reference_end or read.reference_start + 1
    return read.reference_name == window.contig and read.reference_start <= window.end and end >= window.start


class RecordIndex:
    """A source's records by position, so a window finds its records without a scan."""

    BIN = 10_000

    def __init__(self, records):
        self.records = list(records)
        self.templates = defaultdict(list)
        self.bins = defaultdict(list)
        for record in self.records:
            self.templates[record.template].append(record)
            read = record.read
            if read.reference_name is None or read.reference_start < 0:
                continue
            end = read.reference_end or read.reference_start + 1
            for b in range(read.reference_start // self.BIN, (end - 1) // self.BIN + 1):
                self.bins[(read.reference_name, b)].append(record)

    def overlapping(self, contig, start, end):
        """Records aligned within [start, end), each once."""
        seen, found = set(), []
        for b in range(start // self.BIN, max(start, end - 1) // self.BIN + 1):
            for record in self.bins.get((contig, b), ()):
                read = record.read
                if (id(record) not in seen and read.reference_start < end
                        and (read.reference_end or read.reference_start + 1) > start):
                    seen.add(id(record))
                    found.append(record)
        return found

    def touching(self, contig, start, end):
        """{template: all of its records} for templates with a record within [start, end)."""
        return {r.template: self.templates[r.template] for r in self.overlapping(contig, start, end)}


def select_allele_balanced(templates, window, *, caps=None, low_quality_alt=2):
    """Pick templates at one target: {template: (class, reason)}.

    templates maps (read group, read name) to that template's FixtureRecords.
    Only records overlapping the window are classified; a template that has
    none isn't considered.
    """
    caps = dict(DEFAULT_CAPS, **(caps or {}))
    by_class, quality = defaultdict(list), {}
    for template, records in templates.items():
        overlapping = [r for r in records if _overlaps(r.read, window)]
        if not overlapping:
            continue
        shown = template_allele([read_allele(r.read, window) for r in overlapping])
        by_class[shown].append(template)
        if shown == "alt":
            quality[template] = min(_window_quality(r.read, window) for r in overlapping)
    chosen = {}
    for name in CLASSES:
        for template in sorted(by_class[name], key=template_order)[:caps[name]]:
            chosen[template] = (name, f"{name} template (hash order)")
    rest = sorted((t for t in by_class["alt"] if t not in chosen), key=lambda t: (quality[t], template_order(t)))
    for template in rest[:low_quality_alt]:
        chosen[template] = ("alt", "alt template (lowest base quality)")
    return chosen, {name: len(by_class[name]) for name in CLASSES}


def match_required(records, lines):
    """Exact record counts for SAM lines, matched against a source's records by SAM text.

    records is a list, or a RecordIndex (then only records near the lines are
    rendered as text). lines may repeat, as duplicate records do. Raises
    IntegrityError when a line isn't there, appears fewer times than required,
    or matches records that differ in binary form.
    """
    if isinstance(records, RecordIndex):
        spans = {span for span in map(sam_span, lines) if span}
        near = {id(r): r for span in spans for r in records.overlapping(*span)}
        names = {line.split("\t", 1)[0] for line in lines if not sam_span(line)}
        near.update((id(r), r) for r in records.records if names and r.read.query_name in names)
        records = list(near.values())
    by_text = defaultdict(Counter)
    for record in records:
        by_text[record.read.to_string()][record.digest] += 1
    counts, missing = Counter(), []
    for line, n in Counter(lines).items():
        found = by_text.get(line)
        if not found or sum(found.values()) < n:
            missing.append(line.split("\t", 1)[0])
            continue
        if len(found) != 1:
            raise IntegrityError(f"SAM text of {line.split(chr(9), 1)[0]} matches records that differ in binary form")
        counts[next(iter(found))] += n
    if missing:
        raise IntegrityError(f"{len(missing)} required records aren't in the source, e.g. {missing[:3]}")
    return counts


def match_named(records, names):
    """Every record of the named reads (read names, any read group), with its multiplicity."""
    counts, found = Counter(), set()
    for record in (records.records if isinstance(records, RecordIndex) else records):
        if record.read.query_name in names:
            counts[record.digest] += 1
            found.add(record.read.query_name)
    if missing := sorted(set(names) - found):
        raise IntegrityError(f"{len(missing)} required reads aren't in the source, e.g. {missing[:3]}")
    return counts


def match_records(records, wanted):
    """Pinned records ({bam-record-v1 checksum: count}) a source must still hold, matched by checksum."""
    available = Counter(r.digest for r in (records.records if isinstance(records, RecordIndex) else records))
    if missing := Counter(wanted) - available:
        raise IntegrityError(f"{len(missing)} pinned records aren't in the source, e.g. {sorted(missing)[:3]}")
    return Counter(wanted)


def required_spans(subset):
    """Where a required subset's records lie: its SAM lines' spans, or the regions it names."""
    if "sam" in subset:
        return [span for span in map(sam_span, subset["sam"]) if span]
    return [tuple(region) for region in subset.get("regions", ())]


def merge_spans(spans):
    """Overlapping or touching (contig, start, end) spans merged, sorted: the same bases, fewer spans."""
    merged = []
    for contig, start, end in sorted(spans):
        if merged and merged[-1][0] == contig and start <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], end)
        else:
            merged.append([contig, start, end])
    return merged


def sam_span(line):
    """The zero-based reference span (contig, start, end) of a SAM line, or None if unplaced."""
    fields = line.split("\t")
    contig, position, cigar = fields[2], int(fields[3]), fields[5]
    if contig == "*" or position == 0:
        return None
    length, number = 0, ""
    for char in cigar:
        if char.isdigit():
            number += char
            continue
        if char in "MDN=X":
            length += int(number)
        number = ""
    return contig, position - 1, position - 1 + max(1, length)


def _context_window(contig, position, ref, alt, assembly, *, cache, reference_length=None):
    """The allele window, with the reference context it was found in."""
    for flank in CONTEXT_FLANKS:
        start = max(0, position - 1 - flank)
        end = position - 1 + len(ref) + flank
        sequence = reference_sequence(contig, start, end, assembly, cache=cache, reference_length=reference_length)
        try:
            window = allele_window(contig, position, ref, alt, sequence, start)
        except CoordinateError as error:
            if "too short" in str(error) or "doesn't cover" in str(error):
                continue
            raise
        context = dict(start=start, end=end, sha256=hashlib.sha256(sequence.encode()).hexdigest())
        return window, context
    raise CoordinateError(f"{contig}:{position}: a repeat runs past {CONTEXT_FLANKS[-1]} reference bases")


def spec_targets(spec, dataset, *, cache=None):
    """The spec's small-variant targets as recipe targets, and their allele windows."""
    cache = cache or dataset.cache
    targets, windows = {}, {}

    def small(name, contig, position, ref, alt, assembly, reference, label, reference_length=None):
        target = dict(kind="small_variant", assembly=assembly, coordinates="one-based", contig=contig,
                      position=position, ref=ref, alt=alt, reference=reference, label=label)
        if reference_length:
            target["reference_length"] = reference_length
        try:
            window, context = _context_window(contig, position, ref, alt, assembly, cache=cache,
                                              reference_length=reference_length)
        except CoordinateError as error:
            target["window"] = dict(unavailable=str(error))
        else:
            windows[name] = window
            target["window"] = dict(start=window.start, end=window.end, ref=window.ref, alt=window.alt,
                                    context=context)
        targets[name] = target

    chosen = spec["targets"].get("ids")
    for variant in (dataset.variants("all", ids=chosen) if chosen is not None
                    else dataset.variants(spec["targets"].get("catalog", "site"))):
        if variant.status != "ready":
            targets[variant.id] = dict(kind="unresolved", label="current",
                                       reason=f"{variant.status} in snapshot {dataset.name}")
            continue
        chrom, position, ref, alt = variant.allele
        small(variant.id, chrom, position, ref, alt, variant.assembly,
              dict(source="osteosarc", snapshot_id=dataset.id, variant_id=variant.id), "current")
    for extra in spec["targets"].get("extra", []):
        if extra["name"] in targets:
            raise SchemaError(f"Extra target {extra['name']!r} repeats a catalog ID; give it its own name")
        small(extra["name"], extra["contig"], extra["position"], extra["ref"], extra["alt"], extra["assembly"],
              dict(source="library", variant_id=extra.get("variant_id"), used_by=extra.get("used_by", []),
                   note=extra.get("reason", "")),
              extra["label"], extra.get("reference_length"))
    return targets, windows


def structural_targets(spec):
    """The spec's fusions and SVs as recipe targets, their breakends, and the
    RNA sources osteosarc's SV catalogue saw each one's junction in."""
    from .fixtures import load_panel
    from .sv_candidates import load_sv_candidates
    targets, breakends, observed_in = {}, {}, defaultdict(set)
    candidates = None
    for entry in spec["targets"].get("structural", []):
        name, origin = entry["name"], entry.get("from", {})
        if "panel" in origin:
            base = load_panel(origin["panel"])[origin["id"]]
            reference = dict(source="osteosarc", panel=origin["panel"], id=origin["id"])
        elif "sv_candidates" in origin:
            candidates = candidates or load_sv_candidates()["targets"]
            base = candidates[origin["sv_candidates"]]
            reference = dict(source="osteosarc", panel="sv-candidates-v1", id=origin["sv_candidates"])
            observed_in[name].update(e["source_url"] for e in base.get("rna_evidence", ())
                                     if e.get("status") == "adjacency_geometry_observed")
        else:
            base, reference = entry, dict(source="library", note=entry.get("reason", ""))
        reference["used_by"] = entry.get("used_by", [])
        ends = [dict(contig=e["contig"], position=e["position"], orientation=e.get("orientation"))
                for e in base.get("breakends", ())]
        if base.get("kind") == "unresolved" or len(ends) < 2:
            targets[name] = dict(kind="unresolved", label=entry["label"], reference=reference,
                                 reason=base.get("reason") or "Fewer than two breakends: no partner window to "
                                                               "select reads by")
            continue
        targets[name] = dict(kind="sv", assembly=base.get("assembly", "GRCh38"), coordinates="zero-based-interbase",
                             breakends=ends, reference=reference, label=entry["label"])
        breakends[name] = ends
    return targets, breakends, observed_in


#: How far a gap's ends may sit from the breakends it joins (alignments of a junction wobble).
JUNCTION_SLACK = 10


def _joins(records, breakends, pad):
    """Whether a template's alignments join every breakend, as a split read, a
    discordant pair or a chimeric long read would.

    Every breakend's window (pad bases either side) must hold aligned bases of
    the template, but no single aligned block may run across all of them. The
    template must then be split (a supplementary alignment or an SA tag), a pair
    the aligner didn't call proper, or a read whose intron or deletion runs from
    one breakend to another. So a proper pair on either side of the breakends, or
    a read spliced from an exon near one to an exon near another, doesn't count.
    """
    reads = [r.read for r in records if not r.read.is_unmapped]
    windows = [(contig, max(0, position - pad), position + pad) for contig, position in breakends]
    blocks = [(r.reference_name, start, end) for r in reads for start, end in r.get_blocks()]

    def touches(block, window):
        return block[0] == window[0] and block[1] < window[2] and block[2] > window[1]
    if not all(any(touches(b, w) for b in blocks) for w in windows):
        return False
    if any(all(touches(b, w) for w in windows) for b in blocks):
        return False
    split = any(r.is_supplementary or r.has_tag("SA") for r in reads)
    discordant = any(r.is_paired and not r.is_proper_pair and not r.mate_is_unmapped for r in reads)
    return split or discordant or any(_junction(r, breakends) for r in reads)


def _junction(read, breakends):
    """Whether one of the read's introns or deletions runs from one breakend to another."""
    here = [p for c, p in breakends if c == read.reference_name]
    position = read.reference_start
    for op, length in read.cigartuples or ():
        if op in (2, 3):  # D, N
            left, right = position, position + length
            if any(i != j and abs(left - a) <= JUNCTION_SLACK and abs(right - b) <= JUNCTION_SLACK
                   for i, a in enumerate(here) for j, b in enumerate(here)):
                return True
        if op in (0, 2, 3, 7, 8):  # M, D, N, =, X consume the reference
            position += length
    return False


def select_breakend_templates(templates, breakends, *, pad=1000, cap=50):
    """Templates whose alignments join every breakend (see _joins): {template: reason}.

    breakends are (contig, position) pairs. Up to cap, in hash order.
    """
    joined = [t for t, records in templates.items() if _joins(records, breakends, pad)]
    chosen = sorted(joined, key=template_order)[:cap]
    return {t: "joins the breakends (hash order)" for t in chosen}, len(joined)


def bundle_fixtures(bundle):
    """A bundle's library fixtures as required subsets pinned by checksum, so the
    next bundle keeps exactly their records once libraries no longer keep their
    own copies: {name: dict(consumer, source, records, regions, description)}.

    Each fixture is planned from the regions its target records; bundles built
    before targets recorded them (openvax-v1) use the spans of their records.
    """
    from .bundles import verify_bundle
    bundle = Path(bundle)
    manifest = verify_bundle(bundle)
    recipe = read_json(bundle / "recipe.json")
    subsets, unplanned = {}, defaultdict(list)
    for name, member in sorted(recipe["members"].items()):
        target = recipe["targets"][member["target"]]
        if target["kind"] != "fixture":
            continue
        url = recipe["sources"][member["source"]]["identity"].get("url")
        if member["policy"]["kind"] != "exact" or not target.get("consumer") or not url:
            raise SchemaError(f"Can't carry {name}: only exact fixtures of a named library, from a source "
                              "with a URL, carry forward")
        subsets[name] = dict(consumer=target["consumer"], source=url, records=dict(manifest["members"][name]["records"]),
                             description=target.get("description", name))
        if "regions" in target:
            subsets[name]["regions"] = target["regions"]
        else:
            unplanned[member["source"]].append(name)
    for sid, names in unplanned.items():
        wanted = {key for name in names for key in subsets[name]["records"]}
        spans = {r.digest: (r.read.reference_name, r.read.reference_start, r.read.reference_end)
                 for r in _source_records(bundle, manifest, sid) if r.digest in wanted and not r.read.is_unmapped}
        for name in names:
            subsets[name]["regions"] = merge_spans(spans[key] for key in subsets[name]["records"] if key in spans)
    return subsets


def merge_required(carried, paths):
    """Required subsets for the next bundle: the carried ones, except that each
    library with a required file (even one listing nothing) gets that file's
    subsets in place of all it had."""
    documents = [read_json(path) for path in paths]
    required = load_required(paths)
    fresh = {document["consumer"] for document in documents}
    kept = {name: subset for name, subset in carried.items() if subset["consumer"] not in fresh}
    if clash := sorted(set(kept) & set(required)):
        raise SchemaError(f"Carried and required subsets share names: {clash[:3]}")
    return {**kept, **required}


def _source_records(bundle, manifest, sid):
    from .bundles import safe_path
    return read_records(safe_path(Path(bundle), manifest["sources"][sid]["bam"]))


def load_required(paths):
    """Each library's required records, as {subset name: dict(consumer, source, sam, names or records)}.

    A subset lists its records' SAM lines ("sam"), names whole reads ("names"),
    or pins records by bam-record-v1 checksum and count ("records", as carried
    from a previous bundle); names and records come with the regions they lie
    in ("regions": [contig, start, end], zero-based).
    """
    subsets = {}
    for path in paths:
        document = read_json(path)
        for name, subset in document["subsets"].items():
            if name in subsets:
                raise SchemaError(f"Required subset {name!r} is named twice")
            kinds = [kind for kind in ("sam", "names", "records") if kind in subset]
            if len(kinds) != 1 or ("names" in subset and not subset.get("regions")) or (
                    "records" in subset and subset["records"] and not subset.get("regions")):
                raise SchemaError(f"Required subset {name!r} needs one of sam, names with regions, "
                                  "or records with regions")
            subsets[name] = dict(subset, consumer=document["consumer"])
    return subsets


def _source_entry(dataset, file, header, regions, label):
    from .reads import assembly_from_header
    assembly = assembly_from_header(header)
    if assembly is None:
        raise IntegrityError(f"Can't tell the assembly of {file.key} from its header")
    samples = file.samples
    return dict(identity=dict(id=file.id, key=file.key, url=file.url, size=file.size, modified=file.modified),
                assembly=assembly, sample=samples[0] if len(samples) == 1 else None,
                library=file.resolved("library"), product=file.key, label=label,
                acquisition=dict(fetch_pairs=True, timeout=EXTRACTION_TIMEOUT), snapshot_id=dataset.id,
                regions=[dict(contig=r.contig, start=r.start, end=r.end, assembly=r.assembly,
                              **({"reference_length": r.reference_length} if r.reference_length else {}))
                         for r in regions])


def _plan_source(url, dataset, targets, windows, breakends, observed_in, mine, everywhere, sv_everywhere, pad):
    """What one source covers, and the regions to extract for it; None if nothing."""
    from .reads import assembly_from_header
    file = dataset.file(url)
    header = dataset.inspect_alignment(file).header
    assembly = assembly_from_header(header)
    lengths = {row["SN"]: row["LN"] for row in header.get("SQ", [])}
    mito = {lengths[c] for c in ("chrM", "MT", "M", "chrMT") if c in lengths}

    def region(contig, start, end, reference_length=None):
        # GRCh37 mitochondria come in two lengths; say which this source has.
        if contig in ("chrM", "MT", "M", "chrMT") and reference_length is None:
            reference_length = next(iter(mito), None)
        return Region(contig, start, end, assembly, reference_length)
    spans = [span for s in mine.values() for span in required_spans(s)]
    covered = []
    for name, window in windows.items():
        if normalize(targets[name]["assembly"]) != assembly:
            continue
        if targets[name].get("reference_length") and targets[name]["reference_length"] not in mito:
            continue  # the other mitochondrial sequence
        if url in everywhere or any(_near(span, window) for span in spans):
            covered.append(name)
    sv_covered = [name for name, ends in breakends.items()
                  if normalize(targets[name]["assembly"]) == assembly
                  and (url in sv_everywhere or url in observed_in.get(name, ())
                       or any(_near(span, _BreakendWindow(e, pad)) for span in spans for e in ends))]
    # A variant's window with the aligned base on each side (an insertion's window can be empty).
    variant_regions = {n: region(windows[n].contig, max(0, windows[n].start - 1), windows[n].end + 1,
                                 targets[n].get("reference_length")) for n in covered}
    breakend_regions = {n: [region(e["contig"], max(0, e["position"] - pad), e["position"] + pad)
                            for e in breakends[n]] for n in sv_covered}
    wanted = list(variant_regions.values())
    wanted += [r for rs in breakend_regions.values() for r in rs]
    wanted += [region(contig, start, end) for contig, start, end in spans]
    if not wanted:
        return None
    label = dataset.short_name(file)  # names members as osteosarc reads --to names files
    return dict(file=file, label=label, mine=mine, covered=covered, sv_covered=sv_covered, pad=pad,
                source=_source_entry(dataset, file, header, resolve_regions(wanted, header), label),
                # Each target's window and breakends, in this source's contig names.
                variant_regions={n: resolve_regions([r], header)[0] for n, r in variant_regions.items()},
                breakend_regions={n: [resolve_regions([r], header)[0] for r in rs]
                                  for n, rs in breakend_regions.items()})


def _select_source(plan, index, windows, targets, caps, low_quality_alt, cap):
    """The members one source contributes, and any required fixtures it lacks."""
    members, fixtures, problems = {}, {}, []
    available = Counter(r.digest for r in index.records)

    def exact(chosen):
        kept, reasons = Counter(), {}
        for template, reason in chosen.items():
            for record in index.templates[template]:
                kept[record.digest] = available[record.digest]
                reasons.setdefault(record.digest, set()).add(reason)
        return dict(version=1, kind="exact", records=dict(sorted(kept.items())),
                    reasons={k: sorted(v) for k, v in sorted(reasons.items())})
    label = plan["label"]
    for name in plan["covered"]:
        where = plan["variant_regions"][name]
        window = windows[name]
        window = type(window)(where.contig, window.start, window.end, window.ref, window.alt)
        candidates = index.touching(where.contig, window.start - 1, window.end + 1)
        chosen, seen = select_allele_balanced(candidates, window, caps=caps, low_quality_alt=low_quality_alt)
        members[f"{label}.{name}"] = dict(target=name, source=label, observed=seen, policy=dict(
            exact({t: why for t, (_, why) in chosen.items()}), reason="allele-balanced selection"))
    for name in plan["sv_covered"]:
        ends = [(r.contig, r.start, r.end) for r in plan["breakend_regions"][name]]
        candidates = index.touching(*ends[0])
        for end in ends[1:]:
            near = index.touching(*end)
            candidates = {t: rs for t, rs in candidates.items() if t in near}
        breakends = [(r.contig, r.end - plan["pad"]) for r in plan["breakend_regions"][name]]
        chosen, joined = select_breakend_templates(candidates, breakends, pad=plan["pad"], cap=cap)
        members[f"{label}.{name}"] = dict(target=name, source=label, observed=dict(joined=joined),
                                          policy=dict(exact(chosen), reason="templates joining the breakends"))
    for subset_name, required_subset in sorted(plan["mine"].items()):
        try:
            if "sam" in required_subset:
                counts = match_required(index, required_subset["sam"])
            elif "records" in required_subset:
                counts = match_records(index, required_subset["records"])
            else:
                counts = match_named(index, set(required_subset["names"]))
        except IntegrityError as error:
            problems.append(f"{subset_name}: {error}")
            continue
        fixtures[subset_name] = (required_subset, counts)
    return members, fixtures, problems


def build_shared_recipe(spec, dataset, *, required=(), log=print, workers=6):
    """Select the spec's records and return its frozen recipe (see the module docs).

    Sources are extracted in parallel (workers at a time). dataset must be
    opened with offline=False the first time, to stream reads; extractions
    are cached, so rebuilding is quick and offline.
    """
    from concurrent.futures import ThreadPoolExecutor, as_completed
    spec = read_json(spec) if isinstance(spec, (str, Path)) else spec
    if spec.get("snapshot", {}).get("id") not in (None, dataset.id):
        raise IntegrityError(f"The spec pins snapshot {spec['snapshot']['id'][:12]}, not {dataset.id[:12]}")
    selection = spec.get("selection", {})
    caps = dict(DEFAULT_CAPS, **selection.get("caps", {}))
    low_quality_alt = selection.get("low_quality_alt", 2)
    targets, windows = spec_targets(spec, dataset)
    sv_targets, breakends, observed_in = structural_targets(spec)
    targets.update(sv_targets)
    structural = dict(dict(window=1000, cap=50), **selection.get("structural", {}))
    subsets = load_required(required) if not isinstance(required, dict) else required
    everywhere = set(spec.get("sources", {}).get("all_targets", []))
    sv_everywhere = set(spec.get("sources", {}).get("structural", []))
    # openvax-v1 also reads each RNA BAM the SV catalogue saw a junction in; a bundle
    # made from given files reads only those.
    observed = {url for found in observed_in.values() for url in found} if spec["sources"].get("observed", True) else set()
    urls = sorted(everywhere | sv_everywhere | {s["source"] for s in subsets.values()} | observed)
    recipe = dict(schema_version=1, id=spec["id"], kind="shared",
                  snapshot=dict(name=dataset.name, id=dataset.id),
                  selection=dict(classifier="osteosarc.alleles v1", caps=caps, low_quality_alt=low_quality_alt,
                                 structural=structural, order="SHA-256 of read group, tab, read name"),
                  aliases=spec["targets"].get("renamed", {}), targets=targets, sources={}, members={},
                  redistribution=spec.get("redistribution", {"license": "unresolved"}))
    plans = []
    for url in urls:
        mine = {n: s for n, s in subsets.items() if s["source"] == url}
        plan = _plan_source(url, dataset, targets, windows, breakends, observed_in, mine, everywhere,
                            sv_everywhere, structural["window"])
        if plan is None:
            if mine:
                raise IntegrityError(f"Required fixtures with no aligned records to find them by: {sorted(mine)}")
            key = dataset.file(url).key
            log(f"skip {key}: " + ("none of the targets is on its genome build" if url in everywhere | sv_everywhere
                                   else "no targets or library fixtures there"))
            continue
        if any(plan["label"] == other["label"] for other in plans):
            raise IntegrityError(f"Two sources share the name {plan['label']}")
        log(f"{plan['file'].key}: {plural(len(plan['covered']), 'variant')}, "
            f"{plural(len(plan['sv_covered']), 'SV or fusion', 'SVs or fusions')}"
            + (f", {plural(len(mine), 'library fixture')}" if mine else ""))
        plans.append(plan)

    # Every member's name, checked before any reads are fetched.
    names = {}
    for plan in plans:
        for name in plan["covered"] + plan["sv_covered"]:
            _add_member(names, f"{plan['label']}.{name}", None)
    for name in subsets:
        _add_member(names, name, None)

    def extract(plan):
        subset = dataset.extract_reads(plan["file"], [Region(**r) for r in plan["source"]["regions"]],
                                       **plan["source"]["acquisition"])
        return plan, subset

    problems, results = [], {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        for future in as_completed([pool.submit(extract, plan) for plan in plans]):
            plan, subset = future.result()
            index = RecordIndex(read_records(subset.path))
            members, fixtures, missing = _select_source(plan, index, windows, targets, caps, low_quality_alt,
                                                        structural["cap"])
            results[plan["label"]] = (plan, members, fixtures)
            problems += missing
            log(f"{plan['file'].key}: read {plural(len(index.records), 'record')}, "
                f"kept {plural(len(members) + len(fixtures), 'member')}")
    if problems:
        # Every source was read (and its extraction cached), so a rerun after fixing these is quick.
        raise IntegrityError(f"{len(problems)} required fixtures couldn't be matched:\n" + "\n".join(problems))
    def add(name, member):
        _add_member(recipe["members"], name, member)
    for label in sorted(results):  # the same recipe whatever order extractions finished in
        plan, members, fixtures = results[label]
        recipe["sources"][label] = plan["source"]
        for name, member in members.items():
            add(name, member)
        for subset_name, (required_subset, counts) in fixtures.items():
            target = "fixture:" + subset_name
            targets[target] = dict(kind="fixture", assembly=plan["source"]["assembly"],
                                   reference=dict(source="library", consumer=required_subset["consumer"]),
                                   consumer=required_subset["consumer"],
                                   description=required_subset.get("description", subset_name),
                                   # So a later bundle that carries this one plans it the same way.
                                   regions=merge_spans(required_spans(required_subset)))
            add(subset_name, dict(
                target=target, source=label,
                policy=dict(version=1, kind="exact", records=dict(sorted(counts.items())),
                            reason=f"required by {required_subset['consumer']}: {subset_name}")))
    return recipe


def _add_member(members, name, member):
    """Add a member under a name that can be a file name, refusing to replace another."""
    from .bundles import safe_path
    if name in members:
        raise IntegrityError(f"Two members share the name {name}")
    safe_path(Path("."), name)  # export writes the member to a file named after it
    members[name] = member


def plural(n, word, words=None):
    return f"{n:,} {word if n == 1 else words or word + 's'}"


def normalize(assembly):
    from .reads import normalize_assembly
    return normalize_assembly(assembly)


class _BreakendWindow:
    """A breakend's window, shaped like an AlleleWindow for _near."""

    def __init__(self, breakend, pad):
        self.contig, self.start, self.end = breakend["contig"], breakend["position"] - pad, breakend["position"] + pad


def _near(span, window, slack=1):
    contig, start, end = span
    bare = contig.removeprefix("chr")
    return (bare == window.contig.removeprefix("chr") or {bare, window.contig.removeprefix("chr")} <= {"M", "MT"}) \
        and start <= window.end + slack and end >= window.start - slack


# ---------------------------------------------------------------------------
# Published bundles: fetch and check
# ---------------------------------------------------------------------------

BUNDLES = Path(__file__).with_name("data") / "bundles"


def published(name):
    """Where a published bundle (such as openvax-v1) is, and its checksums."""
    path = BUNDLES / f"{name}.release.json"
    if not path.is_file():
        known = sorted(p.name.removesuffix(".release.json") for p in BUNDLES.glob("*.release.json"))
        raise KeyError(f"No bundle folder or published bundle {name!r}; published: {', '.join(known) or 'none'}")
    return json.loads(path.read_text())


#: The dataset's license, which bundles of its reads carry.
REDISTRIBUTION = {"license": "CC0-1.0", "source": "https://registry.opendata.aws/sid-osteosarc/"}


def make_bundle(dataset, to, *, variants=(), svs=(), files=(), caps=None, size_budget=64 * 1024 * 1024,
                log=None):
    """Make a bundle of test reads with openvax-v1's selection, and return its folder.

    For each variant, each file gives up to caps templates of each allele class
    (by default 20 alt, 10 ref, 5 other and 2 uncallable, plus the two
    lowest-quality alt templates); for each SV or fusion, up to 50 templates that
    join its breakends. Members are named FILE.TARGET, like the files osteosarc
    reads --to writes. See Dataset.make_bundle.
    """
    from .bundles import generate_bundle
    to = Path(to)
    if to.exists():
        raise FileExistsError(f"{to} already exists; bundles go in a new folder")
    spec = bundle_spec(dataset, to.name, variants=variants, svs=svs, files=files, caps=caps)
    recipe = build_shared_recipe(spec, dataset, required={}, log=log or (lambda text: None))
    generate_bundle(recipe, to, dataset=dataset, size_budget=size_budget)
    return to


def bundle_spec(dataset, name, *, variants=(), svs=(), files=(), caps=None):
    """The spec make_bundle builds from: these variants and SVs, read from these files only."""
    import difflib

    from .models import File, Sample
    ids = [variant if isinstance(variant, str) else variant.id for variant in _items(variants)]
    svs = list(_items(svs))
    if not ids and not svs:
        raise ValueError("Give the variants (IDs, or data.variants(...)) or SVs to take reads for")
    found = {variant.id: variant for variant in dataset.variants("all", ids=ids)} if ids else {}
    for vid in ids:
        if vid not in found:
            close = difflib.get_close_matches(vid, [v.id for v in dataset.variants("all")], n=3)
            raise ValueError(f"No variant {vid} in snapshot {dataset.name}"
                             + (f"; did you mean {', '.join(close)}?" if close else ""))
        if found[vid].status != "ready":
            raise ValueError(f"{vid} has no usable allele ({found[vid].status}); see osteosarc variants {vid}")
    urls = []
    for item in _items(files):
        chosen = ([f for f in item.files.select(kind="alignment") if f.index_urls] if isinstance(item, Sample)
                  else [item if isinstance(item, File) else dataset.file(item)])
        urls += [f.url for f in chosen if f.url not in urls]
    if not urls:
        raise ValueError("Give the indexed BAMs to take reads from: files, keys, or samples")
    return dict(id=name, snapshot=dict(name=dataset.name, id=dataset.id),
                selection=dict(caps=dict(DEFAULT_CAPS, **(caps or {}))),
                sources=dict(all_targets=urls, structural=urls if svs else [], observed=False),
                targets=dict(ids=ids, structural=[_sv_entry(sv) for sv in svs]),
                redistribution=REDISTRIBUTION)


def _items(value):
    """A single item or a collection of them, as a list."""
    from .models import File, Sample
    if value is None:
        return []
    if isinstance(value, (str, File, Sample)) or not hasattr(value, "__iter__"):
        return [value]
    return list(value)


def _sv_entry(sv):
    """A spec's structural entry for an SV candidate or SV regression ID."""
    from .fixtures import load_panel
    from .sv_candidates import load_sv_candidates
    if sv in load_sv_candidates()["targets"]:
        return dict(name=sv, label="sv", **{"from": {"sv_candidates": sv}})
    if sv in load_panel("sv-regressions-v1"):
        return dict(name=sv, label="sv", **{"from": {"panel": "sv-regressions-v1", "id": sv}})
    raise ValueError(f"No SV {sv}: give an SV candidate ID, from load_sv_candidates() (see the SV candidates "
                     "docs), or a regression target from load_panel('sv-regressions-v1')")


def bundle_file(bundle, member, *, format="bam", cache=None, offline=False):
    """One member of a bundle as a local file, for a test to read.

    bundle is a published bundle's name (such as openvax-v1, downloaded the
    first time) or a bundle folder. The member is exported once into the cache
    as an indexed BAM (or format="sam" or "sam.gz"), made read-only, and reused
    afterward, offline. A member already named x.bam is exported as x.bam.

        bam = osteosarc.bundle_file("openvax-v1", "topiary/osteosarc/bulk_star_t0.sam.gz")
    """
    import difflib
    import tempfile

    from .bundles import _write_members, export_target
    from .cache import Cache, file_lock
    cache = cache if isinstance(cache, Cache) else Cache(cache)
    if offline and not cache.offline:
        cache = Cache(cache.root, offline=True)
    folder = bundle_folder(bundle, cache=cache)
    stat = (folder / "manifest.json").stat()
    manifest_sha = _manifest_digest(str(folder), (stat.st_mtime_ns, stat.st_size))
    exports = cache.workspace / "exports" / f"{folder.name}-{manifest_sha[:16]}"
    target = export_target(exports / format, member, format)
    index = Path(str(target) + ".bai")
    if target.is_file() and (format != "bam" or index.is_file()):
        return target
    manifest = _verified_manifest(str(folder), manifest_sha)
    if member not in manifest["members"]:
        close = difflib.get_close_matches(member, manifest["members"], n=3)
        raise KeyError(f"{member} isn't a member of {bundle}" + (f"; did you mean {', '.join(close)}?" if close
                       else f"; osteosarc test-data list {bundle} lists them"))
    status = manifest["members"][member]["status"]
    if status in ("unresolved", "omitted"):
        raise KeyError(f"{member} has no reads in {bundle} ({status})")
    exports.mkdir(parents=True, exist_ok=True)
    with file_lock(exports / f".{format}.lock"):
        if not (target.is_file() and (format != "bam" or index.is_file())):
            with tempfile.TemporaryDirectory(dir=exports, prefix=".export-") as work:
                made = _write_members(folder, manifest, [member], work, format)[member]
                target.parent.mkdir(parents=True, exist_ok=True)
                # The index first, so a BAM in place always has its index beside it.
                for source, destination in ([(Path(str(made) + ".bai"), index)] if format == "bam" else []) + [
                        (made, target)]:
                    os.chmod(source, 0o444)
                    os.replace(source, destination)
    return target


@functools.lru_cache(maxsize=8)
def _manifest_digest(folder, stamp):
    from .cache import digest
    return digest(Path(folder) / "manifest.json")


@functools.lru_cache(maxsize=4)
def _verified_manifest(folder, manifest_sha):
    """A bundle's manifest, verified in full once per process."""
    from .bundles import verify_bundle
    return verify_bundle(folder, sha256=manifest_sha)


def bundle_folder(value, *, cache=None):
    """The folder of a published bundle's name (fetched and verified), or else of a bundle path.

    A published name always means the published bundle; write ./NAME for a
    local folder that shares its name.
    """
    text = str(value)
    if "/" not in text and os.sep not in text and (BUNDLES / f"{text}.release.json").is_file():
        return fetch_bundle(text, cache=cache)
    if Path(value).exists():
        return Path(value)
    return fetch_bundle(text, cache=cache)  # an unknown name raises, listing the published ones


def pack_release(bundle, archive):
    """Write a bundle as a reproducible .tar.gz and return its release record."""
    import tarfile

    from .bundles import verify_bundle
    from .cache import digest
    bundle, archive = Path(bundle), Path(archive)
    verify_bundle(bundle)
    with open(archive, "wb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed, \
            tarfile.open(fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT) as tar:
        for path in sorted(p for p in bundle.rglob("*") if p.is_file()):
            info = tar.gettarinfo(str(path), arcname=path.relative_to(bundle).as_posix())
            info.uid = info.gid = 0
            info.uname = info.gname = ""
            info.mtime = 0
            info.mode = 0o644
            with open(path, "rb") as handle:
                tar.addfile(info, handle)
    return dict(sha256=digest(archive), size_bytes=archive.stat().st_size,
                manifest_sha256=digest(bundle / "manifest.json"))


def fetch_bundle(name, *, cache=None, offline=False):
    """Download a published bundle into the cache, verify it, and return its directory.

    The bundle is verified in full when it arrives, and its files are made
    read-only; later calls reuse it offline, checking only its manifest. With
    offline=True it never downloads, and raises OfflineError if the bundle
    isn't cached yet.
    """
    import tarfile

    from .bundles import _publication, safe_path, verify_bundle
    from .cache import Cache, digest
    release = published(name)
    cache = cache if isinstance(cache, Cache) else Cache(cache)
    if offline and not cache.offline:
        cache = Cache(cache.root, offline=True)
    root = cache.workspace / "bundles" / f"{name}-{release['manifest_sha256'][:16]}"
    if not root.exists():
        try:
            receipt = cache.fetch(release["url"], sha256=release["sha256"], size=release["size_bytes"])
        except OfflineError as error:
            raise OfflineError(f"{name} isn't in the cache at {cache.root} yet: fetch it once with network "
                               "access, or set OSTEOSARC_CACHE to a cache that has it") from error
        try:
            with _publication(root) as work, tarfile.open(cache.path(receipt), "r:gz") as tar:
                for member in tar.getmembers():
                    if not member.isfile():
                        raise IntegrityError(f"Unexpected archive entry {member.name!r}")
                    target = safe_path(work, member.name)
                    target.parent.mkdir(parents=True, exist_ok=True)
                    with tar.extractfile(member) as source, open(target, "wb") as handle:
                        handle.write(source.read())
                verify_bundle(work, sha256=release["manifest_sha256"])
                for path in work.rglob("*"):
                    if path.is_file():
                        os.chmod(path, 0o444)
        except FileExistsError:
            pass  # another process published it first
    if digest(root / "manifest.json") != release["manifest_sha256"]:
        raise IntegrityError(f"{root} has changed since it was downloaded; delete it to download it again")
    return root


def sam_lines(value):
    """Every SAM record line anywhere in a JSON value: in strings (one line or
    several joined by newlines), lists and objects, whatever the field names.

    A line counts when it has at least 11 tab-separated fields with a numeric
    FLAG and POS; other strings (names, versions, URLs) and non-strings are
    skipped, so a pointer may lead to records or to anything that holds them.
    """
    if isinstance(value, str):
        for line in value.split("\n"):
            line = line.rstrip("\r")
            fields = line.split("\t")
            if len(fields) >= 11 and fields[1].isdigit() and fields[3].isdigit():
                yield line
    elif isinstance(value, dict):
        for item in value.values():
            yield from sam_lines(item)
    elif isinstance(value, list):
        for item in value:
            yield from sam_lines(item)


def resolve_pointer(document, pointer):
    """The part of a JSON document a JSON Pointer (RFC 6901) names; "" or "/" is all of it."""
    if pointer in ("", "/"):
        return document
    if not pointer.startswith("/"):
        raise SchemaError(f"JSON pointer {pointer!r} must start with /")
    here = ""
    for token in pointer[1:].split("/"):
        token = token.replace("~1", "/").replace("~0", "~")
        here += "/" + token
        if isinstance(document, list) and token.isdigit() and int(token) < len(document):
            document = document[int(token)]
        elif isinstance(document, dict) and token in document:
            document = document[token]
        else:
            raise SchemaError(f"JSON pointer {pointer!r}: nothing at {here}")
    return document


def _local_sam(value, root=Path(".")):
    """SAM text lines from a local fixture: a BAM/SAM/SAM.gz path, or
    {"json": path, "pointer": "/a/b"} for the SAM lines under a JSON Pointer."""
    import pysam
    if isinstance(value, dict):
        return list(sam_lines(resolve_pointer(read_json(root / value["json"]), value.get("pointer", ""))))
    path = root / value
    if path.suffix == ".gz" and path.name.endswith(".sam.gz"):
        import tempfile
        with tempfile.NamedTemporaryFile(suffix=".sam") as plain:
            plain.write(gzip.decompress(path.read_bytes()))
            plain.flush()
            with pysam.AlignmentFile(plain.name, "r", check_sq=False) as handle:
                return [read.to_string() for read in handle]
    with pysam.AlignmentFile(str(path), check_sq=False) as handle:
        return [read.to_string() for read in handle.fetch(until_eof=True)]


def check_fixtures(bundle, fixtures, *, root=Path("."), cache=None):
    """Compare a library's local fixtures with a bundle's members, as multisets of SAM text.

    bundle is a published bundle's name (such as openvax-v1, fetched with cache)
    or a bundle folder.
    fixtures maps member names to local files (see _local_sam), relative to root.
    Returns {member: dict(missing=n, extra=n)} for every member that differs, and
    {member: dict(error=why)} for one that isn't a member or can't be read.
    """
    from .bundles import verify_bundle
    bundle = bundle_folder(bundle, cache=cache)
    manifest = verify_bundle(bundle)
    by_source = {}
    problems = {}
    for name, value in sorted(fixtures.items()):
        if name not in manifest["members"]:
            problems[name] = dict(error="not a member of this bundle")
            continue
        member = manifest["members"][name]
        sid = member["source"]
        expected = Counter()
        if member["records"] and sid not in by_source:
            by_source[sid] = {record.digest: record.read.to_string()
                              for record in _source_records(bundle, manifest, sid)}
        for key, n in member["records"].items():
            expected[by_source[sid][key]] += n
        try:
            found = Counter(_local_sam(value, Path(root)))
        except (OsteosarcError, OSError, ValueError, KeyError, TypeError) as error:
            problems[name] = dict(error=f"can't read {value}: {error}")
            continue
        if found != expected:
            problems[name] = dict(missing=sum((expected - found).values()), extra=sum((found - expected).values()))
    return problems
