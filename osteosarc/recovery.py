"""Bounded indexed mate/SA recovery; pointers never substitute for observations."""

from __future__ import annotations

import os
import tempfile
from collections import Counter
from dataclasses import asdict, dataclass

from .cache import Cache, digest, file_lock, stable_id, write_json
from .errors import IntegrityError
from .models import Region
from .reads import _cached_subset, extract_reads
from .records import read_records


@dataclass(frozen=True)
class RecoveryPolicy:
    """Limits apply to the seed union and all indexed partner queries together.

    All available seed/context records survive. Added windows contribute only
    actually matched partners. No full-file scan, synthetic sequence, trimming
    or completeness claim is made, even when every SA pointer was resolved.
    """

    mates: bool = True
    supplementary: bool = True
    max_rounds: int = 4
    max_intervals: int = 128
    max_bases: int = 1_000_000
    max_records: int = 100_000

    def __post_init__(self):
        for name in ("max_rounds", "max_intervals", "max_bases", "max_records"):
            if type(getattr(self, name)) is not int or getattr(self, name) < 1:
                raise ValueError(f"{name} must be a positive integer")


def _leads(record, policy):
    r = record.read
    base = dict(rg=record.template[0], qname=r.query_name)
    leads, problems = [], []
    if r.query_sequence is None:
        problems.append(dict(record=record.digest, reason="missing SEQ"))
    if policy.mates and r.is_paired:
        if r.mate_is_unmapped or r.next_reference_id < 0 or r.next_reference_start < 0:
            problems.append(dict(record=record.digest, reason="unmapped or unavailable mate placement"))
        elif record.segment not in (64, 128):
            problems.append(dict(record=record.digest, reason="ambiguous mate segment flags"))
        else:
            leads.append(dict(base, kind="mate", contig=r.next_reference_name,
                              start=r.next_reference_start, segment=192 ^ record.segment,
                              reverse=r.mate_is_reverse))
    if policy.supplementary and r.has_tag("SA"):
        for entry in r.get_tag("SA").split(";"):
            if not entry:
                continue
            try:
                contig, pos, strand, cigar, mapq, nm = entry.split(",")
                start = int(pos) - 1
                if start < 0 or strand not in ("+", "-") or not cigar or not 0 <= int(mapq) <= 255 or int(nm) < 0:
                    raise ValueError
                leads.append(dict(base, kind="SA", contig=contig, start=start,
                                  segment=record.segment, reverse=strand == "-", cigar=cigar,
                                  mapq=int(mapq), nm=int(nm)))
            except (TypeError, ValueError):
                problems.append(dict(record=record.digest, reason="malformed SA entry", entry=entry))
    return leads, problems


def _matches(record, lead):
    r = record.read
    if (record.template != (lead["rg"], lead["qname"]) or record.segment != lead["segment"]
            or r.reference_name != lead["contig"] or r.reference_start != lead["start"]
            or r.is_reverse != lead["reverse"]):
        return False
    if lead["kind"] == "mate":
        return not r.is_secondary and not r.is_supplementary
    return (r.cigarstring == lead["cigar"] and r.mapping_quality == lead["mapq"]
            and (not r.has_tag("NM") or r.get_tag("NM") == lead["nm"]))


def recover_reads(source, regions, *, policy=None, cache=None, **kwargs):
    """Extend extract_reads with bounded, source/RG/segment-scoped partner recovery.

    Returns a ReadSubset with every requested/visited interval, unresolved lead,
    cap, missing sequence and repeated/cyclic pointer in its receipt. A changed
    source/index fails rather than mixing acquisitions. Cache hits resume only
    verified complete intermediates. Seed record overflow fails explicitly.
    """
    import pysam
    policy = policy or RecoveryPolicy()
    if not isinstance(policy, RecoveryPolicy):
        policy = RecoveryPolicy(**policy)
    cache = cache if isinstance(cache, Cache) else Cache(cache)
    if kwargs.pop("fetch_pairs", False):
        raise ValueError("Choose RecoveryPolicy.mates instead of combining fetch_pairs with recovery")
    seed = extract_reads(source, regions, cache=cache, **kwargs)
    request = dict(schema_version=1, operation="recover_reads", seed=seed.receipt["request"],
                   seed_files=seed.receipt["files"], recovery=asdict(policy))
    directory = cache.workspace / "derived" / stable_id(request)
    with file_lock(cache.workspace / "locks" / (directory.name + ".lock")):
        if cached := _cached_subset(directory, request):
            return cached
        with seed.open() as bam:
            header = bam.header.to_dict()
        lengths = {s["SN"]: s["LN"] for s in header["SQ"]}
        seed_records = list(read_records(seed.path))
        counts = Counter(r.digest for r in seed_records)
        records = {r.digest: r for r in seed_records}
        reasons = {key: {"seed/context region"} for key in records}
        visited = [dict(r, round=0) for r in seed.receipt["resolved_regions"]]
        bases = sum(r["end"] - r["start"] for r in visited)
        if len(visited) > policy.max_intervals or bases > policy.max_bases or sum(counts.values()) > policy.max_records:
            raise IntegrityError("Seed acquisition exceeds recovery interval/base/record limit")
        queried = {(r["contig"], r["start"], r["end"]) for r in visited}
        receipts = [seed.receipt]
        requested, unresolved, problems, repeated = {}, {}, {}, set()
        handled_records, handled_leads = set(), set()
        limits = set()
        assembly = seed.receipt["resolved_regions"][0]["assembly"]

        def source_identity(receipt):
            # Header and index content are checked as well as remote HTTP identity.
            return [receipt["remote_identity"], receipt["header_receipt"]["files"],
                    receipt["request"]["source_sha256"], receipt["request"]["index_sha256"],
                    (receipt.get("index_receipt") or {}).get("sha256")]

        for round_number in range(1, policy.max_rounds + 2):
            leads = {}
            for key in sorted(records.keys() - handled_records):
                outgoing, issues = _leads(records[key], policy)
                handled_records.add(key)
                for issue in issues:
                    problems[stable_id(issue)] = issue
                for lead in outgoing:
                    identity = stable_id(lead)
                    requested.setdefault(identity, dict(lead, from_records=[]))["from_records"].append(key)
                    if identity in handled_leads:
                        repeated.add(identity)
                    else:
                        leads[identity] = lead
            if not leads:
                break
            # Already present seed/partner records satisfy pointers without queries.
            pending = {}
            for identity, lead in sorted(leads.items()):
                matches = [r for r in records.values() if _matches(r, lead)]
                if matches:
                    for r in matches:
                        reasons[r.digest].add("paired mate" if lead["kind"] == "mate" else "SA-linked partner")
                    handled_leads.add(identity)
                    if len(matches) > 1:
                        problems[identity] = dict(lead=identity, reason="multiple matching placements")
                else:
                    pending[identity] = lead
            if not pending:
                continue
            if round_number > policy.max_rounds:
                limits.add("max_rounds")
                for identity in pending:
                    unresolved[identity] = "max_rounds"
                break
            intervals = []
            for identity, lead in sorted(pending.items()):
                interval = (lead["contig"], lead["start"], lead["start"] + 1)
                if lead["contig"] not in lengths or lead["start"] >= lengths[lead["contig"]]:
                    unresolved[identity] = "absent contig or out-of-bounds placement"
                elif any(c == interval[0] and start <= interval[1] < end for c, start, end in queried):
                    unresolved[identity] = "no matching record in visited interval"
                elif interval not in intervals:
                    if len(queried) + len(intervals) >= policy.max_intervals:
                        unresolved[identity] = "max_intervals"
                        limits.add("max_intervals")
                    elif bases + len(intervals) >= policy.max_bases:
                        unresolved[identity] = "max_bases"
                        limits.add("max_bases")
                    else:
                        intervals.append(interval)
            handled_leads.update(pending)
            if not intervals:
                continue
            partner_regions = [Region(c, start, end, assembly, lengths[c]) for c, start, end in sorted(intervals)]
            subset = extract_reads(source, partner_regions, cache=cache, **kwargs)
            if source_identity(subset.receipt) != source_identity(seed.receipt):
                raise IntegrityError("Alignment/header/index changed across recovery acquisitions")
            receipts.append(subset.receipt)
            queried.update(intervals)
            bases += len(intervals)
            visited.extend(dict(asdict(r), round=round_number) for r in partner_regions)
            fetched = list(read_records(subset.path))
            if len(fetched) > policy.max_records:
                raise IntegrityError("Partner acquisition exceeds recovery record limit")
            found = Counter()
            found_records = {}
            for identity, lead in pending.items():
                matches = [r for r in fetched if _matches(r, lead)]
                if matches:
                    multiplicities = Counter(r.digest for r in matches)
                    found |= multiplicities  # union, not sum across overlapping pointers
                    found_records.update((r.digest, r) for r in matches)
                    for r in matches:
                        reasons.setdefault(r.digest, set()).add("paired mate" if lead["kind"] == "mate" else "SA-linked partner")
                    unresolved.pop(identity, None)
                    if len(multiplicities) > 1:
                        problems[identity] = dict(lead=identity, reason="multiple matching placements")
                else:
                    unresolved.setdefault(identity, "missing partner or conflicting alignment fields")
            if sum((counts | found).values()) > policy.max_records:
                raise IntegrityError("Retained records exceed recovery record limit")
            counts |= found
            records.update(found_records)
        directory.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=directory.parent, prefix=".recovery-") as temporary:
            from pathlib import Path
            work = Path(temporary)
            output = work / "reads.bam"
            ordered = sorted(records, key=lambda k: (
                records[k].read.reference_id if records[k].read.reference_id >= 0 else len(lengths),
                records[k].read.reference_start, k))
            with pysam.AlignmentFile(output, "wb", header=header) as out:
                for key in ordered:
                    for _ in range(counts[key]):
                        out.write(records[key].read)
            pysam.index(str(output))
            receipt = dict(request=request, files={name: digest(work / name) for name in ("reads.bam", "reads.bam.bai")},
                           records=sum(counts.values()), scope="bounded_mate_SA_context", complete_template=False,
                           status="truncated" if limits else "bounded", acquisition=receipts,
                           requested_leads=requested, visited_intervals=visited, unresolved=unresolved,
                           observations=list(problems.values()), repeated_or_cyclic_leads=sorted(repeated),
                           limits=sorted(limits), reasons={k: sorted(v) for k, v in sorted(reasons.items())},
                           pysam_version=pysam.__version__)
            write_json(work / "receipt.json", receipt)
            os.replace(work, directory)
        return _cached_subset(directory, request)
