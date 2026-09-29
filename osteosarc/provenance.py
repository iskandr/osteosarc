"""Reviewed evidence overlap, distinct from biological sample identity."""

import copy
import json
from collections import Counter
from functools import lru_cache
from pathlib import Path

import pysam

from .errors import IntegrityError


@lru_cache(maxsize=1)
def _registry():
    return json.loads((Path(__file__).with_name("data") / "evidence-overlaps.json").read_text())

def evidence_overlaps(url):
    """Known regional overlaps for this exact processing-product URL.

    No sample ID is rewritten. Empty results mean no reviewed overlap is listed,
    not that this source is biologically independent of every other source.
    """
    return copy.deepcopy([group for group in _registry()["groups"]
                          if any(s["url"] == url for s in group["sources"])])


def overlap_groups(files):
    urls = {file.url for file in files}
    return copy.deepcopy([dict(group, selected_urls=sorted(urls & {s["url"] for s in group["sources"]}))
                          for group in _registry()["groups"]
                          if len(urls & {s["url"] for s in group["sources"]}) > 1])


def require_no_known_overlaps(files):
    """Fail before counting known overlapping products as separate evidence."""
    groups = overlap_groups(files)
    if groups:
        raise IntegrityError("Sources have known overlapping evidence: " + ", ".join(g["id"] for g in groups)
                             + ". Keep their labels and receipts; select one explicitly with "
                             "files.without_known_overlaps(prefer=[URL]). See Osteosarc #100.")
    return files


def without_known_overlaps(files, *, prefer=()):
    """Keep explicitly preferred products from overlap groups; report exclusions.

    Preference selects evidence for an analysis, not the correct biological label.
    Both underlying source objects and all their original metadata remain intact.
    """
    prefer = (prefer,) if isinstance(prefer, str) else prefer
    preferred = {files[name].url for name in prefer}
    excluded = []
    for group in overlap_groups(files):
        chosen = preferred & set(group["selected_urls"])
        if len(chosen) != 1:
            raise ValueError(f"Choose exactly one preferred source in overlap group {group['id']}")
        kept, = chosen
        excluded += [dict(url=url, kept_url=kept, evidence_overlap=group["id"], issue=group["issue"])
                     for url in group["selected_urls"] if url != kept]
    removed = {entry["url"] for entry in excluded}
    result = type(files)((file for file in files if file.url not in removed), source=dict(files.source,
        excluded_evidence_overlaps=[*files.source.get("excluded_evidence_overlaps", []), *excluded]))
    return require_no_known_overlaps(result)


def audit_alignment_overlap(path_a, path_b):
    """Compare regional BAMs by record multiset, including base qualities.

    Ignore all optional tags (including RG), but keep QNAME, flags, mapping
    quality, sequence, qualities, CIGAR, positions, mate positions and TLEN.
    Normalize chr prefixes and M/MT for comparison only; do not rewrite files.
    This describes these BAM subsets, not the complete source libraries.
    """
    def contig(name):
        name = name.removeprefix("chr") if name else None
        return "MT" if name in ("M", "MT") else name

    def read(path):
        records, names = Counter(), set()
        with pysam.AlignmentFile(path) as bam:
            groups = bam.header.to_dict().get("RG", [])
            for r in bam:
                names.add(r.query_name)
                key = (r.query_name, r.flag, contig(r.reference_name), r.reference_start,
                       r.mapping_quality, r.cigarstring, contig(r.next_reference_name),
                       r.next_reference_start, r.template_length, r.query_sequence,
                       tuple(r.query_qualities) if r.query_qualities is not None else None)
                records[key] += 1
        return records, names, groups

    a, an, ag = read(path_a)
    b, bn, bg = read(path_b)
    return dict(records_a=sum(a.values()), records_b=sum(b.values()),
                names_a=len(an), names_b=len(bn), shared_names=len(an & bn), union_names=len(an | bn),
                shared_records_excluding_tags=sum((a & b).values()),
                different_a=sum((a - b).values()), different_b=sum((b - a).values()),
                read_groups_a=ag, read_groups_b=bg)
