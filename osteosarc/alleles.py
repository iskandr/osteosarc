"""Classify reads against a small variant by the bases they show across it.

Each read shows the reference allele, the alternate allele, something else, or
nothing that decides between them. The comparison covers the variant's minimal
edit, widened across any repeat the edit could slide within, so an indel that
an aligner placed elsewhere in a homopolymer still counts as the same allele.
It uses only the read's own alignment and the reference sequence around the
variant, which the caller supplies; this module never fetches anything and
never looks at base or mapping quality.
"""

from __future__ import annotations

from dataclasses import dataclass

from .errors import CoordinateError

#: What a read or template shows at a variant.
CLASSES = ("alt", "ref", "other", "uncallable")

_ALIGNED, _INSERTION, _DELETION, _SKIP, _SOFT_CLIP, _HARD_CLIP, _PAD, _EQUAL, _DIFF = range(9)


@dataclass(frozen=True)
class AlleleWindow:
    """The reference bases a read must span to show one allele or the other.

    start and end are zero-based, half-open, and cover the minimal edit
    widened across repeats. ref is the reference sequence of that span and alt
    the same span with the variant applied. A read must also align the base on
    each side of the span.
    """

    contig: str
    start: int
    end: int
    ref: str
    alt: str


def allele_window(contig, position, ref, alt, context, context_start):
    """The window for a VCF-style allele (one-based position) on a reference.

    context is the reference sequence beginning at zero-based context_start,
    and must reach at least one base past the widened window on each side.
    """
    ref, alt, context = ref.upper(), alt.upper(), context.upper()
    start = position - 1
    if not ref or not alt or ref == alt or any(c not in "ACGTN" for c in ref + alt):
        raise CoordinateError(f"Not a literal small variant: {ref}>{alt}")
    offset = start - context_start
    if offset < 1 or context[offset:offset + len(ref)] != ref:
        raise CoordinateError(f"{contig}:{position} REF {ref} doesn't match the reference context")
    # Trim shared bases to the minimal edit: reference bases [s, e) become `inserted`.
    s, e, inserted = start, start + len(ref), alt
    while e > s and inserted and ref[e - start - 1] == inserted[-1]:
        e, inserted = e - 1, inserted[:-1]
    while e > s and inserted and ref[s - start] == inserted[0]:
        s, inserted = s + 1, inserted[1:]
    lo, hi = s, e

    def base(i):
        j = i - context_start
        if j < 0 or j >= len(context):
            raise CoordinateError(f"{contig}:{position}: the reference context is too short for this repeat")
        return context[j]

    if s == e or not inserted:  # a pure insertion or deletion can slide within a repeat
        moved = context[s - context_start:e - context_start] if not inserted else inserted
        a, b, seq = s, e, moved
        while base(a - 1) == seq[-1]:
            a, b, seq = a - 1, b - 1, seq[-1] + seq[:-1]
        lo = a
        a, b, seq = s, e, moved
        while base(b) == seq[0]:
            a, b, seq = a + 1, b + 1, seq[1:] + seq[0]
        hi = b
    if lo < context_start + 1 or hi + 1 > context_start + len(context):
        raise CoordinateError(f"{contig}:{position}: the reference context doesn't cover the widened window")
    window_ref = context[lo - context_start:hi - context_start]
    window_alt = context[lo - context_start:s - context_start] + inserted + context[e - context_start:hi - context_start]
    return AlleleWindow(contig, lo, hi, window_ref, window_alt)


def read_allele(read, window):
    """What one alignment shows across the window: alt, ref, other or uncallable.

    Uncallable: unmapped, secondary or QC-failed; on another contig; not
    aligning the bases either side of the window; a splice (N) inside it; or an
    N base in what it shows. Otherwise the read's bases across the window,
    including insertions next to it, are compared with each allele.
    """
    if (read.is_unmapped or read.is_secondary or read.is_qcfail or read.reference_name != window.contig
            or read.query_sequence is None or read.cigartuples is None):
        return "uncallable"
    left, right = window.start - 1, window.end
    sequence = read.query_sequence.upper()
    shown, anchors, qpos, rpos = [], set(), 0, read.reference_start
    for op, length in read.cigartuples:
        if op in (_ALIGNED, _EQUAL, _DIFF):
            for i in range(length):
                r = rpos + i
                if r in (left, right):
                    anchors.add(r)
                elif left < r < right:
                    shown.append(sequence[qpos + i])
            qpos, rpos = qpos + length, rpos + length
        elif op == _INSERTION:
            if left <= rpos - 1 and rpos <= right:  # between the anchors
                shown.append(sequence[qpos:qpos + length])
            qpos += length
        elif op == _DELETION:
            rpos += length
        elif op == _SKIP:
            if rpos < right and rpos + length > left + 1:
                return "uncallable"
            rpos += length
        elif op == _SOFT_CLIP:
            qpos += length
    if anchors != {left, right}:
        return "uncallable"
    shown = "".join(shown)
    if "N" in shown:
        return "uncallable"
    return "alt" if shown == window.alt else "ref" if shown == window.ref else "other"


def template_allele(classes):
    """Combine the classes of one template's alignments.

    The template shows an allele when every alignment that shows something
    agrees; mates that disagree make it other.
    """
    shown = {c for c in classes if c != "uncallable"}
    if not shown:
        return "uncallable"
    return next(iter(shown)) if len(shown) == 1 else "other"
