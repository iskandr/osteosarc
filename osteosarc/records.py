"""Lossless BAM record identities, independent of compression and reference IDs.

The v1 encoding retains the stored CIGAR, SEQ, QUAL and auxiliary bytes (including
floating-point payloads and integer widths). It ignores only the derived bin and
numeric reference IDs. Tag order is deliberately significant. See SAMv1 §4.2.
"""

from __future__ import annotations

import gzip
import hashlib
import struct
from collections import Counter
from dataclasses import dataclass

from .cache import stable_id
from .errors import IntegrityError

RECORD_ENCODING = "bam-record-v1"


def _read(handle, size):
    if size < 0:
        raise IntegrityError("Negative BAM block length")
    value = handle.read(size)
    if len(value) != size:
        raise IntegrityError("Truncated BAM record/header")
    return value


def _integer(handle):
    return struct.unpack("<i", _read(handle, 4))[0]


def _stored_records(path):
    """Each record's stored bytes, and a function giving its v1 identity."""
    with gzip.open(path, "rb") as handle:
        if _read(handle, 4) != b"BAM\x01":
            raise IntegrityError("Lossless identity requires BAM input")
        _read(handle, _integer(handle))
        references = []
        for _ in range(_integer(handle)):
            references.append(_read(handle, _integer(handle))[:-1].decode("utf-8"))
            _integer(handle)
        prefixes = {}  # by reference IDs: records share a few

        def identity(block):
            tid, pos, bin_mq_nl, flag_nc, length, mate_tid, mate_pos, tlen = struct.unpack("<iiIIiiii", block[:32])
            prefix = prefixes.get((tid, mate_tid))
            if prefix is None:
                if any(t < -1 or t >= len(references) for t in (tid, mate_tid)):
                    raise IntegrityError("Invalid BAM reference ID")
                names = [references[t] if t >= 0 else None for t in (tid, mate_tid)]
                prefix = prefixes[tid, mate_tid] = bytes.fromhex(stable_id([RECORD_ENCODING, names]))
            result = hashlib.sha256(prefix)
            result.update(struct.pack("<iIIiii", pos, bin_mq_nl & 0xffff, flag_nc, length, mate_pos, tlen))
            result.update(memoryview(block)[32:])  # no copies
            return result.hexdigest()
        while size := handle.read(4):
            if len(size) != 4:
                raise IntegrityError("Truncated BAM block size")
            block = _read(handle, struct.unpack("<i", size)[0])
            if len(block) < 32:
                raise IntegrityError("Invalid BAM core")
            yield block, identity


def bam_record_digests(path):
    """Yield SHA256 identities from original BAM bytes, without SAM conversion."""
    for block, identity in _stored_records(path):
        yield identity(block)


def record_multiset(path):
    """Return the exact stored record multiplicities of a local BAM."""
    return Counter(bam_record_digests(path))


def read_template(read):
    """A pysam read's template: its read group (or None) and name."""
    return read.get_tag("RG") if read.has_tag("RG") else None, read.query_name


@dataclass(frozen=True)
class FixtureRecord:
    read: object
    digest: str

    @property
    def template(self):
        return read_template(self.read)

    @property
    def segment(self):
        # Retain both bits for malformed/unusual source flags rather than repairing.
        return self.read.flag & 0xc0


def read_records(path, *, keep=None, digests=None):
    """Read local BAM records with identities computed from their stored bytes.

    keep (a test of a pysam read) and digests (a set of identities) choose which
    to read: those either chooses, or all if neither is given. Others aren't
    kept, and are hashed only to look for digests.
    """
    import pysam
    stored = _stored_records(path)
    with pysam.AlignmentFile(path, "rb") as bam:
        for read in bam:
            block, identity = next(stored, (None, None))
            if block is None:
                raise IntegrityError("BAM readers disagree on record count")
            if keep is None and digests is None:
                yield FixtureRecord(read, identity(block))
            elif keep is not None and keep(read):
                yield FixtureRecord(read, identity(block))
            elif digests is not None and (digest := identity(block)) in digests:
                yield FixtureRecord(read, digest)
        if next(stored, None) is not None:
            raise IntegrityError("BAM readers disagree on record count")
