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


def bam_record_digests(path):
    """Yield SHA256 identities from original BAM bytes, without SAM conversion."""
    with gzip.open(path, "rb") as handle:
        if _read(handle, 4) != b"BAM\x01":
            raise IntegrityError("Lossless identity requires BAM input")
        _read(handle, _integer(handle))
        references = []
        for _ in range(_integer(handle)):
            references.append(_read(handle, _integer(handle))[:-1].decode("utf-8"))
            _integer(handle)
        while size := handle.read(4):
            if len(size) != 4:
                raise IntegrityError("Truncated BAM block size")
            block = _read(handle, struct.unpack("<i", size)[0])
            if len(block) < 32:
                raise IntegrityError("Invalid BAM core")
            tid, pos, bin_mq_nl, flag_nc, length, mate_tid, mate_pos, tlen = struct.unpack("<iiIIiiii", block[:32])
            if any(t < -1 or t >= len(references) for t in (tid, mate_tid)):
                raise IntegrityError("Invalid BAM reference ID")
            names = [references[t] if t >= 0 else None for t in (tid, mate_tid)]
            core = struct.pack("<iIIiii", pos, bin_mq_nl & 0xffff, flag_nc, length, mate_pos, tlen)
            yield hashlib.sha256(bytes.fromhex(stable_id([RECORD_ENCODING, names])) + core + block[32:]).hexdigest()


def record_multiset(path):
    """Return the exact stored record multiplicities of a local BAM."""
    return Counter(bam_record_digests(path))


@dataclass(frozen=True)
class FixtureRecord:
    read: object
    digest: str

    @property
    def template(self):
        return (self.read.get_tag("RG") if self.read.has_tag("RG") else None, self.read.query_name)

    @property
    def segment(self):
        # Retain both bits for malformed/unusual source flags rather than repairing.
        return self.read.flag & 0xc0


def read_records(path):
    """Read local BAM records with identities computed from their stored bytes."""
    import pysam
    identities = iter(bam_record_digests(path))
    with pysam.AlignmentFile(path, "rb") as bam:
        for read in bam:
            yield FixtureRecord(read, next(identities))
        if next(identities, None) is not None:
            raise IntegrityError("BAM readers disagree on record count")
