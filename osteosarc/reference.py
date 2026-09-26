"""Short stretches of reference sequence, for checking alleles against the genome.

Sequence comes from the UCSC Genome Browser's public API and is cached like any
other download, so a stretch fetched once is available offline. GRCh38 is hg38
and GRCh37 is hg19. Mitochondrial names follow the length the caller states:
16569 is the revised Cambridge sequence (hg38's chrM, and GRCh37's MT) and
16571 is hg19's chrM.
"""

from __future__ import annotations

import json
from urllib.parse import quote

from .cache import Cache
from .errors import CoordinateError

UCSC = "https://api.genome.ucsc.edu/getData/sequence"
GENOMES = {"GRCh38": "hg38", "GRCh37": "hg19"}
MITOCHONDRIA = ("chrM", "MT", "M", "chrMT")


def ucsc_location(contig, assembly, *, reference_length=None):
    """The UCSC genome and chromosome name for a contig of an assembly."""
    from .reads import normalize_assembly
    assembly = normalize_assembly(assembly)
    if assembly not in GENOMES:
        raise CoordinateError(f"No UCSC genome for assembly {assembly!r}")
    if contig in MITOCHONDRIA:
        if reference_length == 16569 or (reference_length is None and assembly == "GRCh38"):
            return "hg38", "chrM"  # the revised Cambridge sequence, in either assembly
        if reference_length == 16571 and assembly == "GRCh37":
            return "hg19", "chrM"
        raise CoordinateError("Give the mitochondrial reference_length (16569 or 16571) to choose its sequence")
    return GENOMES[assembly], contig if contig.startswith("chr") else "chr" + contig


def reference_sequence(contig, start, end, assembly, *, cache=None, reference_length=None):
    """Reference bases [start, end), zero-based, as upper-case text."""
    if not 0 <= start < end or end - start > 100_000:
        raise CoordinateError("Ask for 1 to 100,000 reference bases at a time")
    genome, chrom = ucsc_location(contig, assembly, reference_length=reference_length)
    cache = cache if isinstance(cache, Cache) else Cache(cache)
    url = f"{UCSC}?genome={genome};chrom={quote(chrom)};start={start};end={end}"
    document = json.loads(cache.path(cache.fetch(url, max_bytes=1_000_000)).read_text())
    sequence = document.get("dna", "").upper()
    if len(sequence) != end - start:
        raise CoordinateError(f"UCSC returned {len(sequence)} bases for {genome} {chrom}:{start}-{end}")
    return sequence
