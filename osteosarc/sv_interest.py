"""Shared, source-pinned SV nominations and their scoped evidence snapshot."""

import json
from pathlib import Path


def load_sv_interest():
    """Load the complete versioned SV interest catalogue without network I/O.

    Returns
    -------
    dict
        ``targets`` maps nomination IDs to fixture-compatible targets, original
        calls, annotation, gene-expression and per-source RNA evidence.
        ``sources`` retains URLs and SHA-256 pins. Every call returns a new
        copy. Single breakends and unresolved geometries remain explicit
        ``kind='unresolved'`` targets, rather than disappearing.

    Notes
    -----
    Nominations include aliases and components of complex events. The
    ``adjacency_group_id`` groups identical retained-flank geometry, not alleles
    or independent mutations. Do not sum its RNA counts. T2 gene TPM is not
    mutant-transcript or protein abundance; incomplete assays remain missing.
    Protein hypotheses and ranking are produced by Isovar/Topiary, not inferred
    from membership in this catalogue.
    """
    path = Path(__file__).with_name("data") / "sv_interest.json"
    return json.loads(path.read_text())
