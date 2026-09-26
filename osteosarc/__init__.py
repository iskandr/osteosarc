"""Browse and download the public osteosarc.com dataset, and make test data from it.

    from osteosarc import Dataset
    data = Dataset.sync()      # once: the website's metadata, about 57 MB
    data = Dataset.open()      # later: the newest snapshot, offline
    data                       # what's there, and how to get it

Test data:

    data.extract_reads(key, variants=[...], to="tests/data")      # a quick BAM of the reads there
    data.make_bundle("tests/data/mine", variants=[...], files=[...])  # a bundle anyone can rebuild
    osteosarc.bundle_file("openvax-v1", member)                   # one member, as a local BAM

Importing performs no network I/O.
"""

from .bundles import (  # noqa: F401
    export_bundle,
    generate_bundle,
    list_bundle,
    verify_bundle,
    verify_digest,
    verify_gzip_digests,
    verify_manifest_files,
)
from .cache import (
    Cache,
    digest,  # noqa: F401
)
from .curation import CORRECTIONS, Change, Correction, CurationWarning, glob
from .dataset import Dataset
from .errors import (
    CoordinateError,
    IntegrityError,
    NoSnapshotsError,
    OfflineError,
    OsteosarcError,
    SchemaError,
)
from .fixtures import load_panel, validate_recipe
from .models import File, Files, Region, Sample, Samples, Variant, Variants
from .parsing import (
    Table,
    parse_variants,  # noqa: F401
)
from .reads import (  # noqa: F401
    ReadFilter,
    ReadSubset,
    assembly_from_header,
    extract_reads,
    inspect_alignment,
)
from .recovery import RecoveryPolicy
from .shared import bundle_file, check_fixtures, fetch_bundle
from .sv_candidates import load_sv_candidates
from .timeline import Event, Timeline

__version__ = "0.12.0"

# The public API. A few other helpers are importable from here because the OpenVax
# libraries use them (digest, parse_variants, inspect_alignment, ...); everything
# else lives in its module: osteosarc.bundles, osteosarc.reads, osteosarc.records ...
__all__ = [
    # Browsing and downloading
    "Dataset", "Sample", "Samples", "File", "Files", "Variant", "Variants", "Event", "Timeline", "Table", "Region",
    "Cache",
    # Reads and test data
    "extract_reads", "ReadFilter", "ReadSubset", "RecoveryPolicy",
    "bundle_file", "fetch_bundle", "list_bundle", "export_bundle", "check_fixtures", "verify_bundle",
    "generate_bundle", "validate_recipe",
    "load_panel",
    "load_sv_candidates",
    # Corrections
    "CORRECTIONS", "Change", "Correction", "CurationWarning", "glob",
    # Errors
    "OsteosarcError", "CoordinateError", "IntegrityError", "NoSnapshotsError", "OfflineError", "SchemaError",
]
