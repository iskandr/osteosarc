"""Browse and download the public osteosarc.com dataset.

    from osteosarc import Dataset
    data = Dataset.sync()      # once: the website's metadata, about 57 MB
    data = Dataset.open()      # later: the newest snapshot, offline
    data                       # what's there, and how to get it

Importing performs no network I/O.
"""

from .bundles import (
    compact_header,
    export_bundle,
    generate_bundle,
    list_bundle,
    pack_bundle,
    safe_path,
    verify_bundle,
    verify_digest,
    verify_gzip_digests,
    verify_manifest_files,
)
from .cache import Cache, Receipt, digest
from .catalog import SNAPSHOT_SOURCES, TABLE_SOURCES, TIMELINE_SOURCES
from .curation import CORRECTIONS, Change, Correction, CurationWarning, glob
from .dataset import Dataset
from .discovery import list_bucket
from .errors import (
    CoordinateError,
    IntegrityError,
    NoSnapshotsError,
    OfflineError,
    OsteosarcError,
    SchemaError,
)
from .fixtures import (
    FixtureSelection,
    load_panel,
    select_fixture_records,
    select_fixtures,
    validate_recipe,
)
from .models import File, Files, Region, Sample, SampleClaim, Samples, Variant, Variants
from .parsing import Table, parse_file, parse_table, parse_variant_index, parse_variants
from .reads import (
    AlignmentInfo,
    ReadFilter,
    ReadSubset,
    assembly_from_header,
    extract_reads,
    inspect_alignment,
    normalize_assembly,
    resolve_regions,
    subset_templates,
)
from .records import RECORD_ENCODING, bam_record_digests, read_records, record_multiset
from .recovery import RecoveryPolicy, recover_reads
from .sv_candidates import load_sv_candidates
from .timeline import Event, Timeline

__version__ = "0.9.0"

__all__ = [
    "load_sv_candidates",
    "FixtureSelection", "load_panel", "select_fixture_records", "select_fixtures", "validate_recipe",
    "RECORD_ENCODING", "bam_record_digests", "read_records", "record_multiset",
    "RecoveryPolicy", "recover_reads",
    "compact_header", "export_bundle", "generate_bundle", "list_bundle", "pack_bundle",
    "safe_path", "verify_bundle", "verify_digest", "verify_gzip_digests", "verify_manifest_files",
    "AlignmentInfo", "File", "Files", "CORRECTIONS", "Cache", "Change", "CoordinateError", "Correction",
    "CurationWarning", "Dataset", "Event", "IntegrityError", "NoSnapshotsError", "TIMELINE_SOURCES", "Timeline", "glob",
    "OfflineError", "OsteosarcError", "ReadFilter", "ReadSubset", "Receipt", "Region",
    "SNAPSHOT_SOURCES", "Sample", "SampleClaim", "Samples", "SchemaError", "TABLE_SOURCES", "Table", "Variant",
    "Variants", "assembly_from_header", "digest", "extract_reads", "inspect_alignment", "list_bucket",
    "normalize_assembly", "parse_file", "parse_table", "parse_variant_index", "parse_variants",
    "resolve_regions", "subset_templates",
]
