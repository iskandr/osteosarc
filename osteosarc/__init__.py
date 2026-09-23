"""Reproducible access to osteosarc.com. Importing performs no network I/O."""

from .cache import Cache, Receipt, digest
from .catalog import SNAPSHOT_SOURCES, TABLE_SOURCES, TIMELINE_SOURCES
from .curation import CORRECTIONS, Change, Correction, CurationWarning, glob
from .dataset import Dataset
from .discovery import list_bucket
from .errors import CoordinateError, IntegrityError, OfflineError, OsteosarcError, SchemaError
from .fixtures import (
    FixtureSelection,
    load_panel,
    select_fixture_records,
    select_fixtures,
    validate_recipe,
)
from .models import Asset, Assets, Region, SampleClaim, Variant, Variants
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
from .timeline import Event, Timeline

__version__ = "0.2.0"

__all__ = [
    "FixtureSelection", "load_panel", "select_fixture_records", "select_fixtures", "validate_recipe",
    "RECORD_ENCODING", "bam_record_digests", "read_records", "record_multiset",
    "AlignmentInfo", "Asset", "Assets", "CORRECTIONS", "Cache", "Change", "CoordinateError", "Correction",
    "CurationWarning", "Dataset", "Event", "IntegrityError", "TIMELINE_SOURCES", "Timeline", "glob",
    "OfflineError", "OsteosarcError", "ReadFilter", "ReadSubset", "Receipt", "Region",
    "SNAPSHOT_SOURCES", "SampleClaim", "SchemaError", "TABLE_SOURCES", "Table", "Variant",
    "Variants", "assembly_from_header", "digest", "extract_reads", "inspect_alignment", "list_bucket",
    "normalize_assembly", "parse_file", "parse_table", "parse_variant_index", "parse_variants",
    "resolve_regions", "subset_templates",
]
