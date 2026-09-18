"""Reproducible access to osteosarc.com. Importing performs no network I/O."""

from .cache import Cache, Receipt, digest
from .catalog import SNAPSHOT_SOURCES, TABLE_SOURCES
from .dataset import Dataset
from .discovery import list_bucket
from .errors import CoordinateError, IntegrityError, OfflineError, OsteosarcError, SchemaError
from .models import Asset, Assets, Region, SampleClaim, Variant, Variants
from .parsing import Table, parse_file, parse_table, parse_variant_index, parse_variants
from .reads import (
    ReadFilter,
    ReadSubset,
    assembly_from_header,
    extract_reads,
    normalize_assembly,
    resolve_regions,
    subset_templates,
)

__version__ = "0.1.0"

__all__ = [
    "Asset", "Assets", "Cache", "CoordinateError", "Dataset", "IntegrityError",
    "OfflineError", "OsteosarcError", "ReadFilter", "ReadSubset", "Receipt", "Region",
    "SNAPSHOT_SOURCES", "SampleClaim", "SchemaError", "TABLE_SOURCES", "Table", "Variant",
    "Variants", "assembly_from_header", "digest", "extract_reads", "list_bucket",
    "normalize_assembly", "parse_file", "parse_table", "parse_variant_index", "parse_variants",
    "resolve_regions", "subset_templates",
]
