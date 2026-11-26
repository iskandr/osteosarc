"""Browse and download the public osteosarc.com dataset, and make test data from it.

    from osteosarc import Dataset
    data = Dataset.sync()      # once: the website's metadata, about 57 MB
    data = Dataset.open()      # later: the newest snapshot, offline
    data                       # what's there, and how to get it

Test data:

    data.extract_reads(key, variants=[...], to="tests/data")      # a quick BAM of the reads there
    data.make_bundle("tests/data/mine", variants=[...], files=[...])  # a bundle anyone can rebuild
    osteosarc.bundle_file("openvax-v2", member)                   # one member, as a local BAM

Importing performs no network I/O and resolves submodules lazily via PEP 562
so ``osteosarc --help`` and tab completion stay fast and import-light.
"""

from importlib import import_module

from ._version import __version__  # noqa: F401

_MODULE_EXPORTS = {
    ".bundles": (
        "export_bundle",
        "generate_bundle",
        "list_bundle",
        "subset_bundle",
        "verify_bundle",
        "verify_digest",
        "verify_gzip_digests",
        "verify_manifest_files",
    ),
    ".cache": (
        "Cache",
        "digest",
    ),
    ".curation": (
        "CORRECTIONS",
        "Change",
        "Correction",
        "CurationWarning",
        "glob",
    ),
    ".dataset": ("Dataset",),
    ".errors": (
        "CoordinateError",
        "IntegrityError",
        "NoSnapshotsError",
        "OfflineError",
        "OsteosarcError",
        "RecordLimitError",
        "SchemaError",
    ),
    ".fixtures": (
        "load_panel",
        "validate_recipe",
    ),
    ".models": (
        "File",
        "Files",
        "Region",
        "Sample",
        "Samples",
        "Variant",
        "Variants",
    ),
    ".parsing": (
        "Table",
        "parse_variants",
    ),
    ".read_receipts": ("read_receipt_files",),
    ".reads": (
        "ReadCount",
        "ReadFilter",
        "ReadSubset",
        "assembly_from_header",
        "count_reads",
        "extract_reads",
        "inspect_alignment",
    ),
    ".recovery": ("RecoveryPolicy",),
    ".shared": (
        "bundle_file",
        "check_fixtures",
        "fetch_bundle",
    ),
    ".sv_candidates": ("load_sv_candidates",),
    ".timeline": (
        "Event",
        "Timeline",
    ),
}

_EXPORT_TO_MODULE = {
    name: module_name
    for module_name, names in _MODULE_EXPORTS.items()
    for name in names
}

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
    "generate_bundle", "subset_bundle", "validate_recipe",
    "load_panel",
    "load_sv_candidates",
    # Corrections
    "CORRECTIONS", "Change", "Correction", "CurationWarning", "glob",
    # Errors
    "OsteosarcError", "CoordinateError", "IntegrityError", "NoSnapshotsError", "OfflineError", "RecordLimitError", "SchemaError",
]


def __getattr__(name):
    module_name = _EXPORT_TO_MODULE.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module = import_module(module_name, __name__)
    value = getattr(module, name)
    globals()[name] = value
    return value


def __dir__():
    return sorted(set(globals()) | set(__all__) | set(_EXPORT_TO_MODULE) | {"__version__"})
