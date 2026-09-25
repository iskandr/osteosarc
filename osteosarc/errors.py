"""Errors callers can distinguish without matching messages."""


class OsteosarcError(Exception):
    """Base exception for dataset operations."""


class IntegrityError(OsteosarcError, ValueError):
    """Cached bytes or source claims do not match the recorded receipt."""


class OfflineError(OsteosarcError):
    """An operation needs bytes that are not present in the local cache."""


class SchemaError(OsteosarcError, ValueError):
    """An upstream schema is missing required fields or is ambiguous."""


class CoordinateError(OsteosarcError, ValueError):
    """A region cannot be safely interpreted on an alignment's reference."""


class NoSnapshotsError(OsteosarcError, FileNotFoundError):
    """No snapshot has been saved in this cache yet."""

    def __init__(self, root=None):
        self.root = root
        where = f" in {root}" if root is not None else ""
        super().__init__(f"No saved snapshots{where}; run Dataset.sync() to download the website's metadata")
