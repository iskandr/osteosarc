"""Errors callers can distinguish without matching messages."""


class OsteosarcError(Exception):
    """Base exception for dataset operations."""


class IntegrityError(OsteosarcError):
    """Cached bytes or source claims do not match the recorded receipt."""


class OfflineError(OsteosarcError):
    """An operation needs bytes that are not present in the local cache."""


class SchemaError(OsteosarcError, ValueError):
    """An upstream schema is missing required fields or is ambiguous."""


class CoordinateError(OsteosarcError, ValueError):
    """A region cannot be safely interpreted on an alignment's reference."""
