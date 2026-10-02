"""Lossless disk storage for receipts with repeated query-name filters."""

import json
import os
import re
import tempfile
from hashlib import sha256
from pathlib import Path

from .cache import Cache, digest, file_lock, share, write_json
from .errors import IntegrityError

SCHEMA = "osteosarc.read-receipt.v1"


def query_names_asset(names, workspace):
    """Store the exact canonical SAM query-name input once per cache."""
    if not names or any(not isinstance(n, str) or not re.fullmatch(r"[!-?A-~]{1,254}", n) for n in names):
        raise ValueError("Invalid SAM query-name list")
    data = ("\n".join(names) + "\n").encode("ascii")
    checksum = sha256(data).hexdigest()
    path = Path(workspace) / "query-names" / (checksum + ".txt")
    with file_lock(Path(workspace) / "locks" / (checksum + ".names.lock")):
        if path.exists():
            if digest(path) != checksum:
                raise IntegrityError("Cached query-name asset was modified")
        else:
            path.parent.mkdir(parents=True, exist_ok=True)
            fd, name = tempfile.mkstemp(dir=path.parent, prefix=".names-")
            try:
                with os.fdopen(fd, "wb") as handle:
                    handle.write(data)
                share(name)
                os.replace(name, path)
            finally:
                Path(name).unlink(missing_ok=True)
    return path


def write_read_receipt(path, receipt, workspace):
    """Store large filters by hash; request identities and the public API stay expanded."""
    assets = {}
    seen = {}

    def encode(value):
        if isinstance(value, dict):
            result = {}
            for key, item in value.items():
                if key == "query_names" and isinstance(item, (list, tuple)) and item:
                    names = tuple(item)
                    if names not in seen:
                        seen[names] = query_names_asset(names, workspace)
                    asset = seen[names]
                    assets[asset.stem] = len(item)
                    result[key] = {"sha256": asset.stem}
                else:
                    result[key] = encode(item)
            return result
        if isinstance(value, (list, tuple)):
            return [encode(item) for item in value]
        return value

    encoded = encode(receipt)
    write_json(path, dict(storage_schema=SCHEMA, query_name_assets=assets, receipt=encoded)
               if assets else receipt)


def _stored_receipt(path, workspace):
    value = json.loads(Path(path).read_text())
    if "storage_schema" not in value:
        return value, {}
    if value["storage_schema"] != SCHEMA:
        raise IntegrityError("Unrecognized read receipt storage schema")
    assets = {}
    for checksum, count in value["query_name_assets"].items():
        if not re.fullmatch("[0-9a-f]{64}", checksum) or type(count) is not int or count < 1:
            raise IntegrityError("Invalid query-name asset identity")
        asset = Path(workspace) / "query-names" / (checksum + ".txt")
        if not asset.is_file() or digest(asset) != checksum:
            raise IntegrityError("Missing or modified query-name asset: " + checksum)
        data = asset.read_bytes()
        try:
            names = data.decode("ascii").splitlines()
        except UnicodeDecodeError as error:
            raise IntegrityError("Invalid query-name asset encoding") from error
        if (len(names) != count or not data.endswith(b"\n")
                or any(not re.fullmatch(r"[!-?A-~]{1,254}", name) for name in names)):
            raise IntegrityError("Invalid query-name asset content")
        assets[checksum] = (asset, names)
    used = set()

    def validate(item):
        if isinstance(item, dict):
            for key, child in item.items():
                if key == "query_names" and isinstance(child, dict):
                    checksum = child.get("sha256")
                    if set(child) != {"sha256"} or checksum not in assets:
                        raise IntegrityError("Unresolved query-name asset reference")
                    used.add(checksum)
                else:
                    validate(child)
        elif isinstance(item, list):
            for child in item:
                validate(child)

    validate(value["receipt"])
    if used != set(assets):
        raise IntegrityError("Unreferenced query-name assets in receipt")
    return value["receipt"], assets


def read_read_receipt(path, workspace):
    """Read either storage format, verifying every referenced asset offline."""
    receipt, assets = _stored_receipt(path, workspace)
    if not assets:
        return receipt

    def expand(value):
        if isinstance(value, dict):
            return {key: list(assets[item["sha256"]][1])
                    if key == "query_names" and isinstance(item, dict)
                    else expand(item) for key, item in value.items()}
        if isinstance(value, list):
            return [expand(item) for item in value]
        return value

    return expand(receipt)


def read_receipt_files(path, workspace):
    """Return verified receipt/asset paths and hashes for downstream provenance pins.

    This does not verify the BAM/index; callers retain their original file pins.
    It avoids expanding a recovery receipt's repeated query-name lists.
    """
    _, assets = _stored_receipt(path, workspace)
    return {Path(path): digest(path), **{asset: checksum for checksum, (asset, _) in assets.items()}}


def compact_read_cache(cache):
    """Opt in to lossless metadata compaction of existing read derivatives.

    BAMs, indexes, semantic requests and derivative paths remain unchanged.
    Receipt file bytes change: do this before pinning their on-disk hashes in
    downstream audits. Old recorded command-input paths remain valid hard links.
    Each receipt is round-trip checked before atomic replacement; interruption
    leaves old or complete new metadata. No network or read acquisition occurs.
    """
    cache = cache if isinstance(cache, Cache) else Cache(cache, offline=True)
    derived = cache.workspace / "derived"
    paths = sorted(derived.glob("*/receipt.json")) + sorted(derived.glob("*.over-limit.json"))
    count = 0
    for path in paths:
        key = path.parent.name if path.name == "receipt.json" else path.name.removesuffix(".over-limit.json")
        with file_lock(cache.workspace / "locks" / (key + ".lock")):
            receipt = read_read_receipt(path, cache.workspace)
            fd, name = tempfile.mkstemp(dir=path.parent, prefix=".compact-")
            os.close(fd)
            temporary = Path(name)
            try:
                write_read_receipt(temporary, receipt, cache.workspace)
                if read_read_receipt(temporary, cache.workspace) != receipt:
                    raise IntegrityError("Receipt compaction changed its contents")
                os.replace(temporary, path)
                if path.name == "receipt.json":
                    old_names = path.parent / "query-names.txt"
                    if old_names.exists():
                        names = receipt.get("request", {}).get("filters", {}).get("query_names")
                        if not names:
                            raise IntegrityError("Query-name input has no recorded filter")
                        asset = query_names_asset(names, cache.workspace)
                        if digest(old_names) != asset.stem:
                            raise IntegrityError("Query-name command input differs from its receipt")
                        os.link(asset, temporary)
                        os.replace(temporary, old_names)
                count += 1
            finally:
                temporary.unlink(missing_ok=True)
    return count
