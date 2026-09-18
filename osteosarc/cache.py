"""Verified, shared, content-addressed downloads (curl; Linux and macOS)."""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote, urlsplit

from .errors import IntegrityError, OfflineError


def digest(path, algorithm="sha256"):
    """Hash a file incrementally, without loading large alignments into RAM."""
    result = hashlib.new(algorithm)
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            result.update(block)
    return result.hexdigest()


def stable_id(value):
    """Return a deterministic SHA256 identity for a JSON-serializable value."""
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def write_json(path, value):
    """Atomically replace a JSON file on the same filesystem."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, name = tempfile.mkstemp(dir=path.parent, prefix=".json-")
    try:
        with os.fdopen(fd, "w") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(name, path)
    finally:
        Path(name).unlink(missing_ok=True)


@contextmanager
def file_lock(path):
    """Serialize writers; OS releases the advisory lock after interruption."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a") as handle:
        fcntl.flock(handle, fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(handle, fcntl.LOCK_UN)


@dataclass(frozen=True)
class Receipt:
    url: str
    sha256: str
    size: int
    filename: str
    retrieved_at: str
    etag: str | None = None
    last_modified: str | None = None
    md5: str | None = None

    def to_dict(self):
        return asdict(self)


class Cache:
    """A shared cache. Construction does no I/O; fetching is always explicit.

    Defaults to OSTEOSARC_CACHE, or $XDG_CACHE_HOME/osteosarc (otherwise
    ~/.cache/osteosarc). Existing objects are verified on every fetch. Refresh
    changes the URL's latest receipt but never removes old snapshot objects.
    """

    def __init__(self, root=None, *, offline=False, timeout=600):
        default = Path(os.environ.get("XDG_CACHE_HOME", Path.home() / ".cache")) / "osteosarc"
        self.root = Path(root or os.environ.get("OSTEOSARC_CACHE", default)).expanduser().resolve()
        self.offline = offline
        self.timeout = timeout

    def path(self, receipt, *, verify=True):
        """Resolve a receipt, detecting missing or modified bytes."""
        if not isinstance(receipt, Receipt):
            receipt = Receipt(**receipt)
        if (len(receipt.sha256) != 64 or any(c not in "0123456789abcdef" for c in receipt.sha256)
                or Path(receipt.filename).name != receipt.filename or receipt.filename in ("", ".", "..")):
            raise IntegrityError("Invalid cache receipt path")
        path = self.root / "objects" / receipt.sha256 / receipt.filename
        if not path.is_file():
            raise OfflineError(f"Cached object is missing: {receipt.url}")
        if verify and (path.stat().st_size != receipt.size or digest(path) != receipt.sha256):
            raise IntegrityError(f"Cached object was modified: {path}")
        return path

    def fetch(self, url, *, refresh=False, sha256=None, md5=None, size=None, max_bytes=None):
        """Download or verify cached bytes, returning their immutable receipt.

        Optional checksums are upstream claims checked against actual bytes;
        ETags are retained as opaque HTTP identities, never treated as MD5s.
        Interrupted downloads are discarded and retried on the next call.
        """
        if urlsplit(url).scheme not in ("https", "http"):
            raise ValueError("Downloads require an HTTP(S) URL; use Cache.import_file for local data")
        if refresh and self.offline:
            raise OfflineError("Cannot refresh in offline mode")
        pointer = self.root / "urls" / (stable_id(url) + ".json")
        with file_lock(self.root / "locks" / (stable_id(url) + ".lock")):
            if pointer.exists() and not refresh:
                receipt = Receipt(**json.loads(pointer.read_text()))
                if receipt.url != url:
                    raise IntegrityError("Cache URL mismatch")
                path = self.path(receipt)
                self._validate(path, sha256=sha256, md5=md5, size=size, max_bytes=max_bytes)
                return receipt
            if self.offline:
                raise OfflineError(f"Not cached: {url}")
            staging = self.root / "staging"
            staging.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory(dir=staging) as temporary:
                path = Path(temporary) / "download"
                headers = Path(temporary) / "headers"
                command = ["curl", "--fail", "--location", "--silent", "--show-error",
                           "--proto", "=http,https", "--proto-redir", "=http,https",
                           "--connect-timeout", "30", "--max-time", str(self.timeout),
                           "--retry", "2", "--dump-header", str(headers), "--output", str(path)]
                if max_bytes is not None:
                    command += ["--max-filesize", str(max_bytes)]
                subprocess.run(command + [url], check=True, capture_output=True,
                               timeout=self.timeout * 3 + 120)
                self._validate(path, sha256=sha256, md5=md5, size=size, max_bytes=max_bytes)
                # Keep the final HTTP response, rather than a redirect's headers.
                final = headers.read_text().strip().split("\n\n")[-1]
                fields = {k.lower(): v.strip() for line in final.splitlines()
                          if ":" in line for k, v in [line.split(":", 1)]}
                receipt = self._store(path, url, fields, md5, move=True)
            write_json(pointer, receipt.to_dict())
            return receipt

    def import_file(self, path, url, *, sha256=None, md5=None, size=None):
        """Adopt an already downloaded file with its original source URL.

        Does not assert when the upstream object was downloaded. retrieved_at
        records the local import time. Useful for existing project caches.
        """
        with file_lock(self.root / "locks" / (stable_id(url) + ".lock")):
            path = Path(path)
            self._validate(path, sha256=sha256, md5=md5, size=size)
            receipt = self._store(path, url, {}, md5)
            write_json(self.root / "urls" / (stable_id(url) + ".json"), receipt.to_dict())
            return receipt

    @staticmethod
    def _validate(path, *, sha256=None, md5=None, size=None, max_bytes=None):
        if size is not None and path.stat().st_size != size:
            raise IntegrityError(f"Size mismatch: {path}")
        if max_bytes is not None and path.stat().st_size > max_bytes:
            raise IntegrityError(f"Download exceeds {max_bytes} bytes")
        for algorithm, expected in (("sha256", sha256), ("md5", md5)):
            if expected is not None and digest(path, algorithm) != expected.lower():
                raise IntegrityError(f"{algorithm} mismatch: {path}")

    def _store(self, path, url, headers, md5, *, move=False):
        checksum = digest(path)
        filename = Path(unquote(urlsplit(url).path)).name or "index"
        if filename in (".", ".."):
            filename = "index"
        receipt = Receipt(url, checksum, path.stat().st_size, filename,
                          datetime.now(timezone.utc).isoformat(), headers.get("etag"),
                          headers.get("last-modified"), md5)
        output = self.root / "objects" / checksum / filename
        output.parent.mkdir(parents=True, exist_ok=True)
        if move:
            os.replace(path, output)
            return receipt
        fd, temporary = tempfile.mkstemp(dir=output.parent, prefix=".object-")
        try:
            with os.fdopen(fd, "wb") as handle, path.open("rb") as source:
                shutil.copyfileobj(source, handle)
            os.replace(temporary, output)
        finally:
            Path(temporary).unlink(missing_ok=True)
        return receipt
