"""Verified, shared, content-addressed downloads (curl; Linux and macOS).

Files live in the OpenVax shared cache, the same layout vaxrank and other
OpenVax tools use, so identical bytes are stored once:

    <root>/objects/sha256/<sha256><original suffixes>   shared content
    <root>/osteosarc/...                                 osteosarc's own records

<root> is OSTEOSARC_CACHE (an isolated cache), else OPENVAX_DATA_CACHE, else
the platform cache directory for "openvax" (as appdirs/datacache choose it).
"""

from __future__ import annotations

import fcntl
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
from contextlib import contextmanager
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import unquote, urlsplit

from .errors import IntegrityError, OfflineError


def digest(path, algorithm="sha256"):
    """Hash a file incrementally, without loading large alignments into RAM."""
    return digests(path, (algorithm,))[algorithm]


def digests(path, algorithms=("sha256",)):
    """Several hashes of a file in one read."""
    results = {name: hashlib.new(name) for name in algorithms}
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            for result in results.values():
                result.update(block)
    return {name: result.hexdigest() for name, result in results.items()}


def file_identity(path):
    """Size, modification time and inode: changes whenever a file is rewritten."""
    status = Path(path).stat()
    return [str(Path(path).resolve()), status.st_size, status.st_mtime_ns, status.st_ino, status.st_dev]


_UMASK = os.umask(0)
os.umask(_UMASK)


def share(path):
    """Give files the permissions the process umask allows (mkstemp/mkdtemp default to owner-only)."""
    os.chmod(path, (0o777 if Path(path).is_dir() else 0o666) & ~_UMASK)


def default_root():
    """The shared OpenVax cache root, resolved exactly as vaxrank/datacache resolve it."""
    if os.environ.get("OPENVAX_DATA_CACHE"):
        return Path(os.environ["OPENVAX_DATA_CACHE"])
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Caches" / "openvax"
    return Path(os.environ.get("XDG_CACHE_HOME") or Path.home() / ".cache") / "openvax"


def object_name(sha256, filename):
    """OpenVax content name: the SHA256 followed by the original file's suffixes."""
    return sha256 + "".join(Path(filename).suffixes)


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
        share(name)
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

    root defaults to OSTEOSARC_CACHE, else the OpenVax shared cache (see the
    module docstring). Content is shared under objects/sha256; receipts,
    snapshots and derived reads are kept under osteosarc/. Existing objects are
    verified on every fetch. Refresh changes the URL's latest receipt but never
    removes old snapshot objects.
    """

    def __init__(self, root=None, *, offline=False, timeout=600):
        self.root = Path(root or os.environ.get("OSTEOSARC_CACHE") or default_root()).expanduser().resolve()
        self.objects = self.root / "objects" / "sha256"
        self.workspace = self.root / "osteosarc"
        self.offline = offline
        self.timeout = timeout
        self._verified = {}

    def path(self, receipt, *, verify=True):
        """Resolve a receipt, detecting missing or modified bytes.

        Each object is hashed once per Cache instance; a rewritten file (new
        size, mtime or inode) is hashed again.
        """
        if not isinstance(receipt, Receipt):
            receipt = Receipt(**receipt)
        if (len(receipt.sha256) != 64 or any(c not in "0123456789abcdef" for c in receipt.sha256)
                or Path(receipt.filename).name != receipt.filename or receipt.filename in ("", ".", "..")):
            raise IntegrityError("Invalid cache receipt path")
        path = self.objects / object_name(receipt.sha256, receipt.filename)
        if not path.is_file():
            raise OfflineError(f"Cached object is missing: {receipt.url}")
        if verify:
            identity = json.dumps(file_identity(path))
            if self._verified.get(identity) != receipt.sha256:
                if path.stat().st_size != receipt.size or digest(path) != receipt.sha256:
                    raise IntegrityError(f"Cached object was modified: {path}")
                self._verified[identity] = receipt.sha256
        return path

    def file_digest(self, path):
        """SHA256 of a local file, remembered on disk by file identity.

        Local BAMs are large; an unchanged file (same size, mtime and inode) is
        not re-read. Any rewrite changes its identity and forces a new hash.
        """
        identity = file_identity(path)
        memo = self.workspace / "digests" / (stable_id(identity) + ".json")
        if memo.is_file():
            record = json.loads(memo.read_text())
            if record.get("identity") == identity:
                return record["sha256"]
        value = digest(path)
        if file_identity(path) == identity:
            write_json(memo, dict(identity=identity, sha256=value))
        return value

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
        pointer = self.workspace / "urls" / (stable_id(url) + ".json")
        with file_lock(self.workspace / "locks" / (stable_id(url) + ".lock")):
            if pointer.exists() and not refresh:
                receipt = Receipt(**json.loads(pointer.read_text()))
                if receipt.url != url:
                    raise IntegrityError("Cache URL mismatch")
                try:
                    path = self.path(receipt)
                except OfflineError:
                    if self.offline:
                        raise
                    path = None  # the object was pruned; download it again below
                if path is not None:
                    # path() verified sha256; an md5 already checked when the object was stored is not re-read.
                    known_md5 = md5 and receipt.md5 and md5.lower() == receipt.md5.lower()
                    self._validate(path, sha256=None if sha256 == receipt.sha256 else sha256,
                                   md5=None if known_md5 else md5, size=size, max_bytes=max_bytes)
                    return receipt
            if self.offline:
                raise OfflineError(f"Not cached: {url}")
            staging = self.workspace / "staging"
            staging.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory(dir=staging) as temporary:
                path = Path(temporary) / "download"
                headers = Path(temporary) / "headers"
                command = ["curl", "--fail", "--location", "--silent", "--show-error",
                           "--proto", "=http,https", "--proto-redir", "=http,https",
                           "--connect-timeout", "30", "--max-time", str(self.timeout),
                           "--retry", "2", "--retry-all-errors",
                           "--dump-header", str(headers), "--output", str(path)]
                if max_bytes is not None:
                    command += ["--max-filesize", str(max_bytes)]
                subprocess.run(command + [url], check=True, capture_output=True,
                               timeout=self.timeout * 3 + 120)
                hashes = digests(path, ("sha256", "md5") if md5 else ("sha256",))
                self._validate(path, sha256=sha256, md5=md5, size=size, max_bytes=max_bytes, hashes=hashes)
                # Keep the final HTTP response, rather than a redirect's headers.
                final = headers.read_text().strip().split("\n\n")[-1]
                fields = {k.lower(): v.strip() for line in final.splitlines()
                          if ":" in line for k, v in [line.split(":", 1)]}
                receipt = self._store(path, url, fields, md5, move=True, checksum=hashes["sha256"])
            write_json(pointer, receipt.to_dict())
            return receipt

    def import_file(self, path, url, *, sha256=None, md5=None, size=None):
        """Adopt an already downloaded file with its original source URL.

        Does not assert when the upstream object was downloaded. retrieved_at
        records the local import time. Useful for existing project caches.
        """
        with file_lock(self.workspace / "locks" / (stable_id(url) + ".lock")):
            path = Path(path)
            self._validate(path, sha256=sha256, md5=md5, size=size)
            receipt = self._store(path, url, {}, md5)
            write_json(self.workspace / "urls" / (stable_id(url) + ".json"), receipt.to_dict())
            return receipt

    @staticmethod
    def _validate(path, *, sha256=None, md5=None, size=None, max_bytes=None, hashes=None):
        if size is not None and path.stat().st_size != size:
            raise IntegrityError(f"Size mismatch: {path}")
        if max_bytes is not None and path.stat().st_size > max_bytes:
            raise IntegrityError(f"Download exceeds {max_bytes} bytes")
        wanted = {name: value for name, value in (("sha256", sha256), ("md5", md5)) if value is not None}
        if wanted:
            hashes = dict(hashes or {})
            missing = tuple(name for name in wanted if name not in hashes)
            hashes.update(digests(path, missing) if missing else {})
            for algorithm, expected in wanted.items():
                if hashes[algorithm] != expected.lower():
                    raise IntegrityError(f"{algorithm} mismatch: {path}")

    def _store(self, path, url, headers, md5, *, move=False, checksum=None):
        checksum = checksum or digest(path)
        filename = Path(unquote(urlsplit(url).path)).name or "index"
        if filename in (".", ".."):
            filename = "index"
        receipt = Receipt(url, checksum, path.stat().st_size, filename,
                          datetime.now(timezone.utc).isoformat(), headers.get("etag"),
                          headers.get("last-modified"), md5)
        output = self.objects / object_name(checksum, filename)
        output.parent.mkdir(parents=True, exist_ok=True)
        # The object may already be shared by another OpenVax tool; never rewrite valid bytes.
        if output.is_file() and output.stat().st_size == receipt.size and digest(output) == checksum:
            self._verified[json.dumps(file_identity(output))] = checksum
            return receipt
        if move:
            share(path)
            os.replace(path, output)
            self._verified[json.dumps(file_identity(output))] = checksum
            return receipt
        fd, temporary = tempfile.mkstemp(dir=output.parent, prefix=".object-")
        try:
            with os.fdopen(fd, "wb") as handle, path.open("rb") as source:
                shutil.copyfileobj(source, handle)
            self._validate(Path(temporary), sha256=checksum, size=receipt.size)
            share(temporary)
            os.replace(temporary, output)
        finally:
            Path(temporary).unlink(missing_ok=True)
        return receipt
