"""Catalogues built from a snapshot once and kept in the cache.

A snapshot never changes, but building its file catalogue (some 400,000 bucket
objects) takes seconds, so each catalogue is saved the first time and loaded after.
A saved copy serves one snapshot, one version of osteosarc's code and one set of
corrections.

Copies are pickles, which can run code when loaded, so each user keeps their own
and loads only files they own that no one else can write. Copies unused for 30
days are removed, and each catalogue keeps its few most recent versions, so two
installed versions of osteosarc sharing a cache don't undo each other's.
"""

import functools
import gzip
import hashlib
import os
import pickle
import tempfile
import time
from pathlib import Path

from .cache import read_own

#: Copies of one catalogue kept (for different versions of osteosarc's code).
KEEP = 3
#: Copies unused for this long are removed.
UNUSED_DAYS = 30


def load_or_build(folder, stem, build):
    """The saved catalogue named stem in folder, or build() it, save it and return it."""
    code = code_key()
    if code is None:  # no source to fingerprint: build every time
        return build()
    path = Path(folder) / f"{stem}-{code}-u{os.getuid()}.pickle.gz"
    saved = read_own(path)
    if saved is not None:
        try:
            # Decompressed whole, not streamed: 40% faster, for a brief ~150 MB of bytes.
            value = pickle.loads(gzip.decompress(saved))
        except Exception:  # damaged: build it again
            pass
        else:
            try:
                os.utime(path)  # recently used
            except OSError:
                pass
            return value
    value = build()
    _save(path, value)
    return value


def _save(path, value):
    temporary = None
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        descriptor, temporary = tempfile.mkstemp(dir=path.parent, prefix=".catalog-")
        with os.fdopen(descriptor, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", compresslevel=1,
                                                                 mtime=0) as handle:
            pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)
        os.replace(temporary, path)
        temporary = None
        _prune(path)
    except Exception:
        pass  # a read-only or full cache still works; the catalogue is built each time
    finally:
        if temporary is not None:
            Path(temporary).unlink(missing_ok=True)


def _prune(path):
    """Remove this user's copies unused for UNUSED_DAYS, and all but the KEEP most
    recently used copies of this catalogue."""
    suffix = f"-u{os.getuid()}.pickle.gz"
    stem = path.name[:-len(suffix)].rsplit("-", 1)[0]
    cutoff = time.time() - UNUSED_DAYS * 86400
    mine = [(p.stat().st_mtime, p) for p in path.parent.glob(f"*{suffix}")]
    for mtime, old in mine:
        if mtime < cutoff and old != path:
            old.unlink(missing_ok=True)
    versions = sorted(((m, p) for m, p in mine if p.name.startswith(stem + "-") and p.exists()), reverse=True)
    for _, old in versions[KEEP:]:
        if old != path:
            old.unlink(missing_ok=True)


@functools.lru_cache(maxsize=1)
def code_key():
    """A fingerprint of what builds the catalogues: osteosarc's code and the HTML parser.

    None when the code can't be read (a zipped or compiled-only install)."""
    from importlib.metadata import PackageNotFoundError, version

    from . import __version__
    sources = sorted(Path(__file__).parent.glob("*.py"))
    if not sources:
        return None
    digest = hashlib.sha256(__version__.encode())
    try:
        digest.update(version("beautifulsoup4").encode())
    except PackageNotFoundError:
        pass
    for path in sources:
        digest.update(path.name.encode() + b"\0" + path.read_bytes())
    return digest.hexdigest()[:12]
