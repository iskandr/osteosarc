"""Catalogues built from a snapshot once and kept in the cache.

A snapshot never changes, but building its file catalogue (some 400,000 bucket
objects) takes seconds, so each catalogue is saved the first time and loaded after.
A saved copy serves one snapshot, one version of osteosarc's code and one set of
corrections.

Copies are pickles, which can run code when loaded, so each user keeps their own,
in a folder only they can write in, and loads only files they own that no one else
can write, each checked to be the catalogue it's named as. Where no folder can be
private (a drive without Unix permissions), nothing is saved. Each kind of catalogue
keeps its three most recently used copies (for other snapshots or versions of
osteosarc), and copies unused for 30 days are removed.
"""

import functools
import gzip
import hashlib
import os
import pickle
import time
from pathlib import Path

from .cache import private_folder, read_own, write_own

#: Copies of each kind of catalogue kept.
KEEP = 3
#: Copies unused for this long are removed.
UNUSED_DAYS = 30


def load_or_build(folder, name, key, build):
    """The saved catalogue called name for key (what it's built from), or build()
    it, save it and return it. Nothing is saved when key is None."""
    code = code_key()
    private = private_folder(folder) if code is not None and key is not None else None
    if private is None:
        return build()
    stem = f"{name}-{key}-{code}"
    path = private / f"{stem}.pickle.gz"
    saved = read_own(path)
    if saved is not None:
        try:
            # Decompressed whole, not streamed: 40% faster, for a brief ~150 MB of bytes.
            saved_stem, value = pickle.loads(gzip.decompress(saved))
        except Exception:  # damaged: build it again
            saved_stem = None
        if saved_stem == stem:  # not another catalogue under this one's name
            try:
                os.utime(path)  # recently used
            except OSError:
                pass
            return value
    value = build()
    try:
        write_own(path, gzip.compress(pickle.dumps((stem, value), protocol=pickle.HIGHEST_PROTOCOL),
                                      compresslevel=1, mtime=0))
        _prune(private, name, path)
    except Exception:
        pass  # a read-only or full cache still works; the catalogue is built each time
    return value


def _prune(folder, name, current):
    """Remove copies unused for UNUSED_DAYS, all but the KEEP most recently used of
    this kind of catalogue, and files left by interrupted saves."""
    now = time.time()
    copies = sorted(((p.stat().st_mtime, p) for p in folder.glob("*.pickle.gz")), reverse=True)
    extra = [p for _, p in [c for c in copies if c[1].name.startswith(name + "-")][KEEP:]]
    for mtime, path in copies:
        if path != current and (mtime < now - UNUSED_DAYS * 86400 or path in extra):
            path.unlink(missing_ok=True)
    for leftover in folder.glob(".own-*"):
        if leftover.stat().st_mtime < now - 86400:
            leftover.unlink(missing_ok=True)


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
