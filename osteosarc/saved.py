"""Catalogues built from a snapshot once and kept in the cache.

A snapshot never changes, but building its file catalogue (some 400,000 bucket
objects) takes seconds, so each catalogue is saved the first time and loaded after.
A saved copy serves one snapshot, one version of osteosarc's code and one set of
corrections.

Copies are pickles, loaded so that they can hold only osteosarc's catalogue classes
and plain values: a pickle could otherwise run code. Each user keeps their own, in a
folder only they can write in, and loads only files they own that no one else can
write, each checked to be the catalogue it's named as. Where no folder can be
private (a drive without Unix permissions), nothing is saved. Each kind of catalogue
keeps its three most recently used copies (for other snapshots or versions of
osteosarc), and copies unused for 30 days are removed.
"""

import functools
import gzip
import hashlib
import io
import os
import pickle
import sys
from pathlib import Path

from .cache import private_folder, prune_own, read_own, write_own

#: The classes a saved catalogue may hold; a copy naming any other fails to load.
CLASSES = {("osteosarc.models", name) for name in ("File", "Files", "SampleClaim", "Variant", "Variants")}
#: Copies of each kind of catalogue kept.
KEEP = 3
#: Copies unused for this long are removed.
UNUSED_DAYS = 30


def load_or_build(folder, name, key, build):
    """The saved catalogue called name for key (what it's built from), or build()
    it, save it and return it. Nothing is saved when key is None."""
    code = code_key()
    if code is None or key is None:
        return build()
    stem = f"{name}-{key}-{code}"
    filename = f"{stem}.pickle.gz"
    with private_folder(folder) as private:
        if private is None:
            return build()
        saved = read_own(private, filename)
        if saved is not None:
            try:
                # Decompressed whole, not streamed: 40% faster, for a brief ~150 MB of bytes.
                saved_stem, value = _Unpickler(io.BytesIO(gzip.decompress(saved))).load()
            except Exception:  # damaged: build it again
                saved_stem = None
            if saved_stem == stem:  # not another catalogue under this one's name
                try:
                    os.utime(filename, dir_fd=private)  # recently used
                except OSError:
                    pass
                return value
        value = build()
        try:
            write_own(private, filename, lambda handle: _dump(handle, (stem, value)))
            _prune(private, name, filename)
        except Exception:
            pass  # a read-only or full cache still works; the catalogue is built each time
        return value


class _Unpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if (module, name) not in CLASSES:
            raise pickle.UnpicklingError(f"{module}.{name} isn't a catalogue class")
        return super().find_class(module, name)


def _dump(handle, value):
    with gzip.GzipFile(fileobj=handle, mode="wb", compresslevel=1, mtime=0) as compressed:
        pickle.dump(value, compressed, protocol=pickle.HIGHEST_PROTOCOL)


def _prune(folder, name, current):
    """Remove copies unused for UNUSED_DAYS, all but the KEEP most recently used of
    this kind of catalogue, and leftovers of interrupted saves."""
    kept = prune_own(folder, ".pickle.gz", UNUSED_DAYS)
    same_kind = sorted(((mtime, n) for n, mtime in kept.items() if n.startswith(name + "-")), reverse=True)
    for _, old in same_kind[KEEP:]:
        if old != current:
            try:
                os.unlink(old, dir_fd=folder)
            except OSError:
                pass  # already gone (another process pruning), or not a file


def _fingerprint():
    """What builds the catalogues: osteosarc's code, and Python. (The only data they
    use is the corrections', which are in each catalogue's name.) None when the
    code can't be read (a zipped or compiled-only install)."""
    sources = sorted(p for p in Path(__file__).parent.glob("*.py") if not p.name.startswith("."))
    digest = hashlib.sha256(sys.version.encode())
    try:
        for path in sources:
            digest.update(path.name.encode() + b"\0" + path.read_bytes())
    except OSError:  # e.g. an editor's lock file, or a file being replaced
        return None
    return digest.hexdigest() if sources else None


# Taken when osteosarc is imported (dataset imports this module), so it's the code
# that's running, even if the files change later.
_CODE = _fingerprint()


@functools.lru_cache(maxsize=1)
def code_key():
    """_CODE, and the version of the HTML parser the file catalogue uses."""
    from importlib.metadata import PackageNotFoundError, version
    if _CODE is None:
        return None
    try:
        parser = version("beautifulsoup4")
    except PackageNotFoundError:
        parser = ""
    return hashlib.sha256(f"{_CODE} {parser}".encode()).hexdigest()[:12]
