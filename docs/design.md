# Snapshots and cache

## Save and reopen a snapshot

```python
from osteosarc import Dataset

data = Dataset.sync()
print(data.name, data.downloaded, data.id)
data = Dataset.open()
```

A snapshot records metadata URLs and SHA-256 receipts. `sync()` downloads the
website's current metadata into a snapshot named by the UTC date, such as
`2026-09-24`; running it again the same day reopens that snapshot. `open()` reopens
the most recently downloaded snapshot, verifies its saved files, and defaults to
offline operation.

To permit additional downloads, use `Dataset.open(offline=False)`.
An uncached request while offline raises `OfflineError`; changed cached bytes
raise `IntegrityError`.

## List and choose snapshots

```python
for row in Dataset.snapshots():
    print(row["name"], row["downloaded"], row["id"][:12])
same = Dataset.open(data.name)
```

`Dataset.snapshots()` lists snapshots newest first by download time.

| Python | CLI | Opens |
| --- | --- | --- |
| `Dataset.open()` | (nothing) | The most recently downloaded snapshot |
| `Dataset.open(date="2026-09")` | `--snapshot 2026-09` | The newest snapshot downloaded in that UTC year, month or day |
| `Dataset.open("2026-09-24.2")` | `--snapshot 2026-09-24.2` | The snapshot with that exact name |
| `Dataset.open("4b07fd4b")` | `--snapshot 4b07fd4b` | The snapshot whose ID starts with six or more given characters |

A name or ID pins one snapshot exactly; a missing name is an error, never read as
a date. On the command line, a value shaped like a date always means a download
date. `osteosarc snapshots` prints the same list and names the default. A snapshot's
`downloaded` time is its latest source download, which is earlier than its creation
if it reused cached sources.

## Choose a cache directory

```python
from osteosarc import Cache

cache = Cache(".cache/my-analysis")
```

Pass `cache=cache` to `Dataset.sync` or `Dataset.open`. Otherwise the cache
location is chosen in this order:

| Setting | Location |
| --- | --- |
| `OSTEOSARC_CACHE` | An isolated cache for this package |
| `OPENVAX_DATA_CACHE` | A shared OpenVax cache |
| macOS default | `~/Library/Caches/openvax` |
| Linux default | `$XDG_CACHE_HOME/openvax`, or `~/.cache/openvax` |

Downloaded objects live at `objects/sha256/<sha256><original suffixes>`, so
other OpenVax tools can reuse them. Osteosarc's snapshots, receipts, and
extracted reads live under `osteosarc/` within that root.

## Refresh metadata

```python
new_data = Dataset.sync(refresh=True)
print(new_data.name, new_data.receipts()["bucket"].sha256)
```

`refresh=True` downloads the website's metadata again, even if a snapshot from
today exists; a second snapshot on one day is named `2026-09-24.2`. Existing
snapshots are never overwritten. A file's first full download is bound to its
snapshot, so refreshing the same URL elsewhere cannot replace those bytes. A
snapshot cannot recover historical contents of a file that was never downloaded.

To label a snapshot for a project, pass a name, as in `Dataset.sync("paper-2026")`
or `osteosarc sync paper-2026`. It is created once and reopened on later runs. A
named snapshot reuses source bytes already in the cache unless `refresh=True`.

The CLI's `sync --source-revision <commit>` pins GitLab source resources to a
full commit. Site-served files have no equivalent versioning.

## Import a file you already downloaded

```python
from osteosarc import SNAPSHOT_SOURCES, digest

# This example uses an existing snapshot's table; substitute your file and URL.
old_path = data.source_path("vafs")
receipt = Cache().import_file(
    old_path, SNAPSHOT_SOURCES["vafs"], sha256=digest(old_path)
)
print(receipt.sha256, receipt.size)
```

An import records the import time. Retain your old manifest if the original
acquisition date matters.

## Download behavior

Datacache downloads and verifies files, retries transient failures, and
publishes complete files atomically. Osteosarc adds per-URL locks, immutable
snapshot receipts, and MD5 checks where the source supplies them. Existing
0.1.0 snapshots and cached reads remain usable without conversion.

HTTP headers are checked before and after downloads when the server supports
HEAD. Regional BAM/CRAM access uses SAMtools; its indexes use datacache and
its extracted BAMs have separate checksum receipts.

Interrupted whole-file downloads restart from the
beginning; byte-range resume is not supported. Cached objects keep their
original suffixes so format-specific readers can recognize them.

[Read extraction](reads.md) has a separate cache keyed by its complete request.
Importing the Python package performs no downloads.
