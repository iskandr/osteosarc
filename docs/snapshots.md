# Snapshots and cache

## Save and reopen a snapshot

```python
from osteosarc import Dataset

data = Dataset.sync()
print(data.name, data.downloaded, data.id)
data = Dataset.open()
```

`sync()` downloads the website's current metadata into a snapshot named by the
UTC date, such as `2026-09-24`, with a checksum for every file. Running it again
the same day reopens that snapshot. `open()` reopens the most recent snapshot,
checks its files, and stays offline.

Use `Dataset.open(offline=False)` to allow new downloads. Offline, a request for a
file that isn't cached raises `OfflineError`; a cached file that changed raises
`IntegrityError`.

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

A name or ID always means one exact snapshot, and a missing name is an error. On
the command line, a value that looks like a date always means a download date.
`osteosarc snapshots` prints the same list and names the default. A snapshot's
`downloaded` time is when its files were downloaded, which can be earlier than
when the snapshot was made if it reused cached files.

## Choose a cache directory

```python
from osteosarc import Cache

cache = Cache(".cache/my-analysis")
```

Pass `cache=cache` to `Dataset.sync` or `Dataset.open`. Otherwise the first of
these applies:

| Setting | Location |
| --- | --- |
| `OSTEOSARC_CACHE` | An isolated cache for this package |
| `OPENVAX_DATA_CACHE` | A shared OpenVax cache |
| macOS default | `~/Library/Caches/openvax` |
| Linux default | `$XDG_CACHE_HOME/openvax`, or `~/.cache/openvax` |

Downloads are stored by checksum, at `objects/sha256/<sha256><original suffixes>`,
so other OpenVax tools can reuse them. Snapshots, receipts and extracted reads go
under `osteosarc/`.

## Find what's downloaded

```python
print(data.local_path("vafs"))  # a site table every snapshot holds
print(data.downloads())
```

Cached files are named by checksum, so ask for them by key instead.
`data.local_path(file)` gives a file's downloaded copy, or `None`, without using
the network. `data.downloads()` lists every downloaded file and every read extract
with its local path; `osteosarc downloads` prints the same, and
`osteosarc files --downloaded` lists the downloaded files among the others.

To put a file in a folder under its own name, use `data.download(file, to="data")`
or `osteosarc download KEY --to data`. A BAM's or VCF's index comes too. Each is a
read-only hard link to the cached copy when the folder is on the same disk, so it
takes no more space and can't be changed by accident; otherwise it's a copy.

## Refresh metadata

```python
new_data = Dataset.sync(refresh=True)
print(new_data.name, new_data.receipts()["bucket"].sha256)
```

`refresh=True` downloads the metadata again even if you have a snapshot from
today; a second snapshot on one day is named `2026-09-24.2`. Snapshots are never
overwritten, and a file downloaded through a snapshot stays tied to it. A snapshot
can't recover an old version of a file it never downloaded.

To name a snapshot for a project, use `Dataset.sync("paper-2026")` or
`osteosarc sync paper-2026`. It's created once and reopened after that. A named
snapshot reuses files already in the cache unless you pass `refresh=True`.

`osteosarc sync NAME --source-revision <commit>` takes the tables that come from
the website's GitLab repository at one commit, given as the full 40-character SHA.
Files served by osteosarc.com itself have no versions.

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

The receipt records when you imported the file, not when you first downloaded it;
keep your own notes if that date matters.

## Download behavior

Files are downloaded with datacache, which checks them, retries on network
errors, and only stores complete files. Osteosarc adds locks so two processes
don't fetch the same file at once, checks MD5s where the source publishes them,
and compares a file's HTTP headers before and after downloading it. An interrupted
download starts over from the beginning.

Reading regions of a BAM uses SAMtools, and each result gets its own receipt (see
[Extract reads](reads.md)). Importing the package never downloads anything.
