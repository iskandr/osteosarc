# Snapshots and cache

The website changes, so osteosarc keeps dated copies of its metadata: snapshots.
A snapshot never changes, works offline, and records a checksum for every file it
holds. Sequencing files stay in the bucket until you ask for them.

## Save and reopen

```python
from osteosarc import Dataset

data = Dataset.sync()
print(data.name, data.downloaded, data.id[:12])
data = Dataset.open()
```

`Dataset.sync()` downloads the website's metadata, about 57 MB, into a snapshot
named by the UTC date, such as 2026-09-25. Running it again the same day reopens
that snapshot. `Dataset.open()` reopens the newest snapshot and stays offline; pass
offline=False to let it download files or read regions of BAMs.

## Choose a snapshot

```python
for row in Dataset.snapshots():
    print(row["name"], row["downloaded"])
```

| Python | Command line | Opens |
| --- | --- | --- |
| `Dataset.open()` | (the default) | The newest snapshot |
| `Dataset.open(date="2026-09")` | `--snapshot 2026-09` | The newest snapshot downloaded in that year, month or day |
| `Dataset.open("2026-09-25.2")` | `--snapshot 2026-09-25.2` | The snapshot with that name, or an ID starting with those characters |

To get today's metadata when you already have a snapshot from today, use
`Dataset.sync(refresh=True)` or `osteosarc sync --refresh`; the new one is named
2026-09-25.2. `Dataset.sync("paper-2026")` names a snapshot for a project; it's
made once and reopened after that. Snapshots are never overwritten.

## Where things are kept

Osteosarc shares a cache with the other OpenVax tools, in the first of these:

| Setting | Location |
| --- | --- |
| OSTEOSARC_CACHE | A cache for osteosarc only |
| OPENVAX_DATA_CACHE | A shared OpenVax cache |
| macOS | ~/Library/Caches/openvax |
| Linux | $XDG_CACHE_HOME/openvax, or ~/.cache/openvax |

`Cache("some/folder")`, passed as cache= to sync and open, or `--cache` on the
command line, uses another folder. Downloads are stored by checksum, so a file
downloaded once is never stored twice, even by another tool.

## Find what's downloaded

```python
print(data.local_path("vafs"))
print(data.downloads())
```

Files in the cache are named by checksum, so ask for them by key.
`data.local_path(file)` gives a file's copy, or None, without using the network,
and `data.downloads()` lists everything downloaded and every read extract, with
local paths. `osteosarc downloads` prints the same.

To put a file in a folder under its own name, use `data.download(file, to="data")`
or `osteosarc download KEY --to data`. It's a read-only hard link to the cached
copy where it can be, so it takes no extra space and can't be changed by accident.

## Import a file you already have

```python
from osteosarc import Cache

# This reuses a table from the snapshot; substitute your own file and its URL.
table = data.file("vafs")
receipt = Cache().import_file(data.local_path(table), table.url)
print(receipt.sha256, receipt.size)
```

`Cache().import_file(path, url)` adds a file you already downloaded, so osteosarc
never downloads it again. The receipt records when you imported it, not when you
first downloaded it.

## How downloads work

Downloads are checked against their published size and checksum, retried on
network errors, and stored only when complete; an interrupted download starts
over. A file downloaded through a snapshot stays tied to it: if the file changes on
the server later, that snapshot keeps the version it saw, and a new snapshot gets
the new one. Importing the package never downloads anything.
