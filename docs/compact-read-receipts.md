# Compact read acquisition receipts

Bounded recovery repeats a seed's query-name list across partner requests. Store
canonical newline-delimited lists once under their SHA-256 in the cache. Both
samtools `-N` inputs and on-disk receipt references use those assets. Keep expanded
in-memory receipts and semantic request/cache identities unchanged. Read legacy
plain JSON receipts offline and reject missing, modified or malformed assets.

Expose the complete receipt-file dependency set so downstream audits can pin
provenance without embedding it. Compaction of existing receipts must round-trip
exactly, preserve every BAM/index hash, and replace old duplicate command inputs
with identical shared files while keeping their recorded paths valid.

Validate offline reuse, multiple partner requests sharing one asset, corruption,
legacy receipts and interrupted writes. Issue: #107.

To compact an existing cache offline before pinning receipt file hashes:

```python
from osteosarc import Cache
from osteosarc.read_receipts import compact_read_cache

compact_read_cache(Cache('/path/to/cache', offline=True))
```

`ReadSubset.receipt` still contains expanded filter lists. To avoid copying those
lists into another report, pin the verified dependency set instead:

```python
from osteosarc import read_receipt_files

pins = read_receipt_files(subset.receipt_path, cache.workspace)
```

The mapping includes the receipt and all referenced name assets. Pin and verify
the BAM and index separately, as before. Legacy receipts remain readable. Old
versions of osteosarc cannot read compact receipts; keep version 0.15.2 or newer
on all clients of a shared cache after compaction.
