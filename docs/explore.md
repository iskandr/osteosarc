# Files, samples, and local caching

Start with an existing snapshot:

```python
from osteosarc import Dataset

data = Dataset.open("baseline")
rna = data.assets.select(kind="alignment", assay="rna-seq", timepoint="T2")
for asset in rna:
    print(asset.key, asset.size, asset.conflicts)
```

## Select sequencing products

```python
ont = data.assets.select(kind="alignment", platform="ont")
dna = data.assets.select(kind="alignment", assay="wgs")
vcfs = data.assets.select(kind="variants", format="vcf")
fastqs = data.assets.select(kind="reads")
pvac = data.assets.select(prefix="neoantigen_prediction/pvactools/", format="tsv")
```

Assays include `rna-seq`, `scrna-seq`, `cite-seq`, `wgs`, and `wes`.
Known platform names are `ont`, `pacbio`, and `illumina`. Unknown metadata is
retained as unknown. Every listed object remains accessible, even when its
format or sample identity cannot be classified.

Use an exact key to choose a processing product; substring searches may find
several versions of the same library:

```python
source = data.asset(
    "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.bam"
)
print(source.index_urls)
print(source.claims)
```

## Examine sample claims

```python
print(data.timepoints.rows[:5])
print(data.samples.columns)
conflicted = data.assets.where(lambda a: bool(a.conflicts))
for asset in conflicted[:3]:
    print(asset.key, asset.values("timepoint"), asset.resolved("timepoint"))
```

An asset is a file or processing product, not an independent biological sample.
`samples` contains source-attributed claims and the associated asset IDs. Dates
retain their published precision; a month and a day within it do not conflict.
Conflicting claims do not match metadata filters unless `include_conflicts=True`;
path inferences require `include_inferred=True`. Review those claims before
combining samples. For biological specimens (T0_tumor, blood_2025-06-26, ...)
linked to their files, use `data.specimens` (see [timelines](timeline.md)).

Label vocabulary lives in `osteosarc/curation.py`. For example, "Normal" and
"Blood" are the same tissue (`blood`), because every normal here is a blood
normal. A `CITE` viewer label means `cite-seq`. With the default corrections,
the stale viewer labels for BG009368 and SARC0277 are fixed, and no
sample-metadata conflicts remain. `Dataset.open(..., corrections=False)` shows
the three genuine source disagreements. Unrecognized source labels are listed
in `data.unrecognized`.

## Download exactly what you need

```python
data = Dataset.open("baseline", offline=False)
path = data.download("snv_top")
table = data.table("snv_top")
print(path.name, table.columns)
```

`download` returns the shared cache object's `pathlib.Path`, named by its
SHA256 plus the original suffixes (for example `…9c2e.genes.results`), so
format-specific readers still recognize it. The original name is kept in the
receipt (`data.cache.fetch(url).filename`). `table` accepts any CSV/TSV asset or exact
key, including RSEM `.genes.results` and `.isoforms.results` files.

```python
reports = data.assets.select(contains=".genes.results", format="tsv")
if reports:
    expression = data.table(reports[0])
    print(expression.rows[:2])
```

Tables preserve original strings: `"0"`, `"NA"`, and `""` are distinct.
`table.to_dataframe()` needs pandas and performs no automatic coercion.
`data.parse(asset)` also supports JSON and FASTA. For large tables/FASTA,
use `download` with a streaming reader; built-in parsers materialize their input.

## Reproducibility and refresh

Metadata snapshots pin URL, SHA256, size, and acquisition receipts. A full
object's first download is bound to that snapshot. Reuse verifies its checksum,
and later URL refreshes cannot replace bound bytes. Missing offline data raises
`OfflineError`; modified cached bytes raise `IntegrityError`.

```python
new_data = Dataset.sync("follow-up", refresh=True)
print(new_data.receipts()["bucket"].sha256)
```

Choose a new name to refresh metadata. Data objects are still acquired lazily:
a snapshot cannot recover historical bytes of an object never downloaded.
Interrupted whole-file downloads restart; byte-range resume is not implemented.

Adopt a previously downloaded object using its original URL and receipt:

```python
from osteosarc import Cache, digest

# Stand-ins for an entry in your existing manifest: the table downloaded above.
old_path, original_url = path, data.asset("snv_top").url
cache = Cache()
receipt = cache.import_file(old_path, original_url, sha256=digest(old_path))
print(receipt.url, receipt.sha256)
```

## Explicit live discovery

```python
from osteosarc import Cache, list_bucket

listing = list_bucket(Cache(), "ONT/", refresh=True)
print(len(listing["files"]))
```

This lists every page of the requested public S3 prefix and keeps page receipts.
It returns a separate inventory and does not mutate the dated Dataset snapshot.
