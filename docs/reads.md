# Extract reads

Fetch the reads around selected variants or regions from a remote BAM into a
cached, indexed local BAM, without downloading the whole file. These examples
open your most recent snapshot; see [Get started](index.md#get-started).

## Requirements

Install `osteosarc` and put `samtools` on PATH. SAMtools 1.21 has been
tested. Header inspection requires `view --no-PG`; extraction also requires
`-M` and `-X`. Paired-mate extraction needs `--fetch-pairs`, barcode filtering
needs `-D`, and query-name filtering and bounded partner recovery need `-N`.
These capabilities are checked before acquisition, with an error explaining how
to upgrade. Cached reads can be reopened without SAMtools.

## Fetch reads around a variant

```python
from osteosarc import Dataset

data = Dataset.open(offline=False)
source = data.asset(
    "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam"
)
targets = data.variants(ids=["DYNC1H1-chr14-101980529"], status="ready")
subset = data.extract_reads(source, variants=targets, padding=100)
print(subset.path, subset.receipt["records"])
```

The result is a local indexed BAM. The extractor checks the source assembly,
downloads the index, and retrieves the requested regions without a full-BAM scan.
Choose a source with [`assets_for_sample`](explore.md#bulk-and-single-cell-data)
or `data.assets.select(...)`.

The command line takes the same variant IDs and returns the same cached BAM:

```sh
osteosarc reads rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam --variant DYNC1H1-chr14-101980529 --padding 100
```

## Specify coordinates

```python
from osteosarc import Region

regions = [Region("chr14", 101980528, 101980530, "GRCh38")]
subset = data.extract_reads(source, regions)
with subset.open() as bam:
    print(bam.count("chr14", 101980528, 101980530))
```

Python regions are **zero-based, half-open**. CLI regions are **one-based,
inclusive**, like SAMtools. To convert a SAMtools region:

```python
region = Region.from_samtools("chr14:101980529-101980530", assembly="GRCh38")
```

```sh
osteosarc reads rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam chr14:101980529-101980530 --assembly GRCh38
```

Overlapping intervals are queried as a union. Original duplicate records,
flags, qualities, and tags are retained. No quality or allele filter is applied
by default. An empty result is valid; an empty region list is an error.

## Filter reads

```python
from osteosarc import ReadFilter

filtered = data.extract_reads(
    source, regions,
    filters=ReadFilter(min_mapq=20, exclude_flags=0x100 | 0x400),
)
```

Here `0x100` excludes secondary alignments and `0x400` excludes duplicate-marked
reads. `require_flags` keeps only reads with the given flags. To select cells,
use `ReadFilter(barcodes=("cell-barcode",), barcode_tag="CB")`. To select
templates by name, use `ReadFilter(query_names=("read-a", "read-b"))`; an empty
tuple leaves names unrestricted. Filters are part of the cache identity, so
each combination is cached separately.

## Recover mates and split alignments

`fetch_pairs=True` also retrieves paired mates outside the intervals. It does
not recover supplementary (split) alignments:

```python
paired = data.extract_reads(source, regions, fetch_pairs=True)
```

To follow mate and supplementary-alignment (`SA`) links, pass a bounded
`RecoveryPolicy` instead. It replaces `fetch_pairs`; combining them is an error:

```python
from osteosarc import RecoveryPolicy

linked = data.extract_reads(source, regions, recovery=RecoveryPolicy(max_rounds=2))
print(linked.receipt["status"], linked.receipt["records"])
```

`mates` and `supplementary` choose which links to follow; both default to true.
The limits default to `max_rounds=4`, `max_intervals=128`, `max_bases=1_000_000`
and `max_records=100_000`. The CLI equivalent is `reads --recover-linked`.

Every available seed record is kept. Partner windows keep only records that
match the seed's source, read group and segment. An `SA` partner must also match
position, strand, CIGAR and MAPQ, plus NM when present. Records keep their hard
clipping, flags, tags and any absent SEQ/QUAL. No query synthesizes a missing
partner or a reverse-complement read. Duplicate counts use the maximum occurrence
count across indexed queries, never a sum of repeated retrievals.

Seed interval and base limits are checked against the resolved union before
extraction. The record budget covers every acquired candidate, including records
that match no partner pointer. A later link into an already visited window reuses
the records fetched there without another query.

The receipt records visited intervals, requested leads, missing or conflicting
partners, malformed `SA` tags, missing sequence, ambiguous placement, repeated or
cyclic leads, and the limits reached. `complete_template` is always false; even a
fully resolved `SA` graph may omit unreported alignments. Reaching the interval or
round limit produces an explicit `truncated` receipt. Exceeding the record budget
fails instead of publishing a partial result. Missing indexes fail without falling
back to a full scan, and a source that changes between steps fails. An
interrupted run reuses only verified intermediate extractions.

At dense loci, partner windows first select the seed query names with
[SAMtools `view -N`](https://www.htslib.org/doc/samtools-view.html) and only then
apply the record cap, so unrelated reads can't exhaust it. The cap still counts
matching names that later read-group or `SA` checks reject. A shared query name
alone never establishes template identity. `recover_reads(source, regions,
policy=...)` is the same operation for local files and URLs.

## Check the alignment's reference

```python
info = data.inspect_alignment(source)
print(info.assembly)
print(info.header["SQ"][:2])
```

Inspection uses reference lengths in the header and does not require an index.
An unresolved assembly is `None`. Extraction rejects assembly conflicts,
ambiguous contigs, missing indexes, out-of-bounds regions, and headers too sparse
to establish an assembly. Coordinates are never lifted over automatically.

GRCh37 mitochondrial queries need `Region(..., reference_length=...)` because
hg19 and hs37d5 differ there. For CRAM, pass `reference="local-reference.fa"`
with an existing `.fai` index.

## Reuse the result offline

```python
offline = Dataset.open()
assert offline.extract_reads(source, regions).path == subset.path
```

Results are cached by the request and source identity. Receipts record intervals,
filters, checksums, source headers, tool versions, and record counts. New remote
extractions check HTTP identity before and after the query; this is not a
checksum of the entire remote BAM.

## Use a local BAM or sample a small fixture

```python
from osteosarc import extract_reads, subset_templates

local = str(subset.path)  # Any local indexed BAM works here
regional = extract_reads(local, regions)
fixture = subset_templates(regional, count=48, seed="fixture-v1")
print(fixture.path, fixture.receipt["records"])
```

Sampling chooses templates by `(read group, query name)` without using alleles
or quality, and keeps their available regional records. The receipt marks the
result as sampled. Use unsampled reads to estimate VAF. For versioned fixtures
with explicit witnesses and controls, use [fixture recipes](fixtures.md).

## Generate a panel for every sample

This writes one indexed BAM per registry-linked RNA-seq BAM product for every
sample on GRCh38. Each BAM contains the union of the nominated loci; products
from the same specimen stay separate. It queries every RNA-seq BAM, so expect it
to take an hour or more; each extraction accepts a longer `timeout=` in seconds
(default 600) for slow remote queries.

<!-- docs-check: skip (batch job over every RNA-seq BAM; takes over an hour) -->
```python
import json
import shutil
from pathlib import Path

from osteosarc import Dataset

data = Dataset.sync("panel-v1")  # Reuses this pinned snapshot on later runs
variants = [v for v in data.variants(status="ready")
            if v.gene in {"NTF3", "MAP2"} and v.assembly == "GRCh38"]
regions = [v.region(padding=500) for v in variants]
assert regions, "No eligible variants in this snapshot"

for sample in data.specimens:
    sources = data.assets_for_sample(
        sample["sample_id"], kind="alignment", format="bam", assay="rna-seq",
    )
    for source in sources:
        if not source.index_urls:
            print("Skipping alignment published without an index:", source.key)
            continue
        if data.inspect_alignment(source).assembly != "GRCh38":
            print("Skipping incompatible or unresolved assembly:", source.key)
            continue
        subset = data.extract_reads(source, regions, fetch_pairs=True)
        output = Path("panel") / sample["sample_id"] / source.id
        output.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(subset.path, output / "reads.bam")
        shutil.copyfile(subset.index_path, output / "reads.bam.bai")
        (output / "receipt.json").write_text(json.dumps(subset.receipt, indent=2))
        print(sample["sample_id"], source.key, output)
```

For explicit loci, replace the `variants`/`regions` lines with:

```python
from osteosarc import Region

regions = [Region.from_samtools(locus, assembly="GRCh38") for locus in [
    "chr2:165658600-165659700",
    "chr12:5439800-5440900",
]]
```

Remove `assay="rna-seq"` to include the other BAM assays. GRCh37 products need
their own verified GRCh37 coordinates. Some vendor BAMs are published without an
index; extraction never falls back to downloading them whole, so the loop skips
them. Samples without matching BAM products have no output. Extraction is indexed and cached; paired mates can lie outside
the loci. No template sampling or allele filtering is applied.
