# Extract reads

Copy just the reads around some variants or regions out of a remote BAM into a
small local BAM, without downloading the whole file. The examples open your most
recent snapshot; see [Get started](index.md#get-started).

## Requirements

You need `samtools` on your PATH; version 1.21 is tested. Osteosarc checks that
your version supports what a request needs before downloading anything, and says
how to upgrade if it doesn't. Reopening reads you already fetched doesn't need
SAMtools.

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

You get a local indexed BAM. Osteosarc checks that the remote BAM uses the
same genome build as the variants, downloads its index, and reads only the parts
it needs. Pick a BAM with [`assets_for_sample`](explore.md#bulk-and-single-cell-data)
or `data.assets.select(...)`.

The command line does the same and reuses the same cached result:

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

Overlapping regions are merged. Reads come back exactly as stored, including
duplicates, flags, qualities and tags, with no filtering unless you ask for it.
A region with no reads gives an empty BAM.

## Filter reads

```python
from osteosarc import ReadFilter

filtered = data.extract_reads(
    source, regions,
    filters=ReadFilter(min_mapq=20, exclude_flags=0x100 | 0x400),
)
```

`0x100` drops secondary alignments and `0x400` drops reads marked as duplicates.
`require_flags` keeps only reads with the given flags. To pick single cells, use
`ReadFilter(barcodes=("cell-barcode",), barcode_tag="CB")`; to pick reads by
name, `ReadFilter(query_names=("read-a", "read-b"))`. Each combination of
filters is cached separately.

## Recover mates and split alignments

`fetch_pairs=True` also fetches each read's mate, even outside your regions. It
doesn't fetch the other pieces of split reads:

```python
paired = data.extract_reads(source, regions, fetch_pairs=True)
```

To follow both mates and split-read (`SA` tag) links, pass a `RecoveryPolicy`
instead of `fetch_pairs`:

```python
from osteosarc import RecoveryPolicy

linked = data.extract_reads(source, regions, recovery=RecoveryPolicy(max_rounds=2))
print(linked.receipt["status"], linked.receipt["records"])
```

`mates` and `supplementary` choose which links to follow (both on by default).
The limits default to `max_rounds=4`, `max_intervals=128`, `max_bases=1_000_000`
and `max_records=100_000`. On the command line, use `reads --recover-linked`.

How it works:

- All reads in your regions are kept. Each link is then fetched from the same
  BAM, and a partner is kept only if its read group, read name, which mate it is,
  position and strand all match the link. A mate must be a primary alignment; a
  split-read piece must also match the link's CIGAR and MAPQ, and its NM tag
  when it has one.
- Nothing is invented: a partner that isn't found stays missing, and reads keep
  their original clipping, flags and tags.
- Hitting the round or interval limit gives a result marked `truncated`.
  Exceeding `max_records` is an error, so you never get a silently partial BAM.
- The receipt lists every region visited, every link followed, and every
  partner that was missing, conflicting or out of bounds. `complete_template`
  is always false: the BAM can't prove that no other pieces exist.
- At busy loci, partner queries first filter by the names of your reads, so
  unrelated reads don't use up `max_records`.

A remote partner query can time out after your own reads were already fetched.
By default the whole request then fails. With
`RecoveryPolicy(on_timeout="incomplete")` you get your reads and any partners
already found instead, marked `status="incomplete"`. The receipt's
`failed_queries` names the query that timed out, and its links are listed as
unresolved. Calling again retries that query and reuses everything already
fetched. A BAM that changed between queries is still an error.

`recover_reads(source, regions, policy=...)` does the same for local files and
URLs.

## Check the alignment's reference

```python
info = data.inspect_alignment(source)
print(info.assembly)
print(info.header["SQ"][:2])
```

Osteosarc works out the genome build from the chromosome lengths in the BAM
header, so it doesn't need the index for this. If it can't tell, `assembly` is
`None`. A request fails, rather than guessing, when the builds don't match, a
region falls off the end of a chromosome, or the BAM has no index. Coordinates
are never converted between builds.

GRCh37 mitochondrial regions need `Region(..., reference_length=...)`, because
hg19 and hs37d5 disagree on its length. For CRAM files, pass
`reference="local-reference.fa"` with a `.fai` index next to it.

## Reuse the result offline

```python
offline = Dataset.open()
assert offline.extract_reads(source, regions).path == subset.path
```

Each result is cached by its request, and its receipt records the regions,
filters, checksums, BAM header, tool versions and read count. For a remote BAM,
Osteosarc checks the file's HTTP headers before and after reading, to catch a file
that changes mid-download.

## Use a local BAM or sample a small fixture

```python
from osteosarc import extract_reads, subset_templates

local = str(subset.path)  # Any local indexed BAM works here
regional = extract_reads(local, regions)
fixture = subset_templates(regional, count=48, seed="fixture-v1")
print(fixture.path, fixture.receipt["records"])
```

This picks 48 reads at random (reproducibly, from `seed`), ignoring alleles and
quality, and keeps every record in the input that shares their read group and name,
such as a mate that's also in the input. It never fetches missing mates; start from
a `fetch_pairs=True` result if you want them. The receipt marks the result as
sampled, so don't use it to estimate VAF. For test data with specific reads and controls, use
[fixture recipes](fixtures.md).

## Generate a panel for every sample

This writes one indexed BAM, with mates, per RNA-seq BAM of every sample,
covering the chosen loci. It queries every RNA-seq BAM, so expect it to take an hour
or more; pass a longer `timeout=` (in seconds, default 600) for slow remote queries.

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

To use your own loci, replace the `variants` and `regions` lines with:

```python
from osteosarc import Region

regions = [Region.from_samtools(locus, assembly="GRCh38") for locus in [
    "chr2:165658600-165659700",
    "chr12:5439800-5440900",
]]
```

Remove `assay="rna-seq"` to include DNA and single-cell BAMs. GRCh37 BAMs need
GRCh37 coordinates. Some vendor BAMs have no index, so the loop skips them rather
than download them whole.
