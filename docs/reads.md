# Extract reads

## Requirements

Install `osteosarc` and put `samtools` on PATH. SAMtools 1.21 has been
tested. Header inspection requires `view --no-PG`; extraction also requires
`-M` and `-X`. Paired-mate extraction needs `--fetch-pairs`, and barcode
filtering needs `-D`; query-name filtering and bounded partner recovery need `-N`. These capabilities are checked before acquisition, with
an error explaining how to upgrade. Cached reads can be reopened without
SAMtools.

These examples use the `baseline` snapshot from [Get started](index.md).

## Fetch reads around a variant

```python
from osteosarc import Dataset

data = Dataset.open("baseline", offline=False)
source = data.asset(
    "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.bam"
)
targets = data.variants(ids=["DYNC1H1-chr14-101980529"], status="ready")
subset = data.extract_reads(source, variants=targets, padding=100)
print(subset.path, subset.receipt["records"])
```

The result is a local indexed BAM. The extractor checks the source assembly,
downloads the index, and retrieves the requested regions without a full-BAM scan.

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

Overlapping intervals are queried as a union. Original duplicate records,
flags, qualities, and tags are retained. No quality or allele filter is applied
by default. An empty result is valid; an empty region list is an error.

## Generate a panel for every sample

This writes one indexed BAM per registry-linked RNA-seq BAM product for every
sample on GRCh38. Each BAM contains the union of the nominated loci; products
from the same specimen stay separate. Install SAMtools first.

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
their own verified GRCh37 coordinates. Samples without matching BAM products
have no output. Extraction is indexed and cached; paired mates can lie outside
the loci. No template sampling or allele filtering is applied. For compact,
versioned regression fixtures with explicit witnesses, use [fixture recipes](fixtures.md).

## Filter reads or recover mates

```python
from osteosarc import ReadFilter

filtered = data.extract_reads(
    source, regions,
    filters=ReadFilter(min_mapq=20, exclude_flags=0x100 | 0x400),
)
paired = data.extract_reads(source, regions, fetch_pairs=True)
```

Here `0x100` excludes secondary alignments and `0x400` excludes duplicate-marked
reads. To select cells, use
`ReadFilter(barcodes=("cell-barcode",), barcode_tag="CB")`.

`fetch_pairs=True` retrieves paired mates outside the intervals. It does not
recover every supplementary alignment for a template.

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
offline = Dataset.open("baseline")
assert offline.extract_reads(source, regions).path == subset.path
```

Results are cached by the request and source identity. Receipts record intervals,
filters, checksums, source headers, tool versions, and record counts. New remote
extractions check HTTP identity before and after the query; this is not a
checksum of the entire remote BAM.

## Use a local BAM or make a test fixture

```python
from osteosarc import extract_reads, subset_templates

local = str(subset.path)  # Any local indexed BAM works here
regional = extract_reads(local, regions)
fixture = subset_templates(regional, count=48, seed="fixture-v1")
print(fixture.path, fixture.receipt["records"])
```

Sampling chooses templates by `(read group, query name)` without using alleles
or quality, and keeps their available regional records. The receipt marks the
result as sampled. Use unsampled reads to estimate VAF.


## Bound partner acquisition at dense loci

`recover_reads` keeps every seed/context record allowed by the requested filters.
For subsequent indexed partner windows it first selects the seed query names
with [SAMtools `view -N`](https://www.htslib.org/doc/samtools-view.html), then
applies the record cap. Unrelated names at dense loci cannot exhaust that cap.
All seed names are selected in every round so a later lead can reuse records
from an already visited window. Source, read group, segment, strand and exact
mate/SA placement checks still determine which acquired records are retained.
A shared query name alone never establishes template identity.

Explicit acquisition can also use `ReadFilter(query_names=("read-a", "read-b"))`.
An empty tuple leaves names unrestricted. Other filters still apply, and the
canonical name set is part of the cache receipt. Caps still count acquired
records with matching names, including records rejected by later RG/SA checks.
