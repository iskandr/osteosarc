# Inspect and extract reads

Install the `reads` extra and `samtools`. The same API accepts Dataset assets
or, through standalone functions, local BAM/CRAM paths and HTTP(S) URLs.

## Inspect assembly before choosing regions

```python
from osteosarc import Dataset

data = Dataset.open("baseline", offline=False)
source = data.asset(
    "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.bam"
)
info = data.inspect_alignment(source)
print(info.assembly)  # GRCh38, GRCh37, or None when unresolved
print(info.header["SQ"][:2])
```

Header inspection does not need an index. It caches the original SAM header
with a checksum and source evidence. Assembly requires at least two
distinguishing canonical contigs and no conflicting lengths. A viewer-wide
hg38 label is insufficient. Liftover remains an explicit consumer operation.

## Query a region union

```python
from osteosarc import Region

regions = [Region("chr14", 101980528, 101980530, "GRCh38")]
subset = data.extract_reads(source, regions)
print(subset.path, subset.receipt["records"])
with subset.open() as bam:
    for read in bam.fetch("chr14", 101980528, 101980530):
        print(read.query_name, read.cigarstring)
```

Python regions are **zero-based, half-open**. For SAMtools coordinates:

```python
region = Region.from_samtools("chr14:101980529-101980530", assembly="GRCh38")
```

The extractor requires an index and queries the union of overlapping intervals.
It retains original record multiplicity, flags, qualities, barcodes, UMIs, and
tags. Defaults impose no quality, allele, duplicate, or template-count filter.
An empty region list is an error; a valid query with no reads is a valid result.

## Explicit filters and mate recovery

```python
from osteosarc import ReadFilter

filtered = data.extract_reads(
    source, regions,
    filters=ReadFilter(min_mapq=20, exclude_flags=0x100 | 0x400),
)
paired = data.extract_reads(source, regions, fetch_pairs=True)
```

`fetch_pairs=True` retrieves paired mates outside the intervals, not every
supplementary alignment. Barcode selection uses
`ReadFilter(barcodes=("cell-barcode",), barcode_tag="CB")`.

Assembly conflicts, ambiguous contig aliases, missing indexes, and intervals
beyond reference bounds raise errors. GRCh37 mitochondrial regions additionally
require `reference_length` to distinguish hg19/hs37d5. CRAM extraction requires
`reference="local-reference.fa"` with an existing `.fai`.

## Cache and provenance

Each result contains an indexed BAM and a receipt with the full request,
resolved intervals, input header, source/index identities, checksums, tool
versions, and record count. Repeating the same request reuses verified bytes
offline. New remote queries check HTTP identity before and after extraction
and reject changes since the cached header. This records HTTP evidence; it
does not claim a whole-BAM checksum for an object only accessed regionally.

```python
offline = Dataset.open("baseline")
assert offline.extract_reads(source, regions).path == subset.path
```

## Local files and small regression fixtures

```python
from osteosarc import extract_reads, inspect_alignment, subset_templates

# Any local indexed BAM works; here, the subset extracted above.
local = str(subset.path)
info = inspect_alignment(local)
regional = extract_reads(local, regions)
fixture = subset_templates(regional, count=48, seed="fixture-v1")
print(info.assembly, regional.receipt["records"], fixture.receipt["records"])
```

The sampler selects templates by `(read group, query name)` independently of
alleles/quality, keeping all available regional records for selected templates.
Its receipt marks the result as a sampled fixture unsuitable for VAF estimation.
Consumer stress-fixture selection policies remain in those repositories.
