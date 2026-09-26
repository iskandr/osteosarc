# Reads

Copy just the reads you need out of the remote BAMs into small local BAMs: test
data for a unit test in seconds, without downloading any BAM whole. You need
SAMtools on your PATH (1.21 is tested); reopening reads you already fetched
doesn't.

## Make test data

From a terminal, name a sample and a variant:

```sh
osteosarc reads T0_tumor --assay rna-seq --variant DYNC1H1-chr14-101980529 --padding 100 --to tests/data
```

A sample stands for each of its indexed BAMs, here its bulk RNA-seq ones.
Osteosarc streams only the reads within 100 bases of the variant, saves each result
in tests/data as a small indexed BAM named for its source and the variant, and
prints the paths. It skips a BAM whose genome build doesn't match the variant's,
saying so. Give a file's key instead of a sample to read one BAM, and repeat
`--variant` for several variants.

In Python:

```python
from osteosarc import Dataset

data = Dataset.open(offline=False)
source = data.file("rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam")
targets = data.variants(ids=["DYNC1H1-chr14-101980529"], status="ready")
subset = data.extract_reads(source, variants=targets, padding=100, to="tests/data")
print(subset.path)  # tests/data/BG003082.DYNC1H1-chr14-101980529.bam
```

Without a folder to put it in, the BAM stays in the cache. Asking again reuses it,
even offline. To find a BAM, see [samples and files](samples.md);
`data.samples["T0_tumor"]` lists a sample's BAMs.

## Regions instead of variants

```python
from osteosarc import Region

regions = [Region("chr14", 101980528, 101980530, "GRCh38")]
subset = data.extract_reads(source, regions)
with subset.open() as bam:
    print(bam.count("chr14", 101980528, 101980530))
```

In Python, regions are zero-based and half-open. On the command line they're
one-based and inclusive, like SAMtools, and need the genome build:

```sh
osteosarc reads rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam chr14:101980529-101980530 --assembly GRCh38
```

`Region.from_samtools("chr14:101980529-101980530", assembly="GRCh38")` converts one.
Overlapping regions are merged. Reads come back exactly as stored, duplicates and
tags included; a region with no reads gives an empty BAM.

## Filter reads

```python
from osteosarc import ReadFilter

filtered = data.extract_reads(source, regions, filters=ReadFilter(min_mapq=20, exclude_flags=0x500))
```

Flag 0x100 marks secondary alignments and 0x400 duplicates. A filter can also keep
only some cell barcodes, as in `ReadFilter(barcodes=("AAACCTG-1",))`, or some read
names. On the command line, use `--min-mapq` and `--exclude-flags`.

## Mates and split reads

`fetch_pairs=True` (or `--fetch-pairs`) also fetches each read's mate, even
outside your regions. To follow split reads as well, pass a recovery policy
(or `--recover-linked`):

```python
from osteosarc import RecoveryPolicy

linked = data.extract_reads(source, regions, recovery=RecoveryPolicy(max_rounds=2))
print(linked.receipt["status"], linked.receipt["records"])
```

Recovery fetches each mate and split-read piece from the same BAM and keeps a
partner only if its name, position and strand match the link. It never invents
one: a partner it can't find stays missing, and the receipt lists it. It stops at
limits on rounds, intervals, bases and records. Running out of rounds marks the
result truncated, and too many records is an error. With
`RecoveryPolicy(on_timeout="incomplete")`, a query that times out gives what was
found so far, marked incomplete, and the next call retries it.

## Genome builds

```python
info = data.inspect_alignment(source)
print(info.assembly)
```

Osteosarc reads the genome build from the chromosome lengths in the BAM header,
and refuses a request whose regions or variants come from another build; it never
converts coordinates. GRCh37 mitochondrial regions also need the chromosome's
length, as in `Region("MT", 12990, 13000, "GRCh37", reference_length=16569)`,
because the two GRCh37 references differ there. CRAM files need a local reference
FASTA, given with `reference=`.

## What's cached

Each result is cached by its request, with a receipt of its regions, filters,
checksums, header, tool versions and read count. Osteosarc checks a remote BAM's
HTTP headers before and after reading, to catch a file that changed mid-read.
`osteosarc downloads` lists every extract with its local path.

## Local BAMs and random samples

The same extraction works on any local indexed BAM, and a reproducible random
sample of its reads makes a quick fixture:

```python
from osteosarc import extract_reads
from osteosarc.reads import subset_templates

regional = extract_reads(str(subset.path), regions)
sample = subset_templates(regional, count=48, seed="fixture-v1")
print(sample.path, sample.receipt["records"])
```

A random sample ignores alleles and quality, so don't use it to estimate allele
fractions. For test data with chosen reads and controls, and the shared test data
the OpenVax libraries use, see [test data](test-data.md).
