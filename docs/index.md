# Get started

Osteosarc is a Python library and command-line tool for the public
[osteosarc.com](https://osteosarc.com/data/) dataset. Use it to find files,
select variants and vaccine peptides, or fetch sequencing reads for an analysis.

## Install

```sh
python -m pip install osteosarc
```

You need Python 3.10+ and Linux or macOS. Read extraction also requires
`samtools` on PATH; see [requirements](reads.md#requirements).

## Save the metadata

```python
from osteosarc import Dataset

data = Dataset.sync("baseline")
```

This downloads about 57 MB of metadata into a local cache. It leaves sequencing
files remote. `baseline` is your name for this snapshot.

## Find files and variants

```python
print(data.describe_samples())
rna = data.assets_for_sample("T0_tumor", kind="alignment", assay="rna-seq")
for asset in rna:
    print(asset.key, asset.size)

for variant in data.variants("vaccine", status="ready"):
    print(variant.gene, variant.allele)
```

A `ready` variant has one literal allele with consistent coordinates. Check
[source corrections](curation.md) before comparing against published results;
they are applied by default.

## Fetch reads for selected variants

```python
targets = data.variants(ids=["DYNC1H1-chr14-101980529"], status="ready")
source = rna["rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam"]
reads = data.extract_reads(source, variants=targets, padding=100)
print(reads.path)
```

This produces an indexed BAM containing the overlapping reads. Downloads go
through datacache; repeating the same request reuses the cached result. Pass
the BAM and variants to [Isovar](consumers.md#isovar) for read-support analysis.

## Pick up where you left off

```python
data = Dataset.open("baseline")  # Opens offline
counts = data.table("vafs").select(gene="SMC5")
print(counts.rows[:2])
```

Use `offline=False` when you want to download additional files or extract reads.
See [snapshots and cache](design.md) to change the cache directory or refresh data.

## Choose an example

| I want to… | Guide |
| --- | --- |
| Find RNA, DNA, or nanopore files and read a table | [Find files](explore.md) |
| Get alleles, read counts, or vaccine peptide sequences | [Select variants](variants.md) |
| Fetch a region, filter reads, or make a small BAM fixture | [Extract reads](reads.md) |
| Browse treatments, specimens, and measurements | [Browse the timeline](timeline.md) |
| Pass data to Varcode, Isovar, Topiary, or Vaxrank | [Use other libraries](consumers.md) |
| Compare a corrected allele with the original reads | [Worked example: MAP2](tour.md) |

## Prefer the terminal?

```sh
osteosarc sync baseline
osteosarc samples baseline
osteosarc assets baseline --assay rna-seq --timepoint T2
osteosarc variants baseline --gene SMC5
osteosarc explore baseline
```

Type `help` in the explorer for commands and `quit` to exit.
