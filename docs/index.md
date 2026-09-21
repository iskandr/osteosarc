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

## Choose a sample and assay

```python
print(data.describe_samples())
```

`T0_tumor` names the primary tumor specimen collected at **T0** (2022-12-16).
`T0_blood` is blood from the same timepoint. These are sample IDs; `T0` alone
is a collection timepoint. Sample type is the `tissue` field (`tumor`, `blood`,
`organoid`); the sequencing assay is independent of both.

`T0_tumor` has **bulk RNA-seq, whole-exome DNA (WES), and whole-genome DNA
(WGS)**. Choose `rna-seq`, `wes`, or `wgs` respectively:

```python
rna = data.assets_for_sample("T0_tumor", kind="alignment", assay="rna-seq")
for asset in rna:
    print(asset.key)
```

`rna-seq` selects bulk RNA; `scrna-seq` selects single-cell RNA. `T1_tumor`
has both. Single-cell data can also be selected by platform, such as
`platform="ont"` or `platform="pacbio"`. See [sample IDs and sequencing types](explore.md)
for examples and the source-label vocabulary.

## Select variants

```python
targets = data.variants(gene="DYNC1H1", status="ready")
for variant in targets:
    print(variant.id, variant.allele)
```

Variant `status` describes genomic-allele usability. `ready` means one
consistent chromosome, position, REF and ALT, with literal DNA bases. It does
not establish read support, somatic status or a protein effect. Omit the
filter to include unresolved entries; see [all statuses](variants.md#variant-status).
The allele tuple is `(chromosome, one-based position, REF, ALT)`.
[Source corrections](curation.md) are applied by default.

## Fetch reads for those variants

```python
source = rna["rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam"]
reads = data.extract_reads(source, variants=targets, padding=100)
print(reads.path)
```

This checks the alignment's assembly and produces an indexed BAM of overlapping
reads. [Isovar](consumers.md#isovar) can classify which reads support the reference
or alternate allele. Downloads use datacache; repeating the same extraction
request reuses the cached result.

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
