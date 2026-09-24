# Osteosarc

Osteosarc is a Python library and command-line tool for the public
[osteosarc.com](https://osteosarc.com/data/) dataset. The dataset shares one
patient's osteosarcoma sequencing, variant calls, cancer vaccine designs and
clinical history. Osteosarc pins the metadata you use and corrects documented
errors in it. It fetches only the sequencing reads you ask for.

| You can… | Guide |
| --- | --- |
| Search nearly 400,000 files by sample, timepoint and assay without downloading them | [Find samples and files](explore.md) |
| Select catalogue variants, read counts, vaccine peptides and ELISPOT results | [Select variants](variants.md) |
| Fetch reads around variants or regions from remote BAMs into a cached local BAM | [Extract reads](reads.md) |
| Browse treatments, specimens, MRD and lab results on a timeline | [Browse the timeline](timeline.md) |
| Do all of this from a terminal or an interactive explorer | [Command line](cli.md) |
| Pass data to Varcode, Isovar, Topiary or Vaxrank | [Use other libraries](consumers.md) |
| Build small, verifiable test BAMs from pinned recipes | [Read fixtures](fixtures.md) |
| Explore 637 candidate structural variants with their evidence | [SV interest catalogue](sv-interest.md) |

!!! note "Corrections are on by default"
    Osteosarc applies 32 documented, evidence-backed [corrections](curation.md)
    to the published data. For example, it replaces the MAP2 vaccine target's
    allele with the complex event that Tempus and CeGaT report and that the
    [tumor reads support](tour.md). Each load checks every correction against
    the snapshot's source records. Pass `corrections=False` to
    `Dataset.open`, or use `osteosarc --no-corrections`, to see the published values.

New to the dataset? Read [Key concepts](concepts.md) for sample IDs, variant
statuses and coordinate conventions.

## Get started

### Install

```sh
python -m pip install osteosarc
```

You need Python 3.9+ on Linux or macOS. Read extraction also requires
`samtools` on PATH; see [requirements](reads.md#requirements).

### 1. Save the metadata

```python
from osteosarc import Dataset

data = Dataset.sync()
```

This downloads about 57 MB of the website's metadata into a local cache, as a
snapshot named by today's date (UTC), and leaves the sequencing files remote.
Running it again the same day reopens that snapshot.

### 2. Choose a sample and assay

```python
print(data.describe_samples())
```

`T0_tumor` names the primary tumor specimen collected at **T0** (2022-12-16).
`T0_blood` is blood from the same timepoint. These are sample IDs; `T0` alone
is a collection timepoint. Sample type is the `tissue` field (`tumor`, `blood`,
`organoid`). The sequencing assay is independent of both.

`T0_tumor` has **bulk RNA-seq, whole-exome DNA (WES), and whole-genome DNA
(WGS)**. Choose `rna-seq`, `wes`, or `wgs` respectively:

```python
rna = data.assets_for_sample("T0_tumor", kind="alignment", assay="rna-seq")
for asset in rna:
    print(asset.key)
```

`rna-seq` selects bulk RNA; `scrna-seq` selects single-cell RNA. `T1_tumor`
has both. Single-cell data can also be selected by platform, such as
`platform="ont"` or `platform="pacbio"`. See [Find samples and files](explore.md)
for more filters and the full assay vocabulary.

### 3. Select variants

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

### 4. Fetch reads for those variants

```python
source = rna["rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam"]
reads = data.extract_reads(source, variants=targets, padding=100)
print(reads.path)
```

This checks the alignment's assembly and produces an indexed BAM of overlapping
reads, without downloading the whole file. [Isovar](consumers.md#isovar) can
classify which reads support the reference or alternate allele. Repeating the
same request reuses the cached result.

### 5. Pick up where you left off

```python
data = Dataset.open()  # Your most recent snapshot, offline
counts = data.table("vafs").select(gene="SMC5")
print(counts.rows[:2])
```

The website changes over time; your snapshot doesn't. Run `Dataset.sync(refresh=True)`
on a later day to save a new one, and `Dataset.snapshots()` to list them. Pass a
name, download date or month, such as `Dataset.open("2026-09")`, to reopen an older one.
Use `offline=False` when you want to download additional files or extract new reads.
See [snapshots and cache](design.md) for details.

## Prefer the terminal?

```sh
osteosarc sync
osteosarc samples
osteosarc assets --sample T0_tumor --kind alignment --assay rna-seq
osteosarc variants --gene SMC5
osteosarc explore
```

Every command uses your most recent snapshot; `osteosarc snapshots` lists them, and
`--snapshot 2026-09` picks another. Type `help` in the explorer for commands and
`quit` to exit. The [command-line guide](cli.md) lists every command.
