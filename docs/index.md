# Osteosarc

Osteosarc is a Python library and command-line tool for the public
[osteosarc.com](https://osteosarc.com/data/) dataset: one patient's osteosarcoma
sequencing, variant calls, cancer vaccines and clinical history. It keeps a fixed
copy of the site's metadata so your results don't change under you, fixes known
errors in it, and downloads only the reads you ask for.

| You can… | Guide |
| --- | --- |
| Search nearly 400,000 files by sample, timepoint and assay without downloading them | [Find samples and files](explore.md) |
| Get variants, read counts, vaccine peptides and ELISPOT results | [Select variants](variants.md) |
| Copy the reads around variants out of remote BAMs, without downloading them whole | [Extract reads](reads.md) |
| Browse treatments, specimens, MRD and lab results on a timeline | [Browse the timeline](timeline.md) |
| Do all of this from a terminal or an interactive explorer | [Command line](cli.md) |
| Pass data to Varcode, Isovar, Topiary or Vaxrank | [Use other libraries](consumers.md) |
| Build small, reproducible test BAMs | [Read fixtures](fixtures.md) |
| Look through 637 candidate structural variants | [SV catalogue](sv-interest.md) |

!!! note "Corrections are on by default"
    Osteosarc fixes 34 known problems in the published data, each with its
    evidence ([corrections](curation.md)). For example, it replaces the MAP2 vaccine target's
    allele with the complex event that Tempus and CeGaT report and that the
    [tumor reads support](tour.md). Each load checks every correction against
    the snapshot's source records. Pass `corrections=False` to
    `Dataset.open`, or use `osteosarc --no-corrections`, to see the published values.

New to the dataset? [Key concepts](concepts.md) explains where the data comes from,
sample IDs, variant statuses and coordinates.

## Get started

### Install

```sh
python -m pip install osteosarc
```

You need Python 3.9+ on Linux or macOS, and `samtools` on your PATH to fetch
reads.

### 1. Save the metadata

```python
from osteosarc import Dataset

data = Dataset.sync()
```

This downloads about 57 MB of the website's metadata, saved as a snapshot named
by today's date (UTC). Sequencing files stay remote. Running it again the same day
reuses the snapshot.

### 2. Choose a sample and assay

```python
print(data.describe_samples())
```

`T0_tumor` is the primary tumor, collected at timepoint **T0** (2022-12-16), and
`T0_blood` is blood from the same visit. The sample type is its `tissue` (`tumor`,
`blood` or `organoid`), and what was sequenced is its `assay`.

`T0_tumor` has **bulk RNA-seq, whole-exome (WES) and whole-genome (WGS)** data.
Its bulk RNA alignments:

```python
rna = data.assets_for_sample("T0_tumor", kind="alignment", assay="rna-seq")
for asset in rna:
    print(asset.key)
```

`rna-seq` means bulk RNA and `scrna-seq` single-cell RNA; `T1_tumor` has both.
You can also pick by platform, such as `platform="ont"` or `platform="pacbio"`. See
[Find samples and files](explore.md)
for more filters and the full assay vocabulary.

### 3. Select variants

```python
targets = data.variants(gene="DYNC1H1", status="ready")
for variant in targets:
    print(variant.id, variant.allele)
```

`status="ready"` keeps variants with one usable allele, given as
`(chromosome, one-based position, REF, ALT)`. It doesn't mean the variant has
reads, is somatic or changes the protein. Leave it out to see every entry; see
[all statuses](variants.md#variant-status).

### 4. Fetch reads for those variants

```python
source = rna["rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam"]
reads = data.extract_reads(source, variants=targets, padding=100)
print(reads.path)
```

This copies the reads around those variants into a small local BAM, after checking
that the BAM and the variants use the same genome build. Asking again reuses the
cached copy. [Isovar](consumers.md#isovar) can tell you which reads carry the variant.

### 5. Pick up where you left off

```python
data = Dataset.open()  # Your most recent snapshot, offline
counts = data.table("vafs").select(gene="SMC5")
print(counts.rows[:2])
```

The website changes over time; your snapshot doesn't. Run `Dataset.sync()` on a
later day to save a new one, and `Dataset.snapshots()` to list them.
`Dataset.open(date="2026-09")` reopens the newest from that month. Pass
`offline=False` to download more files or fetch new reads. See
[snapshots and cache](design.md) for more.

## Prefer the terminal?

```sh
osteosarc sync
osteosarc samples
osteosarc assets --sample T0_tumor --kind alignment --assay rna-seq
osteosarc variants --gene SMC5
osteosarc explore
```

Every command uses your most recent snapshot; `--snapshot 2026-09` picks the newest
from that month. In the explorer, type `help` for commands and `quit` to leave. The
[command-line guide](cli.md) lists every command.
