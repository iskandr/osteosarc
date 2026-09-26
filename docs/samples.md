# Samples and files

- **Samples** are what was collected: tumor biopsies and resections at time points T0
  to T3, an organoid grown from the T1 tumor, and blood draws.
- **Files** are what's in the public S3 bucket: aligned reads (BAMs), raw reads
  (FASTQs), variant calls, expression tables, scans and more, nearly 400,000 in all.

Each sample has BAMs and FASTQ folders. Most files belong to no sample: imaging,
pathology slides, analyses and references. The examples open your most recent
snapshot; see [get started](index.md#get-started).

## Samples

```sh
osteosarc samples                                  # every sample, and what was sequenced
osteosarc samples T1_tumor                         # one sample's files, with commands to get them
osteosarc samples --files                          # every sample's BAMs and FASTQ folders
osteosarc samples --tissue blood --assay cite-seq  # filter by time point, tissue, assay or platform
```

The same in Python, where each line shows a table in a notebook or REPL:

```python
from osteosarc import Dataset

data = Dataset.open()
data.samples
data.samples.select(timepoint="T0", tissue="tumor")
data.samples.select(assay="scrna-seq", platform="ont")
```

The single-cell Oxford Nanopore samples in the 2026-09-25 snapshot:

```text
3 tumor samples

sample    timepoint  date        where  sequencing                               BAMs  FASTQ folders
--------  ---------  ----------  -----  ---------------------------------------  ----  -------------
T1_tumor  T1         2024-06-06  UCLA   rna-seq, wes, wgs, scrna-seq (ont,       14    14
                                        pacbio)
T2_tumor  T2         2025-01-28  UCLA   rna-seq, wgs, scrna-seq (ont), cite-seq  8     7
T3_tumor  T3         2025-04-17  MSKCC  scrna-seq (ont), cite-seq                2     5
```

One sequencing run can have several BAMs, such as a vendor's and a reprocessed one.
`data.samples["T1_tumor"]` shows one sample: where it came from, each BAM with its
assay, provider, size and whether you've downloaded it, each FASTQ folder with its
library, file count and size, and the calls that fetch them.

### What the fields mean

| Field | Example | Meaning |
| --- | --- | --- |
| Sample ID | T0_tumor | A sample, or a fraction of one |
| Time point | T0 | When it was collected; tumor and blood from one visit share it |
| Tissue | tumor | tumor, blood or organoid |
| Assay | rna-seq | What was sequenced, and whether bulk or single-cell |
| Platform | ont | The sequencing technology, when the site says |
| Provider | BostonGene | Who produced the data |

T0_tumor is the primary tumor resection on 2022-12-16, and T0_blood is blood from
the same visit. T1_tumor and T1_organoid share a time point but not a tissue.
T3_tumor_CD45neg is a CD45-negative fraction of T3_tumor. Later blood draws are
named by date, such as blood_2025-06-26, and have no time point. A snapshot's date
is when you downloaded the website's metadata, not a collection date.

In Python, a sample's fields are attributes:

```python
t1 = data.samples["T1_tumor"]
print(t1.id, t1.timepoint, t1.date, t1.tissue, t1.site, t1.description)
print(t1.assays, t1.providers)
print(t1.bams[:2], t1.fastq_folders[:2])
```

A sample's disagreements list dates or sites that other sources give differently,
and its corrections name the [fixes](corrections.md) that changed it.

### Where sequencing comes from

A sample's sequencing combines the site's sample registry with its FASTQ table,
which records each folder's assay. The four 2026 blood draws have nothing in the
registry, but their FASTQ folders show single-cell gene expression, T-cell receptor
and CITE-seq libraries.

The site's labels map to osteosarc's names like this; using a site label as a filter
gives an error that names the right one:

| Site label | Assay | Library | Data |
| --- | --- | --- | --- |
| RNA, bulk RNA | rna-seq | | Bulk RNA sequencing |
| WES | wes | | Whole-exome DNA sequencing |
| WGS | wgs | | Whole-genome DNA sequencing |
| scRNA, scRNA_GEX | scrna-seq | gene expression | Single-cell RNA sequencing |
| scRNA_TCR, scRNA_TCRgd, scRNA_BCR | scrna-seq | αβ TCR, γδ TCR, BCR | T- and B-cell receptor libraries from the same cells |
| scRNA_ONT | scrna-seq, platform ont | long reads | Single-cell Oxford Nanopore RNA sequencing |
| PacBio | scrna-seq, platform pacbio | long reads | Single-cell PacBio RNA sequencing |
| CITE | cite-seq | antibody tags | The antibody-tag library of a CITE-seq capture |

Single cells aren't samples: their barcodes are inside the single-cell files.

## A sample's files

A sample's files are its BAMs and everything in its FASTQ folders:

```python
t0 = data.samples["T0_tumor"]
for file in t0.files.select(kind="alignment", assay="rna-seq"):
    print(file.key, file.size)
```

T1_tumor has both bulk and single-cell RNA, so pick one by assay, and long reads by
platform:

```python
t1 = data.samples["T1_tumor"]
bulk = t1.files.select(kind="alignment", assay="rna-seq")
single_cell = t1.files.select(kind="alignment", assay="scrna-seq")
print(len(bulk), len(single_cell), len(single_cell.select(platform="ont")))
```

Many files have no platform label; the T0 alignments, for instance, don't. Filtering
by platform illumina leaves them out, so filter by assay instead.

## Find any file

`osteosarc files` with no filters summarizes the bucket: how many files of each
kind, their size, and the largest folders. Filters list files:

```sh
osteosarc files
osteosarc files --kind alignment --assay wgs
osteosarc files --sample T1_tumor --kind alignment
osteosarc files --prefix hudson_lab/PBMC_scRNAseq/FASTQ/Jan2026/
```

Each file is listed with its complete key, which is what `osteosarc download` and
`osteosarc reads` take, and marked local if you've downloaded it. In Python:

```python
ont = data.files.select(kind="alignment", platform="ont")
vcfs = data.files.select(kind="variants", format="vcf")
scans = data.files.select(kind="image", format="dcm")
pvac = data.files.select(prefix="neoantigen_prediction/pvactools/", format="tsv")
source = data.file("rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam")
print(source.size, source.samples)
```

The kinds are alignment, reads, variants, expression, annotation, table, reference,
index, image and other. A prefix or text in the key finds files whatever their
labels say.

When the site's sources disagree about a file's sample, filters skip it;
`include_conflicts=True` matches any of their values. A file's claims show what each
source says.

## Get a file

Most of the time you want a region of a BAM, not the whole thing:
[reads](reads.md) streams just those reads into a small local BAM. To download a
whole file:

<!-- docs-check: skip (downloads a 7 GB BAM) -->
```sh
osteosarc download rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam --to .
```

```sh
osteosarc downloads
osteosarc files --downloaded
```

`--to` puts the file, and its index, in a folder under its own name, as a read-only
hard link to the cached copy where it can, so it isn't stored twice. `downloads`
lists what you have and where. In Python, `data.download(file, to=".")`,
`data.local_path(file)` and `data.downloads()` do the same. For a whole FASTQ folder
the AWS command line is faster; `osteosarc samples SAMPLE` prints the command.

## Read a table

The site's read counts, one row per variant and BAM:

```python
counts = data.vafs.select(gene="SMC5")
print(counts.rows[:2])
```

Other tables download when you parse them: the site's own (snv_top, dna_fusions,
rna_fusions) or any CSV or TSV in the bucket.

```python
data = Dataset.open(offline=False)
top = data.parse("snv_top")
print(top.columns, top.rows[:2])
```

Values stay as text, so "0", "NA" and an empty cell stay different; `to_dataframe()`
turns a table into a pandas data frame.

## Look for newer files

A snapshot uses the file listing the website published. To list a folder of the
bucket as it is now, without changing your snapshot:

```python
from osteosarc import Cache
from osteosarc.discovery import list_bucket

listing = list_bucket(Cache(), "ONT/", refresh=True)
print(len(listing["files"]))
```

`osteosarc sync --refresh` saves a new snapshot with the website's latest metadata.
