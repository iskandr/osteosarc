# Samples and files

The dataset has two kinds of things:

- **Samples** are what was collected: tumor biopsies and resections at time points
  T0 to T3, an organoid grown from the T1 tumor, and blood draws.
- **Files** are what's in the public S3 bucket: aligned reads (BAMs), raw reads
  (FASTQs), variant calls, expression tables, scans and more, nearly 400,000 in all.

Each sample has BAMs and FASTQ folders. Most files belong to no sample: imaging,
pathology slides, analyses and references. The examples open your most recent
snapshot; see [Get started](index.md#get-started).

## Samples

```python
from osteosarc import Dataset

data = Dataset.open()
data.samples
data.samples.select(timepoint="T0", tissue="tumor")
data.samples.select(assay="scrna-seq", platform="ont")
```

In a notebook or REPL each of these shows a table. The single-cell Oxford Nanopore
samples in the 2026-09-25 snapshot:

```text
3 tumor samples

sample    timepoint  date        where  sequencing                               BAMs  FASTQ folders
--------  ---------  ----------  -----  ---------------------------------------  ----  -------------
T1_tumor  T1         2024-06-06  UCLA   rna-seq, wes, wgs, scrna-seq (ont,       14    14
                                        pacbio)
T2_tumor  T2         2025-01-28  UCLA   rna-seq, wgs, scrna-seq (ont), cite-seq  8     7
T3_tumor  T3         2025-04-17  MSKCC  scrna-seq (ont), cite-seq                2     5
```

The sequencing column uses the same names as the file filters below. One
sequencing run can have several BAMs, such as a vendor's and a reprocessed one.

`data.samples["T1_tumor"]` shows one sample: where it came from, each BAM (with its
assay, platform, provider, size and whether you've downloaded it) and each FASTQ
folder (with its assay, platform, library, file count and size), followed by the
calls that fetch them. Assays and platforms use the filter names below everywhere;
the library says which of a single-cell capture's libraries a folder holds.

From the command line:

```sh
osteosarc samples                                  # the table above, for every sample
osteosarc samples T1_tumor                         # one sample, with commands to get its files
osteosarc samples --files                          # every sample's BAMs and FASTQ folders
osteosarc samples --tissue blood --assay cite-seq  # filter by timepoint, tissue, assay or platform
```

### Sample ID, time point and sample type

These fields describe different things:

| Field | Example | Meaning |
| --- | --- | --- |
| Snapshot | `2026-09-25` | When you downloaded the website's metadata, not a collection date |
| Sample ID | `T0_tumor` | A sample, or a fraction of one |
| `timepoint` | `T0` | When it was collected; tumor and blood from one visit share it |
| `tissue` | `tumor` | Sample type: `tumor`, `blood`, or `organoid` |
| `assay` | `rna-seq` | What was sequenced, and whether bulk or single-cell |
| `platform` | `ont` | Sequencing technology, when the site says |
| `provider` | `BostonGene` | Who produced the data |

For example, `T0_tumor` is the primary tumor resection on 2022-12-16, while
`T0_blood` is blood from the same time point. `T1_tumor` and `T1_organoid` also
share a time point but have different sample types. `T3_tumor_CD45neg` is a
CD45-negative enriched fraction; its `tissue` is still `tumor`. Later blood draws
are named by date, such as `blood_2025-06-26`, and have no time point.

A `Sample` has these fields, among others:

```python
t1 = data.samples["T1_tumor"]
print(t1.id, t1.timepoint, t1.date, t1.tissue, t1.site, t1.description)
print(t1.assays)          # ('rna-seq', 'wes', 'wgs', 'scrna-seq')
print(t1.providers)
print(t1.bams[:2])        # file keys
print(t1.fastq_folders[:2])
```

`disagreements` lists dates or sites that other sources give differently, and
`corrections` names the [corrections](corrections.md) that changed the sample.

### Where sequencing comes from

A sample's sequencing combines the site's sample registry with its FASTQ table, which
records each folder's assay. The four 2026 blood draws have nothing in the registry,
but their FASTQ folders say they're single-cell gene expression, T-cell receptor and
CITE-seq libraries.

The site's labels map to the filter names like this. Using a site label as a filter
gives an error that names the right one:

| Site label | `assay` | Library | Data |
| --- | --- | --- | --- |
| `RNA`, `bulk RNA` | `rna-seq` | | Bulk RNA sequencing |
| `WES` | `wes` | | Bulk whole-exome DNA sequencing |
| `WGS` | `wgs` | | Bulk whole-genome DNA sequencing |
| `scRNA`, `scRNA_GEX` | `scrna-seq` | gene expression | Single-cell RNA sequencing |
| `scRNA_TCR`, `scRNA_TCRgd`, `scRNA_BCR` | `scrna-seq` | αβ TCR, γδ TCR, BCR | T- and B-cell receptor libraries from the same cells |
| `scRNA_ONT` | `scrna-seq`, with `platform="ont"` | long reads | Single-cell Oxford Nanopore RNA sequencing |
| `PacBio` | `scrna-seq`, with `platform="pacbio"` | long reads | Single-cell PacBio RNA sequencing in this dataset |
| `CITE` | `cite-seq` | antibody tags | The antibody-tag library of a CITE-seq capture |

Single cells aren't samples: cell barcodes are inside the single-cell files.

## A sample's files

`sample.files` holds its BAMs and every file in its FASTQ folders:

```python
t0 = data.samples["T0_tumor"]
rna = t0.files.select(kind="alignment", assay="rna-seq")
for file in rna:
    print(file.key, file.size)
```

`T1_tumor` has both bulk and single-cell RNA, so pick one with `assay`:

```python
t1 = data.samples["T1_tumor"]
bulk = t1.files.select(kind="alignment", assay="rna-seq")
single_cell = t1.files.select(kind="alignment", assay="scrna-seq")
print(len(bulk), len(single_cell), len(single_cell.select(platform="ont")))
```

Many files have no platform label; the T0 alignments, for instance, don't. A
`platform="illumina"` filter leaves them out, so filter by assay instead, and check a
BAM's header if you need its instrument.

`data.files.select(sample="T1_tumor")` gives the same files, and `file.samples`
names the samples a file belongs to.

## Find any file

`osteosarc files` with no filters summarizes the bucket: how many files of each
kind, how big they are, and the largest folders. Add filters to list files:

```sh
osteosarc files
osteosarc files --kind alignment --assay wgs
osteosarc files --sample T1_tumor --kind alignment
osteosarc files --prefix hudson_lab/PBMC_scRNAseq/FASTQ/Jan2026/
```

Lists show each file's complete key, which is what `osteosarc download` and
`osteosarc reads` take, and `yes` in the `local` column for files you've downloaded.
Add `--json` for full records. In Python:

```python
ont = data.files.select(kind="alignment", platform="ont")
dna = data.files.select(kind="alignment", assay="wgs")
vcfs = data.files.select(kind="variants", format="vcf")
fastqs = data.files.select(kind="reads")
scans = data.files.select(kind="image", format="dcm")
pvac = data.files.select(prefix="neoantigen_prediction/pvactools/", format="tsv")
```

Kinds are `alignment`, `reads`, `variants`, `expression`, `annotation`, `table`,
`reference`, `index`, `image` and `other`. `prefix` and `contains` search by path,
whatever a file's labels say. To search every sample at one time point, use
`data.files.select(timepoint="T2", assay="rna-seq")`.

Use an exact key for one file:

```python
source = data.file("rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam")
print(source.size, source.samples, source.index_urls)
```

## Get a file

Most of the time you want a region of a BAM, not the whole thing:
[extract reads](reads.md) streams just those reads into a small local BAM. To
download a whole file:

<!-- docs-check: skip (downloads a 1 GB BAM) -->
```python
data = Dataset.open(offline=False)
path = data.download(source, to=".")   # the BAM, and its index, in this folder
print(data.local_path(source))         # the cached copy, or None; never downloads
print(data.downloads())                # everything downloaded, with local paths
```

```sh
osteosarc download rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam --to .
osteosarc downloads
osteosarc files --downloaded
```

Downloads are kept in the cache under names derived from their contents; `to` puts a
file in a folder of your choice under its own name, as a read-only hard link where it
can, so the file isn't stored twice. For a whole FASTQ folder, the AWS command line is
faster; `osteosarc samples SAMPLE` prints the command, such as:

```sh
aws s3 cp --recursive --no-sign-request s3://sid-sijbrandij-osteosarc-dataset/hudson_lab/PBMC_scRNAseq/FASTQ/Jan2026/capture2/CITE/ CITE/
```

## Check file metadata

```python
print(source.claims)
print(source.values("timepoint"))
print(source.conflicts)
print(data.claims.rows[:2])
```

Each claim is what one of the site's sources says about a file's sample. When they
disagree, filters skip the file. Pass `include_conflicts=True` to match any of the
values, or `include_inferred=True` to also use values guessed from the file's path.
`data.claims` lists every claim with the files it's made for.

## Read a table

```python
data = Dataset.open(offline=False)
path = data.download("snv_top")
table = data.table("snv_top")
print(path, table.columns)
print(table.rows[:2])

counts = data.table("vafs").select(gene="SMC5")
print(counts.rows[:2])
```

Named tables include `vafs`, `vaf_columns`, `snv_top`, `dna_fusions`, and
`rna_fusions`; `osteosarc table NAME` prints one as TSV. You can also pass a file or
an exact key:

```python
reports = data.files.select(contains=".genes.results", format="tsv")
if reports:
    expression = data.table(reports[0])
    print(expression.rows[:2])
```

Values stay as text, so `"0"`, `"NA"` and `""` stay different.
`table.to_dataframe()` needs pandas. `data.parse(file)` also reads JSON and FASTA.
These read the whole file into memory; for big files, open the downloaded path with
a streaming reader.

## Look for newer files

A snapshot uses the file listing the website published. To list the bucket directly:

```python
from osteosarc import Cache, list_bucket

listing = list_bucket(Cache(), "ONT/", refresh=True)
print(len(listing["files"]))
```

`osteosarc discover ONT/` does the same from the command line. This lists everything
under the prefix without changing your snapshot. See [snapshots and cache](snapshots.md)
to refresh metadata or import files you already have.
