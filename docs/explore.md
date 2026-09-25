# Find samples, files and tables

Browse the samples, find sequencing files by sample, assay, platform or path
without downloading them, and read the site's tables.
The examples open your most recent snapshot; see [Get started](index.md#get-started).

## Browse samples and sequencing types

```python
from osteosarc import Dataset

data = Dataset.open()
print(data.describe_samples())
print(data.describe_samples(timepoint="T0", tissue="tumor"))
print(data.describe_samples(assay="scrna-seq", platform="ont"))
```

The overview of single-cell Oxford Nanopore samples in the 2026-09-24 snapshot is:

```text
sample    date        sequencing                                  BAMs  FASTQ_folders
--------  ----------  ------------------------------------------  ----  -------------
T1_tumor  2024-06-06  rna-seq; wes; wgs; scrna-seq (ont, pacbio)  14    14
T2_tumor  2025-01-28  rna-seq; wgs; scrna-seq (ont)               8     7
T3_tumor  2025-04-17  scrna-seq (ont)                             2     5
```

The sequencing column uses the same names as the file filters. The BAM and
FASTQ-folder columns count files, and one sequencing run can have several (for
example, a vendor BAM and a reprocessed one). From the command line, use
`osteosarc samples` with the same filters.

## Sample ID, timepoint, and sample type

These fields describe different things:

| Field | Example | Meaning |
| --- | --- | --- |
| Snapshot | `2026-09-24` | When you downloaded the website's metadata, not a collection date |
| Sample ID | `T0_tumor` | A specimen, or a fraction of one |
| `timepoint` | `T0` | When it was collected; tumor and blood from one visit share it |
| `tissue` | `tumor` | Sample type: `tumor`, `blood`, or `organoid` |
| `assay` | `rna-seq` | What was sequenced, and whether bulk or single-cell |
| `platform` | `ont` | Sequencing technology, when the site says |
| `provider` | `BostonGene` | Who produced the data |

For example, `T0_tumor` is the primary tumor resection on 2022-12-16, while
`T0_blood` is blood from the same timepoint. `T1_tumor` and `T1_organoid`
also share a timepoint but have different sample types. `T3_tumor_CD45neg`
is a CD45-negative enriched fraction; its `tissue` is still `tumor`.

Take IDs from `describe_samples()` or `data.specimens`; not every timepoint has
every sample type. IDs and filter values are case-sensitive.

```python
for row in data.specimens:
    print(row["sample_id"], row["timepoint"], row["tissue"], row["assays"])
```

`data.specimens` lists the specimens. `data.samples` lists what each of the site's
pages says about each file, which is what the file filters use.

## Bulk and single-cell data

`T0_tumor` has **bulk RNA-seq, WES and WGS**, and no single-cell data. Its bulk
RNA alignments:

```python
rna = data.assets_for_sample("T0_tumor", kind="alignment", assay="rna-seq")
for asset in rna:
    print(asset.key, asset.size)
```

`T1_tumor` has both bulk and single-cell RNA, so pick one with `assay`:

```python
bulk = data.assets_for_sample("T1_tumor", kind="alignment", assay="rna-seq")
single_cell = data.assets_for_sample("T1_tumor", kind="alignment", assay="scrna-seq")
ont_single_cell = single_cell.select(platform="ont")
pacbio_single_cell = single_cell.select(platform="pacbio")
print(len(bulk), len(single_cell), len(ont_single_cell), len(pacbio_single_cell))
```

The site labels sequencing differently from the filters. Its labels map like this,
and using one as a filter gives an error that names the right filter:

| Registry label | API `assay` | Data |
| --- | --- | --- |
| `RNA` | `rna-seq` | Bulk RNA sequencing |
| `WES` | `wes` | Bulk whole-exome DNA sequencing |
| `WGS` | `wgs` | Bulk whole-genome DNA sequencing |
| `scRNA` | `scrna-seq` | Single-cell RNA sequencing |
| `scRNA_ONT` | `scrna-seq`, with `platform="ont"` | Single-cell Oxford Nanopore RNA sequencing |
| `PacBio` | `scrna-seq`, with `platform="pacbio"` | Single-cell PacBio RNA sequencing in this dataset |
| `CITE` | `cite-seq` | Single-cell RNA and antibody-tag profiling |

Single cells aren't samples: cell barcodes are inside the single-cell files.

Many files have no platform label; the T0 alignments, for instance, don't. A
`platform="illumina"` filter leaves them out, so filter by assay instead, and check
a BAM's header if you need its instrument.

`assets_for_sample` finds a sample's files through the specimen registry and its
FASTQ folders. Files not linked to any sample are still in `data.assets`. To
search every sample at one timepoint, use
`data.assets.select(timepoint="T2", assay="rna-seq")`.

The command line takes the same filters and prints a table with each file's
complete key, which is what `osteosarc reads` takes. Add `--json` for full records:

```sh
osteosarc assets --sample T1_tumor --kind alignment --assay scrna-seq --platform ont
osteosarc assets --timepoint T2 --assay rna-seq --limit 5
```

## Find other file types

```python
ont = data.assets.select(kind="alignment", platform="ont")
dna = data.assets.select(kind="alignment", assay="wgs")
vcfs = data.assets.select(kind="variants", format="vcf")
fastqs = data.assets.select(kind="reads")
pvac = data.assets.select(prefix="neoantigen_prediction/pvactools/", format="tsv")
```

You can also search by path with `prefix` or `contains`, which finds files
whatever their sample labels.

Use an exact key when choosing an alignment for [read extraction](reads.md):

```python
source = data.asset(
    "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam"
)
print(source.index_urls)
```

## Check sample metadata

```python
print(source.claims)
print(source.values("timepoint"))
print(source.conflicts)
print(data.samples.rows[:2])
```

When the site's pages disagree about a file, filters skip it. Pass
`include_conflicts=True` to match any of the values, or `include_inferred=True` to
also use values guessed from the file's path. Use [specimens](timeline.md#find-a-specimens-files)
for biological samples and their associated files.

## Download and read a table

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
`rna_fusions`. You can also pass an asset or exact key:

```python
reports = data.assets.select(contains=".genes.results", format="tsv")
if reports:
    expression = data.table(reports[0])
    print(expression.rows[:2])
```

Values stay as text, so `"0"`, `"NA"` and `""` stay different.
`table.to_dataframe()` needs pandas. `data.parse(asset)` also reads JSON and FASTA.
These read the whole file into memory; for big files, open the downloaded path
with a streaming reader.

## Look for newer files

A snapshot uses the file listing the website published. To list the bucket
directly:

```python
from osteosarc import Cache, list_bucket

listing = list_bucket(Cache(), "ONT/", refresh=True)
print(len(listing["files"]))
```

This lists everything under the prefix without changing your snapshot. See
[snapshots and cache](design.md) to refresh metadata or import files you already
have.
