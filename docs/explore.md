# Find files and read tables

Create a snapshot with `Dataset.sync("baseline")` first, as shown in
[Get started](index.md). The examples below use that saved snapshot.

## Browse samples and sequencing types

```python
from osteosarc import Dataset

data = Dataset.open("baseline")
print(data.describe_samples())
print(data.describe_samples(timepoint="T0", tissue="tumor"))
```

The filtered overview in the checked 2026-09-18 snapshot is:

```text
sample    date        sequencing     BAMs  FASTQ_folders
--------  ----------  -------------  ----  -------------
T0_tumor  2022-12-16  RNA; WES; WGS  10    9
```

The sequencing column uses the registry's labels. BAM and FASTQ-folder counts
count file products; several files can come from the same sequencing library.

## Sample ID, timepoint, and sample type

These fields describe different things:

| Field | Example | Meaning |
| --- | --- | --- |
| Snapshot name | `baseline` | Your local name for saved metadata |
| Sample ID | `T0_tumor` | A biological specimen or fraction in the registry |
| `timepoint` | `T0` | Collection timepoint; shared by tumor and blood specimens |
| `tissue` | `tumor` | Sample type: `tumor`, `blood`, or `organoid` |
| `assay` | `rna-seq` | What was sequenced and whether it was bulk or single-cell |
| `platform` | `ont` | Sequencing technology, when the metadata specifies it |
| `provider` | `BostonGene` | Data provider |

For example, `T0_tumor` is the primary tumor resection on 2022-12-16, while
`T0_blood` is blood from the same timepoint. `T1_tumor` and `T1_organoid`
also share a timepoint but have different sample types. `T3_tumor_CD45neg`
is a CD45-negative enriched fraction; its `tissue` is still `tumor`.

Use IDs from `describe_samples()` or `data.specimens`; not every timepoint/type
combination has a registry entry. IDs and filter values are case-sensitive.

```python
for row in data.specimens:
    print(row["sample_id"], row["timepoint"], row["tissue"], row["assays"])
```

`data.specimens` is the biological registry. `data.samples` contains the
per-file metadata claims used for asset filtering.

## Bulk and single-cell data

`T0_tumor` has **bulk RNA-seq, WES and WGS** in the checked snapshot. Its
registry lists no single-cell assay. Select its bulk RNA alignments with:

```python
rna = data.assets_for_sample("T0_tumor", kind="alignment", assay="rna-seq")
for asset in rna:
    print(asset.key, asset.size)
```

`T1_tumor` has both bulk and single-cell RNA. The sample ID alone does not
choose between them:

```python
bulk = data.assets_for_sample("T1_tumor", kind="alignment", assay="rna-seq")
single_cell = data.assets_for_sample("T1_tumor", kind="alignment", assay="scrna-seq")
ont_single_cell = single_cell.select(platform="ont")
pacbio_single_cell = single_cell.select(platform="pacbio")
print(len(bulk), len(single_cell), len(ont_single_cell), len(pacbio_single_cell))
```

| Registry label | API `assay` | Data |
| --- | --- | --- |
| `RNA` | `rna-seq` | Bulk RNA sequencing |
| `WES` | `wes` | Bulk whole-exome DNA sequencing |
| `WGS` | `wgs` | Bulk whole-genome DNA sequencing |
| `scRNA` | `scrna-seq` | Single-cell RNA sequencing |
| `scRNA_ONT` | `scrna-seq`, with `platform="ont"` | Single-cell Oxford Nanopore RNA sequencing |
| `PacBio` | `scrna-seq`, with `platform="pacbio"` | Single-cell PacBio RNA sequencing in this dataset |
| `CITE` | `cite-seq` | Single-cell RNA and antibody-tag profiling |

A sample ID identifies a specimen or fraction; individual cell barcodes live
in the single-cell data. Each asset is one file, and reprocessed alignments
remain separate products.

Platforms are `ont`, `pacbio`, and `illumina`. Platform metadata is incomplete:
the T0 alignments have assay labels but no platform labels in this snapshot.
An `illumina` filter therefore excludes them. Use the assay filter above and
inspect an alignment's header if you need its instrument platform.

`assets_for_sample` uses registry links and FASTQ folders. Unassigned files
remain available through `data.assets`. To search all samples at a timepoint,
use `data.assets.select(timepoint="T2", assay="rna-seq")`.

## Find other file types

```python
ont = data.assets.select(kind="alignment", platform="ont")
dna = data.assets.select(kind="alignment", assay="wgs")
vcfs = data.assets.select(kind="variants", format="vcf")
fastqs = data.assets.select(kind="reads")
pvac = data.assets.select(prefix="neoantigen_prediction/pvactools/", format="tsv")
```

You can also search by path with `prefix` or
`contains`, including files whose sample metadata is unknown.

Use an exact key when choosing an alignment for [read extraction](reads.md):

```python
source = data.asset(
    "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.bam"
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

Metadata filters exclude conflicting values by default. Pass
`include_conflicts=True` to match any published claim, or `include_inferred=True`
to include values inferred from paths. Use [specimens](timeline.md#find-a-specimens-files)
for biological samples and their associated files.

## Download and read a table

```python
data = Dataset.open("baseline", offline=False)
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

Tables preserve strings: `"0"`, `"NA"`, and `""` stay distinct.
`table.to_dataframe()` requires pandas. `data.parse(asset)` also reads JSON and
FASTA. These parsers load their input into memory; for large files, pass the
downloaded path to a streaming reader.

## Look for newer files

The snapshot uses the website's dated bucket listing. To query S3 directly:

```python
from osteosarc import Cache, list_bucket

listing = list_bucket(Cache(), "ONT/", refresh=True)
print(len(listing["files"]))
```

This returns all pages for the prefix, with receipts, without changing your
snapshot. See [snapshots and cache](design.md) to refresh metadata or import
files you already have.
