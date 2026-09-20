# Find files and read tables

Create a snapshot with `Dataset.sync("baseline")` first, as shown in
[Get started](index.md). The examples below use that saved snapshot.

## Browse samples and sequencing types

```python
from osteosarc import Dataset

data = Dataset.open("baseline")
print(data.describe_samples())
print(data.describe_samples(timepoint="T2", tissue="tumor"))
```

The overview lists the registry's sequencing types, BAM counts, and FASTQ
folder counts. Counts are file products, not independent samples.

## Find a sample's RNA alignments

```python
rna = data.assets_for_sample("T0_tumor", kind="alignment", assay="rna-seq")
for asset in rna:
    print(asset.key, asset.size)
```

Each asset is a file. Different processing products from the same sample have
separate entries.

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

Assays are `rna-seq`, `scrna-seq`, `cite-seq`, `wgs`, and `wes`. Platforms are
`ont`, `pacbio`, and `illumina`. You can also search by path with `prefix` or
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
