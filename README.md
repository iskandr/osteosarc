# osteosarc

Python library and command-line tool for the public [osteosarc.com](https://osteosarc.com/data/)
dataset: one patient's osteosarcoma sequencing, variant calls, cancer vaccines and
clinical history. Find a sample's files, pick variants and vaccine peptides, and
fetch the reads around a variant without downloading a whole BAM.

[Documentation](https://iskandr.github.io/osteosarc/) ·
[Key concepts](https://iskandr.github.io/osteosarc/concepts/) ·
[Command line](https://iskandr.github.io/osteosarc/command-line/) ·
[Python API](https://iskandr.github.io/osteosarc/python-api/) ·
[Changelog](https://iskandr.github.io/osteosarc/changelog/)

## Features

- **Samples and their files.** Every tumor, organoid and blood sample, what was
  sequenced, and its BAMs and FASTQ folders, with the commands that fetch them.
- **Browse without downloading.** Search nearly 400,000 files by sample, kind,
  assay and folder, and see which you've already downloaded.
- **Variants and vaccine peptides.** The site's variants with checked alleles, read
  counts, which pipelines found them, vaccine peptides and ELISPOT results.
- **Reads around a variant.** Copy just the reads you need out of a remote BAM into a
  small local one.
- **Corrected by default.** 35 fixes to known problems in the published data, each
  with its evidence, such as the MAP2 vaccine target's allele. Every load checks
  them against the snapshot's sources, and you can turn them off.
- **Clinical timeline.** Every treatment, procedure, scan and MRD result on one
  chart, with a row per drug.
- **Reproducible.** The site's metadata is saved as dated snapshots that reopen
  offline. The website changes; your results don't, until you sync again.
- **OpenVax integration.** Works with Varcode, Isovar, Topiary and Vaxrank, builds
  reproducible test BAMs, and includes a list of 637 candidate structural variants.

## Install

```sh
python -m pip install osteosarc
```

Needs Python 3.9+ on Linux or macOS, and [SAMtools](https://www.htslib.org/) on
your PATH to fetch reads.

## Look around

```sh
osteosarc sync                # Once: about 57 MB of the website's metadata
osteosarc                     # Which snapshot you're using, and every command
osteosarc samples             # Samples and what was sequenced
osteosarc samples T1_tumor    # One sample's files, and commands to get them
osteosarc files               # What's in the bucket, by kind and folder
osteosarc variants --gene MAP2
osteosarc timeline            # Treatments, procedures, scans and MRD
```

Every command prints text for people, and `--json` for scripts. `osteosarc repl`
opens Python with the data loaded as `data`, and in Python or a notebook everything
shows a readable preview:

```python
from osteosarc import Dataset

data = Dataset.sync()   # or Dataset.open() to reopen your newest snapshot offline
data                    # What's here, and how to get it
data.samples            # Samples and their sequencing
data.samples["T1_tumor"]
```

## From samples to reads

```python
from osteosarc import Dataset

data = Dataset.sync()  # Save today's website metadata (about 57 MB)

rna = data.samples["T0_tumor"].files.select(kind="alignment", assay="rna-seq")
targets = data.variants(gene="DYNC1H1", status="ready")

source = rna["rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam"]
reads = data.extract_reads(source, variants=targets, padding=100)
print(reads.path)  # Local indexed BAM
```

The BAM key is the file's path in the dataset's public S3 bucket. Only the reads
near the variants are downloaded, into a small indexed BAM in your local cache.
`data.download(key, to=".")` fetches a whole file instead, and `data.downloads()`
lists what you have.

Later, `Dataset.open()` reopens your most recent snapshot without a network
connection. The website changes over time, so snapshots are saved by download date:
`osteosarc snapshots` lists them, and `Dataset.open(date="2026-09")` or
`--snapshot 2026-09` picks the newest from that month.
[Get started](https://iskandr.github.io/osteosarc/#get-started) explains each step.

The same workflow from the terminal:

```sh
osteosarc sync
osteosarc samples T0_tumor
osteosarc files --sample T0_tumor --kind alignment --assay rna-seq
osteosarc variants --gene DYNC1H1 --status ready
osteosarc reads rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam --variant DYNC1H1-chr14-101980529 --variant DYNC1H1-chr14-102030200 --padding 100
osteosarc downloads
```

## Guides

| I want to… | Read |
| --- | --- |
| Understand samples, files, variant status, coordinates and corrections | [Key concepts](https://iskandr.github.io/osteosarc/concepts/) |
| Find a sample's RNA, DNA, single-cell or long-read files, and download them | [Samples and files](https://iskandr.github.io/osteosarc/samples/) |
| Get alleles, read counts or vaccine peptides | [Select variants](https://iskandr.github.io/osteosarc/variants/) |
| Fetch, filter or pair reads by variant or region | [Extract reads](https://iskandr.github.io/osteosarc/reads/) |
| Browse treatments, MRD and labs | [Browse the timeline](https://iskandr.github.io/osteosarc/timeline/) |
| Pass data to Varcode, Isovar, Topiary or Vaxrank | [Use other libraries](https://iskandr.github.io/osteosarc/openvax/) |
| Build small, verifiable test BAMs | [Read fixtures](https://iskandr.github.io/osteosarc/test-data/) |
| Explore candidate structural variants | [SV candidates](https://iskandr.github.io/osteosarc/sv-candidates/) |

## Data, license and citation

Osteosarc applies [source corrections](https://iskandr.github.io/osteosarc/corrections/)
by default. Use `Dataset.open(corrections=False)` or
`osteosarc --no-corrections` to see the published values.

Code is Apache-2.0. The dataset is listed as CC0-1.0 in the
[AWS Open Data Registry](https://registry.opendata.aws/sid-osteosarc/).
Cite the dataset and your access date when using it.

## Development

```sh
python -m pip install -e '.[test]'
ruff check osteosarc tests scripts
python -m pytest -q
```

See [testing](https://iskandr.github.io/osteosarc/testing/) for documentation
builds and live-example checks.
