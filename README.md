# osteosarc

Python library and command-line tool for the public [osteosarc.com](https://osteosarc.com/data/)
dataset: one patient's osteosarcoma sequencing, variant calls, cancer vaccines and
clinical history. Find files, select variants and vaccine peptides, and fetch reads
around a variant without downloading an entire BAM.

[Documentation](https://iskandr.github.io/osteosarc/) ·
[Key concepts](https://iskandr.github.io/osteosarc/concepts/) ·
[Command line](https://iskandr.github.io/osteosarc/cli/) ·
[Python API](https://iskandr.github.io/osteosarc/api/) ·
[Changelog](https://iskandr.github.io/osteosarc/changelog/)

## Features

- **Browse without downloading.** Search nearly 400,000 files by sample, timepoint and
  assay: bulk and single-cell RNA, exome, genome, Oxford Nanopore and PacBio.
- **Variants and vaccine peptides.** Catalogue variants with checked genomic alleles,
  read counts, pipeline detections, vaccine peptides and ELISPOT results.
- **Reads around a variant.** Indexed queries copy only the reads you need from a
  remote BAM into a cached, indexed local BAM.
- **Corrected by default.** 32 documented, evidence-backed fixes to the published
  data, such as the MAP2 vaccine target's allele. Every load checks them against
  the snapshot's sources, and you can turn them off.
- **Clinical timeline.** Treatments, procedures, imaging, MRD and lab results as a
  text chart or in an interactive terminal explorer.
- **Reproducible.** Named metadata snapshots with SHA-256 receipts reopen offline.
- **OpenVax integration.** Adapters for Varcode, Isovar, Topiary and Vaxrank,
  versioned read-fixture recipes, and a catalogue of 637 structural-variant candidates.

## Install

```sh
python -m pip install osteosarc
```

Requires Python 3.9+ on Linux or macOS. Fetching reads also needs
[SAMtools](https://www.htslib.org/) on PATH.

## Quickstart

```python
from osteosarc import Dataset

data = Dataset.sync("baseline")  # Save about 57 MB of metadata as "baseline"
print(data.describe_samples())

rna = data.assets_for_sample("T0_tumor", kind="alignment", assay="rna-seq")
targets = data.variants(gene="DYNC1H1", status="ready")

source = rna["rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam"]
reads = data.extract_reads(source, variants=targets, padding=100)
print(reads.path)  # Local indexed BAM
```

Later, `Dataset.open("baseline")` reopens the snapshot without a network connection.
[Get started](https://iskandr.github.io/osteosarc/#get-started) explains each step.

The same workflow from the terminal:

```sh
osteosarc sync baseline
osteosarc samples baseline
osteosarc assets baseline --sample T0_tumor --kind alignment --assay rna-seq
osteosarc variants baseline --gene DYNC1H1 --status ready
osteosarc reads baseline rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam --variant DYNC1H1-chr14-101980529 --variant DYNC1H1-chr14-102030200 --padding 100
osteosarc explore baseline
```

`osteosarc explore` opens an interactive browser for specimens, files, variants and
the timeline. Type `help` for commands and `quit` to leave.

## Guides

| I want to… | Read |
| --- | --- |
| Understand sample IDs, variant status, coordinates and corrections | [Key concepts](https://iskandr.github.io/osteosarc/concepts/) |
| Find RNA, DNA, single-cell or long-read files and read tables | [Find samples and files](https://iskandr.github.io/osteosarc/explore/) |
| Get alleles, read counts or vaccine peptides | [Select variants](https://iskandr.github.io/osteosarc/variants/) |
| Fetch, filter or pair reads by variant or region | [Extract reads](https://iskandr.github.io/osteosarc/reads/) |
| Browse treatments, specimens, MRD and labs | [Browse the timeline](https://iskandr.github.io/osteosarc/timeline/) |
| Pass data to Varcode, Isovar, Topiary or Vaxrank | [Use other libraries](https://iskandr.github.io/osteosarc/consumers/) |
| Build small, verifiable test BAMs | [Read fixtures](https://iskandr.github.io/osteosarc/fixtures/) |
| Explore candidate structural variants | [SV interest catalogue](https://iskandr.github.io/osteosarc/sv-interest/) |

## Data, license and citation

Osteosarc applies [source corrections](https://iskandr.github.io/osteosarc/curation/)
by default. Use `Dataset.open("baseline", corrections=False)` or
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

See [testing](https://iskandr.github.io/osteosarc/validation/) for documentation
builds and live-example checks.
