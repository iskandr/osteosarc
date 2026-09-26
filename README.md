# osteosarc

Python library and command-line tool for the public [osteosarc.com](https://osteosarc.com/data/)
dataset: one patient's osteosarcoma sequencing, variant calls, cancer vaccines and
clinical history.

- **Explore it:** every sample and what was sequenced, every file in the bucket,
  the variant catalogue with its vaccines, and the clinical timeline.
- **Make test data from it:** pull the reads around a few variants out of a remote
  BAM into a small local one, without downloading the BAM.

[Documentation](https://iskandr.github.io/osteosarc/) ·
[Concepts](https://iskandr.github.io/osteosarc/concepts/) ·
[Command line](https://iskandr.github.io/osteosarc/command-line/) ·
[Python API](https://iskandr.github.io/osteosarc/python-api/) ·
[Changelog](https://iskandr.github.io/osteosarc/changelog/)

## Install

```sh
python -m pip install osteosarc
```

Needs Python 3.9+ on Linux or macOS, and [SAMtools](https://www.htslib.org/) on
your PATH to fetch reads.

## Explore

```sh
osteosarc sync                # once: about 57 MB of the website's metadata
osteosarc                     # the snapshot in use, and every command
osteosarc samples             # samples and what was sequenced
osteosarc samples T1_tumor    # one sample's files, and commands to get them
osteosarc files               # what's in the bucket, by kind and folder
osteosarc variants --gene MAP2
osteosarc timeline            # treatments, procedures, scans and MRD
```

Every command prints text for people, or JSON with `--json`. In Python or a
notebook, everything shows a readable preview; `osteosarc repl` opens Python with
the data loaded:

```python
from osteosarc import Dataset

data = Dataset.sync()   # later: Dataset.open(), offline
data                    # what's here, and how to get it
data.samples["T1_tumor"]
```

## Make test data

The reads around a variant in each of a sample's RNA-seq BAMs, saved in a folder:

```sh
osteosarc reads T0_tumor --assay rna-seq --variant DYNC1H1-chr14-101980529 --padding 100 --to tests/data
```

Or in Python, for one BAM:

```python
source = data.file("rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam")
dync1h1 = data.variants(gene="DYNC1H1", status="ready")
reads = data.extract_reads(source, variants=dync1h1, padding=100, to="tests/data")
print(reads.path)
```

Only the reads near the variants are fetched, into a small indexed BAM; asking again
reuses it, even offline.

The website changes over time, so its metadata is saved as dated snapshots:
`osteosarc snapshots` lists them, and `--snapshot 2026-09` picks the newest from that
month. Osteosarc also fixes [known errors](https://iskandr.github.io/osteosarc/corrections/)
in the website's data; `osteosarc --no-corrections` shows the published values.

## Guides

| To | Read |
| --- | --- |
| Understand samples, variants, coordinates and corrections | [Concepts](https://iskandr.github.io/osteosarc/concepts/) |
| Find a sample's files, and download them | [Samples and files](https://iskandr.github.io/osteosarc/samples/) |
| Get alleles, read counts and vaccine peptides | [Variants and vaccines](https://iskandr.github.io/osteosarc/variants/) |
| Fetch reads by variant or region | [Reads](https://iskandr.github.io/osteosarc/reads/) |
| Chart treatments, scans and MRD | [Timeline](https://iskandr.github.io/osteosarc/timeline/) |
| Pass data to Varcode, Isovar, Topiary or Vaxrank | [OpenVax libraries](https://iskandr.github.io/osteosarc/openvax/) |
| Look through 637 candidate structural variants | [SV candidates](https://iskandr.github.io/osteosarc/sv-candidates/) |

## Data, license and citation

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
