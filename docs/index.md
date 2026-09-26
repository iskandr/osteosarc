# Osteosarc

Osteosarc is a Python library and command-line tool for the public
[osteosarc.com](https://osteosarc.com/data/) dataset: one patient's osteosarcoma
sequencing, variant calls, cancer vaccines and clinical history. It does two things
well:

- **Exploring.** See every sample and what was sequenced, every file in the bucket,
  the variant catalogue and the clinical timeline, from a terminal or Python.
- **Making test data.** Pull just the reads around a few variants out of a remote
  BAM into a small local one, ready for a unit test, without downloading the BAM.

It keeps a dated copy of the website's metadata, so results don't change under
you, and fixes [known errors](corrections.md) in it.

## Get started

```sh
python -m pip install osteosarc
osteosarc sync                # once: about 57 MB of the website's metadata
osteosarc                     # the snapshot in use, and every command
osteosarc samples             # samples and what was sequenced
osteosarc samples T1_tumor    # one sample's files, and commands to get them
osteosarc variants --gene MAP2
osteosarc timeline            # treatments, procedures, scans and MRD
```

Every command prints text for people, or JSON with `--json`. In Python or a
notebook, everything shows a readable preview, and the dataset itself lists what's
there:

```python
from osteosarc import Dataset

data = Dataset.sync()   # later: Dataset.open(), offline
data
data.samples["T1_tumor"]
print(data.variants(gene="MAP2"))
```

`osteosarc repl` opens Python with the data already loaded as `data`.

## Make test data

The reads around a variant in each of a sample's RNA-seq BAMs, saved in a folder
with readable names:

```sh
osteosarc reads T0_tumor --assay rna-seq --variant DYNC1H1-chr14-101980529 --padding 100 --to tests/data
```

Only those reads are fetched, with an index, so each file is a few hundred
kilobytes. A BAM on another genome build is skipped with a note. The same in
Python, for one BAM:

```python
rna = data.samples["T0_tumor"].files.select(kind="alignment", assay="rna-seq")
source = rna["rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam"]
dync1h1 = data.variants(gene="DYNC1H1", status="ready")
reads = data.extract_reads(source, variants=dync1h1, padding=100, to="tests/data")
print(reads.path)
```

Asking again reuses the cached result, even offline. [Reads](reads.md) covers
regions, filters and mates; [test data](test-data.md) covers openvax-v1, the test reads
the OpenVax libraries share, and bundles of your own.

## Where to go next

| To | Read |
| --- | --- |
| Understand samples, variants, coordinates and corrections | [Concepts](concepts.md) |
| Find a sample's files, and download them | [Samples and files](samples.md) |
| Get alleles, read counts and vaccine peptides | [Variants and vaccines](variants.md) |
| Fetch reads by variant or region | [Reads](reads.md) |
| Chart treatments, scans and MRD | [Timeline](timeline.md) |
| Use every command | [Command line](command-line.md) |
| Pass data to Varcode, Isovar, Topiary or Vaxrank | [OpenVax libraries](openvax.md) |

You need Python 3.9+ on Linux or macOS, and SAMtools on your PATH to fetch reads.
The website changes over time; a snapshot doesn't. Run `osteosarc sync` on a later
day for a new one; see [snapshots and cache](snapshots.md).
