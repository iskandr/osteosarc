# Osteosarc

Osteosarc is a Python library and command-line tool for the public
[osteosarc.com](https://osteosarc.com/data/) dataset: one patient's osteosarcoma
sequencing, variant calls, cancer vaccines and clinical history. It does two things
well:

- **Exploring.** See every sample and what was sequenced, every file in the bucket,
  the variant catalogue and the clinical timeline, from a terminal or Python.
- **Making test data.** Pull just the reads around a few variants out of a remote
  BAM into a small local one, or a bundle of them, ready for unit tests, without
  downloading any BAM whole.

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

`osteosarc repl` opens Python with the data already loaded as `data`. You need
Python 3.9+ on Linux or macOS, and SAMtools on your PATH to fetch reads.

## Make test data

The reads around a variant in each of a sample's RNA-seq BAMs, saved in a folder
with readable names:

```sh
osteosarc reads T0_tumor --assay rna-seq --variant DYNC1H1-chr14-101980529 --padding 100 --to tests/data
```

Only those reads are fetched, with an index, so each file is a few hundred
kilobytes, and asking again reuses them, even offline. For a bundle that anyone can
rebuild, with a balanced set of reads at each variant and every record pinned, use
`osteosarc test-data make`; see [test data](test-data.md), and [reads](reads.md) for
regions, filters and mates.

## How the data fits together

**Where it comes from.** Everything is public: osteosarc.com's pages and tables, the
dataset's S3 bucket (sid-sijbrandij-osteosarc-dataset), which holds every sequencing
file, and the [website's source](https://gitlab.com/slowkow/osteosarc.com) for tables
the site doesn't offer as downloads. Reference genomes (Ensembl, NCBI, UCSC) check
alleles.

**Snapshots.** The website changes; a snapshot, its metadata as downloaded on one
day, doesn't. Commands use the newest; `--snapshot 2026-09` picks another, and every
read extraction records the snapshot it came from. See
[snapshots and cache](snapshots.md).

**Samples and files.** A sample is something collected from the patient: a tumor, an
organoid grown from one, or a blood draw, with a time point such as T0. A file is
anything in the bucket; samples have BAMs and FASTQ folders, while scans, slides and
analyses belong to none. Assays are rna-seq (bulk), scrna-seq (single-cell), wes, wgs
and cite-seq; platforms are illumina, ont and pacbio, when the site says. See
[samples and files](samples.md).

**Variants.** Each has an ID such as DYNC1H1-chr14-101980529; a ready variant has one
usable allele, which doesn't make it somatic, covered by reads or protein-changing.
See [variants and vaccines](variants.md).

**Coordinates.** Alleles are one-based, VCF-style. Regions are zero-based and
half-open in Python, and one-based and inclusive on the command line, like SAMtools:
`Region("chr14", 101980528, 101980530, "GRCh38")` and chr14:101980529-101980530 are
the same two bases. Osteosarc never converts between genome builds; it checks each
BAM's build and refuses regions from another.

**Corrections.** Osteosarc fixes 35 known problems in the website's data by default,
such as misplaced alleles, mislabeled samples and wrong dates, checking each when a
snapshot opens; `--no-corrections` shows the published values. See
[corrections](corrections.md).

**Missing is not zero.** Values stay as published: NA, 0 and an empty cell stay
different, an ELISPOT that wasn't run isn't negative, and no reads in a limited fetch
doesn't mean there are none. Every download and extraction keeps a receipt with
checksums, and a cached file that no longer matches it is an error.

## Where to go next

| To | Read |
| --- | --- |
| Find a sample's files, and download them | [Samples and files](samples.md) |
| Get alleles, read counts and vaccine peptides | [Variants and vaccines](variants.md) |
| Fetch reads by variant or region | [Reads](reads.md) |
| Make test data, or use the OpenVax libraries' | [Test data](test-data.md) |
| Chart treatments, scans and MRD | [Timeline](timeline.md) |
| Use every command | [Command line](command-line.md) |
| Pass data to Varcode, Isovar, Topiary or Vaxrank | [OpenVax libraries](openvax.md) |
