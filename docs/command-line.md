# Command line

The `osteosarc` command does what the Python API does and shares its cache. Run
it on its own to see which snapshot it uses and every command:

```text
Browse:
  samples [SAMPLE]   Tumor, organoid and blood samples; one sample's files and how to get them
  files              Files in the S3 bucket: an overview, or a list with --kind, --sample or --prefix
  variants [ID]      The variant catalogue, with alleles, vaccines and read counts
  vaccines           Vaccine targets and ELISPOT results
  timeline           Treatments, procedures, scans and MRD over time (--around DATE for one week)
  corrections [ID]   Known problems in the website's data, and the fixes applied

Get data:
  reads FILE|SAMPLE  The reads around variants or in regions, as a small BAM: test data in seconds
  download FILE      A whole file (--to DIR puts it, and its index, in DIR)
  downloads          What's already on this computer, and where

Test data for libraries:
  test-data ...      Shared test reads (openvax-v1): list, export, and check your copies

Snapshots of the website's metadata:
  sync               Download the current metadata (about 57 MB); commands use the newest
  snapshots          Saved snapshots, by download date
  repl               Python with the newest snapshot loaded as `data`
```

Every command prints text for people; add `--json` for records to use in scripts.
Each command's `--help` lists its options, with examples.

## Make test data from real reads

```sh
osteosarc reads T0_tumor --assay rna-seq --variant DYNC1H1-chr14-101980529 --padding 100 --to test-reads
```

This streams only the reads within 100 bases of the variant from each of T0_tumor's
RNA-seq BAMs, never a whole file, and saves each one as a small indexed BAM in the
test-reads folder, named for its source and the variant. BAMs on a different genome
build, or without an index, are skipped with a note. Asking again reuses the
result, even offline.

Give a file's key instead of a sample to read one BAM, and repeat `--variant` for
several variants. For a region instead, give it one-based and inclusive, like
SAMtools, with its assembly:

```sh
osteosarc reads rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam chr14:101980429-101980629 --assembly GRCh38
```

| Option | Effect |
| --- | --- |
| `--assay`, `--platform` | With a sample, only its BAMs of that assay or platform |
| `--to DIR` | Also save each BAM and its index in DIR |
| `--min-mapq N`, `--exclude-flags 0x500` | Leave out low-quality or flagged reads |
| `--fetch-pairs` | Also fetch mates that fall outside the regions |
| `--recover-linked` | Follow mates and split reads, within limits |
| `--reference FASTA` | A local indexed reference, needed for CRAM |

Without `--to`, the command prints the path of each BAM in the cache, so you can
write `samtools view $(osteosarc reads ...)`. See [Reads](reads.md).

## Browse

```sh
osteosarc samples
osteosarc samples T1_tumor
osteosarc files --sample T0_tumor --kind alignment --assay rna-seq
osteosarc variants --gene MAP2
osteosarc variants MAP2-chr2-209694768
osteosarc vaccines
osteosarc timeline --since 2024-05 --until 2024-09
osteosarc timeline --around 2025-01-28
osteosarc corrections
```

- **samples** lists every sample with its sequencing and how many BAMs and FASTQ
  folders it has. With a sample ID, it shows that sample's files with their size and
  whether they're downloaded, followed by the commands that fetch them; `--files`
  does that for every sample.
- **files** with no filters summarizes the bucket by kind and folder. With filters it
  lists files, BAMs first, each with the key that download and reads take. Filter by
  sample, kind, format, path prefix, time point, assay, platform, tissue or
  provider; `--downloaded` keeps the files already on this computer.
- **variants** lists the variants on the site with their alleles, protein changes,
  vaccines and the pipelines that found them. With an ID, it shows one variant with
  its read counts.
- **vaccines** lists the vaccine targets and their ELISPOT results.
- **timeline** charts treatments, procedures, scans and MRD results by month;
  `--around DATE` lists everything within a week of a date, and `--list` every event.
- **corrections** lists the fixes osteosarc applies to the website's data; with an
  ID, it shows one with its evidence. `--strict` exits with an error if a fix no
  longer matches the data.

Assay and platform values are the same everywhere: rna-seq, wes, wgs, scrna-seq and
cite-seq; illumina, ont and pacbio. A value no file uses is an error that lists the
valid ones. See [Samples and files](samples.md) and [Variants](variants.md).

## Download files

```sh
osteosarc download snv_top
osteosarc downloads
```

download saves one whole file in the cache and prints its path; with `--to DIR` it
also puts the file, and a BAM's or VCF's index, in DIR under its own name.
downloads lists what you've downloaded and the reads you've extracted, with their
local paths. BAMs are large, so for tests read a region instead.

## Test data for libraries

<!-- docs-check: skip (downloads openvax-v1) -->
```sh
osteosarc test-data list openvax-v1
osteosarc test-data export openvax-v1 tests/data --member MEMBER
osteosarc test-data check openvax-v1 fixtures.json
```

openvax-v1 is the bundle of reads the OpenVax libraries share; list downloads it the
first time, export writes members as indexed BAMs, and check compares your own copies
with it. The same commands take a bundle folder, and generate builds one from your
own recipe. See [Test data](test-data.md).

## Snapshots

```sh
osteosarc sync
osteosarc snapshots
```

sync saves the website's metadata as a snapshot named by today's UTC date; running
it again that day reuses it, and `--refresh` makes a new one. Commands use the
newest snapshot; `--snapshot` picks another by download date (a year, month or
day, as in `--snapshot 2026-09`) or by name. repl opens Python with the snapshot loaded as `data`, and offers
to sync first if you haven't. See [Snapshots and cache](snapshots.md).

## Global options

These go before the command, as in `osteosarc --offline variants`.

| Option | Effect |
| --- | --- |
| `--cache DIR` | Use this cache instead of the shared OpenVax one |
| `--offline` | Never use the network |
| `--no-corrections` | Show the website's values unchanged |
| `--version` | Print the installed version |

Only sync, download, reads, repl, test-data generate, and the first use of openvax-v1
use the network.
