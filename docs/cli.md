# Command line

The `osteosarc` command does the same things as the Python API and shares its
cache. Run it on its own to see which snapshot it's using and every command,
grouped by what you're doing:

```text
Browse:
  samples [SAMPLE]   Tumor, organoid and blood samples; one sample's files and how to get them
  files              Files in the S3 bucket: an overview, or a list with --kind, --sample or --prefix
  variants [ID]      The variant catalogue, with alleles, vaccines and read counts
  vaccines           Vaccine targets and ELISPOT results
  timeline           Treatments, procedures, scans and MRD over time
  on DATE            Everything within a week of a date
  corrections [ID]   Known problems in the website's data, and the fixes applied

Get data:
  download FILE      Download a whole file (--to DIR puts it, and its index, in DIR)
  reads FILE REGION  Stream the reads in a region, or around a variant, into a small local BAM
  downloads          What's already on this computer, and where
  table NAME         Print one of the site's tables, or a bucket CSV/TSV, as TSV

Snapshots of the website's metadata:
  sync               Download the current metadata (about 57 MB); commands use the newest
  snapshots          Saved snapshots, by download date

More:
  repl               Python with the newest snapshot loaded as `data`
  discover PREFIX    List a bucket folder live, without a snapshot
  fixtures ...       Build and check read fixtures for tests
```

Every command prints text for people; add `--json` for records to use in scripts.
`osteosarc COMMAND --help` lists a command's options, with examples.

## Global options

Global options go before the command, as in `osteosarc --offline variants`.

| Option | Effect |
| --- | --- |
| `--cache DIR` | Use this cache root instead of `OSTEOSARC_CACHE` or the shared OpenVax cache |
| `--offline` | Never use the network |
| `--no-corrections` | Show the website's values unchanged; see [corrections](curation.md) |
| `--version` | Print the installed version |

Only `sync`, `download`, `table`, `reads`, `discover`, `repl` (when you haven't
synced yet) and `fixtures generate` use the network.

## Snapshots

```sh
osteosarc sync
osteosarc snapshots
```

`sync` saves the website's metadata as a snapshot named by today's UTC date, such
as `2026-09-25`; running it again that day reuses it. `--refresh` downloads a new
one (`2026-09-25.2`), and `sync NAME` gives it your own name. `sync NAME
--source-revision COMMIT` takes the tables that come from the website's GitLab
repository at one full, 40-character commit. `snapshots` lists what you have,
newest first.

Commands use the newest snapshot. To use another, pass `--snapshot` after the
command. A date, month or year picks the newest snapshot downloaded then
(`--snapshot 2026-09` is the newest from September 2026); a name or ID prefix
picks exactly one:

```sh
osteosarc variants --gene MAP2 --snapshot 2026-09
```

## Samples and files

```sh
osteosarc samples
osteosarc samples T1_tumor
osteosarc samples --files
osteosarc samples --assay scrna-seq --platform ont
osteosarc files
osteosarc files --sample T0_tumor --kind alignment --assay rna-seq
osteosarc files --prefix vendor/tempus/ --limit 20
```

`samples` lists every sample with its sequencing and how many BAMs and FASTQ
folders it has. `samples SAMPLE` shows one: each BAM with its assay, provider, size
and whether it's downloaded, each FASTQ folder with its file count and size, and
the commands that fetch them, written out. `samples --files` lists every sample's
files.

`files` with no filters summarizes the bucket by kind and folder. With filters it
lists files, BAMs first, each with its complete key, 50 at a time (`--limit`).
Filter with `--sample`, `--kind`, `--format`, `--prefix`, `--contains`,
`--timepoint`, `--assay`, `--platform`, `--tissue`, `--provider` and `--library`;
`--downloaded` keeps the files already on this computer.

The sequencing column of `samples` uses the filter names. `rna-seq, wes, wgs,
scrna-seq (ont, pacbio)` means bulk RNA, exome and genome data, plus single-cell
RNA, some of it Oxford Nanopore or PacBio. `--assay scrna-seq` picks all of that
single-cell data and `--platform ont` narrows it. A value no file uses is an error
that lists the valid ones. See [Samples and files](samples.md).

## Variants, vaccines and corrections

```sh
osteosarc variants --gene MAP2
osteosarc variants MAP2-chr2-209694768
osteosarc variants --set vaccine --status ready
osteosarc vaccines
osteosarc corrections
osteosarc corrections allele-MAP2-chr2-209694768
osteosarc corrections --strict
```

`variants` lists variants with their alleles, protein changes, vaccines and the
pipelines that found them. `--set` is `site` (default), `all` or `vaccine`; filter
with `--gene`, `--vaccine`, `--pipeline` and `--status`. `variants ID` shows one
variant with its corrections and the read counts the site published.
`vaccines` lists the vaccine targets and their ELISPOT results. See
[Select variants](variants.md).

`corrections` lists the fixes Osteosarc applies to the website's data, and
`corrections ID` shows one with the records it checks and its evidence.
`--strict` exits with an error if a correction no longer matches the snapshot's
data, or the data uses a label Osteosarc doesn't know. It checks the snapshot you
have; run `osteosarc sync --refresh` first to check the live website. See
[Source corrections](curation.md).

## Get data

<!-- docs-check: skip (downloads a 1 GB BAM) -->
```sh
osteosarc download rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam --to .
```

```sh
osteosarc download snv_top
osteosarc downloads
osteosarc table vaf_columns
```

`download` saves one whole file in the cache and prints its path. With `--to DIR`
it also puts the file, and a BAM's or VCF's index, in `DIR` under their own names
(as hard links where possible, so they aren't stored twice). `downloads` lists
what you've downloaded and the reads you've extracted, with their local paths.

`table` prints a named table (`vafs`, `vaf_columns`, `snv_top`, `dna_fusions`,
`rna_fusions`) or any CSV or TSV file in the bucket, as TSV, or as JSON with
`--json`.

## Fetch reads

```sh
osteosarc reads rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam --variant DYNC1H1-chr14-101980529 --padding 100
osteosarc reads rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam chr14:101980429-101980629 --assembly GRCh38 --min-mapq 20 --exclude-flags 0x500
```

Give variant IDs with `--variant` (repeat it for several) and optional
`--padding`, or `contig:start-end` regions with `--assembly`. The two commands above
fetch the same bases. Regions are **one-based and inclusive**, like SAMtools. Only
those reads are downloaded, into a small indexed BAM in the cache; the command
prints its path, so you can write `samtools view $(osteosarc reads ...)`. Add
`--json` for its index and receipt.

| Option | Effect |
| --- | --- |
| `--min-mapq N`, `--exclude-flags FLAGS` | Filter reads; flags can be hex, such as `0x500` |
| `--fetch-pairs` | Also fetch mates outside the regions |
| `--recover-linked` | Follow mates and split reads, [within limits](reads.md#recover-mates-and-split-alignments) |
| `--reference FASTA` | Local indexed reference, needed for CRAM |
| `--index PATH_OR_URL` | The index, if the file listing doesn't have one |
| `--reference-length N` | Expected chromosome length, for mitochondrial regions |

See [Extract reads](reads.md) for requirements and caching.

## Browse the timeline

```sh
osteosarc timeline
osteosarc timeline --since 2024-05 --until 2024-09
osteosarc timeline --lane MRD --since 2025 --list
osteosarc on 2025-01-28 --days 5
```

`timeline` draws a chart with a row per treatment. Filter with `--since`,
`--until`, `--lane` and `--contains`; `--all` also charts lab draws, DICOM studies
and other frequent records. `--list` prints one line per event, and `--json` the
event records. `on` lists events within a few days of a date. See
[Browse the timeline](timeline.md).

## Python

```sh
osteosarc repl
```

`repl` opens Python (IPython, if it's installed) with the newest snapshot loaded as
`data`, and shows what it holds. If you haven't synced yet, it offers to download
the website's metadata first. See the [Python API](api.md).

## Test fixtures

```sh
osteosarc fixtures panel vaccine-loci-v1
```

| Command | Effect |
| --- | --- |
| `fixtures panel NAME` | Print a shipped target panel: `vaccine-loci-v1`, `sv-regressions-v1` or `sv-interest-v1` |
| `fixtures select RECIPE --source ID=BAM` | Show which reads a recipe picks, and why |
| `fixtures generate RECIPE OUTPUT` | Fetch the reads, apply the recipe, and write a bundle |
| `fixtures pack RECIPE OUTPUT --source ID=BAM` | Write a bundle from local BAMs |
| `fixtures verify BUNDLE [--sha256 HASH]` | Check a bundle offline |
| `fixtures list BUNDLE` | Each member's status, read count and reasons |
| `fixtures export BUNDLE OUTPUT --member NAME` | Write a member as an indexed BAM (or `--format sam`) |

These print JSON. See [Read fixtures and bundles](fixtures.md).

## List newer bucket files

```sh
osteosarc discover neoantigen_prediction/pvactools/
```

`discover` lists the bucket under a prefix, without changing any snapshot. It
reuses a listing it already fetched; add `--refresh` for the bucket's current
contents. See [Look for newer files](samples.md#look-for-newer-files).
