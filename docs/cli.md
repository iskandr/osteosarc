# Command line

The `osteosarc` command does the same things as the Python API and shares its
cache. Every command uses your most recent snapshot unless you pass `--snapshot`.
`samples`, `specimens`, `assets`, `snapshots`, `timeline` and `on` print text (the
first five also have `--json`); the others print JSON.
`osteosarc --help` and `osteosarc COMMAND --help` list every option.

## Global options

Global options go before the command, as in `osteosarc --offline variants`.

| Option | Effect |
| --- | --- |
| `--cache DIR` | Use this cache root instead of `OSTEOSARC_CACHE` or the shared OpenVax cache |
| `--offline` | Never use the network |
| `--no-corrections` | Show the website's values unchanged; see [corrections](curation.md) |
| `--version` | Print the installed version |

Only `sync`, `download`, `table`, `reads`, `discover` and `fixtures generate` use
the network.

## Snapshots

```sh
osteosarc sync
osteosarc snapshots
osteosarc curation --strict
```

`sync` saves the website's metadata as a snapshot named by today's UTC date, such
as `2026-09-24`; running it again that day reuses it. `--refresh` downloads a new
one (`2026-09-24.2`), and `sync NAME` gives it your own name. `sync NAME
--source-revision COMMIT` takes the tables that come from the website's GitLab
repository at one full, 40-character commit. `snapshots` lists what you have,
newest first.

To use an older snapshot, pass `--snapshot`. A date, month or year picks the newest
snapshot downloaded then (`--snapshot 2026-09` is the newest from September 2026);
a name or ID prefix picks exactly one:

```sh
osteosarc variants --gene MAP2 --snapshot 2026-09
```

`curation --strict` fails if a correction no longer matches the snapshot's data,
or the data uses a label Osteosarc doesn't know. It checks the snapshot you have;
run `osteosarc sync --refresh` first to check the live website.

## Browse samples, files and variants

```sh
osteosarc samples
osteosarc samples --assay scrna-seq --platform ont
osteosarc specimens T2_tumor
osteosarc assets --sample T0_tumor --kind alignment --assay rna-seq
osteosarc assets --timepoint T2 --assay rna-seq --limit 5
osteosarc assets --sample T0_tumor --kind alignment --json
osteosarc variants --gene MAP2
osteosarc variants --set vaccine --status ready
osteosarc vaccines
osteosarc timepoints
```

| Command | Prints |
| --- | --- |
| `samples` | Samples with their sequencing and file counts. Filter with `--timepoint`, `--tissue`, `--assay` and `--platform` |
| `specimens [SAMPLE_ID]` | All specimens, or one specimen's files, corrections and nearby events |
| `assets` | Matching files with their full keys, 50 at a time (`--limit`). Filter with `--sample`, `--kind`, `--format`, `--prefix`, `--contains`, `--timepoint`, `--assay`, `--platform`, `--tissue`, `--provider` and `--library` |
| `variants` | Variants (JSON). `--set` is `site` (default), `all` or `vaccine`; filter with `--gene`, `--vaccine`, `--pipeline`, `--status` and `--vaccine-source` |
| `vaccines` | Vaccine peptides with ELISPOT results (JSON) |
| `timepoints` | Timepoints and their dates (JSON) |

The sequencing column of `samples` uses the filter names. `rna-seq; wes; wgs;
scrna-seq (ont, pacbio)` means bulk RNA, exome and genome data, plus single-cell
RNA, some of it Oxford Nanopore or PacBio. `--assay scrna-seq` picks all of that
single-cell data and `--platform ont` narrows it. A value no file uses is an error
that lists the valid ones.

See [Find samples and files](explore.md) and [Select variants](variants.md) for
what the filters mean.

## Download files and tables

```sh
osteosarc download snv_top
osteosarc table vaf_columns
```

`download` saves one whole file to the cache and prints its path. `table` prints
the rows of a named table (`vafs`, `vaf_columns`, `snv_top`, `dna_fusions`,
`rna_fusions`) or of any CSV or TSV file in the bucket, as JSON.

## Fetch reads

```sh
osteosarc reads rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam --variant DYNC1H1-chr14-101980529 --padding 100
osteosarc reads rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam chr14:101980429-101980629 --assembly GRCh38 --min-mapq 20 --exclude-flags 0x500
```

Give variant IDs with `--variant` (repeat it for several) and optional
`--padding`, or `contig:start-end` regions with `--assembly`. The two commands above
fetch the same bases. Regions are **one-based and inclusive**, like SAMtools. The
command prints the local BAM path, its index and the receipt.

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
osteosarc timeline --since 2024-05 --until 2024-09
osteosarc timeline --lane MRD --since 2025 --list
osteosarc on 2025-01-28 --days 5
```

`timeline` draws a text chart; `--list` prints one line per event, and `--json`
the event records. Filter with `--since`, `--until`, `--lane` and `--contains`.
`on` lists events within a few days of a date. See [Browse the timeline](timeline.md).

## Interactive explorer

```sh
osteosarc explore
```

The explorer opens with a snapshot summary. Its commands:

| Command | Shows |
| --- | --- |
| `summary` | Counts of corrections, events, variants and files |
| `samples [timepoint=T2] [assay=scrna-seq]` | Samples and their sequencing |
| `specimens`, `specimen SAMPLE_ID` | The registry, or one specimen's files, disagreements and nearby events |
| `assets [key=value ...]` | Files, for example `assets sample=T1_tumor kind=alignment` |
| `variants [GENE] [status=ready] [vaccine=mRNA]` | Variants |
| `corrections [ID]` | Every correction's status, or one correction's changes and evidence |
| `timeline [SINCE [UNTIL]] [lane=TEXT]` | The timeline chart |
| `zoom SINCE [UNTIL]`, `only TEXT`, `reset` | Set or clear the date window and lane filter |
| `lanes`, `events [TEXT]`, `on DATE [DAYS]` | Lane names, matching events, and events near a date |
| `help [COMMAND]`, `quit` | Help, and leave |

`zoom` and `only` persist until `reset`.

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

See [Read fixtures and bundles](fixtures.md).

## List newer bucket files

```sh
osteosarc discover neoantigen_prediction/pvactools/
```

`discover` lists the bucket under a prefix, without changing any snapshot. It
reuses a listing it already fetched; add `--refresh` for the bucket's current
contents.
See [Look for newer files](explore.md#look-for-newer-files).
