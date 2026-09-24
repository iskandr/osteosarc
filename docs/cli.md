# Command line

The `osteosarc` command calls the same `Dataset` API as Python, so both produce
the same selections and share one cache. Commands that list data print JSON;
`samples`, `specimens`, `timeline` and `on` print readable text. Run
`osteosarc --help` or `osteosarc COMMAND --help` for every option.
Every command uses your most recent snapshot unless you pass `--snapshot`.
See the [changelog](changelog.md) for when each option arrived.

## Global options

Global options go before the command, as in `osteosarc --offline variants`.

| Option | Effect |
| --- | --- |
| `--cache DIR` | Use this cache root instead of `OSTEOSARC_CACHE` or the shared OpenVax cache |
| `--offline` | Forbid all network access |
| `--no-corrections` | Use the published sources unchanged; see [corrections](curation.md) |
| `--version` | Print the installed version |

Listing commands never download data. `download`, `table` and `reads` fetch what
they need unless `--offline` is set.

## Snapshots

```sh
osteosarc sync
osteosarc snapshots
osteosarc curation --strict
```

`sync` saves the website's current metadata as a snapshot named by today's UTC
date, such as `2026-09-24`; repeating it the same day reopens that snapshot.
`--refresh` downloads a new one (`2026-09-24.2`), and `sync NAME` saves a snapshot
under a name you choose. `--source-revision COMMIT` pins the site's GitLab sources
to one commit. `snapshots` lists what you have, newest first, and names the default.

Every other command uses the most recent snapshot. Pick another with `--snapshot`:
an exact name, a UTC download year, month or day (`--snapshot 2026-09` opens the
newest from September 2026), or an ID prefix:

```sh
osteosarc variants --gene MAP2 --snapshot 2026-09
```

`curation --strict` exits nonzero if a correction no longer matches its source or
a source label is unrecognized. The older form with the snapshot name first, as in
`osteosarc samples baseline`, still works but is deprecated.

## Browse samples, files and variants

```sh
osteosarc samples --timepoint T0
osteosarc specimens T2_tumor
osteosarc assets --sample T0_tumor --kind alignment --assay rna-seq
osteosarc assets --timepoint T2 --assay rna-seq --limit 5
osteosarc variants --gene MAP2
osteosarc variants --set vaccine --status ready
osteosarc vaccines
osteosarc timepoints
```

| Command | Prints |
| --- | --- |
| `samples` | Specimens with their sequencing types and file counts; `--json` prints the original per-file sample claims |
| `specimens [SAMPLE_ID]` | The specimen registry, or one specimen's files, corrections and nearby events |
| `assets` | The total and the first `--limit` (default 50) matching files. Filter with `--sample`, `--kind`, `--format`, `--prefix`, `--contains`, `--timepoint`, `--assay`, `--platform`, `--tissue`, `--provider` and `--library` |
| `variants` | Catalogue entries. Choose `--set` `site` (default), `all` or `vaccine`, and filter with `--gene`, `--vaccine`, `--pipeline`, `--status` and `--vaccine-source` |
| `vaccines` | Vaccine-overlap rows with ELISPOT results |
| `timepoints` | Published timepoint and date pairs |

See [Find samples and files](explore.md) and [Select variants](variants.md) for
what the filters mean.

## Download files and tables

```sh
osteosarc download snv_top
osteosarc table vaf_columns
```

`download` fetches one complete file into the cache and prints its path. `table`
parses a named table (`vafs`, `vaf_columns`, `snv_top`, `dna_fusions`,
`rna_fusions`) or any CSV/TSV asset key, and prints its rows.

## Fetch reads

```sh
osteosarc reads rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam --variant DYNC1H1-chr14-101980529 --padding 100
osteosarc reads rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam chr14:101980429-101980629 --assembly GRCh38 --min-mapq 20 --exclude-flags 0x500
```

Give catalogue variants with `--variant` (repeat it for several) and optional
`--padding`, or give `contig:start-end` regions with `--assembly`. The two commands
above request the same bases. Regions are
**one-based and inclusive**, like SAMtools. The command prints the local BAM path,
its index and the receipt.

| Option | Effect |
| --- | --- |
| `--min-mapq N`, `--exclude-flags FLAGS` | Filter reads; flags accept hex such as `0x500` |
| `--fetch-pairs` | Also retrieve paired mates outside the regions |
| `--recover-linked` | Follow mate and `SA` links within [bounded limits](reads.md#recover-mates-and-split-alignments) |
| `--reference FASTA` | Local indexed reference, required for CRAM |
| `--index PATH_OR_URL` | Explicit index when the catalogue lists none |
| `--reference-length N` | Expected contig length for mitochondrial regions |

See [Extract reads](reads.md) for requirements and caching.

## Browse the timeline

```sh
osteosarc timeline --since 2024-05 --until 2024-09
osteosarc timeline --lane MRD --since 2025 --list
osteosarc on 2025-01-28 --days 5
```

`timeline` draws a text chart; `--list` prints one line per event and `--json`
prints event records. Filter with `--since`, `--until`, `--lane` and `--contains`.
`on` lists events within some days of a date. See [Browse the timeline](timeline.md).

## Interactive explorer

```sh
osteosarc explore
```

The explorer opens with a snapshot summary. Its commands:

| Command | Shows |
| --- | --- |
| `summary` | Snapshot, correction, timeline, variant and asset counts |
| `samples [timepoint=T2] [tissue=tumor]` | Specimens and sequencing types |
| `specimens`, `specimen SAMPLE_ID` | The registry, or one specimen's files, disagreements and nearby events |
| `assets [key=value ...]` | Files, for example `assets kind=alignment timepoint=T1 assay=rna-seq` |
| `variants [GENE] [status=ready] [vaccine=mRNA]` | Catalogue entries |
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
| `fixtures select RECIPE --source ID=BAM` | Print a recipe's record membership and reasons |
| `fixtures generate RECIPE OUTPUT` | Acquire declared inputs, select, and write a bundle |
| `fixtures pack RECIPE OUTPUT --source ID=BAM` | Write a bundle from local inputs |
| `fixtures verify BUNDLE [--sha256 HASH]` | Check a bundle offline |
| `fixtures list BUNDLE` | Member status, record counts and reasons |
| `fixtures export BUNDLE OUTPUT --member NAME` | Write named, indexed BAM (or `--format sam`) exports |

See [Read fixtures and bundles](fixtures.md).

## List newer bucket files

```sh
osteosarc discover neoantigen_prediction/pvactools/
```

`discover` lists a live S3 prefix with receipts, without changing any snapshot.
See [Look for newer files](explore.md#look-for-newer-files).
