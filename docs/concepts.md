# Key concepts

This page covers the ideas that the rest of the documentation assumes. Each section
links to its full guide.

## Snapshots

The osteosarc.com website changes: files are added, variants are renamed and
records are corrected. A snapshot is a local copy of its metadata as downloaded on
one date: the file listing, variants, samples, vaccines and timeline sources, with
SHA-256 receipts. A snapshot never changes, so your results don't shift under you,
and it works offline. Sequencing files stay remote until you request them.

- `Dataset.sync()` saves today's metadata as a snapshot named by the UTC date, such
  as `2026-09-24`. Running it again that day reopens it. `refresh=True` downloads a
  new one.
- `Dataset.open()` reopens the most recently downloaded snapshot, offline. Pass
  `offline=False` to allow new downloads.
- `Dataset.snapshots()` lists them. `Dataset.open("2026-09")` picks the newest
  snapshot downloaded in a month; a day, a year, an exact name or an ID prefix also
  work.

On the command line, `osteosarc sync` and `osteosarc snapshots` do the same, and every
other command takes `--snapshot`. Each snapshot's content ID (`data.id`) is recorded in
every read extraction and fixture bundle made from it. Files are stored in a shared
OpenVax cache, so other tools can reuse them. See [Snapshots and cache](design.md).

## Samples, timepoints and assays

Several separate fields describe each sample, and each file's sequencing:

| Field | Example | Meaning |
| --- | --- | --- |
| Sample ID | `T0_tumor` | A biological specimen or fraction |
| `timepoint` | `T0` | Collection timepoint, shared by tumor and blood specimens |
| `tissue` | `tumor` | Sample type: `tumor`, `blood`, or `organoid` |
| `assay` | `rna-seq` | What was sequenced, and whether bulk or single-cell |
| `platform` | `ont` | Sequencing technology, when the metadata states it |

`rna-seq` means **bulk** RNA; `scrna-seq` means **single-cell** RNA. The other
assays are `wes` (whole exome), `wgs` (whole genome) and `cite-seq`. Platforms are
`illumina`, `ont` (Oxford Nanopore) and `pacbio`, but many files have no platform
label. `describe_samples()` lists every sample ID. See
[Find samples and files](explore.md).

## Assets

Every object in the public bucket is an `Asset`. It has a `key` (its path), `url`,
`kind` (`alignment`, `reads`, `variants`, `table`, `expression` and others),
`format`, and the metadata claims that sources make about it. Filters match only
unambiguous published claims unless you opt in to conflicting or path-inferred
values. See [Check sample metadata](explore.md#check-sample-metadata).

## Variants and their status

Each catalogue entry has a stable ID such as `DYNC1H1-chr14-101980529`.
`data.variants()` selects the website's entries, `"vaccine"` selects vaccinated
entries, and `"all"` adds entries found only in the count export or source JSON.

`status="ready"` means the entry has one literal genomic allele with consistent
coordinates. It doesn't mean the allele is somatic, supported by reads, or
changes the protein. Other statuses explain why an allele is unusable.
See [Variant status](variants.md#variant-status).

## Coordinates

| Where | Convention | Example |
| --- | --- | --- |
| `variant.allele` | One-based position with a VCF-style anchored REF/ALT | `("chr14", 101980529, "G", "A")` |
| `Region` in Python | Zero-based, half-open | `Region("chr14", 101980528, 101980530, "GRCh38")` |
| CLI regions and `Region.from_samtools` | One-based, inclusive | `chr14:101980529-101980530` |
| Fixture SV breakends | Zero-based interbase | `{"contig": "chr9", "position": 22510295}` |

The two regions above describe the same two bases. Osteosarc never lifts coordinates
over between assemblies. Read extraction checks the alignment header and rejects a
region whose assembly doesn't match. GRCh37 alignments need GRCh37 coordinates.
See [Specify coordinates](reads.md#specify-coordinates).

## Corrections

Osteosarc applies 32 documented corrections to the published records by default.
They fix misplaced alleles, relabeled samples, wrong specimen dates and wrongly
mapped read counts. Each load checks that a correction still matches its source.
A correction the source has since fixed is reported as `fixed_upstream`. An
unexpected change is reported as `stale`, and that correction is skipped.
Corrected records carry the correction IDs. Open a snapshot with
`corrections=False` to see the published values. See [Source corrections](curation.md).

## Missing is not zero

Osteosarc keeps source values as strings and marks what wasn't measured:

- Count rows cleared by a correction are `""` (unmeasured), not zero.
- `NA`, `0` and an empty cell stay distinct in tables.
- An untested ELISPOT assay is not a negative result.
- A `below_loq` measurement records the reported limit, not a concentration.
- A bounded read extraction, or an incomplete RNA assay, is not evidence of zero support.

## Receipts and caching

Every download and read extraction records a receipt with checksums, sizes and
source identity. Repeating a request reuses the verified cached result, including
offline. If cached bytes don't match their receipt, you get an `IntegrityError`.
If a request needs bytes that aren't cached while offline, you get an
`OfflineError`. See [Snapshots and cache](design.md) and
[Reuse the result offline](reads.md#reuse-the-result-offline).
