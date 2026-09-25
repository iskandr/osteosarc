# Key concepts

The ideas the rest of the documentation assumes, each linking to its full guide.

## Where the data comes from

Everything Osteosarc reads is public:

- **osteosarc.com**: the variants page, read-count tables, file listing, vaccine
  and timeline data, and analysis tables.
- **The dataset's S3 bucket**, `sid-sijbrandij-osteosarc-dataset`, which holds every
  sequencing file listed on osteosarc.com.
- **The website's source code** on GitLab
  ([slowkow/osteosarc.com](https://gitlab.com/slowkow/osteosarc.com)), for eight
  tables the site is built from but doesn't offer as downloads: the variant data
  file, the specimen registry, sample and file metadata, the timeline sheet, the
  timepoint summary, and the imaging and pathology indexes.

Beyond these, Osteosarc only uses public reference genomes and gene annotation
(Ensembl, NCBI, UCSC) to check alleles and name genes.

## Snapshots

The website changes: files get added, variants renamed and records fixed. A
snapshot is a copy of the site's metadata as downloaded on one day: the file
listing, variants, samples, vaccines and timeline, with checksums. A snapshot never
changes, so your results don't shift under you, and it works offline. Sequencing
files stay remote until you ask for them.

- `Dataset.sync()` saves today's metadata as a snapshot named by the UTC date, such
  as `2026-09-24`. Running it again that day reopens it. `refresh=True` downloads a
  new one.
- `Dataset.open()` reopens the most recently downloaded snapshot, offline. Pass
  `offline=False` to allow new downloads.
- `Dataset.snapshots()` lists them. `Dataset.open(date="2026-09")` picks the
  newest snapshot downloaded in a month, day or year. `Dataset.open(name)` opens
  one exactly, by name or ID prefix.

On the command line, `osteosarc sync` and `osteosarc snapshots` do the same, and every
other command takes `--snapshot`. Every read extraction and fixture bundle records the
ID of the snapshot it came from (`data.id`). Downloads go in a shared OpenVax cache
that other tools can reuse. See [Snapshots and cache](design.md).

## Samples, timepoints and assays

Separate fields describe each sample and each file's sequencing:

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

Every file in the bucket is an `Asset`, with a `key` (its path in the bucket), a
`url`, a `kind` (`alignment`, `reads`, `variants`, `table`, `expression` and
others), a `format`, and what the site's pages say about its sample. When the site's
pages disagree about a file, filters skip it unless you ask for conflicting values.
See [Check sample metadata](explore.md#check-sample-metadata).

## Variants and their status

Each variant has an ID such as `DYNC1H1-chr14-101980529`. `data.variants()` gives
the variants page, `"vaccine"` the ones in a vaccine, and `"all"` adds entries that
appear only in the site's other files.

`status="ready"` means the entry has one usable allele. It doesn't mean the variant
is somatic, has reads or changes the protein. See [Variant status](variants.md#variant-status).

## Coordinates

| Where | Convention | Example |
| --- | --- | --- |
| `variant.allele` | One-based position with a VCF-style anchored REF/ALT | `("chr14", 101980529, "G", "A")` |
| `Region` in Python | Zero-based, half-open | `Region("chr14", 101980528, 101980530, "GRCh38")` |
| CLI regions and `Region.from_samtools` | One-based, inclusive | `chr14:101980529-101980530` |
| Fixture SV breakends | Zero-based interbase | `{"contig": "chr9", "position": 22510295}` |

The two regions above are the same two bases. Osteosarc never converts coordinates
between genome builds: it checks each BAM's build and refuses regions from another,
so GRCh37 BAMs need GRCh37 coordinates. See [Specify coordinates](reads.md#specify-coordinates).

## Corrections

Osteosarc applies 32 documented corrections to the published records by default.
They fix misplaced alleles, relabeled samples, wrong specimen dates and wrongly
mapped read counts. Each load checks that a correction still matches its source.
A correction the source has since fixed is reported as `fixed_upstream`. An
unexpected change is reported as `stale`, and that correction is skipped.
Corrected records carry the correction IDs. Open a snapshot with
`corrections=False` to see the published values. See [Source corrections](curation.md).

## Missing is not zero

Osteosarc keeps values exactly as published and never fills in a gap with zero:

- Counts cleared by a correction are `""` (not measured).
- `NA`, `0` and an empty cell stay different.
- An ELISPOT that wasn't run isn't a negative result.
- A `below_loq` value is the lab's reporting limit, not a measurement.
- Finding no reads in a limited or incomplete fetch doesn't mean there are none.

## Receipts and caching

Every download and read extraction keeps a receipt with checksums and sizes, and
repeating a request reuses the cached result, even offline. If a cached file no
longer matches its receipt you get an `IntegrityError`; if you're offline and a
file isn't cached, an `OfflineError`. See [Snapshots and cache](design.md) and
[Reuse the result offline](reads.md#reuse-the-result-offline).
