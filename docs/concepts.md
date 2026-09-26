# Concepts

The ideas the rest of the documentation assumes.

## Where the data comes from

Everything osteosarc reads is public:

- **osteosarc.com**: the variants page, read counts, file listing, vaccines,
  timeline and analysis tables.
- **The dataset's S3 bucket** (sid-sijbrandij-osteosarc-dataset), which holds every
  sequencing file the site lists.
- **The website's source** on [GitLab](https://gitlab.com/slowkow/osteosarc.com), for
  eight tables the site is built from but doesn't offer as downloads, such as the
  sample registry and the timeline sheet.

Osteosarc also uses public reference genomes (Ensembl, NCBI, UCSC) to check alleles.

## Snapshots

The website changes: files get added, variants renamed, records fixed. A snapshot is
a copy of its metadata as downloaded on one day. It never changes, works offline,
and leaves the sequencing files remote until you ask for one.

`osteosarc sync` (or `Dataset.sync()`) saves one, named by the date; commands and
`Dataset.open()` use the newest. `osteosarc snapshots` lists them, and
`--snapshot 2026-09` picks the newest from a month. Every read extraction records the
snapshot it came from. See [snapshots and cache](snapshots.md).

## Samples and files

A sample is something collected from the patient: a tumor biopsy or resection, an
organoid grown from one, or a blood draw. A file is anything in the bucket. Samples
have files (BAMs of aligned reads and folders of FASTQ raw reads), but most files,
such as scans, slides and analyses, belong to no sample.

| Field | Example | Meaning |
| --- | --- | --- |
| Sample ID | T0_tumor | A sample, or a fraction of one |
| Time point | T0 | When it was collected; tumor and blood from one visit share it |
| Tissue | tumor | tumor, blood or organoid |
| Assay | rna-seq | What was sequenced, and whether bulk or single-cell |
| Platform | ont | The sequencing technology, when the site says |

rna-seq means bulk RNA and scrna-seq single-cell RNA; the other assays are wes
(whole exome), wgs (whole genome) and cite-seq. Platforms are illumina, ont (Oxford
Nanopore) and pacbio, but many files have no platform label.

`osteosarc samples` and `data.samples` list the samples; `osteosarc files` and
`data.files` list the files. See [samples and files](samples.md).

## Variants

Each variant has an ID such as DYNC1H1-chr14-101980529. `data.variants()` gives the
variants page; a variant whose status is ready has one usable allele. That doesn't
mean it's somatic, has reads or changes a protein. See
[variants and vaccines](variants.md).

## Coordinates

| Where | Convention | Example |
| --- | --- | --- |
| A variant's allele | One-based, VCF-style REF and ALT | chr14, 101980529, G, A |
| Regions in Python | Zero-based, half-open | `Region("chr14", 101980528, 101980530, "GRCh38")` |
| Regions on the command line | One-based, inclusive, like SAMtools | `chr14:101980529-101980530` |
| SV breakends | Zero-based, between bases | chr9 at 22510295 |

The two regions above are the same two bases. Osteosarc never converts coordinates
between genome builds: it checks each BAM's build and refuses regions from another,
so GRCh37 BAMs need GRCh37 coordinates. See [reads](reads.md).

## Corrections

Osteosarc fixes 35 known problems in the website's data by default: misplaced or
missing alleles, mislabeled samples, wrong sites and dates, read counts measured in
the wrong place, and accession typos. Each fix is checked when a snapshot opens; one
the website has since made is marked fixed upstream, and one whose data changed
unexpectedly is skipped as stale. `osteosarc --no-corrections` shows the published
values. See [corrections](corrections.md).

## Missing is not zero

Osteosarc keeps values as published and never fills a gap with zero. A count a
correction cleared is empty, not 0; NA, 0 and an empty cell stay different; an
ELISPOT that wasn't run isn't a negative result; a value below the limit of
quantification is the lab's reporting limit, not a measurement; and no reads in a
limited fetch doesn't mean there are none.

## Receipts and caching

Every download and read extraction keeps a receipt with checksums and sizes, and
asking again reuses the cached result, even offline. A cached file that no longer
matches its receipt raises an IntegrityError; a file that isn't cached, when you're
offline, raises an OfflineError. See [snapshots and cache](snapshots.md).
