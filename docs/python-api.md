# Python API

Everything here imports from `osteosarc`, and importing never downloads anything.
In a notebook or REPL, the dataset, samples, files, variants, tables and timeline
all show a readable preview, and `data` itself lists what's there.
`osteosarc repl` opens Python with the newest snapshot loaded as `data`.

```python
from osteosarc import Dataset

data = Dataset.open()
data.samples["T1_tumor"]
```

## Open a snapshot

| Call | What it does |
| --- | --- |
| `Dataset.sync()` | Download the website's metadata as a snapshot named by today's UTC date; the same day, reopen it |
| `Dataset.open()` | Reopen the newest snapshot, offline; pass offline=False to allow downloads |
| `Dataset.open(date="2026-09")` | The newest snapshot downloaded in a year, month or day; or give a name |
| `Dataset.snapshots()` | Every saved snapshot, newest first |
| `data.name`, `data.id`, `data.downloaded` | The snapshot's name, checksum ID and download time |
| `data.summary()` | Counts of everything, and what to try next |

A snapshot applies osteosarc's [corrections](corrections.md) unless you open it with
corrections=False. See [Snapshots and cache](snapshots.md).

## Explore

| Call | What it does |
| --- | --- |
| `data.samples` | Every tumor, organoid and blood sample, as a table |
| `data.samples["T1_tumor"]` | One sample: where and when it was collected, and its BAMs and FASTQ folders |
| `data.samples.select(tissue="blood")` | Samples by time point, tissue, assay or platform |
| `sample.files` | A sample's BAMs and the files in its FASTQ folders |
| `data.files` | Every file in the bucket, and the site's tables |
| `data.files.select(kind="alignment", assay="rna-seq")` | Files by kind, format, path prefix, sample, time point, assay, platform, tissue or provider |
| `data.file(key)` | One file, by its key in the bucket, its URL or a table's name |
| `data.variants(gene="MAP2")` | The variants page, with alleles; also by ID, vaccine, pipeline or status |
| `data.variants("all")`, `data.variants("vaccine")` | Every entry, or just the vaccine targets |
| `data.vaccines`, `data.vaccine_peptides("mRNA")` | Vaccine targets with ELISPOT results; the vaccines' peptides |
| `data.timeline` | Every dated event: treatments, procedures, scans, MRD and labs |
| `data.timeline.render()`, `.listing()`, `.around("2025-01-28")` | A chart by month, one line per event, or a week around a date |
| `data.corrections` | Each fix osteosarc makes to the website's data, and whether it applied |

Collections filter with select, and also with where and a function; index them by
ID or key, and turn them into dictionaries with to_records. See
[Samples and files](samples.md), [Variants](variants.md) and [Timeline](timeline.md).

## Get data

| Call | What it does |
| --- | --- |
| `data.extract_reads(file, variants=data.variants(gene="SMC5"), padding=100)` | Stream just those reads into a small indexed BAM; asking again reuses it |
| `data.extract_reads(..., to="tests/data")` | Also save the BAM and its index in a folder, named for the source and the variants |
| `data.download(file, to=".")` | Download a whole file; to also puts it, and its index, in a folder |
| `data.local_path(file)` | The file's downloaded copy, or None, without using the network |
| `data.downloads()` | Everything downloaded and extracted, with local paths |
| `data.inspect_alignment(file)` | A BAM's header and genome build, without its reads |
| `data.parse(file)` | A table or JSON file's contents; values stay as text |
| `data.vafs` | The site's read counts for every variant |
| `data.measurements` | MRD, lab and cytometry values, as published |

Regions are zero-based and half-open, as in `Region("chr14", 101980428, 101980630, "GRCh38")`;
`Region.from_samtools("chr14:101980429-101980630", assembly="GRCh38")` reads the
one-based form. Read filters, mates and split reads are options of extract_reads;
see [Reads](reads.md).

## Test data

| Call | What it does |
| --- | --- |
| `data.make_bundle("dir", variants=[...], files=[...])` | Make a bundle of test reads: a balanced set at each variant (or SV) in each BAM, every record pinned |
| `bundle_file("dir", member)` | One member as a local indexed BAM (or SAM), exported once into the cache and reused offline |
| `list_bundle("dir")` | Each member, with its records and why they were kept |
| `export_bundle("dir", "tests/data", members=[...])` | Write members into a folder as indexed BAMs (or SAM), named after them |
| `check_fixtures("openvax-v1", {"member": "file.bam"}, root="tests/data")` | Compare your own files with a bundle's members; returns those that differ |
| `verify_bundle("dir")` | Check every file, record and index |
| `fetch_bundle("openvax-v1")` | Download a published bundle once, and return its folder (`offline=True` never downloads) |
| `generate_bundle(recipe, "dir", dataset=data)` | Make a bundle from a [recipe](test-data.md#recipes) |
| `validate_recipe(recipe)` | Check a recipe before fetching anything |
| `load_panel("vaccine-loci-v1")`, `load_sv_candidates()` | Shipped lists of targets, and the [SV candidates](sv-candidates.md) |

Each of these takes a bundle's folder or the name of a published bundle, such as
openvax-v1, the reads the OpenVax libraries share.
`osteosarc.shared.published("openvax-v1")` gives its release record, with the
manifest checksum to record alongside your results. See [Test data](test-data.md).

## Corrections

| Call | What it does |
| --- | --- |
| `CORRECTIONS` | The built-in fixes |
| `Correction(...)`, `Change(...)`, `glob(...)` | Write your own fix |
| `Dataset.open(corrections=[*CORRECTIONS, mine])` | Open a snapshot with your fixes added |

A fix that no longer matches the data is skipped with a CurationWarning. See
[Corrections](corrections.md).

## Errors

Osteosarc's errors all derive from OsteosarcError:

- **OfflineError:** something isn't downloaded and the snapshot is offline.
- **IntegrityError:** a checksum or a file's identity doesn't match.
- **CoordinateError:** a region or variant can't be used, for example on another genome build.
- **SchemaError:** a source or recipe isn't in the expected shape.
- **NoSnapshotsError:** there's no snapshot yet; run `Dataset.sync()`.

Unknown sample, file and variant IDs raise KeyError. Osteosarc never falls back to
downloading a whole BAM, converts between genome builds, or ranks peptides.

## Everything else

The lower-level pieces live in their modules: osteosarc.reads (reading regions of
alignments), osteosarc.bundles (bundle files), osteosarc.records (record checksums)
and osteosarc.discovery (listing the bucket live). `help()` on any of them describes
it.
