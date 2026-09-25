# Python API

Everything below imports from `osteosarc`. Importing never downloads anything or
needs the OpenVax libraries. The guides show each call in context; this page lists
them by task.

[Snapshots and cache](#snapshots-and-cache) ·
[Files and tables](#files-and-tables) ·
[Variants and vaccines](#variants-and-vaccines) ·
[Reads](#reads) ·
[Corrections, timelines and specimens](#corrections-timelines-and-specimens) ·
[Fixtures and bundles](#fixtures-and-bundles) ·
[Errors and limitations](#errors-and-limitations)

## Snapshots and cache

| Call | Returns / behavior |
| --- | --- |
| `Dataset.sync(name=None, cache=None, refresh=False, sources=None, corrections=True)` | Save the website's metadata as a snapshot named by UTC date (or `name`); reopens today's or an existing named one unless `refresh=True` |
| `Dataset.open(name=None, date=None, cache=None, offline=True, corrections=True)` | Verify and reopen the most recent snapshot, one with an exact name or ID prefix, or the newest downloaded on a UTC `date` (`2026-09-24`, `2026-09`, `2026`); `corrections=False` or a list of `Correction`s |
| `Dataset.snapshots(cache=None)` | Table of saved snapshots, newest download first: `name`, `downloaded`, `created`, `id` |
| `data.id`, `data.name`, `data.downloaded`, `data.receipts()` | Snapshot content identity, name, latest source download time, and source receipts |
| `data.source_path(name)` | Verified local metadata path |
| `Cache(root=None, offline=False, timeout=600)` | Shared local object cache |
| `cache.fetch(url, refresh=False, sha256=None, md5=None, size=None, max_bytes=None)` | Download and check a file; returns a `Receipt` |
| `cache.path(receipt)` | Check a receipt and return the file's local path |
| `cache.import_file(path, url, sha256=None, md5=None, size=None)` | Add a file you already have, without downloading it |
| `list_bucket(cache, prefix, refresh=False, max_pages=1000)` | List everything in the bucket under a prefix; `refresh=True` for the current contents |
| `digest(path, algorithm="sha256")` | A file's checksum |
| `SNAPSHOT_SOURCES`, `TABLE_SOURCES`, `TIMELINE_SOURCES` | The URLs a snapshot downloads |

See [Snapshots and cache](design.md).

## Files and tables

| Call / property | Returns / behavior |
| --- | --- |
| `data.assets`, `data.samples`, `data.timepoints` | Every file; what the site says about each file's sample; timepoint dates |
| `data.describe_samples(timepoint=None, tissue=None, assay=None, platform=None, width=None)` | Readable sample and sequencing overview, with sequencing shown as filter names |
| `data.assets_for_sample(sample_id, **filters)` | A sample's BAMs and the files in its FASTQ folders |
| `data.asset(key_or_id)` | One file, by key, URL, ID or table name |
| `data.download(asset, refresh=False, verify_size=True)` | Download one whole file and return its path |
| `data.table(asset)`, `data.parse(asset)` | Parse a CSV/TSV table, or a supported JSON/FASTA resource |
| `data.open_variants(asset)` | Open a VCF or BCF with pysam |
| `assets.select(kind=..., format=..., prefix=..., contains=...)` | Filter files by type or path |
| `assets.select(timepoint=..., assay=..., platform=..., tissue=..., provider=..., library=...)` | Filter files by sample and sequencing |
| `asset.values(field)`, `asset.resolved(field)`, `asset.conflicts` | What the site says about a file, and where it disagrees |
| `collection.where(predicate)`, `collection[:n]`, `collection.to_records()` | Filter with a function, slice, or convert to dictionaries |
| `table.select(**fields)`, `table.where(predicate)` | Filter rows; values stay as text |
| `table.rows`, `table.columns`, `table.source`, `table.to_dataframe()` | The rows, columns, source, or a pandas DataFrame |
| `table.diagnostics` | Rows with too few or too many fields: their `row`, `line` and `fields` |
| `parse_table(text, delimiter="\t", strict=True)` | Parse CSV or TSV text; `strict=False` keeps broken rows |
| `parse_file(path, format=None)` | Parse a local JSON, CSV, TSV or FASTA file |

Asset selections accept `include_conflicts` and `include_inferred` opt-ins. An
assay, platform or tissue that no file uses raises `ValueError` listing the valid
values; a registry label such as `scRNA_ONT` names the filters to use instead.
The named tables are `vafs`, `vaf_columns`, `snv_top`, `dna_fusions`, and
`rna_fusions`. Other tables use exact asset keys or Asset instances.
See [Find samples and files](explore.md).

## Variants and vaccines

| Call / property | Returns / behavior |
| --- | --- |
| `data.variants(set="site", **filters)` | Variants: `site`, `all` or `vaccine` |
| `variants.select(gene=..., ids=..., vaccine=..., pipeline=..., status=...)` | Filter variants |
| `variant.allele`, `variant.region(padding=0)` | The allele, and the region it covers; error if not `ready` |
| `variants.regions(padding=0)` | The regions of several variants |
| `variants.to_varcode(genome=..., assembly=None, ...)` | Varcode variants, with the original entries in their metadata; see [naming options](consumers.md#varcode) |
| `data.annotations`, `data.vafs` | The site's variant records and read counts |
| `data.vaccines`, `data.vaccine_names`, `data.vaccine_peptides(vaccine=None)` | Vaccine rows, vaccine names, and the peptides with their ELISPOT experiments |
| `data.pipeline_names` | The pipelines that detected variants |
| `parse_variants(index, vafs)` | Build variants from the variants page and read-count table; broken rows go in `annotations["parse_errors"]` |
| `parse_variant_index(html)` | Read the variants page |

Variant filters additionally include `vaccinated`, `on_site`, and
`vaccine_source="overlap"` or `"source_variants"`.
See [Select variants](variants.md).

## Reads

```text
Region(contig, start, end, assembly, reference_length=None)
ReadFilter(min_mapq=0, exclude_flags=0, require_flags=0,
           barcodes=(), barcode_tag="CB", query_names=())
RecoveryPolicy(mates=True, supplementary=True, max_rounds=4,
               max_intervals=128, max_bases=1_000_000, max_records=100_000,
               on_timeout="fail")
```

| Call | Returns / behavior |
| --- | --- |
| `data.inspect_alignment(asset)` | `AlignmentInfo`: original `.header`, `.assembly`, `.path`, `.receipt` |
| `data.extract_reads(asset, regions=None, variants=None, padding=0, **options)` | Supply regions or selected variants; returns `ReadSubset` with `.path`, `.index_path`, `.receipt`, `.open()` |
| `inspect_alignment(source, cache=None, snapshot_id=None, timeout=600)` | Same inspection for local files or HTTP(S) URLs |
| `extract_reads(source, regions, cache=None, index=None, filters=None, reference=None, fetch_pairs=False, snapshot_id=None, timeout=600, recovery=None, max_records=None)` | Explicit indexed region union |
| `recover_reads(source, regions, policy=None, cache=None, **options)` | Bounded mate/SA partner recovery; the same as `extract_reads(..., recovery=policy)` |
| `subset_templates(source, count, cache=None, seed="0")` | Deterministic local fixture sampling |
| `Region.from_samtools(text, assembly=..., reference_length=None)` | Convert a one-based inclusive `contig:start-end` region |
| `assembly_from_header(header)`, `resolve_regions(regions, header)` | Assembly evidence and validated/merged intervals |
| `normalize_assembly(value)` | Canonical assembly family name, such as `GRCh38` |

`data.inspect_alignment` accepts `timeout`; `data.extract_reads` supplies
its own cache and snapshot ID and accepts the remaining extraction options.
`ReadSubset.open()` is a pysam context manager. `max_records` stops acquisition
on overflow and discards the partial output. See [Extract reads](reads.md).

## Corrections, timelines and specimens

| Call / property | Returns / behavior |
| --- | --- |
| `data.corrections` | Table: each correction's `status` (`applied`, `fixed_upstream`, `stale`, `unavailable`, `disabled`), `changes`, `summary`, `evidence` |
| `data.unrecognized` | Table of source labels outside the vocabulary (upstream drift) |
| `CORRECTIONS`, `Correction(id, summary, changes, evidence, verified, alternatives=())`, `Change(source, match, expect, set, absent=False)`, `glob(pattern)` | The built-in registry and its building blocks ([details](curation.md)) |
| `CurationWarning` | Warning emitted when a correction no longer matches its source |
| `variant.annotations["corrections"]`, `["count_corrections"]` | Corrections to the variant itself, and to only some of its count rows |
| `data.timeline` | `Timeline` of `Event`s from every dated source |
| `timeline.select(lane=, category=, kind=, source=, track=, timepoint=, contains=, since=, until=)` | Filtered Timeline; dates are `YYYY`, `YYYY-MM`, or `YYYY-MM-DD` |
| `timeline.around(date, days=7)`, `timeline.lanes()` | Neighborhood of a date; ordered lane names |
| `timeline.render(width=None, since=None, until=None, legend=True)`, `timeline.listing()` | ASCII chart; one line per event |
| `event.date`, `.end`, `.precision`, `.open_end`, `.timepoint`, `.value`, `.source`, `.corrections`, `.details` | Published precision, correction IDs, and the original record |
| `timeline.source["undated"]` | Source rows left off the timeline because their dates could not be read |
| `data.specimens` | Table: registry rows with `assets`, `fastq_folders`, `disagreements`, `corrections` |
| `data.measurements` | Table: MRD, lab, and cytometry values with raw strings and `kind` |
| `osteosarc.explore.Explorer(data)` | Interactive shell (`osteosarc explore`) |
| `osteosarc.explore.specimen_view`, `specimens_view`, `assets_view`, `variants_view`, `corrections_view`, `summary_view` | The shell's views as strings |

A snapshot made before the timeline sources existed raises `SchemaError` from
`timeline`, `specimens`, and `measurements`; everything else works.
See [Source corrections](curation.md) and [Browse the timeline](timeline.md).

## Fixtures and bundles

| Call | Returns / behavior |
| --- | --- |
| `validate_recipe(recipe)` | Validated copy of a v1 recipe; raises `SchemaError` |
| `select_fixtures(recipe, sources)`, `data.select_fixtures(recipe, sources)` | `FixtureSelection` from local BAMs or `ReadSubset`s; `.manifest` holds membership and reasons |
| `select_fixture_records(records, policy, regions=(), context_regions=())` | Execute one policy over records: counts, reasons and status |
| `load_panel(name)` | Shipped target panel: `vaccine-loci-v1`, `sv-regressions-v1` or `sv-interest-v1` |
| `load_sv_interest()` | The SV interest catalogue with targets and source pins ([details](sv-interest.md)) |
| `generate_bundle(recipe, destination, sources=None, cache=None, **options)`, `data.generate_bundle(recipe, destination, **options)` | Acquire, select and publish a bundle; returns its manifest |
| `pack_bundle(selection, destination, header_policy="full", size_budget=64 MiB)` | Publish a bundle from an existing selection |
| `verify_bundle(directory, sha256=None)`, `list_bundle(directory)` | Offline verification; member status, counts and reasons |
| `export_bundle(directory, destination, members=None, format="bam")` | Named indexed BAM, SAM or SAM.gz exports |
| `compact_header(header, records)` | Header reduced to assembly, read-group and proven program lineage |
| `record_multiset(path)`, `bam_record_digests(path)`, `read_records(path)` | Lossless `RECORD_ENCODING` (`bam-record-v1`) identities of a local BAM's records |
| `safe_path(root, name)`, `verify_digest(path, expected)`, `verify_gzip_digests(...)`, `verify_manifest_files(root, files)` | Path and checksum helpers for consumer fixture formats |

The historical consumer formats live in their own modules:
`osteosarc.legacy_fixtures` (Isovar exact-SAM fixtures), `osteosarc.regional_corpus`
(Topiary regional corpora), `osteosarc.cohort_bundle` (Vaxrank cohorts) and
`osteosarc.fixture_assets` (Topiary asset manifests). See
[Read fixtures and bundles](fixtures.md) and [OpenVax fixture adoption](fixture-migration.md).

## Errors and limitations

Osteosarc's own errors derive from `OsteosarcError`. `OfflineError` means requested
bytes have not been acquired. `IntegrityError` means receipt, checksum, or source
identity disagrees. `SchemaError` signals an incompatible source schema or recipe.
`CoordinateError` rejects unresolved or incompatible coordinate requests. The last
three are also `ValueError`s. Unknown keys, samples and variant IDs raise
`KeyError`, other invalid arguments raise `ValueError`, and reusing a snapshot name
with `refresh=True` raises `FileExistsError`.
Failed external commands raise `subprocess.CalledProcessError` with captured
diagnostic output; callers can record it in acquisition manifests.

No implicit full-BAM fallback, liftover, biological sample deduplication,
reference installation, variant-effect calculation, or peptide ranking occurs.
See [snapshots and cache](design.md) for storage and refresh behavior.
