# API reference

Import public objects from `osteosarc`. Imports perform no downloads and do not
require the consumer packages. Install optional libraries only for the features you use.

## Dataset and cache

| Call | Returns / behavior |
| --- | --- |
| `Dataset.sync(name, cache=None, refresh=False, sources=None, corrections=True)` | Create a named metadata snapshot; repeated unchanged names reopen it |
| `Dataset.open(name, cache=None, offline=True, corrections=True)` | Verify and reopen a saved snapshot; `corrections=False` or a list of `Correction`s |
| `data.id`, `data.receipts()` | Snapshot identity and source receipts |
| `data.source_path(name)` | Verified local metadata path |
| `data.assets`, `data.samples`, `data.timepoints` | All assets, source-attributed sample claims, published dates |
| `data.describe_samples(timepoint=None, tissue=None, width=None)` | Readable sample and sequencing overview |
| `data.assets_for_sample(sample_id, **filters)` | Registry-linked alignments and files in the sample's FASTQ folders |
| `data.asset(key_or_id)` | Resolve exactly one key, URL, asset ID, or named resource |
| `data.download(asset)` | Fetch one complete object, bind bytes to the snapshot, return Path |
| `data.table(asset)`, `data.parse(asset)` | Parse a table or a supported JSON/FASTA resource |
| `data.open_variants(asset)` | Context-managed native pysam VCF/BCF reader |
| `Cache(root=None, offline=False, timeout=600)` | Shared local object cache |
| `cache.fetch(url, refresh=False, sha256=None, md5=None, size=None, max_bytes=None)` | Verified download receipt |
| `cache.path(receipt)` | Verify receipt and return its local Path |
| `cache.import_file(path, url, sha256=None, md5=None, size=None)` | Adopt existing bytes without downloading |
| `list_bucket(cache, prefix, refresh=False)` | Separate complete paginated public S3 inventory with receipts |

The named tables are `vafs`, `vaf_columns`, `snv_top`, `dna_fusions`, and
`rna_fusions`. Other tables use exact asset keys or Asset instances.

## Corrections, timelines, and specimens

| Call / property | Returns / behavior |
| --- | --- |
| `data.corrections` | Table: each correction's `status` (`applied`, `fixed_upstream`, `stale`, `unavailable`, `disabled`), `changes`, `summary`, `evidence` |
| `data.unrecognized` | Table of source labels outside the vocabulary (upstream drift) |
| `CORRECTIONS`, `Correction(id, summary, changes, evidence, verified)`, `Change(source, match, expect, set)`, `glob(pattern)` | The built-in registry and its building blocks ([details](curation.md)) |
| `data.timeline` | `Timeline` of `Event`s from every dated source |
| `timeline.select(lane=, category=, kind=, source=, track=, timepoint=, contains=, since=, until=)` | Filtered Timeline; dates are `YYYY`, `YYYY-MM`, or `YYYY-MM-DD` |
| `timeline.around(date, days=7)`, `timeline.lanes()` | Neighbourhood of a date; ordered lane names |
| `timeline.render(width=None, since=None, until=None)`, `timeline.listing()` | ASCII chart; one line per event |
| `event.date`, `.end`, `.precision`, `.open_end`, `.timepoint`, `.value`, `.source`, `.corrections`, `.details` | Published precision, correction IDs, and the original record |
| `timeline.source["undated"]` | Source rows left off the timeline because their dates could not be read |
| `data.specimens` | Table: registry rows with `assets`, `fastq_folders`, `disagreements`, `corrections` |
| `variant.annotations["corrections"]`, `["count_corrections"]` | Corrections to the variant itself, and to only some of its count rows |
| `data.measurements` | Table: MRD, lab, and cytometry values with raw strings and `kind` |
| `osteosarc.explore.Explorer(data)` | Interactive shell (`osteosarc explore`) |
| `osteosarc.explore.specimen_view`, `specimens_view`, `assets_view`, `variants_view`, `corrections_view`, `summary_view` | The shell's views as strings |

A snapshot made before the timeline sources existed raises `SchemaError` from
`timeline`, `specimens`, and `measurements`; everything else works.

## Selection and source interpretation

| Call / property | Meaning |
| --- | --- |
| `assets.select(kind=..., format=..., prefix=..., contains=...)` | Select file products |
| `assets.select(timepoint=..., assay=..., platform=..., tissue=..., provider=..., library=...)` | Match unambiguous published metadata by default |
| `asset.values(field)`, `asset.resolved(field)`, `asset.conflicts` | Inspect metadata assertions |
| `data.variants(set="site", **filters)` | Select `site`, `all`, or `vaccine` entries |
| `variants.select(gene=..., ids=..., vaccine=..., pipeline=..., status=...)` | Compose exact variant filters |
| `variant.allele`, `variant.region(padding=0)` | Unique literal allele / anchored reference span; unresolved entries raise |
| `variants.regions(padding=0)` | Explicit Region tuple for selected alleles |
| `variants.to_varcode(genome=..., assembly=None, ...)` | Native collection with original entries in metadata; [contig naming options](consumers.md#varcode) |
| `data.annotations`, `data.vafs` | Original source annotations and count rows |
| `data.vaccines`, `data.vaccine_names`, `data.vaccine_peptides(vaccine=None)` | Overlap rows, vaccine names, exact published peptides/experiments |
| `data.pipeline_names` | Available source-JSON pipeline labels |
| `collection.where(predicate)`, `collection[:n]`, `collection.to_records()` | Predicate filtering, slicing, dictionaries |
| `table.select(**fields)`, `table.where(predicate)` | Original-row filtering; no type coercion |
| `table.rows`, `table.columns`, `table.source`, `table.to_dataframe()` | Raw table values and optional pandas conversion |
| `parse_variants(index, vafs)` | Join index rows/HTML to VAF TSV text or a Table; preserve malformed entries with `annotations["parse_errors"]` |
| `parse_table(text, strict=True)` | Parse CSV/TSV; opt into `strict=False` to retain ragged rows |
| `table.diagnostics` | Ragged-row diagnostics with zero-based `row`, source `line`, and original `fields`; filtering keeps the associated diagnostics |

Asset selections accept `include_conflicts` and `include_inferred` opt-ins.
Variant filters additionally include `vaccinated`, `on_site`, and
`vaccine_source="overlap"` or `"source_variants"`.

## Alignments and fixtures

```text
Region(contig, start, end, assembly, reference_length=None)
ReadFilter(min_mapq=0, exclude_flags=0, require_flags=0,
           barcodes=(), barcode_tag="CB")
```

| Call | Returns / behavior |
| --- | --- |
| `data.inspect_alignment(asset)` | `AlignmentInfo`: original `.header`, `.assembly`, `.path`, `.receipt` |
| `data.extract_reads(asset, regions=None, variants=None, padding=0, **options)` | Supply regions or selected variants; returns `ReadSubset` with `.path`, `.index_path`, `.receipt`, `.open()` |
| `inspect_alignment(source, cache=None, snapshot_id=None, timeout=600)` | Same inspection for local files or HTTP(S) URLs |
| `extract_reads(source, regions, cache=None, index=None, filters=None, reference=None, fetch_pairs=False, snapshot_id=None, timeout=600)` | Explicit indexed region union |
| `subset_templates(source, count, cache=None, seed="0")` | Deterministic local fixture sampling |
| `assembly_from_header(header)`, `resolve_regions(regions, header)` | Assembly evidence and validated/merged intervals |

`data.inspect_alignment` accepts `timeout`; `data.extract_reads` supplies
its own cache and snapshot ID and accepts the remaining extraction options.
`ReadSubset.open()` is a pysam context manager.

## Errors and limitations

`OfflineError` means requested bytes have not been acquired. `IntegrityError`
means receipt, checksum, or source identity disagrees. `SchemaError` signals an
incompatible source schema. `CoordinateError` rejects unresolved or incompatible
coordinate requests. Failed external commands raise `subprocess.CalledProcessError`
with captured diagnostic output; callers can record it in acquisition manifests.

No implicit full-BAM fallback, liftover, biological sample deduplication,
reference installation, variant-effect calculation, or peptide ranking occurs.
See [snapshots and cache](design.md) for storage and refresh behavior.
