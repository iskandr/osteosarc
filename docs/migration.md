# Replace existing download helpers

Use the [library examples](consumers.md) for the API calls. Migrate acquisition
separately from changes to reference releases, allele selection, or analysis
settings.

## Replace one operation at a time

Install `osteosarc` normally: pysam and datacache are standard dependencies.
The older `[reads]` extra remains accepted for compatibility. The CLI now
prints a sample overview by default; use `samples --json` for its original
source-attributed claims.

| Existing code | Replacement |
| --- | --- |
| Download website metadata and save checksums | `Dataset.sync(name)` |
| Reopen pinned metadata | `Dataset.open(name)` |
| Parse the variant page and join count-export alleles | `data.variants()` |
| Find files by assay, timepoint, or path | `data.assets.select(...)` |
| Download and verify a whole file | `data.download(asset)` |
| Survey an alignment's reference | `data.inspect_alignment(asset)` |
| Fetch indexed regions and optional paired mates | `data.extract_reads(asset, regions, ...)` |
| Construct native Varcode alleles | `variants.to_varcode(genome=...)` |

## Preserve the analysis

Keep reference releases, transcript selection, structural-variant definitions,
read filters, and scoring settings in the downstream project. Declare the
assembly for custom-named references with `to_varcode(..., assembly="GRCh38")`.

Before retiring an old extractor, compare complete SAM records, including tags
and the number of occurrences of each record. Equal read counts alone do not
show that the records match.

`subset_templates` samples without using allele support or quality. It does not
reproduce Isovar's alternate-read enrichment or another project's fixture
selection policy.

## Review corrections

Osteosarc applies [source corrections](curation.md) by default, including the
MAP2 allele change and five relocated Tempus alleles. Open a snapshot with
`corrections=False` when comparing against original published inputs. Review
changes to biological expectations separately from the acquisition migration.

## Reuse old downloads

Use the shared `OPENVAX_DATA_CACHE` directory where possible. Import existing
files with `Cache.import_file(path, original_url, sha256=...)`; see
[snapshots and cache](design.md#import-a-file-you-already-downloaded).
Keep the old manifests for their acquisition dates and provenance.
