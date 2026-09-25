# Replace existing download helpers

For projects that download osteosarc.com data with their own code. Switch the
downloading first, and change references, allele selection or analysis settings in
a separate step. The [library examples](consumers.md) show the calls.

## Replace one operation at a time

`pip install osteosarc` brings pysam and datacache with it.

| Your code does | Use |
| --- | --- |
| Download website metadata and save checksums | `Dataset.sync()` |
| Reopen pinned metadata | `Dataset.open()`, or `Dataset.open(name_or_id)` for one exact snapshot |
| Parse the variants page and read-count table | `data.variants()` |
| Find files by assay, timepoint, or path | `data.assets.select(...)` |
| Download and verify a whole file | `data.download(asset)` |
| Check a BAM's genome build | `data.inspect_alignment(asset)` |
| Fetch reads in regions, with mates | `data.extract_reads(asset, regions, ...)` |
| Make Varcode variants | `variants.to_varcode(genome=...)` |

## Preserve the analysis

Your project keeps its reference releases, transcripts, SV definitions, read
filters and scoring. For a reference with a custom name, pass
`to_varcode(..., assembly="GRCh38")`.

Before dropping an old read extractor, compare whole records, tags included, and
how many times each appears; matching read counts aren't enough.

`subset_templates` picks reads at random, ignoring alleles and quality, so it won't
reproduce hand-picked test reads. Use [fixture recipes](fixtures.md) for those.

## Review corrections

Osteosarc applies [corrections](curation.md) by default, which changes some
alleles. Open a snapshot with `corrections=False` to compare with results from the
published values, and review any changed results separately.

## Reuse old downloads

Point `OPENVAX_DATA_CACHE` at a shared directory, and add files you already have
with `Cache.import_file(path, original_url, sha256=...)` (see
[snapshots and cache](design.md#import-a-file-you-already-downloaded)). Keep your old
notes if you need the original download dates.
