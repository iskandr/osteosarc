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
| List a sample's BAMs and FASTQs | `data.samples["T1_tumor"].files` |
| Find files by assay, timepoint, or path | `data.files.select(...)` |
| Download and verify a whole file | `data.download(file)`, or `data.download(file, to=DIR)` for a copy under its own name |
| Check whether a file is already downloaded | `data.local_path(file)` |
| Check a BAM's genome build | `data.inspect_alignment(file)` |
| Fetch reads in regions, with mates | `data.extract_reads(file, regions, ...)` |
| Make Varcode variants | `variants.to_varcode(genome=...)` |

## Names changed in 0.9

| Before 0.9 | Now |
| --- | --- |
| `data.assets`, `data.asset(key)`, `Asset`, `Assets` | `data.files`, `data.file(key)`, `File`, `Files` |
| `data.specimens`, `data.describe_samples()` | `data.samples` (`Sample` objects) |
| `data.assets_for_sample("T1_tumor", ...)` | `data.samples["T1_tumor"].files.select(...)` or `data.files.select(sample="T1_tumor", ...)` |
| `data.samples` (source claims) | `data.claims`, with `file_ids` |
| `data.explore()`, `osteosarc explore` | `osteosarc repl` |
| `osteosarc assets`, `specimens`, `curation` | `osteosarc files`, `samples`, `corrections` |
| JSON output from `variants`, `vaccines`, `corrections`, `sync`, `table`, `discover` | Text; add `--json` |

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
