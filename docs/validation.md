# Testing

## Run the local checks

From a checkout:

```sh
python -m pip install -e '.[test,docs]'
ruff check osteosarc tests scripts
python -m pytest -q
python -m build
python -m mkdocs build --strict
```

Tests run offline with public metadata excerpts and synthetic indexed BAMs.
They cover cache integrity, source conflicts, corrections, coordinate checks,
record and tag preservation, filtering, and offline reuse. Excerpt checksums
are recorded in `tests/data/provenance.json`.

The consumer CI job tests Varcode 7.0.0 and 9.3.7 with Isovar 1.17.0 and
Topiary 5.55.1. It checks native variant conversion, mitochondrial and custom
contig annotation, source-name preservation, Isovar read evidence, Topiary report
loading, and the paired-read acquisition pattern used by Vaxrank. Full protein
reconstruction and vaccine ranking remain downstream integration tests.

## Run the documentation examples

```sh
python scripts/check_docs.py
python scripts/check_docs.py docs/reads.md
```

The script executes Python blocks in page order and runs `osteosarc` commands.
Install and development commands are listed without execution. It uses a
temporary cache by default; `--cache DIR` reuses an existing one.

Live examples need network access and SAMtools. The library examples also need
the consumer packages and an indexed Ensembl 95 human reference. Start with
`README.md` or `docs/index.md` to create the `baseline` snapshot when running
selected pages against an empty cache.

The weekly `drift` workflow checks live sources for stale corrections and
unrecognized labels, then executes the documentation examples.

## Repeat the live read check

```python
from osteosarc import Dataset, Region

data = Dataset.sync("smoke")
source = data.asset(
    "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.bam"
)
regions = [Region("chr14", 101980528, 101980530, "GRCh38")]
subset = data.extract_reads(source, regions, timeout=180)
print(subset.receipt["records"])
assert Dataset.open("smoke").extract_reads(source, regions).path == subset.path
```

On 2026-09-18, this returned 3,788 records in a 144,646-byte indexed BAM.
The source header and index were checked, and the same result reopened offline.

## Recorded checks: 2026-09-18

| Check | Result |
| --- | --- |
| Catalog | 395,541 assets, including 843 alignments and 323 VCF/BCF files |
| Site variants after corrections | 179 ready, 2 nonliteral, 1 missing literal allele (rechecked 2026-09-20) |
| Corrections | 32 applied, none stale (rechecked 2026-09-20) |
| Unrecognized source labels | 0 |
| Sample metadata conflicts | 0 corrected; 3 with corrections disabled |
| Timeline | 787 events, 21 specimens |

See the [MAP2 example](tour.md) for a read-level check of an allele correction.
These checks did not download every BAM/FASTQ or validate every source's
scientific claims. CRAM extraction is tested with local reference files.
