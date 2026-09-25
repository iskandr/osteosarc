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

The tests run offline, on small excerpts of the website's files and on BAMs made
up for the tests. `tests/data/provenance.json` records where each excerpt came from
and its checksum.

CI also tests the OpenVax libraries against Osteosarc: Varcode 7.0.0 and 9.3.7,
Isovar 1.17.0 and Topiary 5.55.1. Their full protein and ranking tests run in their
own repositories.

## Run the documentation examples

```sh
python scripts/check_docs.py
python scripts/check_docs.py docs/reads.md
```

This runs every example in `README.md` and in each page of the site, against the
live website: each page's Python blocks in order, then its `osteosarc` commands.
Install commands are listed but not run, and a block marked
`<!-- docs-check: skip (reason) -->` is skipped. It uses a fresh temporary cache
unless you pass `--cache DIR`. You need network access, SAMtools, the OpenVax
libraries and the Ensembl 95 annotation. When running single pages on an empty
cache, include `README.md` first so a snapshot exists.

A weekly `drift` workflow checks that every correction still matches the live
website, then runs all the examples.

## Repeat the live read check

```python
from osteosarc import Dataset, Region

data = Dataset.sync()
source = data.asset(
    "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.bam"
)
regions = [Region("chr14", 101980528, 101980530, "GRCh38")]
subset = data.extract_reads(source, regions, timeout=180)
print(subset.receipt["records"])
assert Dataset.open().extract_reads(source, regions).path == subset.path
```

On 2026-09-18 this returned 3,788 reads, and reopened the same result offline.

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

## Source changes checked: 2026-09-21

The current snapshot has 181 entries: 179 ready, one nonliteral and one missing
a literal allele. The site merged two USH2A entries and renamed five variants.
All 32 correction rules pass: 29 apply and three are fixed upstream. There are
no unrecognized labels. The historical snapshot above still passes with its
original IDs and 182 entries. See [source corrections](curation.md) for how the
updated entries and their counts are handled.

A fresh snapshot on 2026-09-24 matched these results: 395,541 assets, 181
entries (179 ready), 29 applied and three fixed-upstream corrections, no
unrecognized labels, and 21 specimens. The timeline had 788 events.
