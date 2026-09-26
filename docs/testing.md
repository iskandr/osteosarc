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

The tests run offline, on small excerpts of the website's files and on BAMs made up
for the tests; tests/data/provenance.json records where each excerpt came from, with
its checksum, and scripts/build_test_fixtures.py remakes the excerpts from a folder of
downloaded website files. CI also runs the OpenVax libraries against osteosarc: Varcode 7.0.0 and
9.3.7, Isovar 1.17.0 and Topiary 5.55.1. Their full protein and ranking tests run in
their own repositories.

## Run the documentation examples

```sh
python scripts/check_docs.py
python scripts/check_docs.py docs/reads.md
```

This runs every example in the README and on each page of the site against the live
website: each page's Python blocks in order, then its osteosarc commands. Install
commands are listed but not run, and a block preceded by
`<!-- docs-check: skip (reason) -->` is skipped. It starts from a fresh cache unless
you pass `--cache DIR`. It needs network access, SAMtools, the OpenVax libraries and the
Ensembl 95 annotation. On an empty cache, check the README first so that a snapshot
exists.

A weekly drift workflow checks that every correction still matches the live website,
then runs all the examples.

## Repeat the live read check

```python
from osteosarc import Dataset, Region

data = Dataset.sync()
source = data.file(
    "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.bam"
)
regions = [Region("chr14", 101980528, 101980530, "GRCh38")]
subset = data.extract_reads(source, regions, timeout=180)
print(subset.receipt["records"])
assert Dataset.open().extract_reads(source, regions).path == subset.path
```

On 2026-09-18 this returned 3,788 reads, and reopened the same result offline.

## Last checked

On a snapshot downloaded 2026-09-25 (UTC):

| Check | Result |
| --- | --- |
| Files | 395,541, including 843 alignments and 323 VCF or BCF files |
| Variants | 181 on the variants page: 179 ready, 1 with a placeholder allele, 1 with none |
| Corrections | 35: 30 applied, 5 already fixed on the site, none stale |
| Unknown labels | None |
| Timeline and samples | 790 events and 21 samples |

The 2026-09-18 snapshot, from before the site renamed five variants and merged two
USH2A entries, still applies all 35 corrections to its 182 entries. Every allele
correction was checked against the reference genome, and none only rewrites an
equivalent form of the published allele. These checks don't download every BAM or test
every scientific claim in the sources; the [MAP2 example](map2.md) checks one correction
against the reads.
