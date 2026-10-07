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
its checksum. scripts/build_test_fixtures.py remakes the excerpts from a folder of the
website's files, downloaded under their own names: variants.html, source-variants.json,
vaccine_overlap.json, bams.json, bucket_listing.json, vafs.tsv, vafs-columns.tsv,
bam-metadata.tsv, data.html, events.json, timeline.csv, mrd.json,
samples-consolidated.tsv, samples.json, fastqs-consolidated.tsv, flow-manifest.json,
dicom-studies.json, pathology-slides.json, lab_results.tsv and cytometry.tsv. CI also runs the OpenVax libraries against osteosarc: Varcode 7.0.0 and
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
and, separately, runs all the examples.

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

On a snapshot downloaded 2026-10-07 (UTC), using the site's bucket inventory
generated 2026-10-01:

| Check | Result |
| --- | --- |
| Files | 397,038 bucket objects; 397,043 catalogue entries including site tables; 913 alignments and 357 VCF or BCF files |
| Variants | 183 on the variants page: 182 ready, 1 with a placeholder allele |
| Corrections | 35: 13 applied, 22 already fixed on the site, none stale |
| Unknown labels | None |
| Timeline and samples | 807 events and 31 samples (24 blood, 5 tumor/fractions, 2 organoids) |

The sample audit checks every registered BAM and FASTQ location, including the
eleven entries represented as filename prefixes or Tempus row aliases rather
than directories. The current tumor dates and collection sites agree with the
timepoint summary. Regression excerpts and original source receipts are in
`tests/data/sample-locations-2026-10-07.json`; older fixtures remain unchanged.
Live S3 listings matched the inventory for all 317 FASTQ locations, and HEAD
requests matched the sizes and modification dates of all 94 sample BAMs. One
published FASTQ-table count differs: `rna-seq/tempus/TL-24-ALMY2X4KMV/RNA` lists
two compressed FASTQs (1,238,829,178 bytes), but also contains `2022_R1.fastq`.
The inventory and sample view retain all three files (4,470,300,364 bytes).

The 2026-09-18 snapshot, from before the site renamed five variants and merged two
USH2A entries, still applies all 35 corrections to its 182 entries. Every allele
correction was checked against the reference genome, and none only rewrites an
equivalent form of the published allele. These checks don't download every BAM or test
every scientific claim in the sources; the [MAP2 example](map2.md) checks one correction
against the reads.
