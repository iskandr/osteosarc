# Validation

Validation uses the same public Dataset API for Python callers and the CLI.
Normal tests have no network dependency. The fixture directory contains small
parsed excerpts from public resources, with original and excerpt SHA256s in
`tests/data/provenance.json`. Indexed BAM fixtures are synthetic and include
exact repeated records, secondary/supplementary/duplicate flags, differing
read groups, barcodes, UMIs, and float-valued tags.

## Offline checks

```sh
ruff check osteosarc tests scripts
python -m pytest -q
python -m build --no-isolation
```

The tests cover snapshot immutability, corruption, missing cache objects,
explicit refresh, concurrent download reuse, filename collisions, failed
download publication, source conflicts, incomplete alleles, ambiguous vaccine
joins, missing versus zero values, complete S3 pagination, assembly and
mitochondrial checks, indexed region unions, real duplicate preservation,
filter effects, offline extraction reuse, and explicit template sampling.

The consumer review additionally exercises custom-named PyEnsembl references
and retained Varcode metadata, native Isovar read evidence, Topiary RSEM/pVAC
imports (including header-only reports), paired-mate recovery outside a panel,
cached header inspection, changed remote identities, and untested vaccine assay
states. The native adapter tests run with Varcode 7.0.0, Isovar 1.17.0, and
Topiary 5.55.1 in a dedicated CI job. No reference downloads or prediction
models are needed. Full consumer scientific pipelines remain their own tests.

Review found and fixed a Linux-specific failure: SAMtools distribution build
flags can contain non-UTF-8 bytes after the version line. Version receipts now
decode only the version line, with a regression test. Parser support is checked
before downloading an unsupported large alignment object.

Documentation builds with `python -m mkdocs build --strict`. PRs validate the
site; merges to main publish it through the GitHub Pages workflow.

## Curation and source drift

`tests/test_curation.py` and `tests/test_timeline.py` exercise every
correction outcome on public excerpts:

* `applied`, `fixed_upstream`, stale because a record changed or disappeared,
  and `unavailable` on older snapshots;
* all-or-nothing application, flags on assets through globs, and the
  `corrections` column;
* `--strict` exit codes and vocabulary drift;
* the registry's completeness (every ID documented, literal replacement alleles).

The fixtures include every record the timeline corrections target, so those
apply in tests exactly as they do live.

Against the full 2026-09-18 sources:

| Check | Result |
| --- | --- |
| Corrections | 29 applied; 0 stale; 0 fixed upstream |
| Unrecognized source labels | 0 |
| Site variants, corrected | 177 ready, 2 nonliteral (USH2A transposition, MUC3A), 3 missing literal alleles |
| Site variants, `corrections=False` | 172 ready, 7 nonliteral, 3 missing |
| Asset sample-metadata conflicts | 0 corrected; 3 genuine source disagreements uncorrected (previously 31, of which 28 were normalization artifacts) |
| Timeline | 787 events from 8 sources (2022-10-06 .. 2026-09-16), 21 specimens |

The corrections were derived from the original public records:

* the Tempus VCFs and the BAM header;
* the site's counting script (`pileup-json`), reproducing 171 of 172 website
  SNV counts at the wrong locus;
* Ensembl liftover and reference sequence;
* NCBI RefSeq status;
* regional reads.

In the T1 tumor WGS, 30 reads support the corrected MAP2 allele and none
support the published deletion (see the [tour](tour.md)).
[Corrections](curation.md) lists the audit verdicts.

The `drift` workflow (`.github/workflows/drift.yml`) runs every Monday and on
demand. It syncs the live sources into a new snapshot and runs
`osteosarc curation --strict`, so stale corrections or new labels fail the job.
It then runs every documentation example.

## Documentation examples

```sh
python scripts/check_docs.py                  # every page, fresh temporary cache
python scripts/check_docs.py docs/tour.md     # one page
```

The script executes each page's Python blocks in order, in one namespace, and
each `osteosarc` shell command. Interactive commands receive `quit`. Install,
clone and development commands are listed but not run. It needs network access
and `samtools`. The Varcode/Isovar recipes need the consumer packages and an
installed Ensembl 95.

On 2026-09-18, in a fresh virtualenv installed with
`pip install -e '.[reads,test,docs]'` plus the three consumer packages, all 94
examples passed: Python blocks and `osteosarc` commands across 13 pages. The 13
install, clone, and development commands were listed but not run.

## Public smoke check performed 2026-09-18

The public website bucket listing was generated **2026-09-11T23:56:05Z**.
It contained 395,536 objects. Reconciliation exposed 395,541 assets including
five named website tables: 843 alignment files, 323 VCF/BCF files, 2,396 raw-read
files, and other catalog objects. These are file counts, not biological samples.

The metadata join retained 182 website entries (172 ready, seven nonliteral,
three missing literal alleles before corrections), 203 total entries including
the count export, and 44 entries with a positive site vaccine count.

The live shared downloader fetched `bams/bams.json` (19,722 bytes,
SHA256 `88518bc3f9b0e9a7e6fffe13f6c5506e8f781a7287962eee20e412259313360e`).
The named table endpoints parsed successfully:

| Resource | Rows |
| --- | ---: |
| `snv_top` | 276 |
| `dna_fusions` | 133 |
| `rna_fusions` | 96 |

The actual regional extractor queried the published original
`rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.bam` at
`Region("chr14", 101980528, 101980530, "GRCh38")`. It produced **3,788 original
records in a 144,646-byte indexed BAM**, inspected the source header, downloaded
the published index, and retained source HTTP identities. Repeating the request
with an offline Dataset returned the same verified artifact.

The local smoke receipt is `.cache/live-smoke.json`; raw downloads and large
metadata snapshots are ignored by Git. The public test excerpts are checked in.

## Reproduce the bounded live check

```python
from osteosarc import Dataset, Region

data = Dataset.sync("smoke")
for name in ("snv_top", "dna_fusions", "rna_fusions"):
    print(name, len(data.table(name)))
source = data.asset(
    "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.bam"
)
regions = [Region("chr14", 101980528, 101980530, "GRCh38")]
subset = data.extract_reads(source, regions, timeout=180)
offline = Dataset.open("smoke")
assert offline.extract_reads(source, regions).path == subset.path
```

Full BAM/FASTQ downloads and exhaustive scientific validation of every source
were not performed. Remote CRAM extraction requires an explicit indexed
reference and is tested locally rather than against every public CRAM.

The regional semantics follow the
[SAMtools view documentation](https://www.htslib.org/doc/samtools-view.html).
The upstream sources are the [data portal](https://osteosarc.com/data/),
[variant index](https://osteosarc.com/variants/),
[vaccine page](https://osteosarc.com/vaccines/), and the
[public site repository](https://gitlab.com/slowkow/osteosarc.com).
