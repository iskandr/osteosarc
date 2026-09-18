# Corrections and source drift

Every hand-written interpretation of the osteosarc sources lives in one file,
[`osteosarc/curation.py`](https://github.com/iskandr/osteosarc/blob/main/osteosarc/curation.py).
It has two layers:

* **Vocabulary (always on).** This layer maps source spellings to query names,
  such as `"Boston Gene"` to `BostonGene`, a `CITE` viewer label to `cite-seq`,
  or `"Normal"` to `blood`. The mapping is lossless: claims keep their original
  labels, and metadata rows keep their original values.
* **Corrections (optional).** These target specific source records known to be
  wrong or misleading. Each correction records:
    * the records it touches;
    * the published values it was written against;
    * what it changes, or that it only flags the records;
    * its evidence.

```python
from osteosarc import Dataset

data = Dataset.open("baseline")                          # verified corrections applied
raw = Dataset.open("baseline", corrections=False)        # published sources, unchanged
for row in data.corrections:
    print(row["status"], row["id"], row["summary"][:70])
```

The command line takes `--no-corrections` before any subcommand:

```sh
osteosarc curation baseline            # every correction's status, plus unrecognized labels
osteosarc --no-corrections variants baseline --gene MAP2
```

## What a correction does to your data

Nothing is corrected silently. A correction either **edits** records or **flags**
them, and every touched object says so:

| Object | Where the correction IDs appear |
| --- | --- |
| `Variant` | `variant.annotations["corrections"]` |
| `Asset` | `asset.metadata["corrections"]` |
| count rows (`data.vafs`, `data.table("vafs")`) | an extra final `corrections` column |
| specimens (`data.specimens`) | the `corrections` field |
| timeline events | `event.corrections` (and a note in `timeline.listing()`) |
| `Variants.source` / Varcode metadata | `source["corrections"]` lists every applied ID |

```python
map2 = data.variants()["MAP2-chr2-209694768"]
print(map2.allele, map2.annotations["corrections"])
print(raw.variants()["MAP2-chr2-209694768"].allele)       # the published 22-bp deletion
tempus = data.vafs.select(bam_file="TL-24-5GQLV9WSXQ_T.sorted.bam")
print(tempus.rows[0]["total_reads"] == "", tempus.rows[0]["corrections"])
```

When a correction replaces an allele, that variant's count rows are cleared to
`""` (unmeasured), because their counts were measured for the wrong allele or
locus. Clearing is never presented as a zero. Variant IDs remain the website's
identifiers, so `CABLES1-chr18-23135500` keeps its ID while its allele moves to
chr18:23135764.

## Anticipating upstream changes

The website is rebuilt from sheets and scripts, so the sources will change.
Every time a Dataset loads, each correction is re-checked against its sources:

| Status | Meaning | Applied? |
| --- | --- | --- |
| `applied` | The published values still match what the correction was written against | yes |
| `fixed_upstream` | The source already contains the corrected values | nothing to do |
| `stale` | The source changed some other way, or the record is gone | **no**, and a `CurationWarning` is issued |
| `unavailable` | The snapshot predates a source the correction needs | no |
| `disabled` | `corrections=False` | no |

A correction applies all or nothing: if any of its changes is stale, none is
made. Some corrections carry a *witness*: a change that only checks a
known-wrong published value. For example, `tempus-grch37-counts` requires
PDZRN4's published 7/7. If the website recomputes those counts, the witness
fails and the correction goes stale instead of blanking corrected data.

Vocabulary drift is reported separately. New assay, tissue, provider,
timepoint, pipeline, or vaccine labels appear in `data.unrecognized`, so that
filters are not silently wrong:

```python
print(len(data.unrecognized))   # 0 for the 2026-09-18 sources
```

A routine for each new sync:

```sh
osteosarc sync 2026-10-01
osteosarc curation 2026-10-01 --strict   # exit 1 on stale corrections or unrecognized labels
```

The `drift` workflow in `.github/workflows/drift.yml` runs this check weekly and
on demand. When it fails, run `osteosarc curation <snapshot>` and read each stale
correction's `changes` entries. Each entry shows the source, the record
selector, the state (`missing`, `changed`, `pending`, `already_correct`), and
the fields that differ. Then either delete the correction (fixed upstream) or
update its `expect` values.

## Your own corrections

Corrections are plain data. Add to the built-in list, remove from it, or replace it:

```python
from osteosarc import CORRECTIONS, Change, Correction, glob

mine = Correction(
    "my-lab-note", "Flag every CITE-seq BAM for a local QC review.",
    (Change("bucket", {"key": glob("kamil/blood/output/*CITE*/outs/*.bam")}),),
    evidence=("internal QC log 2026-09",), verified="2026-09-18")
data = Dataset.open("baseline", corrections=[*CORRECTIONS, mine])
without_counts = [c for c in CORRECTIONS if c.id != "tempus-grch37-counts"]
data = Dataset.open("baseline", corrections=without_counts)
```

`Change(source, match, expect={}, set={})`:

* `source` is one of `vafs`, `bam_metadata`, `source_variants`,
  `vaccine_overlap`, `variant_index`, `bams`, `bucket`, or one of the timeline
  sources (`events`, `events_sheet`, `mrd`, `specimens`, `timepoint_summary`,
  `fastqs`, `flow`, `imaging`, `pathology`, `labs`, `cytometry`).
* `match` selects records by exact field values or `glob(...)`, and always
  matches against the *published* records.
* Dotted names reach nested fields, for example `"detection.pVACtools 2025"`.
* An empty `set` flags the records without editing them.

## Built-in corrections (sources of 2026-09-18)

All 29 apply to a snapshot of the 2026-09-18 sources. Each summary and its
evidence URLs are in the registry and in `data.corrections`.

### Read counts and alleles

| ID | Action | What |
| --- | --- | --- |
| `tempus-grch37-counts` | edit | The website counted the GRCh37 Tempus WES BAM (TL-24-5GQLV9WSXQ) at GRCh38 coordinates. All 200 count rows are cleared. |
| `allele-CABLES1-chr18-23135500` | edit | Literal Tempus allele at chr18:23135764 (T>TGGCGGC); the site had chr18:23135500 and `dup`. |
| `allele-CCDC40-chr17-80058951` | edit | Literal Tempus allele at chr17:80090148; the site had `not_reported`. |
| `allele-DCHS2-chr4-154322488` | edit | Literal Tempus delins at chr4:154323273. |
| `allele-GAPVD1-chr9-125299105` | edit | Literal Tempus delins at chr9:125301980 (TAGTGC>ATTGG). |
| `allele-GOLGA6L2-chr15-23441121` | edit | Literal Tempus 120-bp insertion at chr15:23440197. |
| `allele-MAP2-chr2-209694768` | edit | The curated 22-bp deletion, a vaccine target, is not the observed allele. Tempus and CeGaT report one complex −28 bp event, which the catalogue's own protein sequence matches. |
| `map2-split-representations` | flag | Two other MAP2 entries are pieces of that same event. |
| `muc3a-grch38-placement` | flag | The catalogue locus lies in GRCh38-only sequence; the original call does not place there. |
| `ush2a-transposed-duplicate` | flag | chr1:215560752 is a digit transposition of chr1:215650752. |
| `fam157a-withdrawn-protein` | flag | The annotated protein model has been withdrawn by NCBI. |
| `natera-alleles-unavailable` | flag | No public source gives the COL3A1 or OTUD4 alleles. |
| `transcript-DCHS2` | edit | `NM_1142552` becomes `NM_001142552.1`. |
| `gene-symbol-TRMO` | edit | `TMRO` is a typo for `TRMO`; gene-symbol joins with pVACseq otherwise miss it. |

The five relocations and MAP2 were lifted from the original GRCh37 Tempus
TL-24-ALMY2X4KMV VCFs with Ensembl. Each lift was a single ungapped mapping on
the forward strand, and each GRCh38 REF matched the reference sequence. No
website coordinate falls inside its variant's indel-equivalence span, so
normalization cannot explain the differences. The website positions came from
wrong transcript offsets in the curated sheet. For CABLES1 the offset is exactly
its 264-bp 5′UTR. Each wrong position still lands in an exon of the right gene
with a matching reference base, so simple consistency checks pass. **`ready`
therefore does not mean validated.** It means the sources agree on one literal
allele.

### Samples, pipelines, and files

| ID | Action | What |
| --- | --- | --- |
| `viewer-label-BG009368` | edit | The viewer label still says T0 2022-12; the consolidated metadata re-assigned it to T1 2024-06, with VAF evidence. |
| `viewer-label-SARC0277` | edit | Likewise, from T0 BostonGene to T2 UCLA 2025-01. |
| `provider-IPISRC044-T1-rna` | edit | The oncoanalyser T1 RNA was built from BostonGene FASTQs, not UCLA ones. |
| `pvac-2025-detection` | edit | CDC40, PIP5K1A, SMC5 and TECPR1 are in the 2025-04-27 pVACtools runs. |
| `pvac-header-only-filtered-reports` | flag | Seven filtered reports contain only a header. |
| `pvac-rna-fields-na` | flag | RNA depth, VAF and expression are `NA` in all 131,209 report rows. |
| `pvac-extended-run-no-class-i` | flag | The "MHCI.extended" run has only Class-II reports. |

### Timeline and specimens

| ID | Action | What |
| --- | --- | --- |
| `specimen-T1-site` | edit | T1 was a UCLA biopsy, not a UCSF resection. It spans two UCLA biopsies (2024-06-06 and 2024-06-11). |
| `specimen-T2-date-site` | edit | T2 was a UCLA biopsy on 2025-01-28, not UCSF on 2025-01-06 (the delivery folder date). |
| `specimen-T3-site` | edit | T3 was an MSKCC resection, not a UCSF biopsy. |
| `pbmc-capture-dates` | flag | Four PBMC specimen dates are capture dates; the flow-cytometry draws were 2–4 days earlier. |
| `events-duplicate-rows` | flag | SQ3370 and Trabectedin each appear twice. |
| `tempus-timepoint` | flag | The timeline places Tempus at T0, but all public Tempus data are TL-24 (T1 2024-06) accessions. |
| `apheresis-date` | flag | The apheresis is 2024-05-14 in the timeline and 2024-05-15 in the ELISPOT records. |
| `reyagel-end-date` | flag | The sheet's ReyaGel end date is `7/14` with no year, so the site shows a single day. It is the only malformed date among about 28,000. |

## Audit verification (2026-09-18)

An external audit reported six groups of curation problems. Each claim was
re-derived from public sources: the original VCFs, BAM headers, the counting
code, and Ensembl or NCBI.

| Audit claim | Verdict | Notes |
| --- | --- | --- |
| 1. Tempus BAM counted at the wrong build | **Confirmed** | The BAM header is b37 (`human_g1k_v37`). `pileup-json` only toggles the `chr` prefix. 171 of 172 website SNV rows reproduce exactly at the wrong locus. KMT2D 0/1,851, ZNRF3 0/1,287, AKT2 2/1,217, ATRX 1/586 and PDZRN4 7/7→0/0 all reproduce, using primary non-duplicate reads with MAPQ and BQ ≥20. Four "100% VAF" rows (PDZRN4, FHL3, HIC2, BRAT1) are reference reads at the wrong locus. Only this one of the 38 counted BAMs is GRCh37. |
| 2. Five wrong coordinates | **Confirmed** | See above. The calls come from TL-24-ALMY2X4KMV, not the BAM in claim 1. |
| 3. DCHS2 accession, FAM157A, USH2A | **Confirmed** | NM_001145248.1 is suppressed; FAM157A is now lncRNA NR_146164.1. The original USH2A Natera record is not public. |
| 3. MAP2 representations unreconciled | **Partly confirmed** | Tempus and CeGaT agree with each other, and with the site's `CT>AG` plus 28-bp deletion taken together. Only the curated 22-bp deletion is inconsistent. |
| 4. BG009368 and SARC0277 metadata conflicts | **Confirmed** | The consolidated metadata's notes column documents deliberate re-assignment, so the viewer labels are stale. |
| 4. Counting metadata labels 9 BAMs ONT, 2 normals tumour | **Confirmed, no impact** | That file only maps paths to labels; the count export uses the consolidated metadata. Not corrected. |
| 5. Missing pVAC flags, empty RNA fields, incomplete outputs | **Confirmed** | The four missing flags are the complete set. The RNA fields were already `NA` in the pVACseq inputs. |
| 6. MUC3A, COL3A1, OTUD4, LENS, Natera and BostonGene reports | **Confirmed** | No LENS files or clinical reports are in the bucket. The Natera BAMs are GRCh37 (hs37d5). |

The review found further issues beyond the audit. Some are corrected above; the
rest are documented here only:

* The detection key "Tempus 2022" names calls that come from a 2024 accession.
* The consolidated metadata leaves four BAMs with no specimen: three organoid
  DRAGEN BAMs and the T2 UCLA blood DRAGEN BAM.
* Four 2026 blood specimens have no BAMs yet.
* ELISPOT `experiments[].date` is the earliest PBMC sample in an experiment,
  not the assay date.
* The legacy `data/treatment_timeline.json` is still served but stale. It is
  not used here.
* The website's MT-ND5 row reads 0/0 for the Tempus BAM because `chrM` never
  resolves to `MT`. It is covered by `tempus-grch37-counts`.
* The curated sheet contains placeholder strings that went through
  reverse-complement code: `dearoper_aon` is "not_reported" reverse-complemented,
  and `pud` is "dup" reversed.

Three changes in the package itself were also needed:

* CITE-seq libraries are now `cite-seq`.
* A file's own data-page row overrides its directory's row.
* Viewer-label dates and providers are now claims. A month and a day inside
  that month are not a conflict.

Before these fixes, 28 of the snapshot's 31 metadata conflicts were
normalization artifacts. With corrections applied, no sample-metadata conflicts
remain. Without corrections, the three genuine source disagreements are visible.
