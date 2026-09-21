# Source corrections

Osteosarc applies documented corrections to source records by default. You can
inspect each change, disable corrections, or supply your own.

## Compare corrected and published values

```python
from osteosarc import Dataset

data = Dataset.open("baseline")
raw = Dataset.open("baseline", corrections=False)
map2 = data.variants()["MAP2-chr2-209694768"]
print("Corrected:", map2.allele)
print("Published:", raw.variants()[map2.id].allele)
print(map2.annotations["corrections"])
```

Allele corrections clear affected count rows to `""` (unmeasured), because the
original counts describe a different allele or locus. Website IDs are retained
even if coordinates change. See the [MAP2 example](tour.md) to check a
correction against reads.

## Inspect the changes and evidence

```python
for row in data.corrections:
    print(row["id"], row["status"], row["summary"], row["evidence"])
```

Each correction edits or flags records. The affected objects carry its ID:

| Object | Correction IDs |
| --- | --- |
| Variant | `variant.annotations["corrections"]` |
| Some of a variant's count rows | `variant.annotations["count_corrections"]` |
| Asset | `asset.metadata["corrections"]` |
| Count, vaccine, annotation, measurement, and specimen rows | `corrections` field |
| Timeline event | `event.corrections` |
| Variant selection / Varcode metadata | `source["corrections"]` |

## Check for source changes

Every load checks that a correction still matches the published records.

| Status | Meaning |
| --- | --- |
| `applied` | Expected source values match; correction applied |
| `fixed_upstream` | The source already contains the corrected values |
| `stale` | The record changed or disappeared; correction skipped, with a warning |
| `unavailable` | The snapshot lacks a required source |
| `disabled` | Corrections were turned off |

Corrections apply all or nothing. A correction may also check an unchanged
record as evidence that the original problem remains. For example,
`tempus-grch37-counts` checks PDZRN4's published 7/7 count before clearing rows.

```sh
osteosarc curation baseline --strict
osteosarc --no-corrections variants baseline --gene MAP2
```

`--strict` exits nonzero for stale corrections or unrecognized source labels.
Inspect the correction's `changes` field to see which records differ. Check
the new source evidence before revising or removing a correction.

Source labels are normalized separately: for example, `Boston Gene` becomes
`BostonGene`, and `CITE` becomes `cite-seq`. Original labels are retained.
Unknown labels appear in `data.unrecognized`.

## Add a local correction

```python
from osteosarc import CORRECTIONS, Change, Correction, glob

mine = Correction(
    "my-lab-note", "Flag CITE-seq BAMs for QC review.",
    (Change("bucket", {"key": glob("kamil/blood/output/*CITE*/outs/*.bam")}),),
    evidence=("internal QC log 2026-09",), verified="2026-09-18",
)
data = Dataset.open("baseline", corrections=[*CORRECTIONS, mine])
```

`Change(source, match, expect={}, set={})` selects published records by exact
field values or `glob(...)`. Dotted names address nested fields. An empty
`set` flags records without editing them. Pass a filtered list of `CORRECTIONS`
to disable individual corrections.

## Built-in corrections

All 32 applied to the 2026-09-18 snapshot when checked on 2026-09-20. Evidence URLs are available in
`data.corrections` and the [registry source](https://github.com/iskandr/osteosarc/blob/main/osteosarc/curation.py).

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
| `allele-FAM157A-p_W70_Q71ins_14` | edit | Supply the verified 42-base insertion at GRCh38 chr3:198153259; retain the withdrawn protein-model caveat. |
| `allele-COL3A1-Splice` | edit | Supply the public Tempus 737-base deletion matching c.4254+1_4255-1del, anchored at GRCh38 chr2:189010889. |
| `map2-split-representations` | flag | Two other MAP2 entries are pieces of that same event. |
| `muc3a-grch38-placement` | annotate | Record the GRCh37 call and assembly gap; leave GRCh38 placement unresolved. |
| `ush2a-transposed-duplicate` | annotate | Record a possible relationship to chr1:215650752; original Natera identity remains unconfirmed. |
| `fam157a-withdrawn-protein` | flag | The annotated protein model has been withdrawn by NCBI. |
| `natera-alleles-unavailable` | flag | The original Natera report is unavailable. COL3A1 now has independent Tempus evidence. |
| `otud4-source-unavailable` | annotate | No public genomic allele was found; retain the entry and identify the missing source. |
| `transcript-DCHS2` | edit | `NM_1142552` becomes `NM_001142552.1`. |
| `gene-symbol-TRMO` | edit | `TMRO` is a typo for `TRMO`; gene-symbol joins with pVACseq otherwise miss it. |

The relocated alleles and MAP2 were mapped from the original GRCh37 Tempus
TL-24-ALMY2X4KMV VCFs with Ensembl. Each mapped REF matched GRCh38. The
published positions fall outside the indel-equivalence spans, so normalization
does not explain the differences. Several wrong positions still matched a
reference base in the correct gene; `ready` alone cannot validate an allele.

FAM157A and COL3A1 were checked separately against UCSC's GRCh37-to-GRCh38
chain and versioned NCBI RefSeq windows. The complete REF and surrounding
sequence agree across assemblies. Offline fixtures test the source VCFs,
reference bases and equivalent indel representations. See
[the five reviewed entries](variants.md) for usage and remaining limitations.

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

## Other source limitations (2026-09-18)

* The `Tempus 2022` detection label refers to a 2024 accession.
* Three organoid DRAGEN BAMs and the T2 UCLA blood DRAGEN BAM have no specimen
  assignment in the consolidated metadata. Four 2026 blood specimens have no BAMs.
* ELISPOT `experiments[].date` is the earliest PBMC sample date, not the assay date.
* The legacy `data/treatment_timeline.json` is stale and is not used here.
* LENS files and Natera/BostonGene clinical reports were absent from the bucket.

The Tempus count correction was checked against the site's `pileup-json`
script: 171 of 172 SNV rows reproduced at the wrong locus. The BAM uses b37
(`human_g1k_v37`), while the site queried GRCh38 positions. Its MT-ND5 0/0
count also reflects a failed `chrM` to `MT` lookup. These rows are cleared,
not interpreted as negative evidence.
