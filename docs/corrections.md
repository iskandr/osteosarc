# Corrections

Osteosarc fixes 35 known problems in the website's data, and does it by default.
Every time a snapshot is opened, each fix is checked against the data. You can see
what changed, turn the fixes off, or add your own.

## Compare corrected and published values

```python
from osteosarc import Dataset

data = Dataset.open()
raw = Dataset.open(corrections=False)
map2 = data.variants()["MAP2-chr2-209694768"]
print("Corrected:", map2.allele)
print("Published:", raw.variants()[map2.id].allele)
print(map2.annotations["corrections"])
```

When a fix changes a variant's allele, it usually clears that variant's read counts,
because they were measured for a different allele or position. MAP2 is the exception
(see below). Variant IDs are the snapshot's own; newer snapshots may rename them. The
[MAP2 example](map2.md) checks a correction against the reads.

## See every correction

```python
for row in data.corrections:
    print(row["id"], row["status"], row["summary"], row["evidence"])
```

`osteosarc corrections` lists them, and `osteosarc corrections ID` shows one with
the records it checks and its evidence.

Every record a fix touches carries its ID:

| Record | Where the IDs are |
| --- | --- |
| Variant | `variant.annotations["corrections"]` |
| Some of a variant's count rows | `variant.annotations["count_corrections"]` |
| File | `file.metadata["corrections"]` |
| Sample | `sample.corrections` |
| Count, vaccine, annotation and measurement rows | the `corrections` field |
| Timeline event | `event.corrections` |
| Variant selection or Varcode metadata | `source["corrections"]` |

## When the website changes

Each fix records what the data should look like before it's applied, so it notices
when the website changes:

| Status | Meaning |
| --- | --- |
| `applied` | The problem is still in the data, and the fix was applied |
| `fixed_upstream` | The website now has the corrected value, or has removed the bad record |
| `stale` | The data changed in an unexpected way; the fix was skipped, with a warning |
| `unavailable` | The snapshot doesn't include the data the fix needs |
| `disabled` | Corrections are turned off |

A fix applies completely or not at all. Some fixes also check a record they don't
change, to confirm the problem is still there: `tempus-grch37-counts`, for example,
checks PDZRN4's published 7/7 before clearing counts. When the website reorganizes
records, a fix can list each layout it has been checked against; exactly one must
match, or the fix is `stale`.

```sh
osteosarc corrections --strict
osteosarc --no-corrections variants --gene MAP2
```

`--strict` fails when a fix is stale or the data uses a label Osteosarc doesn't
know. It checks the snapshot you have; run `osteosarc sync --refresh` first to check
the live website. A fix's `changes` field shows which records differ. Check the new
data before changing or removing a fix.

Labels are also tidied separately, without being counted as corrections: `Boston
Gene` becomes `BostonGene`, and `CITE` becomes `cite-seq`. The original labels are
kept, and unknown ones are listed in `data.unrecognized`.

## Add your own correction

```python
from osteosarc import CORRECTIONS, Change, Correction, glob

mine = Correction(
    "my-lab-note", "Flag CITE-seq BAMs for QC review.",
    (Change("bucket", {"key": glob("kamil/blood/output/*CITE*/outs/*.bam")}),),
    evidence=("internal QC log 2026-09",), verified="2026-09-18",
)
data = Dataset.open(corrections=[*CORRECTIONS, mine])
```

`Change(source, match, expect={}, set={})` picks records by exact field values or
`glob(...)`; dotted names reach nested fields. A `Change` with no `set` flags records
without editing them. To turn off some built-in fixes, pass a filtered list of
`CORRECTIONS`.

`Correction(..., alternatives=(changes,))` adds another layout the fix accepts, and
`Change(..., absent=True)` requires that no records match. Pair an absence check with
a check on a record that should still be there: a record disappearing doesn't prove
the problem was fixed.

## The corrections

Checked against the website on 2026-09-24: 30 fixes apply and 5 are already fixed on
the site. An older snapshot from 2026-09-18 uses all 35. Each fix's evidence is in
`data.corrections` and in the
[source code](https://github.com/iskandr/osteosarc/blob/main/osteosarc/curation.py).
None of them only rewrites a value in an equivalent form: each fixes a wrong value,
fills in a missing one, or flags something to be careful with.

### Variants and read counts

| ID | Does | What |
| --- | --- | --- |
| `tempus-grch37-counts` | edit | Clear the counts the site measured in a GRCh37 Tempus exome BAM at GRCh38 positions, which describe unrelated places. Newer snapshots drop those rows. |
| `allele-CABLES1-chr18-23135500` | edit | Give the Tempus allele, chr18:23135764 T>TGGCGGC. The site had the call 264 bp away with a `dup` placeholder; it has since moved it, but still without an allele. |
| `allele-CCDC40-chr17-80058951` | edit | Give the Tempus allele at chr17:80090148. The site had it 31 kb away with `not_reported`. |
| `allele-DCHS2-chr4-154322488` | edit | Give the Tempus allele (a complex change) at chr4:154323273. The site had it 785 bp away, with `not_reported` and a wrong REF. |
| `allele-GAPVD1-chr9-125299105` | edit | Give the Tempus allele, chr9:125301980 TAGTGC>ATTGG. The site had it 2.9 kb away, with `not_reported` and a wrong REF. |
| `allele-GOLGA6L2-chr15-23441121` | edit | Give the Tempus allele, a 120-base insertion at chr15:23440197. The site had it 924 bp away, with `dup` and a wrong REF. |
| `allele-MAP2-chr2-209694768` | edit | The vaccine target's 22-bp deletion isn't what's in the tumor. Tempus calls one complex change (c.2599_2630delinsAGGG) and CeGaT calls the same change as three records; the catalogue's own protein sequence and the reads match it. The published counts are kept as an approximation: the site counts any large deletion there, which in the T0 tumor exome is always this change. |
| `allele-FAM157A-p_W70_Q71ins_14` | edit | Fill in the missing allele: a 42-base insertion at chr3:198153259, from a public Tempus call. |
| `allele-COL3A1-Splice` | edit | Fill in the missing allele: a 737-base deletion anchored at chr2:189010889, from a public Tempus call. It removes exactly intron 50, matching the catalogue's cDNA description. |
| `map2-split-representations` | flag | Two other MAP2 entries are pieces of the same change; don't count them separately. |
| `muc3a-grch38-placement` | note | The Tempus call is on GRCh37, in a repeat that differs between genome builds, so its GRCh38 position stays unresolved. |
| `ush2a-transposed-duplicate` | note | USH2A-chr1-215560752 is a typo of USH2A-chr1-215650752: the first position is outside the USH2A gene. The site has merged them; the kept entry's location label and sequence context still came from the typo, and are fixed. |
| `fam157a-withdrawn-protein` | flag | NCBI withdrew the protein model behind FAM157A's insertion. The site has since added the same note. |
| `natera-alleles-unavailable` | flag | COL3A1 and OTUD4 come from a Natera report that isn't public. A public Tempus call now gives COL3A1's allele. |
| `otud4-source-unavailable` | note | No public call gives OTUD4's allele, so the entry stays without one. |
| `transcript-DCHS2` | edit | The accession NM_1142552 is missing two zeros. The fix is NM_001142552, which the site now uses too. |
| `transcript-COL4A2` | edit | The accession `NM_001846.` has a stray dot. |
| `transcript-GTF3C5` | edit | The accession NM_00112283 is missing a digit; it's NM_001122823. |
| `gene-symbol-TRMO` | edit | `TMRO` is a typo for `TRMO`. |

The five moved alleles and MAP2 come from the Tempus TL-24-ALMY2X4KMV calls in the
bucket, converted from GRCh37 to GRCh38 with Ensembl. Each REF matches GRCh38, and no
equivalent way of writing an allele reaches its old position. These fixes make the
catalogue match what Tempus called; several of these calls have few supporting reads
or sit in repeats, so they don't prove a variant is real. FAM157A and COL3A1 were also
checked with UCSC's liftover chain and NCBI's reference sequence.

The GRCh37 counts were wrong because the site's pileup script looked up GRCh38
positions in a BAM aligned to GRCh37 (`human_g1k_v37`). Its MT-ND5 count of 0/0 comes
from the same problem. These counts are cleared rather than read as zero.

### Samples, pipelines and files

| ID | Does | What |
| --- | --- | --- |
| `viewer-label-BG009368` | edit | The site labeled this reprocessed RNA T0 2022-12; its FASTQ and allele fractions show it's T1 2024-06. The site now uses the fixed label. |
| `viewer-label-SARC0277` | edit | Likewise, it's T2 UCLA 2025-01, not T0 BostonGene. The site now uses the fixed label. |
| `provider-IPISRC044-T1-rna` | edit | The oncoanalyser T1 RNA BAM was built from BostonGene's FASTQs, but the site labels it UCLA, in its file list and its read counts. |
| `pvac-2025-detection` | edit | CDC40, PIP5K1A, SMC5 and TECPR1 are in the 2025-04-27 pVACtools runs, which the site's detection flag misses. |
| `pvac-header-only-filtered-reports` | flag | Seven filtered pVACseq reports are empty. An empty filtered report doesn't mean no epitope passed. |
| `pvac-rna-fields-na` | flag | RNA depth, VAF and expression are `NA` in every pVACseq report: RNA was never given to pVACseq. |
| `pvac-extended-run-no-class-i` | flag | The "MHCI.extended" run only has Class II reports. |

### Timeline and samples

| ID | Does | What |
| --- | --- | --- |
| `specimen-T1-site` | edit | T1 was a UCLA biopsy, not a UCSF resection. It spans two UCLA biopsies, on 2024-06-06 and 2024-06-11. |
| `specimen-T2-date-site` | edit | T2 was a UCLA biopsy on 2025-01-28, not UCSF on 2025-01-06. |
| `specimen-T3-site` | edit | T3 was an MSKCC resection, not a UCSF biopsy. |
| `pbmc-capture-dates` | flag | Four blood samples are dated by when their cells were captured; the blood was drawn two to four days earlier. |
| `events-duplicate-rows` | flag | SQ3370 and Trabectedin each appear twice. |
| `tempus-timepoint` | flag | The timeline dates the Tempus tests to T0 (2022), but the site labels the Tempus files T1 (2024-06). The timeline agrees with the data. |
| `tempus-file-labels` | flag | The Tempus files are labeled T1 2024-06, but their variants match the T0 tumor: they carry all three variants seen only at T0 and none of the 35 seen only at T1. The labels are left as published. |
| `apheresis-date` | flag | The apheresis is 2024-05-14 on the timeline and 2024-05-15 in the ELISPOT records. |
| `reyagel-end-date` | flag | The sheet's ReyaGel end date is `7/14`, with no year, so the site shows a single day. It's the only broken date among the sheet's 371. |

## Other things to know

- The Tempus files are labeled T1 but look like the T0 tumor (see `tempus-file-labels`).
- Some DRAGEN BAMs, mostly blood normals and organoid runs, aren't linked to a
  sample; list them with `data.files.select(prefix="kamil/basespace/results/")`.
  The four 2026 blood samples have no BAMs yet, only FASTQs.
- An ELISPOT experiment's date is the date of the earliest blood sample it used, not
  of the assay.
- The site's older `data/treatment_timeline.json` is out of date, and isn't used.
- The site lists LENS pipeline detections, but has no LENS files to download.
