# Corrections

Osteosarc fixes 35 known problems in the website's data, by default. Each time a
snapshot is opened, every fix is checked against the data. You can see what changed,
turn the fixes off, or add your own.

```sh
osteosarc corrections                               # every fix and its status
osteosarc corrections allele-MAP2-chr2-209694768    # one fix, with its evidence
osteosarc --no-corrections variants --gene MAP2     # the website's values
```

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

A fix that changes a variant's allele usually clears its read counts, since they
were measured for a different allele or position; MAP2 is the exception (see below).
Variant IDs are the snapshot's own, and newer snapshots may rename them. The
[MAP2 example](map2.md) checks a correction against the reads.

## See every correction

```python
for row in data.corrections:
    print(row["id"], row["status"], row["summary"])
```

Every record a fix touches names it: a variant in its annotations (under
corrections, or count_corrections when only some of its read counts changed), a file
in its metadata, and a sample, timeline event, or row of counts, vaccines or
measurements in its own corrections field.

## When the website changes

Each fix records what the data should look like before it's applied, so it notices
when the website changes. Its status says what it found:

| Status | Meaning |
| --- | --- |
| applied | The problem is still there, and the fix was applied |
| fixed_upstream | The website now has the right value, or has removed the bad record |
| stale | The data changed unexpectedly; the fix was skipped, with a warning |
| unavailable | The snapshot doesn't have the data the fix needs |
| disabled | Corrections are turned off |

A fix applies completely or not at all. Some also check a record they don't change,
to confirm the problem is still there; tempus-grch37-counts, for example, checks
PDZRN4's published 7/7 before clearing counts. When the website reorganizes records, a
fix can list each layout it accepts; exactly one must match, or the fix is stale.

```sh
osteosarc corrections --strict
```

This exits with an error if a fix is stale or the data uses a label osteosarc doesn't
know. It checks the snapshot you have; run `osteosarc sync --refresh` first to check the
live website, and look at the new data before changing or removing a fix.

Labels are tidied separately, without counting as corrections: Boston Gene becomes
BostonGene, and CITE becomes cite-seq. The original labels are kept.

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

A change picks records by exact field values or glob patterns (dotted names reach
nested fields). Without new values to set, it only flags the records. To turn off some
built-in fixes, pass a shorter list. A correction can accept other record layouts, and
a change can require that no record matches; pair such a check with one on a record
that should still be there, since a record disappearing doesn't prove a fix.

## The corrections

Checked against the website on 2026-09-24: 30 fixes apply and 5 are already fixed on
the site; a snapshot from 2026-09-18 uses all 35. Each fix's evidence is in
`data.corrections` and the
[source](https://github.com/iskandr/osteosarc/blob/main/osteosarc/curation.py). None
only rewrites a value in an equivalent form: each fixes a wrong value, fills in a
missing one, or flags something to be careful with.

### Variants and read counts

| ID | Does | What |
| --- | --- | --- |
| tempus-grch37-counts | edit | Clears counts measured in a GRCh37 Tempus exome BAM at GRCh38 positions, which describe unrelated places. Newer snapshots drop those rows. |
| allele-CABLES1-chr18-23135500 | edit | Gives the Tempus allele, chr18:23135764 T>TGGCGGC; the site had it 264 bp away as a dup placeholder. |
| allele-CCDC40-chr17-80058951 | edit | Gives the Tempus allele at chr17:80090148; the site had it 31 kb away, not reported. |
| allele-DCHS2-chr4-154322488 | edit | Gives the Tempus allele (a complex change) at chr4:154323273; the site had it 785 bp away with a wrong REF. |
| allele-GAPVD1-chr9-125299105 | edit | Gives the Tempus allele, chr9:125301980 TAGTGC>ATTGG; the site had it 2.9 kb away with a wrong REF. |
| allele-GOLGA6L2-chr15-23441121 | edit | Gives the Tempus allele, a 120-base insertion at chr15:23440197; the site had it 924 bp away as a dup. |
| allele-MAP2-chr2-209694768 | edit | The vaccine target's 22-bp deletion isn't what's in the tumor. Tempus calls one complex change (c.2599_2630delinsAGGG), CeGaT the same as three records, and the catalogue's protein and the reads match it. Its counts are kept: the site counts any large deletion there, which in the T0 exome is always this one. |
| allele-FAM157A-p_W70_Q71ins_14 | edit | Fills in the missing allele: a 42-base insertion at chr3:198153259, from a public Tempus call. |
| allele-COL3A1-Splice | edit | Fills in the missing allele: a 737-base deletion at chr2:189010889 that removes exactly intron 50, from a public Tempus call. |
| map2-split-representations | flag | Two other MAP2 entries are pieces of the same change; don't count them separately. |
| muc3a-grch38-placement | note | The Tempus call is on GRCh37, in a repeat that differs between builds, so its GRCh38 position stays unresolved. |
| ush2a-transposed-duplicate | note | USH2A-chr1-215560752 is a typo of USH2A-chr1-215650752. The site has merged them; the kept entry's location label and sequence context are fixed. |
| fam157a-withdrawn-protein | flag | NCBI withdrew the protein model behind FAM157A's insertion; the site now says so too. |
| natera-alleles-unavailable | flag | COL3A1 and OTUD4 come from a Natera report that isn't public; a Tempus call now gives COL3A1's allele. |
| otud4-source-unavailable | note | No public call gives OTUD4's allele, so the entry has none. |
| transcript-DCHS2 | edit | NM_1142552 is missing two zeros: NM_001142552, which the site now uses. |
| transcript-COL4A2 | edit | NM_001846. has a stray dot. |
| transcript-GTF3C5 | edit | NM_00112283 is missing a digit: NM_001122823. |
| gene-symbol-TRMO | edit | TMRO is a typo for TRMO. |

The five moved alleles and MAP2's come from the Tempus TL-24-ALMY2X4KMV calls in the
bucket, converted from GRCh37 to GRCh38 with Ensembl. Each REF matches GRCh38, and no
equivalent way of writing an allele reaches its old position. These fixes make the
catalogue match what Tempus called; several calls have few reads or sit in repeats, so
they don't prove a variant is real. FAM157A and COL3A1 were also checked with UCSC's
liftover and NCBI's reference sequence.

The GRCh37 counts were wrong because the site's pileup script looked up GRCh38
positions in a BAM aligned to GRCh37; MT-ND5's 0/0 has the same cause. These counts are
cleared rather than read as zero.

### Samples, pipelines and files

| ID | Does | What |
| --- | --- | --- |
| viewer-label-BG009368 | edit | This reprocessed RNA is T1 2024-06, not T0 2022-12, by its FASTQ and allele fractions. The site now agrees. |
| viewer-label-SARC0277 | edit | Likewise, it's T2 UCLA 2025-01, not T0 BostonGene. The site now agrees. |
| provider-IPISRC044-T1-rna | edit | The oncoanalyser T1 RNA BAM was built from BostonGene's FASTQs, not UCLA's. |
| pvac-2025-detection | edit | CDC40, PIP5K1A, SMC5 and TECPR1 are in the 2025-04-27 pVACtools runs, which the site's detection flag misses. |
| pvac-header-only-filtered-reports | flag | Seven filtered pVACseq reports are empty, which doesn't mean no epitope passed. |
| pvac-rna-fields-na | flag | RNA depth, VAF and expression are NA in every pVACseq report: RNA was never given to pVACseq. |
| pvac-extended-run-no-class-i | flag | The "MHCI.extended" run has only Class II reports. |

### Timeline and samples

| ID | Does | What |
| --- | --- | --- |
| specimen-T1-site | edit | T1 was a UCLA biopsy, not a UCSF resection; it spans two UCLA biopsies, on 2024-06-06 and 2024-06-11. |
| specimen-T2-date-site | edit | T2 was a UCLA biopsy on 2025-01-28, not UCSF on 2025-01-06. |
| specimen-T3-site | edit | T3 was an MSKCC resection, not a UCSF biopsy. |
| pbmc-capture-dates | flag | Four blood samples are dated by cell capture; the blood was drawn two to four days earlier. |
| events-duplicate-rows | flag | SQ3370 and Trabectedin each appear twice. |
| tempus-timepoint | flag | The timeline dates the Tempus tests to T0 (2022), which the data supports, but the files are labeled T1 (2024-06). |
| tempus-file-labels | flag | The Tempus files are labeled T1, but carry all three variants seen only at T0 and none of the 35 seen only at T1. The labels are left as published. |
| apheresis-date | flag | The apheresis is 2024-05-14 on the timeline and 2024-05-15 in the ELISPOT records. |
| reyagel-end-date | flag | ReyaGel's end date is 7/14, with no year, so it shows as one day; the only broken date of the sheet's 371. |

## Other things to know

- Some DRAGEN BAMs, mostly blood normals and organoid runs, aren't linked to a sample;
  `osteosarc files --prefix kamil/basespace/results/` lists them. The four 2026 blood
  samples have only FASTQs so far.
- An ELISPOT experiment's date is that of the earliest blood sample it used.
- The site's older treatment timeline file is out of date, and isn't used.
- The site lists LENS pipeline detections, but has no LENS files to download.
