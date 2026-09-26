# Variants and vaccines

The website lists the patient's tumor variants, which pipelines found them, how
many reads support them in each BAM, and which went into the cancer vaccines.
The examples open your most recent snapshot; see [Get started](index.md#get-started).

## Select variants

```python
from osteosarc import Dataset

data = Dataset.open()
site = data.variants()
dynein = site.select(gene="DYNC1H1", status="ready")
for variant in dynein:
    print(variant.id, variant.allele)
```

`data.variants()` gives every variant on the website's variants page.
`data.variants("vaccine")` gives the ones in at least one vaccine, and
`data.variants("all")` adds entries found only in the site's read-count table.
Select by gene, vaccine, pipeline or status.

From a terminal:

```sh
osteosarc variants --gene DYNC1H1 --status ready
osteosarc variants --vaccine mRNA
osteosarc variants DYNC1H1-chr14-101980529
```

The last shows one variant with its allele, effect, vaccines, corrections and the
site's read counts.

## Alleles

```python
variant = site["DYNC1H1-chr14-101980529"]
print(variant.allele)               # (chromosome, one-based position, REF, ALT)
print(variant.region(padding=100))  # zero-based, half-open
```

Pass variants straight to [read extraction](reads.md) to get the reads around them.

## Variant status

A variant's status says whether it has one usable allele. It says nothing about
samples or read support.

| Status | Meaning |
| --- | --- |
| ready | One allele, written in DNA bases |
| missing_literal_allele | No allele given |
| non_literal_allele | REF or ALT is a placeholder, such as dup or not_reported |
| ambiguous_literal_allele | More than one allele given |
| conflicting_coordinates | The site's pages disagree on the position or allele |
| malformed_source_row | A read-count row is broken; the other entries are unaffected |

Only ready variants have an allele and a region. Ready doesn't mean the variant is
somatic, expressed or changes the protein. [Corrections](corrections.md), on by
default, fix some alleles and statuses; open the snapshot with
`corrections=False` to see the published ones.

## Read counts

```python
counts = data.vafs.select(variant_id=variant.id)
print(counts.rows[:2])
```

Counts are kept as published, one row per variant and BAM. A missing count isn't
zero: counts measured at the wrong allele or position are cleared by a correction,
so they read as not measured.

## Vaccines

```python
for row in data.vaccine_peptides("mRNA"):
    print(row["variant_id"], row["sequence"])

for row in data.vaccines.select(gene="SMC5"):
    print(row["vaccines"], row["elispot_status"])
```

The vaccines table has one row per vaccine target: which vaccines include it, and
its ELISPOT result. An ELISPOT that wasn't run, or has no recorded response, isn't
a negative result. `osteosarc vaccines` prints the same table.

Which vaccines contain a variant comes from the site's vaccine-overlap file. Its
variant data file flags vaccines separately, and the two can disagree; select
with `vaccine_source="source_variants"` to use the flags.

## Reviewed alleles

Five entries needed their alleles reviewed
([issue #5](https://github.com/iskandr/osteosarc/issues/5)); each variant's
annotations hold the outcome and evidence:

```python
fam157a = site["FAM157A-p_W70_Q71ins_14"]
print(fam157a.allele)
print(fam157a.annotations["allele_resolution"]["summary"])
```

| Entry | Outcome |
| --- | --- |
| FAM157A | A 42-base insertion at chr3:198153259, from a public Tempus call. NCBI withdrew the protein model, so it isn't known to make a neoantigen. |
| COL3A1 | A 737-base deletion at chr2:189010889 that removes intron 50, as the catalogue's cDNA says. |
| MUC3A | A GRCh37 duplication whose GRCh38 position can't be settled: it sits in a repeat that differs between builds. |
| OTUD4 | No public call gives its allele. |
| USH2A-chr1-215560752 | A typo of USH2A-chr1-215650752. The site has since merged the two. |

## VCFs and your own parsing

The bucket's VCFs are ordinary files: find them with
`data.files.select(kind="variants")`, download one with `data.download(file)` and
open it with pysam. For Varcode variants, see
[OpenVax libraries](openvax.md#varcode). To parse the site's variant files yourself,
use `osteosarc.parsing.parse_variants`.
