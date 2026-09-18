# Variants, annotations, and vaccines

```python
from osteosarc import Dataset

data = Dataset.open("baseline")
site = data.variants()
all_entries = data.variants("all")
vaccine_targets = data.variants("vaccine")
```

`site` includes unresolved website entries. `all` also includes additional
entries from the count export and source annotation JSON. `vaccine` selects
site entries with a positive website vaccine count.

## Select exact alleles

```python
ready = site.select(status="ready")
variant = site["DYNC1H1-chr14-101980529"]
print(variant.allele)  # original (chrom, one-based position, REF, ALT)
regions = ready.select(gene="DYNC1H1").regions(padding=100)
unresolved = site.where(lambda v: v.status != "ready")
print([(v.id, v.status) for v in unresolved])
```

`ready` means the sources provide one internally consistent literal allele.
It does not establish biological truth or validate REF against a genome.
Unresolved entries remain inspectable. Their `.allele` and `.region()` raise;
select explicitly before extracting reads. Indel anchors are retained.

## Inspect source annotations

```python
print(data.pipeline_names)
detected = data.variants(pipeline="oncoanalyser")
print(variant.annotations["source_record"])
counts = data.vafs.select(variant_id=variant.id)
print(counts.rows[:2])
```

`data.annotations` preserves the full source records. Caller VCF/BCF files use
their native parser so headers, genotypes, multiallelic records, and symbolic
alleles survive:

```python
data = Dataset.open("baseline", offline=False)
calls = data.assets.select(kind="variants", format="vcf")
if calls:
    with data.open_variants(calls[0]) as vcf:
        print(vcf.header)
        for record in vcf:
            print(record.contig, record.pos, record.ref, record.alts)
            break
```

`open_variants` downloads the selected full file and its listed index, when
present. Use `.fetch(contig, start, end)` for indexed subsetting afterward.

## Vaccine membership and experiments

```python
print(data.vaccine_names)
mrna = data.variants(vaccine="mRNA")
source_flags = data.variants(vaccine="mRNA", vaccine_source="source_variants")
for row in data.vaccine_peptides("mRNA"):
    print(row["variant_id"], row["sequence"], row["experiments"])
```

Default membership uses the vaccine-overlap JSON, joined only at an
unambiguous gene/locus. `source_variants` requests the separate source JSON
flags. The site count, membership, peptide inclusion, and experimental results
are distinct source assertions; disagreements remain in annotations.

```python
for row in data.vaccines.select(gene="SMC5"):
    print(row["elispot_status"], row["elispot_response"])
```

`not_tested` and a missing response do not mean a negative assay. Published
peptides are comparator data, not recommendations for vaccine design.
