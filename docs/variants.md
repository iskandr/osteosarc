# Select variants and vaccine peptides

These examples use the `baseline` snapshot from [Get started](index.md).

## Select variants

```python
from osteosarc import Dataset

data = Dataset.open("baseline")
site = data.variants()
vaccine_targets = data.variants("vaccine", status="ready")
dynein = site.select(gene="DYNC1H1", status="ready")
for variant in dynein:
    print(variant.id, variant.allele)
```

`site` includes all website entries, including unresolved alleles. The `vaccine`
set selects entries with a positive site vaccine count. Use `data.variants("all")`
to include additional count-export and source-JSON entries.

## Get an allele and its read-extraction region

```python
variant = site["DYNC1H1-chr14-101980529"]
print(variant.allele)  # (chromosome, one-based position, REF, ALT)
print(variant.region(padding=100))
regions = dynein.regions(padding=100)
```

`ready` means one literal allele with internally consistent coordinates. It does
not validate REF against a genome or establish somatic status. Unresolved
entries remain visible, but `.allele` and `.region()` raise errors:

```python
unresolved = site.where(lambda v: v.status != "ready")
print([(v.id, v.status) for v in unresolved])
```

[Source corrections](curation.md) are on by default. Use `corrections=False`
when opening the snapshot to inspect the published alleles.

A malformed count-export row marks its entry `malformed_source_row`; other
entries remain available. Inspect `variant.annotations["parse_errors"]` for
the source values and error. For standalone parsing, pass the index and TSV
text to `parse_variants(index, vaf_tsv)`. Missing or duplicate required headers
still raise `SchemaError`.

## Read counts and annotations

```python
counts = data.vafs.select(variant_id=variant.id)
print(counts.rows[:2])
print(variant.annotations["source_record"])
print(data.pipeline_names)
detected = data.variants(pipeline="oncoanalyser")
```

Counts retain raw source values. Missing counts differ from zero; counts for
corrected alleles or wrongly mapped alignments are cleared to unmeasured.
`data.vafs.diagnostics` retains ragged rows' original fields and line numbers;
missing trailing cells are `None`.

## Get vaccine peptides

```python
print(data.vaccine_names)
mrna = data.variants(vaccine="mRNA")
for row in data.vaccine_peptides("mRNA"):
    print(row["variant_id"], row["sequence"], row["experiments"])
```

Vaccine membership comes from the overlap JSON. For the separate flags in the
source variant JSON, use `vaccine_source="source_variants"`. These sources can
disagree; the annotations retain those disagreements.

```python
for row in data.vaccines.select(gene="SMC5"):
    print(row["elispot_status"], row["elispot_response"])
```

An untested assay or missing response is not a negative result.

## Open a VCF

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

This downloads the full VCF and its listed index.
The returned pysam reader preserves headers, genotypes, and multiallelic records.
For native Varcode objects, see [Use other libraries](consumers.md#varcode).
