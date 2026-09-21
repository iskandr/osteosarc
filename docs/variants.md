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

## Variant status

`status` describes whether a catalogue entry has a usable genomic allele.
It is independent of the sample, sequencing assay and measured read support.

| `variant.status` | Meaning |
| --- | --- |
| `ready` | One allele with literal DNA bases and consistent chromosome/position fields |
| `missing_literal_allele` | No usable genomic allele was supplied |
| `non_literal_allele` | REF or ALT contains a placeholder such as `dup` or `not_reported` |
| `ambiguous_literal_allele` | Multiple candidate alleles were supplied for the same entry |
| `conflicting_coordinates` | The index, count export or source JSON disagree on the locus or allele |
| `malformed_source_row` | A count-export row has an invalid position or missing/extra fields |

`ready` permits `.allele`, `.region()` and Varcode conversion. It does not
validate REF against a genome, establish somatic status, guarantee RNA support
or establish a protein effect. These need separate evidence. The source-backed
checks for individual corrected alleles are recorded in their annotations.

Omit the `status` filter to keep all entries. Unresolved entries remain visible,
but `.allele` and `.region()` raise errors:

```python
unresolved = site.where(lambda v: v.status != "ready")
print([(v.id, v.status) for v in unresolved])
```

[Source corrections](curation.md) are on by default. Use `corrections=False`
when opening the snapshot to inspect the published alleles.

The `status` in `data.corrections` describes whether a correction applied.
The nested `allele_resolution["status"]` describes the evidence review's outcome.
Neither replaces the `variant.status` values above.

For the five entries reviewed in [issue #5](https://github.com/iskandr/osteosarc/issues/5),
`allele_resolution` records the outcome and source evidence:

```python
for v in site:
    resolution = v.annotations.get("allele_resolution")
    if resolution:
        print(v.id, resolution["status"], resolution["summary"])

fam157a = site["FAM157A-p_W70_Q71ins_14"]
print(fam157a.allele)
print(fam157a.annotations["allele_resolution"]["protein_interpretation"])
```

| Entry | Outcome |
| --- | --- |
| FAM157A | Verified 42-base insertion at GRCh38 chr3:198153259. The protein model was withdrawn; a usable genomic allele does not establish a neoantigen. |
| COL3A1 | Verified 737-base deletion anchored at GRCh38 chr2:189010889, matching the catalogue's cDNA annotation. |
| MUC3A | Public GRCh37 duplication; GRCh38 placement remains unresolved across an assembly gap. |
| OTUD4 | Source unavailable: the protein label alone does not identify a genomic allele. |
| USH2A-chr1-215560752 | Possible duplicate of USH2A-chr1-215650752. The original report is needed to confirm identity; both entries remain separate. |

The two resolved alleles come from a public Tempus VCF, with versioned RefSeq
checks and source checksums in the annotation. Neither has published count
rows; use [read extraction](reads.md) to examine support.

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
