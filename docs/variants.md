# Select variants and vaccine peptides

Pick variants from the website's catalogue by gene, vaccine, pipeline or status,
then get their alleles, read counts and vaccine peptides. The examples open your most recent snapshot;
see [Get started](index.md#get-started) to save one.

## Select variants

```python
from osteosarc import Dataset

data = Dataset.open()
site = data.variants()
vaccine_targets = data.variants("vaccine", status="ready")
dynein = site.select(gene="DYNC1H1", status="ready")
for variant in dynein:
    print(variant.id, variant.allele)
```

`data.variants()` gives every variant on the website's variants page, including
ones without a usable allele. `"vaccine"` gives the ones in at least one vaccine,
and `"all"` adds entries that appear only in the site's read-count table or
variant data file.

From the command line:

```sh
osteosarc variants --gene DYNC1H1 --status ready
osteosarc variants --set vaccine --vaccine mRNA
```

## Get an allele and its read-extraction region

```python
variant = site["DYNC1H1-chr14-101980529"]
print(variant.allele)  # (chromosome, one-based position, REF, ALT)
print(variant.region(padding=100))
regions = dynein.regions(padding=100)
```

`region()` covers the REF allele, in zero-based, half-open coordinates. To fetch
reads, pass the variants straight to [read extraction](reads.md).

## Variant status

`status` says whether an entry has a usable genomic allele. It says nothing
about samples or read support.

| `variant.status` | Meaning |
| --- | --- |
| `ready` | One allele, written in DNA bases, at a consistent position |
| `missing_literal_allele` | No allele given |
| `non_literal_allele` | REF or ALT is a placeholder such as `dup` or `not_reported` |
| `ambiguous_literal_allele` | More than one allele given for the entry |
| `conflicting_coordinates` | The site's pages and files disagree on the position or allele |
| `malformed_source_row` | A read-count row has a broken position or the wrong number of fields |

Only `ready` entries have `.allele`, `.region()` and Varcode conversion. `ready`
doesn't mean the variant is checked against the genome, somatic, expressed or
changes the protein.

Without the `status` filter you get every entry; `.allele` and `.region()` raise
an error on the unusable ones:

```python
unresolved = site.where(lambda v: v.status != "ready")
print([(v.id, v.status) for v in unresolved])
```

[Corrections](curation.md) are on by default and can change an entry's allele
and status. Open the snapshot with `corrections=False` to see the published ones.

## Read counts and annotations

```python
counts = data.vafs.select(variant_id=variant.id)
print(counts.rows[:2])
print(variant.annotations["source_record"])
print(data.pipeline_names)
detected = data.variants(pipeline="oncoanalyser")
```

Counts are kept exactly as published. A missing count isn't zero. Counts that
were measured for the wrong allele or position are cleared by a correction, so
they read as not measured. `data.vafs.diagnostics` lists rows with too few or too
many fields, with their line numbers.

## Get vaccine peptides

```python
print(data.vaccine_names)
mrna = data.variants(vaccine="mRNA")
for row in data.vaccine_peptides("mRNA"):
    print(row["variant_id"], row["sequence"], row["experiments"])
```

Which vaccines contain a variant comes from the site's vaccine-overlap file. The
site's variant data file flags vaccines separately; use
`vaccine_source="source_variants"` for those flags. The two can disagree, and the
annotations keep both.

```python
for row in data.vaccines.select(gene="SMC5"):
    print(row["elispot_status"], row["elispot_response"])
```

An ELISPOT that wasn't run, or has no recorded response, isn't a negative result.

## Open a VCF

```python
data = Dataset.open(offline=False)
calls = data.assets.select(kind="variants", format="vcf")
if calls:
    with data.open_variants(calls[0]) as vcf:
        print(vcf.header)
        for record in vcf:
            print(record.contig, record.pos, record.ref, record.alts)
            break
```

This downloads the whole VCF and its index, and returns a pysam reader with the
original headers, genotypes and records.
For native Varcode objects, see [Use other libraries](consumers.md#varcode).

## Entries with reviewed alleles

Five entries were reviewed in [issue #5](https://github.com/iskandr/osteosarc/issues/5).
Their `allele_resolution` annotation records the outcome and source evidence:

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
| USH2A-chr1-215560752 | Possible duplicate of USH2A-chr1-215650752. Historical snapshots retain both. The site now merges them; the original report is still needed to confirm identity. |

The two resolved alleles come from a public Tempus VCF, with versioned RefSeq
checks and source checksums in the annotation. Neither has published count
rows; use [read extraction](reads.md) to examine support.

In newer snapshots, the retained USH2A entry keeps its published `C>A` allele
and counts. Its `annotations["source_record"]["upstream_merge"]` records the
retired ID, the upstream commit and the identity caveat. Osteosarc corrects the
outdated location label without creating an allele for the retired entry.

Three different statuses appear here. `variant.status` describes allele
usability, as listed under [Variant status](#variant-status). The `status` column of `data.corrections`
describes whether a correction applied. The nested `allele_resolution["status"]`
describes the evidence review's outcome.

## Malformed source rows

A broken row in the site's read-count table marks its entry
`malformed_source_row` and leaves the other entries alone.
`variant.annotations["parse_errors"]` shows the row and what was wrong. To parse the
files yourself, use `parse_variants(index, vaf_tsv)`. A table missing a required
column raises `SchemaError`.
