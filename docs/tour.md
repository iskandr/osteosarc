# Worked example: check MAP2 against the reads

The website's MAP2 vaccine target is a 22-bp deletion, but Tempus and CeGaT both
call a different, complex change there. This example checks both against T1 tumor
WGS reads, and shows how to check any [correction](curation.md) against the data. It
needs `samtools`.

## Compare the alleles

```python
from osteosarc import Dataset

data = Dataset.sync()
raw = Dataset.open(corrections=False)
variant = data.variants()["MAP2-chr2-209694768"]
print("Published:", raw.variants()[variant.id].allele)
print("Corrected:", variant.allele)
print(variant.annotations["corrections"])
```

In the 2026-09-18 snapshot:

```text
Published: ('chr2', 209694768, 'CCTGGGCTACTGTGTGTTCAATA', 'C')
Corrected: ('chr2', 209694768, 'CCTGGGCTACTGTGTGTTCAATAAGTACACAGT', 'CAGGG')
```

The variant keeps its website ID. See [the corrections](curation.md#variants-and-read-counts)
for the evidence.

## Fetch the surrounding reads

```python
source = data.file(
    "kamil/oncoanalyser/IPISRC044_T1_ucla/alignments/dna/IPISRC044_tumor_T1_ucla.redux.bam"
)
print(data.inspect_alignment(source).assembly)
subset = data.extract_reads(source, [variant.region(padding=30)])
print(subset.path, subset.receipt["records"])
```

This fetches just that region of the BAM, not the whole file.

## Count exact sequence matches

The sequences below include a few bases on each side of the change. Count exact
matches among primary, non-duplicate reads with mapping quality 20 or more:

```python
from collections import Counter

haplotypes = {
    "corrected allele": "CTCGAAGACAGGGCCCATTGCCA",
    "published allele": "CTCGAAGACAGTACACAGTCCCATTGCCA",
    "reference": "CTCGAAGACCTGGGCTACTGTGTGTTCAATAAGTACACAGTCCCATTGCCA",
}
support = Counter()
with subset.open() as bam:
    for read in bam.fetch("chr2", 209694760, 209694800):
        if (read.is_secondary or read.is_supplementary or read.is_duplicate
                or read.mapping_quality < 20):
            continue
        hits = [name for name, seq in haplotypes.items()
                if seq in (read.query_sequence or "")]
        support[hits[0] if hits else "no exact match"] += 1
print(support)
```

Result from the 2026-09-18 check:

```text
Counter({'no exact match': 55, 'reference': 37, 'corrected allele': 30})
```

Thirty reads match the corrected allele and none match the published one. Reads
with `no exact match` may stop short of the change or have sequencing errors. For
RNA evidence and protein sequences, use [Isovar](consumers.md#isovar).
