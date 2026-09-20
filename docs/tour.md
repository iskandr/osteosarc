# Worked example: check MAP2 against the reads

The website's MAP2 vaccine target differs from the complex allele reported by
Tempus and CeGaT. This example compares both sequences against T1 tumor WGS
reads. It requires `osteosarc` and `samtools`.

## Compare the alleles

```python
from osteosarc import Dataset

data = Dataset.sync("baseline")
raw = Dataset.open("baseline", corrections=False)
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

The website ID stays the same after correction. See the
[correction registry](curation.md#read-counts-and-alleles) for its evidence.

## Fetch the surrounding reads

```python
source = data.asset(
    "kamil/oncoanalyser/IPISRC044_T1_ucla/alignments/dna/IPISRC044_tumor_T1_ucla.redux.bam"
)
print(data.inspect_alignment(source).assembly)
subset = data.extract_reads(source, [variant.region(padding=30)])
print(subset.path, subset.receipt["records"])
```

This requests a region of the indexed BAM. It does not download the whole file.

## Count exact sequence matches

The sequences below include bases on both sides of the event. Count matches
among primary, non-duplicate reads with mapping quality at least 20:

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

Thirty reads match the corrected allele and none match the published allele.
This is an exact-sequence check: reads in `no exact match` may stop short of the
event or contain mismatches. For RNA evidence and reconstruction, use
[Isovar](consumers.md#isovar).
