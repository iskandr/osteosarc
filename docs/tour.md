# A guided tour

This page explores the dataset end to end. Each step answers a question, using
the same objects the other pages describe. Every code block on this site is
executed by `scripts/check_docs.py` (see [validation](validation.md)). The
printed results shown here come from the 2026-09-18 sources.

## 1. Make a snapshot and look around

```python
from osteosarc import Dataset
from osteosarc.explore import summary_view

data = Dataset.sync("baseline")      # ~57 MB of metadata, once; later use Dataset.open
print(summary_view(data))
```

```text
snapshot baseline (4580057f4269)
corrections: 29 applied
timeline: 787 events, 2022-10-06 .. 2026-09-16, 31 lanes
specimens: 21
variants: 182 on site, 177 ready, 3 missing_literal_allele, 2 non_literal_allele
assets: 380,033 other, 8,137 table, 2,396 reads, 1,657 index, 1,043 reference, 843 alignment, ...
```

The "assets" are files in the public bucket. They are processing products,
not biological samples. The 21 specimens are the biological samples. The
snapshot ID identifies this exact set of source receipts, so yours will differ.

## 2. What happened, and when?

```python
print(data.timeline.render(width=110, since="2022-11", until="2023-12"))
print(data.timeline.select(lane="Time points").listing())
```

Four sampling time points anchor the omics data: T0 (primary resection,
2022-12-16), T1 (UCLA biopsies, June 2024), T2 (UCLA biopsy, 2025-01-28), and
T3 (MSKCC resection, 2025-04-17). Everything around a time point:

```python
print(data.timeline.around("2025-04-17", days=2).listing())
```

## 3. Which specimens exist, and what was measured on them?

```python
from osteosarc.explore import specimen_view, specimens_view

print(specimens_view(data, width=120))
print(specimen_view(data, "T1_tumor"))
```

The detail view lists the specimen's BAMs by assay, its FASTQ folders, any
correction applied to its registry row, and the events within a week of
collection.

## 4. Choose files for an analysis

Metadata filters use unambiguous published claims:

```python
t2_rna = data.assets.select(kind="alignment", assay="rna-seq", timepoint="T2")
for asset in t2_rna:
    print(asset.key, asset.values("provider"), asset.size)
```

Before choosing coordinates, check an alignment's actual reference. The
catalogue's global "hg38" label is not a guarantee:

```python
data = Dataset.open("baseline", offline=False)
tempus = data.asset("vendor/tempus/TL-24-5GQLV9WSXQ/DNA/TL-24-5GQLV9WSXQ_T.sorted.bam")
print(data.inspect_alignment(tempus).assembly)   # GRCh37: a b37 BAM
```

This one BAM is why the website's Tempus read counts were wrong: they were
counted at GRCh38 positions. The `tempus-grch37-counts` correction clears them.

## 5. Variants, and what the corrections changed

```python
from collections import Counter

site = data.variants()
print(Counter(v.status for v in site))
print(Counter(c for v in site for c in v.annotations["corrections"]).most_common(4))
print(Counter(c for v in site for c in v.annotations["count_corrections"]))

raw = Dataset.open("baseline", corrections=False).variants()
for v in site:
    if v.alleles != raw[v.id].alleles:
        print(v.id, raw[v.id].alleles[0][1:], "->", v.allele[1:])
```

## 6. Check a correction against the reads

The MAP2 vaccine target was published as a 22-bp deletion. The Tempus and
CeGaT calls describe a −28 bp complex event. The T1 tumor WGS reads settle it:

```text
Counter({'ready': 177, 'missing_literal_allele': 3, 'non_literal_allele': 2})
[('pvac-2025-detection', 4), ('natera-alleles-unavailable', 2), ('allele-CABLES1-chr18-23135500', 1), ('allele-CCDC40-chr17-80058951', 1)]
Counter({'tempus-grch37-counts': 179})
CABLES1-chr18-23135500 (23135500, 'G', 'dup') -> (23135764, 'T', 'TGGCGGC')
...
MAP2-chr2-209694768 (209694768, 'CCTGGGCTACTGTGTGTTCAATA', 'C') -> (209694768, 'CCTGGGCTACTGTGTGTTCAATAAGTACACAGT', 'CAGGG')
```

Seventeen variants carry corrections of their own. The allele edits are the
six that changed. `tempus-grch37-counts` cleared one BAM's count row for 179
variants without changing the variants themselves, so it is listed under
`count_corrections`.

```python
from collections import Counter

variant = site["MAP2-chr2-209694768"]
source = data.asset("kamil/oncoanalyser/IPISRC044_T1_ucla/alignments/dna/IPISRC044_tumor_T1_ucla.redux.bam")
subset = data.extract_reads(source, [variant.region(padding=30)])
haplotypes = {"observed complex event": "CTCGAAGACAGGGCCCATTGCCA",
              "published 22-bp deletion": "CTCGAAGACAGTACACAGTCCCATTGCCA",
              "reference": "CTCGAAGACCTGGGCTACTGTGTGTTCAATAAGTACACAGTCCCATTGCCA"}
support = Counter()
with subset.open() as bam:
    for read in bam.fetch("chr2", 209694760, 209694800):
        if read.is_secondary or read.is_supplementary or read.is_duplicate or read.mapping_quality < 20:
            continue
        hit = [name for name, seq in haplotypes.items() if seq in (read.query_sequence or "")]
        support[hit[0] if hit else "does not span the event"] += 1
print(support)
```

```text
Counter({'does not span the event': 55, 'reference': 37, 'observed complex event': 30})
```

Thirty reads carry the observed event, and none carry the published deletion.
Only about 280 regional records were transferred; the whole BAM was never
downloaded.

## 7. MRD and vaccines over time

```python
doses = data.timeline.select(lane="Cancer vaccines", since="2025-06", until="2025-12")
print(doses.listing())
signatera = data.measurements.select(source="mrd", measurement="Signatera")
print([(r["date"], r["kind"], r["value"]) for r in signatera.rows[-5:]])
print(data.timeline.select(category="MRD").render(width=100, since="2025"))
```

## 8. When the website changes

```python
new = Dataset.sync("after-update", refresh=True)     # fresh metadata, a new snapshot
for row in new.corrections:
    if row["status"] != "applied":
        print(row["id"], row["status"], row["changes"])
print(list(new.unrecognized))
```

Corrections that no longer match their sources are not applied, and are
reported instead. See [corrections and source drift](curation.md).

## Terminal only

Everything above is also available without Python:

```sh
osteosarc timeline baseline --since 2024-05 --until 2024-09
osteosarc specimens baseline T2_tumor
osteosarc variants baseline --gene MAP2
osteosarc curation baseline
```

`osteosarc explore baseline` opens the interactive version.
