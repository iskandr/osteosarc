# SV candidates

637 structural variants (SVs) in the patient's tumor that are worth a closer
look: deletions, duplications, inversions, translocations and RNA fusions. Each
comes with its calls, breakpoints, nearby genes, the genes' expression and the RNA
reads that cross it. The list ships with osteosarc and needs no network:

```python
from osteosarc import load_sv_candidates

candidates = load_sv_candidates()["targets"]
macrod2 = candidates["SV0203"]
print(len(candidates), macrod2["genes"], macrod2["breakends"])
```

## Why these

It collects every SV call in the public data so the OpenVax libraries test the
same set, and so none is dropped without a look:

- 523 SV calls from oncoanalyser (ESVEE, PURPLE and LINX) that pass its filters;
- 32 DRAGEN SV calls;
- 78 long-read RNA fusions from CTAT;
- 4 fusions added by hand: FOXO3–STRADA–CCDC47, PARD3B–CDKN2B, GABBR1–SLC29A1 and
  OTUD7A–FMN1.

Being on the list doesn't make a candidate real, expressed or protein-coding, and
none was dropped for lacking RNA support or a gene. The files it was built from,
as of 2026-09-22, are listed with their checksums in
`load_sv_candidates()["sources"]`.

## Reading an entry

- **Coordinates.** original_coordinates are the caller's one-based positions;
  breakends are zero-based positions between bases. A single breakend or an
  unclear shape is marked unresolved, and keeps the original call.
- **Duplicates.** Entries with the same adjacency group join the same two ends,
  often from different callers. Don't add up their counts.
- **Genes** come from Ensembl 115. A caller's gene name can be a nearby gene the
  break doesn't touch.
- **Expression** is gene-level TPM from the T2 bulk RNA-seq: how much the gene is
  expressed, not the fusion.
- **RNA support** counts reads crossing each breakpoint, per RNA BAM in the bucket.
  A missing count means the BAM wasn't checked.

Isovar and Topiary rebuild and rank the proteins these SVs might make.

## As test targets

The sv-candidates-v1 panel gives the same entries as targets for a
[test-data recipe](test-data.md):

```python
from osteosarc import load_panel

targets = load_panel("sv-candidates-v1")
```

A target is only a place; a recipe says which reads to fetch there.
