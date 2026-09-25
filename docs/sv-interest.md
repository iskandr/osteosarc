# SV catalogue

A list of 637 candidate structural variants (SVs) worth a closer look, with their
calls, breakpoints, nearby genes, expression and RNA support. It ships with the
package and loads without a network connection:

```python
from osteosarc import load_panel, load_sv_interest

catalogue = load_sv_interest()
targets = load_panel("sv-interest-v1")  # The same entries, for fixture recipes
macrod2 = targets["SV0203"]
print(len(targets), macrod2["kind"])
```

## Where the candidates come from

633 candidates come from files on osteosarc.com or in its data bucket, as of
2026-09-22:

- 523 PASS calls from the oncoanalyser SV tables (ESVEE, PURPLE and LINX)
- 32 DRAGEN SV calls
- 78 long-read RNA fusions from the CTAT tables

The other 4 were added by hand, as named fusions: FOXO3–STRADA–CCDC47,
PARD3B–CDKN2B, GABBR1–SLC29A1 and OTUD7A–FMN1. One of the CTAT fusions,
TPST1–CRCP, is also on that hand-made list.

`catalogue["sources"]` gives the URL and SHA-256 checksum of every downloaded input.
Two entries, the candidate list and the RNA read counts, were derived from those
inputs and have only a checksum.

A candidate is on the list because it's worth investigating. That doesn't mean
it's real, expressed or makes a protein, and no candidate was dropped for lacking
RNA support or a gene. The list includes breaks that join an expressed gene to
an intergenic region, such as MACROD2 (SV0203), ITGBL1 (SV0089) and KMT2C (SV0377).

## Reading an entry

- **Coordinates.** `original_coordinates` are the caller's one-based positions.
  `breakends` are zero-based positions between bases. Single breakends and other
  unclear shapes are `kind="unresolved"`, with the original call kept.
- **Duplicates.** Entries with the same `adjacency_group_id` join the same two
  ends, often reported by different callers. Don't add up their counts.
- **Genes.** Genes come from Ensembl 115. A caller's gene name can be a nearby
  gene rather than one the break falls in.
- **Expression.** `gene_expression` is gene-level TPM from the T2 (January 2025)
  bulk RNA-seq. It says how much the gene is expressed, not the fusion transcript.
- **RNA support.** `rna_evidence` counts reads that cross each breakpoint in
  RNA-seq BAMs from the bucket, one entry per BAM. Counts are reads, not molecules,
  and a missing count means that BAM wasn't checked, not that it had no reads.

To rebuild and rank the proteins these SVs might make, use Isovar and Topiary.

## Use it in fixtures

`load_panel("sv-interest-v1")` gives the same entries as fixture targets. Listing a
target doesn't fetch any reads; write a [fixture recipe](fixtures.md) for that.

```sh
osteosarc --offline fixtures panel sv-interest-v1 > targets.json
```
