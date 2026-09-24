# Shared SV interest catalogue

```python
from osteosarc import load_panel, load_sv_interest

catalogue = load_sv_interest()          # offline; includes source pins
targets = load_panel("sv-interest-v1")  # same targets for fixture recipes
macrod2 = targets["SV0203"]
```

The versioned catalogue retains **637 nominations** from the September 22,
2026 source audit. It includes the displayed ESVEE/PURPLE/LINX PASS table,
focused DRAGEN calls, CTAT long-read calls and focal RNA definitions. It is not
an exhaustive scan of every DRAGEN call or every possible RNA junction.
Membership is a nomination for investigation, not vaccine admission or proof
of translation. No entry is removed for missing RNA, annotation or an ORF.

The expressed-gene-to-intergenic set includes MACROD2 (SV0203), ITGBL1 (SV0089),
ZYG11B (SV0085), RERG (SV0078), PHACTR1 (SV0368), TFDP2 (SV0381), TENM1
(SV0281), GABBR1 (SV0111), IMMT (SV0499), KLF15 (SV0172) and KMT2C (SV0377).
The previously investigated FOXO3, PARD3B, OTUD7A and TPST1 rearrangements and
additional chr20/ATP8B5P/GABBR1 junctions are retained too.

## Coordinates, aliases and evidence

Each target preserves its original one-based reported coordinates and original
VCF alleles where available. Resolved targets provide zero-based interbase
boundaries and the retained side of each breakend. Traversal orientation is
not inferred from those sides: the original oriented hypotheses remain in a
separate field. Single breakends and unresolved geometries use
`kind="unresolved"` and retain their original call instead of fabricating a mate.

`adjacency_group_id` groups identical retained-flank geometry, including reverse
traversals. It does not establish identical inserted alleles or independent
mutations. Original calls retain sample-specific alleles and LINX complex-event
membership. Neither aliases nor alternative ORFs should have their counts added.

Gene annotation is pinned to Ensembl 115. Exact breakends may be intergenic
even when a historical caller gives them a nearby gene name. The absence of an
overlapping coding model is not a reason to exclude a target.

`gene_expression` contains January 2025 **T2 UCLA bulk RSEM** measurements
matched by unique stable Ensembl gene ID, preserving the original versioned ID.
It does not establish expression in T0, T1 or the organoid, or identify a mutant
allele. TPM is a gene/transcript abundance measure, not protein abundance
([RSEM documentation](https://deweylab.github.io/RSEM/rsem-calculate-expression.html)).

`rna_evidence` keeps each product's sample, original source URL, completion
status and separate split, splice-gap and deletion-gap counts. Incomplete
products have missing counts, not zero. The eight completed broad regional
extractions are bounded around nominated loci; they do not exhaustively test
cryptic splicing away from a DNA breakpoint. Counts are scoped templates,
not independent molecules or calibrated protein abundance. Exact DNA-allele
checks with 20-base flanks remain separate from adjacency geometry.

Isovar supplies reconstructed sequences, frames and complete-interval support;
Topiary supplies the report/ranking policy. No protein hypothesis or biological
likelihood is inferred by this data-access interface. Read-orientation flags
describe processed input reads; they do not establish biological RNA strand
or justify automatic rejection of reverse-complement-supported candidates.

## Reproducibility

`catalogue["sources"]` records original source URLs, sizes and SHA-256 pins,
including the DNA VCFs, call tables, Ensembl annotation and full gene-expression
file. Every API call returns a fresh object; caller changes cannot mutate the
packaged snapshot. Existing `sv-regressions-v1` and `vaccine-loci-v1` panels
keep their established contents and policies.

```sh
osteosarc --offline fixtures panel sv-interest-v1 > targets.json
```

Selection and bundle generation still require an explicit source-scoped fixture
recipe. Merely listing a target does not create a BAM or claim that its reads
have been acquired.
