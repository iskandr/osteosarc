# Native retrieval coordinate regressions

Four unchanged BAM/index pairs and their original case metadata from the Vaxrank
cohort archive at commit `0565cf15754fa3ac45b82fa5fc7fd5ddfbe06cfa`, path
`vaxrank/data/sid-test-data.zip`, members `osteosarc/shared-v1/*`.
The manifest retains each original asset's SHA-256, source URL, selection and
reference-validation receipt. The public source reads are CC0.

Cases 44, 45 and 48 exercise Osteosarc #99: plain `MT` on a GRCh37 nuclear
reference with rCRS mitochondria, UCSC hg19 chrM with its different length and
position, and a GRCh37 NR2F2 locus whose catalogue ID is GRCh38. Case 00 is an
ordinary GRCh38 control. No read or historical coordinate has been edited.

`osteosarc/data/retrieval-cases.json` separately versions the native metadata for
all 49 cases. It was generated with `annotate_retrieval_manifest` against all
original assets: all checksums and all 49 native indexed queries passed. Original
cases and BAM identities remain unchanged; the new fields supply exact query
contigs, coordinates and reference identities. Do not merge all cases as one BAM.
