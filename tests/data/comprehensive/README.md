# Comprehensive vaccine / RNA regression corpus

Pinned to osteosarc snapshot **2026-09-28**, ID
`efb65d4b683bda162879c86c0ea889afcc64704b5228224c061109f306503b91`.
Dataset license: CC0-1.0. Code: Apache-2.0.

This corpus contains every vaccine-associated small variant found in any of the
site's three membership sources, longitudinal RNA evidence, DNA evidence and
matched-normal controls, frameshifts, and structural-variant regressions.
`openvax-v3` is a local candidate; it has not been published or made the default.

## Contents

| File | Contents |
| --- | --- |
| `catalog.json.gz` | Evidence catalogue; ordinary gzip-compressed JSON |
| `inputs.json.gz` | Frozen corrected inputs, original annotations, and source checksums; sufficient to rebuild the catalogue offline |
| `spec.json` | Reproducible read-selection specification; carries every library fixture from openvax-v2 |
| `count-source-aliases.json` | Explicit CellRanger count-table aliases and public sample paths |
| `historical-rna.json` | 293 pre-sampling RNA observations from the verified openvax-v2 recipe, with source and assembly identities |
| `checksums.json` | SHA-256 of the frozen inputs and derived catalogue |
| `summary.json` | Checked catalogue totals |

The catalogue has **202 small-variant entries**, including **51 vaccine
variants** and **12 annotated frameshifts**. It retains unresolved entries and
additional count-table alleles. A genomic length change modulo three is recorded
separately: that alone does not establish a coding frameshift.

The website's index marks only 44 variants as vaccinated. Seven further mRNA
targets appear in its vaccine-overlap table: ABI3BP, ANKRD17, ASPM, EXD3, PPP1R3F,
PRRC2C and SPG11. All **38 overlap rows** are preserved, including the two TECPR1
rows. Membership assertions from the index, source-variant flags and overlap
table stay separate; the inclusion policy is their union. Peptides and ELISPOT
annotations remain source assertions, with unknowns intact.

`allele_support` is a complete **202 × 37 = 7,474** matrix of variant × reported
BAM. This includes **15 tumor RNA products across T0–T3**, spanning bulk RNA,
short-read single-cell RNA and ONT, as well as tumor, organoid and normal DNA.
The 21-sample registry is included. `source_reported_support` additionally
preserves the site's depth/rounded-VAF summaries, including longitudinal blood
GEX/CITE-seq. These summaries are not independently recounted raw allele counts.
Samples without a measurement remain without a measurement.

`evidence_overlaps` preserves the reviewed Tempus ALMY/KCV relationship, its raw
regional audit and published duplicate-delivery evidence. The affected KCV count
source links to that group and flags its unresolved biological specimen/timepoint.
The published T1 label is retained as a source assertion. ALMY's reprocessed BAM
is not a second independent RNA product in this count matrix.

`structural_candidates` preserves **all 637 candidates**, their original calls,
breakends, adjacency-group IDs, expression annotations and sample-specific RNA
evidence. The read specification covers **28 structural targets**, including
the seven existing fusion regressions and additional MAP4, TENM1, CNOT3/LENG1,
RUNX1, TGFBR2 and TFDP2 cases. Unresolved breakends remain explicit. Ordinary
splicing, nearby geometry and exact DNA-allele evidence are separate fields;
joining reads alone do not validate a rearrangement or neoantigen.

## RNA-negative SNV controls

Four SNVs have strong DNA evidence and zero reported alternate reads in
**all 15 tumor RNA products**. Their RNA coverage is explicitly distinguished:

| Variant | Example tumor DNA alt/depth | Timepoint | RNA control class |
| --- | --- | --- | --- |
| ACE chr17:63497361 | 12/90, BostonGene WES | T1 | Covered, no ALT |
| BTN3A3 chr6:26448446 | 11/71, Personalis WGS | T0 | Covered, no ALT |
| CTSE chr1:206022947 | 47/381, BostonGene WES | T1 | Insufficient coverage |
| DBH chr9:133652279 | 181/603, BostonGene WES | T0 | Insufficient coverage |

Each comparison requires tumor DNA depth ≥50, alternate reads ≥10 and VAF ≥10%;
normal DNA from the same provider, assay and timepoint with depth ≥30 and VAF
≤1%. The covered controls additionally require RNA reference reads ≥20 at that
timepoint. CTSE and DBH have fewer than 20 reference reads in every reported
RNA product; they test the distinction between absent alternate support and
insufficient coverage. The catalogue records every qualifying comparison and
its BAM identities. These criteria do not establish conclusively somatic status.

The negative-control scope is **tumor RNA**, not every tissue. BTN3A3 has a
separate blood CITE-seq summary (December 2025) reporting VAF 0.0023 at depth
435. This is preserved in `nonzero_rna_summary_exceptions`, without inferring
an integer ALT count or treating the small signal as independently validated.
ACE, CTSE and DBH have no nonzero RNA summary exceptions in the frozen inputs.

**NR2F2 chr15:96332299 is a disagreement case, not a negative control.** The
current count table reports zero RNA ALT, but the historical GRCh37 CeGaT
P116686_3 fixture records **7 alternate templates** before sampling. Its original
allele, assembly, source and observed counts remain in `known_rna_disagreements`.
Historical RNA support is checked before admitting a negative control.
The read audit also verifies preservation of every selected historical NR2F2
record and identifies the current library members containing them. Its seven
historical ALT templates use the original, unfiltered classifier; they are not
claimed to pass the independent MAPQ/NH/base-quality policy below.

The read audit independently recounts these loci from **full regional extracts
before fixture sampling**. It collapses RG/QNAME templates, excludes duplicate,
secondary, supplementary, unmapped and QC-failed records, and requires MAPQ ≥20
and base quality ≥20. MAPQ 255 is accepted only with NH=1; missing base qualities
are unknown. Overlapping mates with conflicting bases count as `other`.
Any alternate reads contradicting the published zeros are reported in
`published_zero_rna_discrepancies`, not hidden by changing thresholds.
The packaged README reports these discrepancies when present; it claims zero
high-quality RNA ALT only when the independent recounts confirm that result.
Both audit paths bind every count source and delivered BAM/index checksum to
the original acquisition receipt and snapshot before attributing its reads.

## Meaning of a support cell

| Status | Meaning |
| --- | --- |
| `alt_observed` | At least one alternate read reported; not a confidence threshold |
| `covered_no_alt` | Measured reference support and zero alternate reads |
| `no_coverage` | Measured depth is zero |
| `other_only` | Coverage exists, but neither REF nor ALT was observed |
| `not_measured` | Missing row, missing count, or counts invalidated by a correction |
| `unresolved_allele` | No usable literal allele; never a negative observation |

Do not sum processing products as independent biological replicates. The
catalogue reports ambiguous BAM basenames instead of guessing a public file.
Neither a rounded VAF of zero nor an absent RNA record proves absent expression.
Published counts retain their original counting units and unspecified filtering;
the independently recounted template counts have their own explicit filters.

## Build and verify

From the repository root, with the pinned snapshot in the cache:

```sh
python scripts/shared_test_data/catalog_comprehensive.py
python scripts/shared_test_data/controls_comprehensive.py build/comprehensive
python scripts/shared_test_data/verify_evidence_comprehensive.py build/comprehensive
python scripts/shared_test_data/build_comprehensive.py build/comprehensive
python scripts/shared_test_data/audit_comprehensive.py build/comprehensive
python scripts/shared_test_data/package_comprehensive.py build/comprehensive
python -m pytest tests/test_corpus.py tests/test_shared.py -q
```

The read build streams bounded regions, never whole remote BAMs. It writes
`build/comprehensive/openvax-v3/`, its `.tar.gz`, `recipe.json` and `build.json`.
It does not write a published release record. The audit writes
`read-audit.json.gz` and `audit-summary.json`. Reruns reuse cached extractions;
`--offline` on the build proves cached inputs suffice.
`reference-contexts.json.gz` contains the reference windows and sequence hashes
used for allele validation, including repeat-expanded indel windows.
The audit also writes `portable-recipe.json` with the public index URLs. This
can rebuild the exact selected records without the website snapshot:
`osteosarc.generate_bundle(recipe, destination, cache=osteosarc.Cache(...))`,
using an online cache for first acquisition. The portable recipe has a different
recipe checksum because it adds index metadata; selected record identities stay
the same. Reading the supplied bundle itself requires no source BAMs or snapshot.
The packaging command writes `osteosarc-comprehensive-2026-09-28.tar.gz`,
containing `evidence/`, `reads/`, `audit/`, reproduction source code and a
per-file checksum manifest. Its checksum is in `dataset-archive.json`.
The evidence and unsampled controls can also be verified and packaged separately:
`python scripts/shared_test_data/package_comprehensive.py build/comprehensive --evidence-only`.
That archive is explicitly labeled `-evidence`; it does not contain or certify the
larger balanced read bundle. Its checksum is in `dataset-evidence-archive.json`.
`controls/` contains unsampled, indexed BAMs at the four negative-control
positions and the NR2F2 disagreement position, with acquisition receipts, so the independent recount can be
repeated from the delivered files without downloading any BAMs.

Every selected record is pinned with its binary BAM identity and multiplicity.
All openvax-v2 library members are checked for exact preservation. Vaccine
coverage is checked against every all-target source, including empty members.
Carry-forward acquisition queries known record starts, including placed-unmapped
records, avoiding the long
introns spanned by historical fixture windows; every required binary record and
its multiplicity must still be recovered.
The new recipe gives the GABBR1/SLC29A1, OTUD7A/FMN1 and PARD3B/CDKN2B fusions
explicit retained sides, per issue #96, without changing published panels.
The comprehensive selection permits up to 2,000 supplementary regions per source
(the expanded ONT selection exceeds the old 100-region limit).
Supplementary-alignment retrieval permits up to 1,000,000 records per queried
region; MAPQ-zero/invalid SA entries, absent contigs and over-limit regions are
excluded and reported by the acquisition log.
Depth probes stream counts without storing a BAM. Their receipts, limits and
source-identity checks remain cached, and verified extracts from earlier attempts
can supply the same counts. This reduces disk use without changing selection.

Read fixtures are deliberately balanced and capped, and are **unsuitable for
abundance or VAF estimation**. Use the published support matrix or pre-sampling
observations for their respective measurements. Use fixture BAMs for testing
read handling, allele extraction, translation and downstream ranking behavior.

```python
import gzip
import json
import osteosarc

with gzip.open("tests/data/comprehensive/catalog.json.gz", "rt") as handle:
    catalogue = json.load(handle)
vaccine_ids = [vid for vid, v in catalogue["variants"].items()
               if v["vaccine"]["included"]]
rna = [r for r in catalogue["allele_support"]
       if r["variant_id"] in vaccine_ids and r["assay_type"] == "RNA"]

bam = osteosarc.bundle_file(
    "build/comprehensive/openvax-v3",
    "BG003082.Aligned.sortedByCoord.out.md.DYNC1H1-chr14-101980529",
)
```

Tests prove exact offline catalogue reproduction, all vaccine joins, complete
matrix coverage, separation of missing/zero/other-only evidence, same-timepoint
negative-control comparisons, and failure on conflicting alleles, duplicate
rows, missing RNA counts or newly positive RNA evidence.
