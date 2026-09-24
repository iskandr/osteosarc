# OpenVax fixture adoption

Isovar, Topiary, Vaxrank and Varcode share Osteosarc's fixture recipes, selection,
bundles and verification instead of maintaining their own. The shared APIs arrived
in Osteosarc 0.2.0 to 0.2.3 (PRs [#22](https://github.com/iskandr/osteosarc/pull/22)
to [#25](https://github.com/iskandr/osteosarc/pull/25)). The consumers adopted them
on 2026-09-23: [Isovar #351](https://github.com/openvax/isovar/pull/351),
[Topiary #380](https://github.com/openvax/topiary/pull/380),
[Vaxrank #503](https://github.com/openvax/vaxrank/pull/503) and
[Varcode #489](https://github.com/openvax/varcode/pull/489). The scientific baseline
stays separate from the mechanics migration.

## Inventory and ownership

| Consumer and fixture family | Shared replacement | Consumer-owned content |
| --- | --- | --- |
| Isovar `sid_data` exact-record pool (311 historical fixtures) | `legacy_fixtures` generate/pack/verify/export, original SAM order and multiplicity | Reviewed recipe, historical source snapshot and test assertions |
| Isovar expansion: vaccine variants, mitochondrial controls, reference/alt/other/uncallable and quality stress | `select_assigned_names`, `write_selected_names`; explicit historical QName-only policy | Independent CIGAR classification, frozen preferred-witness assignments and producer version |
| Isovar short/long SV, three fusions, empty controls | `select_window_segments` and indexed extraction; explicit once-only or preserve-duplicates policy | Breakends/exon windows, observed-alignment compatibility, reconstruction/translation |
| Isovar figure comparisons, stress/performance families | Existing exact-record recipe pool, shared source acquisition | Scientific comparison, stress design, benchmark outputs |
| Topiary compact Sid, indels, regional RNA overlays, all-variant audit, shared NTF3 | `regional_corpus.generate_regional_corpus`, `fixture_assets.fixture_paths` | Selection declaration, annotations, translation/prediction/ranking expectations |
| Vaxrank 58 exact cohorts (49 vaccine-locus retrieval cases plus context/ranking/SV) | `cohort_bundle` pinning, selection, header checks, nested manifests, packing and loading | Checked-in cohort recipe, reference/prediction metadata and ranking tests |
| Vaxrank regional corpus and reporting analyses | Existing `extract_reads` acquisition; source-specific selection declarations | Library grouping, independent quality/orientation audit and reports |
| Varcode historical snapshot and source VCF fixtures | Optional current adapter; existing `Dataset`/`Cache` | Original historical export version, alleles, VCF event definitions and protein tests; no BAM dependency added to ordinary tests |

The historical adapters keep distinct reviewed transport formats, with explicit
SAM-text fidelity where appropriate. There are no runtime imports of Isovar,
Topiary or Vaxrank inside Osteosarc. Their compatibility scripts contain defaults
and CLI plumbing; the moved acquisition/selection/packing/verification loops have
one owner. Topiary's nested manifest and Vaxrank's cohort format remain supported
without requiring either checkout to use a data bundle.

## Common recipe conformance

All three read consumers accept `--panel-recipe recipe.json --panel-source
rna=archive.bam --output NEW_DIRECTORY --offline` through their existing builders:

- `python -m isovar.sid_data generate ...`
- `python -m scripts.generate_sid_fixtures ...` (Topiary)
- `python examples/osteosarc_test_data/build.py --cache CACHE ...` (Vaxrank)

The shared recipe includes exact source identity, targets, required witnesses,
controls/context, explicit caps and frozen evidence-producer provenance. New
variants need declarations, not another BAM writer. See [fixture recipes](fixtures.md).

Run the same source/target/policy through all three actual CLIs, with Python
network connections disabled, and compare full membership, reasons, source and
header identities:

```sh
python -m scripts.check_fixture_consumers --isovar /path/to/isovar \
  --topiary /path/to/topiary --vaxrank /path/to/vaxrank
```

This explicit integration check requires the consumer checkouts. Ordinary
bundle use, consumer installation and tests require no sibling checkout.

## Added SV research panel

`osteosarc/data/additional_sv_recipe.json` freezes SV0461, SV0402, SV0055
(ATP8B5P region), SV0175 (GABBR1 duplication), and SV0499 (IMMT region). The sdist
ships original BAM inputs, indexes, receipts, hypotheses and a regeneration command
in `tests/data/additional_svs/README.md`. Four source archives retain 251 records,
with 24 witness/context members. Exact record identities preserve missing
qualities, original-query intervals, direction labels and shared source scope.
These annotations retain the saved research interpretation; no regenerated
prediction becomes its own expected-result oracle.

The named SV panel also pins the existing PARD3B/CDKN2B, GABBR1/SLC29A1 and
OTUD7A/FMN1 RNA event definitions. Source-call disagreements and complex-cluster
relationships remain in original records. Shared reads and repeated candidate
matches never imply independent molecules.

## Releases and distribution

Each consumer pins its Osteosarc version deliberately and bumps its own version when
that pin changes. Recorded historical data versions and source, correction and
reference identities do not change with the pin. Isovar's offline fixture API also
runs on Python 3.9, and Osteosarc's Python 3.9 suite is tested explicitly. Varcode
keeps the current adapter optional and the historical collector's exact
old-version requirement separate.

Before changing a pin, run each consumer's lint and test gates and inspect real
sdists and wheels for fixture membership and size. Isovar and Vaxrank's
built-distribution regressions exercise original-read loading outside their
checkouts. No fixtures are acquired during package installation. Consumer
scientific rebaselines require separate review, even when acquisition or
regeneration succeeds.
