# Repository comparison and migration

Reviewed 2026-09-18 in `~/code/varcode`, `~/code/isovar`, `~/code/topiary`,
and `~/code/vaxrank`. The requested `varode` and `~/vaxrank` paths were absent;
the existing sibling repositories were used. This change creates the shared
package in `~/code/osteosarc`; it does not modify the four consumers.

## What overlaps today

| Repository | Existing entry points | Scope and differences |
| --- | --- | --- |
| Isovar | `tests/data/osteosarc/fetch_sources.py`, `rebuild.py`, `expansion/inventory.py`, `discover.py`, `acquire.py`, `metadata.py` | Small selected fixtures plus the broadest inventory/acquisition implementation. Separate immutable receipts, source/header surveys, RNA-oriented bucket classification, explicit GRCh37 handling, variant-page metadata. Initial fixture selection includes allele/quality-blind template sampling **and** selected low-quality alternate-supporting templates. |
| Topiary | `scripts/osteosarc_variant_audit.py`, `osteosarc_rna_overlay.py`, `tests/osteosarc_helpers.py`, `tests/data/pvacseq/osteosarc/regenerate.py` | Another curl/SHA256 snapshot function and nine-cell HTML variant parser; all site variants, selected RNA products, reference subsets, pVAC row selection, expression overlays. Keeps unresolved catalog entries. |
| Vaxrank | `examples/osteosarc_read_corpus/build.py`, `tests/osteosarc_fixture_support.py`, `tests/osteosarc_helpers.py` | Another downloader with resumable curl transfers, read identities, canonical-allele trimming, sample/path inference, DNA + RNA corpus acquisition and receipt checks. Several test assets are copied from Isovar. |
| Varcode | `tests/test_osteosarc_fusions.py`, `output/osteosarc-phasing/*/scripts/`, `output/osteosarc-phasing/*/inputs/generator-source/` | Curated structural-variant and phasing fixtures; analysis scripts import Isovar's test helpers through a hardcoded repository path or retain copied generator modules. Annotation references vary by task. |

These are concrete duplications, not just similarly named helpers:

* Isovar `expansion/inventory.py:fetch_snapshot`, Topiary
  `osteosarc_variant_audit.py:fetch_snapshot`, and Vaxrank
  `osteosarc_read_corpus/build.py:download` all fetch the same website resources
  with curl and SHA256 receipts, using different file layouts and retry policies.
* Isovar and Topiary duplicate the HTML variant table parser. Their initial
  selection differs: vaccine targets versus all website entries. Varcode's
  analysis scripts import the Isovar version and redo joins to VAF alleles.
* All four construct local reference/fixture subsets and verify hashes, but
  reference versions and biological expected results belong to their individual
  tests. Sharing acquisition does not justify changing those expected results.
* Read acquisition differs in index discovery, interval padding, record identity,
  filtering, failure reporting, and whether downstream code is reading a sampled
  fixture or an uncapped regional source.

## Concepts made explicit

| Concept | Shared representation |
| --- | --- |
| Website/S3 metadata at an acquisition date | Named `Dataset` with exact source receipts |
| One original file or processing product | `Asset`, identified by complete URL |
| Sample, timepoint, assay, provider assertion | `SampleClaim` with source, basis, and retained conflicts |
| Shared local download | `Cache` object bytes + SHA256 receipt; optional published MD5 verification |
| File/table subset | `Assets.select`, `Table.select/where` |
| Website entry vs literal genomic allele | `Variant.id`, `alleles`, `status`; source annotations retained |
| Site/vaccine/pipeline variant selection | `Dataset.variants` + explicit source of vaccine membership |
| Coordinate interval | `Region`, zero-based half-open with required assembly |
| Uncapped regional alignment records | `extract_reads`, default `ReadFilter()` |
| Small deterministic regression fixture | `subset_templates`, separate from acquisition |
| Variant effects and transcript choice | Explicit Varcode/PyEnsembl consumer code |
| RNA assemblies, epitopes, ranking | Existing Isovar, Topiary, Vaxrank APIs |

The full inventory contains raw reads, multiple alignment products, VCFs,
pipeline reports, expression matrices, and other files. Classification is
conservative. An unclassified file remains accessible by key or URL; unknown
metadata is not replaced with a guessed platform or sample identity.

## Consumer changes to make next

1. **Replace acquisition helpers first.** Add an `osteosarc` dependency (or test
   extra where acquisition is test-only), use the shared OpenVax cache
   (`OPENVAX_DATA_CACHE`, the layout vaxrank already uses) and a named
   snapshot, and replace URL constants/download helpers with `Dataset` and
   `Cache`. Keep current fixture bytes and expected biological outputs.
2. **Replace private test imports.** Varcode's analysis scripts should use the
   public package rather than `sys.path` edits into Isovar's `tests` package.
   Use `Dataset.variants` instead of copying the HTML/VAF join.
3. **Adopt one regional extractor.** Pass explicit asset IDs/keys, regions, and
   filters. Compare complete SAM record multisets and tag values against each
   old acquisition before deleting it. Matching read counts alone is insufficient.
4. **Separate selection recipes.** Keep Isovar's stress selection and Topiary's
   pVAC feature/row selection as named consumer fixture recipes. The shared
   allele-blind template sampler does not reproduce the existing low-quality
   alternate-read enrichment. Record its policy separately if migrating it.
5. **Keep biological choices in consumers.** Preserve reference release,
   transcript selection, allele normalization, RNA policy, and scoring behavior.
   The current fixture references include Ensembl 87 and 95 for different
   purposes; do not silently select one release for all projects.

Existing cached files can be adopted without redownloading:

See [consumer recipes](consumers.md) for the native API calls and usage-test
coverage. Custom subset genomes keep their unique reference names and use
`selected.to_varcode(genome=genome, assembly="GRCh38")`. Header surveys use
`data.inspect_alignment(asset)` before the consumer chooses coordinates.

```python
from osteosarc import SNAPSHOT_SOURCES, Cache, Dataset, digest

# Stand-in for a file in an older project cache: the snapshot's own copy.
old_path = Dataset.open("baseline").source_path("vafs")
original_url = SNAPSHOT_SOURCES["vafs"]
receipt = Cache().import_file(old_path, original_url, sha256=digest(old_path))
print(receipt.sha256, receipt.size)
```

Imports record the local import time, not a fabricated original download time.
Keep original manifests when historical acquisition dates matter.

## Known source limitations retained in the API

The reviewed source metadata disagree about some BostonGene/UCLA library
timepoints, and the viewer's global hg38 label covers files that require
individual reference checks. Natera WGS tumor/normal assignments are described
as inferred in the public data page. Vaccine-overlap and variant-source JSON
also differ in their membership flags. The package retains these distinctions.
Verified corrections are optional, centralized in `osteosarc/curation.py`, and
re-checked against every snapshot (see [corrections](curation.md)); they do not
declare a validated truth set.

`ready` checks internal literal-allele availability. Independent reference-allele
validation, indel equivalence, liftover, donor demultiplexing, and biological
replicate resolution remain separate analyses. Byte-equivalent aliases are not
automatically collapsed across URLs, and metadata completeness is limited to
what the public sources actually publish.
