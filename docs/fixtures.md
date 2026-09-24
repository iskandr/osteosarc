# Read fixtures and bundles

Fixture recipes turn a few reads from the dataset into small, versioned test BAMs
that any project can regenerate and verify offline. Isovar, Topiary and Vaxrank
build their Osteosarc test data this way.

A recipe declares three things:

- **Targets**: what each fixture is about. A target is a small variant, a
  structural variant (SV), or an explicitly unresolved entry.
- **Sources**: the alignments that reads come from, with their identity,
  assembly, sample, library and processing product.
- **Members**: pairs of one target and one source, each with a selection policy.

Osteosarc executes the policy. It doesn't infer allele, junction or protein
support: that classification stays in the evidence-producing library, and its
versioned assignments become recipe inputs. Selection doesn't use the current
predicted protein, and a fixture is never a full-source abundance or VAF estimate.

## Build a bundle

This recipe keeps up to 25 read templates around DYNC1H1 from a T0 tumor RNA-seq
alignment. `Dataset.generate_bundle` fetches the window with indexed extraction,
selects the records, and publishes a self-contained directory:

```python
import json
from pathlib import Path

from osteosarc import Dataset, verify_bundle

data = Dataset.open(offline=False)
source = data.asset("rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam")
variant = data.variants()["DYNC1H1-chr14-101980529"]
chrom, position, ref, alt = variant.allele
window = variant.region(padding=100)

recipe = {
    "schema_version": 1,
    "id": "dync1h1-rna-example",
    "targets": {
        "DYNC1H1": {
            "kind": "small_variant", "assembly": "GRCh38", "coordinates": "one-based",
            "contig": chrom, "position": position, "ref": ref, "alt": alt,
            "reference": {"source": "osteosarc", "snapshot_id": data.id, "variant_id": variant.id},
        },
    },
    "sources": {
        "rna": {
            "identity": {"key": source.key}, "assembly": "GRCh38",
            "sample": "T0_tumor", "library": "BG003082", "product": source.key,
        },
    },
    "members": {
        "DYNC1H1-rna": {
            "target": "DYNC1H1", "source": "rna",
            "regions": [{"contig": window.contig, "start": window.start,
                         "end": window.end, "assembly": window.assembly}],
            "policy": {"version": 1, "kind": "regional", "cap": 25, "seed": "example"},
        },
    },
}
Path("recipe.json").write_text(json.dumps(recipe, indent=2))
data.generate_bundle(recipe, "dync1h1-bundle")
member = verify_bundle("dync1h1-bundle")["members"]["DYNC1H1-rna"]
print(member["status"], member["record_count"])
```

A regional member reports `truncated` when more templates overlapped the window than
its cap allowed. Every record of a kept template is retained, so the record count
can exceed the template cap. The destination must not exist yet.

Check, list and export the bundle offline:

```sh
osteosarc fixtures verify dync1h1-bundle
osteosarc fixtures list dync1h1-bundle
osteosarc fixtures export dync1h1-bundle dync1h1-exported --member DYNC1H1-rna
```

`select_fixtures(recipe, sources)` runs the same selection on BAMs you already
have, without acquisition. The CLI prints the same membership and reasons as the
Python API:

<!-- docs-check: skip (needs your own local BAM) -->
```sh
osteosarc --offline fixtures select recipe.json --source rna=archive.bam
```

## Recipe reference (v1)

A recipe is a JSON-compatible dictionary with `schema_version: 1`, a stable `id`,
and `targets`, `sources` and `members` mappings. `validate_recipe(recipe)` checks
it before any source is read. `select_fixtures` and `Dataset.select_fixtures` run
the identical implementation. Their `sources` argument maps source IDs to explicit
local BAM paths or acquired `ReadSubset` objects. Remote BAMs must first go
through bounded `extract_reads`.

### Targets

A small variant declares `kind: small_variant`, `assembly`, a pinned `reference`
identity, `coordinates: one-based`, `contig`, `position`, `ref` and `alt`.

An SV declares `kind: sv`, `coordinates: zero-based-interbase`, and at least two
`breakends` with `contig`, `position` and `orientation` (`+`, `-`, or null if
unknown). Retain the inserted sequence and the original source-call and
annotation provenance in the target.

`kind: unresolved` requires a `reason` and does not acquire or select reads.

### Sources

A source declares `identity` (original asset, snapshot and checksum metadata),
`assembly`, `sample`, `library` and `product`. Unknown scope is explicitly null.
Different processing products remain separate sources even if they share a
library. CB/UMI labels never merge templates into inferred molecules. Pin a local
archive with `archive_sha256`.

### Members and selection policies

Each member names a `target` and a `source`, and a `policy` with `version: 1`:

| Policy `kind` | Selection |
| --- | --- |
| `regional` | Candidates overlapping explicit `regions`, optionally template-capped |
| `exact` | Required `records` mapping from record digest to positive multiplicity |
| `witnesses` | Original source/RG/QNAME/segment selectors with pinned reasons |
| `stratified` | Required witnesses plus deterministically sampled optional strata |
| `empty` | Deliberate empty fixture, distinct from unavailable acquisition |
| `omitted` | Explicit omission with `reason` |

`regions` and `context_regions` contain `Region` fields: `contig`, `start`, `end`
and `assembly`, zero-based and half-open. Context records are retained with the
reason `assembly context`. Overlap alone never promotes them to junction support.

### Witnesses, strata and caps

Each witness or stratum `assignment` contains a `selector` (`rg`, `qname`, and an
optional `segment` equal to the original FLAG & 0xc0), a `reason`, and a
`producer` name and version. `required` defaults to true, and a missing required
witness fails. Optional assignments can specify a `stratum`.

`cap` is an optional-template budget **per stratum**; `strata` maps stratum names
to overrides. Required witnesses bypass caps. All records of a sampled template
that are available in the input survive. `seed` defaults to the string `0`.
Hash ordering makes sampling independent of input order. Missing required record
occurrences fail instead of silently weakening a pinned regression.

### Duplicates and results

Duplicate multiplicity is preserved by default. `duplicate_policy:
identical-record-once` explicitly requests the legacy deduplication behavior.
Shared fixture members reference a source record; their counts must not be
summed as independent evidence.

Results distinguish `selected`, `truncated`, `empty`, `unresolved` and `omitted`.
Selection returns a manifest of per-record reasons and multiplicities, along with
the original records in memory for export.

## Record identity

`record_multiset(path)` uses the `bam-record-v1` encoding: stored CIGAR, sequence,
qualities, flags and typed auxiliary bytes, with reference names replacing numeric
IDs and the derived bin omitted. Tag order is significant. Float payloads and
integer widths survive; compression and coordinate-order ties do not affect
equality. The encoding follows [SAM/BAM §4.2](https://samtools.github.io/hts-specs/SAMv1.pdf).

An exact recipe may explicitly use `encoding: sam-text-v1` for historical SAM
checksums. This cannot prove bitwise tag fidelity and rejects a text identity
that ambiguously maps to different binary records. It is not the default.

## Named panels

```sh
osteosarc fixtures panel vaccine-loci-v1
```

`load_panel(name)` returns a shipped panel of targets. A target name alone is
never an acquisition selector.

| Panel | Contents |
| --- | --- |
| `vaccine-loci-v1` | The pinned historical vaccine alleles. Consumers keep their own reference and correction policies. |
| `sv-regressions-v1` | Historical RNA events and the additional 2026-09-23 research panel, with original VCF anchors beside explicitly converted interbase boundaries |
| `sv-interest-v1` | All 637 nominations from the September 22, 2026 SV audit, a broader discovery set than the regression panel. See the [SV interest catalogue](sv-interest.md). |

## Keep mates and split reads

A member can set `retain_partners: true` to keep recovered mates and split
records for its selected templates, carrying the recovery reasons. This needs a
source acquired with a [`RecoveryPolicy`](reads.md#recover-mates-and-split-alignments).
`acquisition_status` is separate from selection status: zero retained records in
bounded, truncated or incomplete input is not evidence of zero support in the
source. A source's `acquisition` can set `{"recovery": {"on_timeout": "incomplete"}}`
so that a partner-query timeout yields an `incomplete` acquisition instead of
failing the bundle.

## Portable bundles

`generate_bundle(recipe, destination, sources=..., cache=...)` acquires declared
inputs, selects, and publishes a self-contained directory. `pack_bundle(selection,
destination)` starts from an existing selection. Destinations must be new.
`Dataset.generate_bundle` also checks source identities against its snapshot.

### Acquisition

A source may declare a small historical `archive` with `url`, `sha256`, and
`size_bytes`, or `identity.url`, `index` and bounded `regions` for live indexed
acquisition. Member regions are the fallback acquisition union. Source
`acquisition` holds explicit filters and recovery settings. Explicit local
`sources` are verified against any declared archive hash. Interrupted extraction
can reuse the existing verified cache derivatives. Normal installation never
acquires data. Indexed acquisition preserves declared inventory size and
modification metadata; the Dataset API also accepts an identity containing only
a snapshot key or ID.

### Contents

The bundle stores a shared source record pool, indexed BAMs, full original headers,
the recipe, acquisition receipts, source/sample/library/product identities, member
multiplicity and reasons, tool versions, parent lineage, and file hashes and sizes.
`redistribution` in the recipe carries source license and citation information;
absent license information remains unresolved. Remote HTTP identity is preserved as
HTTP evidence, separate from archive or full-file SHA-256. Historical receipt paths
are audit strings, never dependencies for offline verification or export.

### Headers

`header_policy="full"` is the default. Compact mode keeps every SQ (including
assembly-identifying contigs), source comments, retained RG/SM/LB metadata, and
required PG ancestry. When producer lineage is unresolved it retains all PGs.
Neither mode invents missing source metadata. Sorting only updates HD sort fields;
the full original header remains archived. This avoids weakening assembly guards
for downstream offline use.

### Verification

`verify_bundle(directory, sha256=pinned_manifest_hash)` checks all files, the
recipe digest, source and member record multisets, and index enumeration. Pin the
manifest hash when consuming an external release; internal consistency checks
alone are not an authenticity signature. `record_multiset` separately supports
lossless equivalence comparisons across compression and tool versions. Exact BAM
byte reproducibility requires the same recorded pysam/HTSlib toolchain.

### Export

`export_bundle(bundle, destination, members=[...])` adds named, coordinate-sorted,
indexed BAM exports, retaining the self-contained source pool and provenance.
`format="sam"` or `"sam.gz"` explicitly requests legacy SAM-text fidelity. Empty
members export valid empty indexed BAMs. Unresolved and omitted members stay
declared without fabricated data. The default size budget is 64 MiB including
metadata; set a smaller `size_budget` for a consumer's package. Publication is
atomic only after verification, and existing destinations are refused.

Exporting an exported bundle is supported: each request replaces the named export
set in the new destination, so changing format leaves no stale export files.
SAM exports use coordinate order, consistent with their headers. Verification
checks member counts and status and SAM field multisets as well as BAM identities.

## Tests as examples

The constructors in `tests/conftest.py`, `tests/test_fixtures.py` and
`tests/test_bundles.py` build recipes and BAMs without network access. They cover
source and read-group collisions, duplicates, required controls, contexts, empty
selections, missing witnesses, float precision, API/CLI conformance, and
generation, packing, export and verification in a fresh offline directory. Corrupt
records, missing duplicates, nested members, recipe changes, swapped indexes,
unsafe paths and size-budget failures are separate regressions.
