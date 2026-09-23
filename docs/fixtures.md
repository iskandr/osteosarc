# Reproducible read fixtures

Fixture recipes freeze a target, source and selection policy. They do not infer
allele/junction/protein support. Classification remains in the evidence-producing
library; its versioned assignments become recipe inputs. Selection does not use
the current predicted protein and is never a full-source abundance/VAF estimate.

## Recipe v1

`select_fixtures(recipe, sources)` accepts a JSON-compatible dictionary with
`schema_version: 1`, a stable `id`, and `targets`, `sources`, `members` mappings.
`Dataset.select_fixtures` runs the identical implementation. `sources` maps source
IDs to explicit local BAM paths or acquired `ReadSubset` objects; remote BAMs
must first go through bounded `extract_reads`. Local archives can be pinned with
`archive_sha256` in the recipe's source declaration.

A small-variant target declares `kind: small_variant`, `assembly`, pinned
`reference` identity, `coordinates: one-based`, `contig`, `position`, `ref`, `alt`.
An SV declares `kind: sv`, `coordinates: zero-based-interbase`, and at least two
`breakends` with `contig`, `position`, `orientation` (`+`, `-`, or null if unknown).
Retain inserted sequence and original source-call/annotation provenance in the
target. `kind: unresolved` requires a reason and does not acquire/select reads.

A source declares `identity` (original asset/snapshot/checksum metadata),
`assembly`, `sample`, `library`, and `product`. Unknown scope is explicitly null.
Different processing products remain separate sources even if they share a
library. CB/UMI labels never merge templates into inferred molecules.

Each member names a `target` and `source`, and a `policy` with `version: 1`:

| Kind | Selection |
| --- | --- |
| `regional` | Candidates overlapping explicit `regions`, optionally template-capped |
| `exact` | Required `records` mapping from record digest to positive multiplicity |
| `witnesses` | Original source/RG/QNAME/segment selectors with pinned reasons |
| `stratified` | Required witnesses plus deterministically sampled optional strata |
| `empty` | Deliberate empty fixture, distinct from unavailable acquisition |
| `omitted` | Explicit omission with `reason` |

`regions` and `context_regions` contain `Region` fields: `contig`, `start`, `end`,
`assembly` (zero-based, half-open). Context records are retained with the reason
`assembly context`, never promoted to junction support by overlap alone.

Each witness/stratum `assignment` contains a `selector` (`rg`, `qname`, optional
`segment` = original FLAG & 0xc0), a `reason`, and `producer` name/version.
`required` defaults true; absence fails. Optional assignments can specify
`stratum`. `cap` is an optional-template budget **per stratum**; `strata` maps
stratum names to overrides. Required witnesses bypass caps. All records of a
sampled template available in the input survive. `seed` defaults to the string
`0`. Hash ordering makes sampling independent of input order. Missing required
record occurrences fail instead of silently weakening a pinned regression.

Duplicate multiplicity is preserved by default. `duplicate_policy:
identical-record-once` explicitly requests the legacy deduplication behavior.
Shared fixture members reference a source record; their counts must not be
summed as independent evidence. Results distinguish selected, truncated, empty,
unresolved and omitted. Selection returns a manifest of per-record reasons and
multiplicities, along with in-memory original records for export.

## Identity and fidelity

`record_multiset(path)` uses `bam-record-v1`: stored CIGAR, sequence, qualities,
flags and typed auxiliary bytes, with reference names replacing numeric IDs and
the derived bin omitted. Tag order is significant. Float payloads and integer
widths survive; compression and coordinate-order ties do not affect equality.
The encoding follows [SAM/BAM §4.2](https://samtools.github.io/hts-specs/SAMv1.pdf).

An exact recipe may explicitly use `encoding: sam-text-v1` for historical SAM
checksums. This cannot prove bitwise tag fidelity and rejects a text identity
that ambiguously maps to different binary records. It is not the default.

## Named panels and CLI

`load_panel("vaccine-loci-v1")` supplies the pinned historical vaccine alleles;
consumers retain their own reference/correction policies. `sv-regressions-v1`
lists shared SV targets with unresolved states until exact reviewed event
selectors are supplied. A target name alone is never an acquisition selector.

```sh
osteosarc fixtures panel vaccine-loci-v1
osteosarc --offline fixtures select recipe.json --source rna=archive.bam
```

The CLI prints the same membership and reasons as the standalone and Dataset
APIs. Constructors in `tests/conftest.py` and `tests/test_fixtures.py` exercise
source/RG collisions, duplicates, required controls, contexts, empty selections,
missing witnesses, float precision and API/CLI conformance without network.

## Bounded partner recovery

```python
from osteosarc import RecoveryPolicy, extract_reads

subset = extract_reads(source, seed_and_context_regions, cache=cache,
                       recovery=RecoveryPolicy(max_rounds=4, max_intervals=128,
                                               max_bases=1_000_000,
                                               max_records=100_000))
```

The seed union preserves every available record. Distant mate/SA windows retain
only matching source/RG/segment observations; SA requires matching position,
strand, CIGAR and MAPQ, plus NM when present. Hard clipping, absent SEQ/QUAL,
secondary/supplementary flags and tags are unchanged. Duplicate counts use the
maximum occurrence count across indexed queries, never a set or a sum of repeated
retrievals. No query synthesizes a missing partner or reverse-complement read.

The receipt records visited intervals, requested leads, missing/conflicting
partners, malformed SA, missing sequence, ambiguous placement, repeated/cyclic
leads, and acquisition limits. `complete_template` is always false. Even a fully
resolved SA graph may omit unreported alignments. Missing indexes fail without a
scan fallback. Interval/round limits produce an explicit truncated receipt;
record overflow fails instead of publishing a partial successful artifact.
Changed sources/headers/indexes fail across steps. Interrupted runs reuse only
verified intermediate extractions.

A member can set `retain_partners: true` to retain recovered mates/split records
for its selected templates, carrying the recovery reasons. `acquisition_status`
is separate from selection status: zero retained records in bounded or truncated
input is not evidence of zero support in the source. CLI extraction also supports
`reads --recover-linked`. The existing `fetch_pairs` behavior remains unchanged;
choose a recovery policy instead of combining those options.
