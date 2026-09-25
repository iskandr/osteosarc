# Read fixtures and bundles

A fixture recipe describes a small test BAM: which reads to take from which files,
and why. Osteosarc turns the recipe into a bundle that anyone can rebuild and
check offline. Isovar, Topiary and Vaxrank build their test data this way.

A recipe has three parts:

- **Targets**: what each fixture is about: a small variant, a structural variant
  (SV), or an entry marked unresolved.
- **Sources**: the BAMs that reads come from, with their sample and assembly.
- **Members**: one target and one source, plus a rule for picking reads.

Osteosarc only follows the rule. Deciding which reads support an allele is up to
the library that uses the fixture, and a fixture is never an estimate of VAF.

## Build a bundle

This recipe keeps up to 25 reads (with their mates) around DYNC1H1 from a T0
tumor RNA-seq BAM:

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

`truncated` means more reads overlapped the window than the cap allowed. Mates
are kept together, so the record count can be larger than the cap. The
destination directory must not exist yet.

Check, list and export the bundle offline:

```sh
osteosarc fixtures verify dync1h1-bundle
osteosarc fixtures list dync1h1-bundle
osteosarc fixtures export dync1h1-bundle dync1h1-exported --member DYNC1H1-rna
```

To run a recipe on BAMs you already have, use `select_fixtures(recipe, sources)`
or the CLI:

<!-- docs-check: skip (needs your own local BAM) -->
```sh
osteosarc --offline fixtures select recipe.json --source rna=archive.bam
```

## Recipe fields

A recipe is a dictionary with `schema_version: 1`, an `id`, and `targets`,
`sources` and `members`. `validate_recipe(recipe)` checks it before any reads are
touched. `select_fixtures` and `fixtures select` only read BAMs you already have
locally; `generate_bundle` and `fixtures generate` fetch remote ones for you, with
the same indexed extraction as [`extract_reads`](reads.md).

**Targets.** A small variant has `kind: small_variant`, `assembly`, `reference`
(where the target came from), `coordinates: one-based`, `contig`, `position`,
`ref` and `alt`. An SV has `kind: sv`, `coordinates: zero-based-interbase`, and two
or more `breakends` with `contig`, `position` and `orientation` (`+`, `-` or null).
An `unresolved` target needs a `reason` and gets no reads.

**Sources.** Each source has an `identity` object that says which file it is,
such as `{"key": "rna-seq/..."}` or `{"url": "https://..."}`, plus `assembly`,
`sample`, `library` and `product`; use null for anything unknown. Two processed
versions of one library are separate sources. `archive_sha256` pins a local file.

**Members.** Each member names a `target` and `source` and has a `policy` with
`version: 1` and one of these kinds:

| `kind` | Picks |
| --- | --- |
| `regional` | Reads overlapping `regions`, optionally capped |
| `exact` | Specific records, by checksum and count |
| `witnesses` | Named reads (read group, read name, mate), each with a reason |
| `stratified` | Named reads plus a reproducible sample of others, by group |
| `empty` | Nothing, on purpose |
| `omitted` | Nothing, with a `reason` |

`regions` use `Region` fields (`contig`, `start`, `end`, `assembly`), zero-based
and half-open. Reads in `context_regions` are kept as context, never as support.

**Named reads and caps.** A witness or stratum `assignment` has a `selector`
(`rg`, `qname`, and optionally `segment`, the read's FLAG & 0xc0), a `reason`,
and the `producer` name and version that chose it. Required reads (the default)
must be present, or selection fails. `cap` limits the optional reads per stratum;
`strata` sets per-stratum caps. Sampling is reproducible from `seed` (default
`"0"`) and doesn't depend on read order.

**Duplicates.** Repeated identical records are kept as repeats unless
`duplicate_policy` is `identical-record-once`. Members that share records
shouldn't be counted as independent evidence.

Each member ends up `selected`, `truncated`, `empty`, `unresolved` or `omitted`,
with the reason every record was picked.

## Record identity

`record_multiset(path)` fingerprints every record of a BAM from its stored bytes
(the `bam-record-v1` encoding), so two BAMs compare equal only if they hold the
same records the same number of times, regardless of compression or order. An
`exact` recipe can instead use `encoding: sam-text-v1` to match older SAM-text
checksums; that's weaker, because SAM text can hide differences in tag types.

## Named panels

```sh
osteosarc fixtures panel vaccine-loci-v1
```

`load_panel(name)` returns a set of shipped targets:

| Panel | Contents |
| --- | --- |
| `vaccine-loci-v1` | The vaccine target alleles used in earlier fixtures |
| `sv-regressions-v1` | RNA fusion events and five candidate SVs used in regression tests |
| `sv-interest-v1` | All 637 entries of the [SV catalogue](sv-interest.md) |

A panel only lists targets; a recipe decides which reads to fetch.

## Keep mates and split reads

A member with `retain_partners: true` also keeps the mates and split alignments of
its reads, if its source was fetched with a
[`RecoveryPolicy`](reads.md#recover-mates-and-split-alignments). A source's
`acquisition` can set `{"recovery": {"on_timeout": "incomplete"}}` to keep going
when a partner query times out. The bundle then records the acquisition as
`incomplete`. Zero reads from a limited or incomplete fetch isn't evidence of zero
support.

## Bundles

`generate_bundle(recipe, destination, sources=..., cache=...)` fetches the reads,
applies the recipe and writes a self-contained directory. `pack_bundle(selection,
destination)` starts from a selection you already made. `Dataset.generate_bundle`
also checks each source against the snapshot.

**Fetching.** A source can point to a small pinned `archive` (`url`, `sha256`,
`size_bytes`), or to a BAM by `identity.url` with an `index`, fetched around the
members' regions. `acquisition` passes read filters and recovery settings.
Installing the package never downloads data.

**Contents.** A bundle holds each source's records once, indexed BAMs, the original
headers, the recipe, fetch receipts, file checksums, and every member's records and
reasons. The recipe's `redistribution` field carries license and citation notes.

**Headers.** `header_policy="full"` (the default) keeps the source headers.
`"compact"` keeps only what's needed: every sequence line, read groups, and the
programs that produced the reads. The full header is archived either way.

**Checking.** `verify_bundle(directory, sha256=...)` checks every file, record and
index offline. Pass the manifest's SHA-256 when using someone else's bundle;
without it, the check shows the bundle is intact but not who made it.

**Exporting.** `export_bundle(bundle, destination, members=[...])` writes one sorted,
indexed BAM per member (or SAM with `format="sam"`). Empty members become valid
empty BAMs. The default size limit is 64 MiB; set `size_budget` for a smaller
package. Existing destinations are never overwritten.

## Tests as examples

`tests/test_fixtures.py` and `tests/test_bundles.py` build recipes and BAMs without a
network connection, and cover each rule above.
