# Test data

Unit tests need small, real sets of reads that anyone can rebuild. Osteosarc makes
them two ways: a quick BAM of the reads around some variants, or a bundle built from
a recipe, which pins every record so rebuilding gives the same bytes.

## A quick test BAM

```sh
osteosarc reads rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam --variant DYNC1H1-chr14-101980529 --padding 100 --to test-data
```

This streams just those reads from the public BAM and saves them, with an index, as
test-data/BG003082.DYNC1H1-chr14-101980529.bam. Give a sample ID instead of a file to
do the same for each of that sample's BAMs:

<!-- docs-check: skip (streams reads from several BAMs) -->
```sh
osteosarc reads T0_tumor --assay rna-seq --variant DYNC1H1-chr14-101980529 --padding 100 --to test-data
```

In Python, `data.extract_reads(file, variants=..., to="test-data")` does the same. See
[reads](reads.md) for regions, filters and mates.

## A bundle from a recipe

A recipe says which reads to take from which BAMs, and why. Osteosarc turns it into a
bundle: indexed BAMs plus a manifest listing every record and the reason it was kept.
Anyone can rebuild the bundle and check it offline. This recipe keeps up to 25 reads,
with their mates, around one variant in a T0 tumor RNA-seq BAM:

```python
import json
from pathlib import Path

from osteosarc import Dataset, generate_bundle, verify_bundle

data = Dataset.open(offline=False)
source = data.file("rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam")
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
generate_bundle(recipe, "dync1h1-bundle", dataset=data)
member = verify_bundle("dync1h1-bundle")["members"]["DYNC1H1-rna"]
print(member["status"], member["record_count"])
```

The status is truncated when more reads overlapped than the cap allowed. Mates stay
together, so there can be more records than the cap. Then, offline:

```sh
osteosarc test-data list dync1h1-bundle
osteosarc test-data verify dync1h1-bundle
osteosarc test-data export dync1h1-bundle dync1h1-exported --member DYNC1H1-rna
```

The same recipe runs from the command line, fetching the reads, or on BAMs you already
have:

<!-- docs-check: skip (needs your own BAM) -->
```sh
osteosarc test-data generate recipe.json bundle
osteosarc --offline test-data generate recipe.json bundle --source rna=archive.bam
```

## Recipes

A recipe has targets (what each piece of test data is about), sources (the BAMs its
reads come from) and members (one target in one source, with a rule for picking
reads). Osteosarc only follows the rule: deciding which reads support an allele is up
to the library that uses them, and a member's reads are never an estimate of VAF.

```text
targets:   small_variant   one-based position, ref and alt, and where it came from
           sv              two or more breakends (zero-based, with orientation + / - / null)
           unresolved      a reason; gets no reads
sources:   identity        which file: {"key": ...} or {"url": ...}
           assembly, sample, library, product   (null when unknown)
members:   target, source, regions (zero-based, half-open), and a policy:
           regional        reads overlapping the regions, optionally capped
           exact           specific records, by checksum and count
           witnesses       named reads, each with a reason
           stratified      named reads, plus a reproducible sample of others by group
           empty, omitted  nothing, on purpose (omitted says why)
```

Named reads give their read group and name (and optionally which mate), a reason, and
the tool that chose them; they must be present unless marked optional. Caps sample the
optional reads reproducibly from the policy's seed, whatever the read order. Identical
records keep their repeats. `validate_recipe(recipe)` checks a recipe before any reads
are touched.

Each member ends up selected, truncated, empty, unresolved or omitted, and the
manifest gives the reason each record was kept.

**Mates and split reads.** A member can also keep the mates and split alignments of its
reads, when its source was fetched with mate and split-read recovery
([reads](reads.md#mates-and-split-reads)). A recovery that timed out is
recorded as incomplete: finding no reads then isn't evidence that there are none.

**Records.** Bundles identify each record by a checksum of its stored bytes, so two BAMs
match only when they hold the same records the same number of times, whatever their
compression or order. Older checksums of SAM text still work in exact members, but can
miss differences in tag types.

## Bundles

A bundle holds each source's records once, as indexed BAMs, with the original headers,
the recipe, what was fetched, file checksums, and every member's records and reasons.
Verifying checks every file, record and index offline; give the manifest's checksum
when using someone else's bundle, since without it the check shows the bundle is
intact but not who made it. Exporting writes one sorted, indexed BAM per member (or SAM),
and empty members become valid empty BAMs. The size limit is 64 MiB unless you set
another, and nothing is ever overwritten. Installing osteosarc never downloads data.

## Target lists

`load_panel(name)` returns a shipped list of targets to build recipes from:

| Name | Holds |
| --- | --- |
| vaccine-loci-v1 | The vaccine target alleles |
| sv-regressions-v1 | Three RNA fusions and five candidate SVs used in regression tests |
| sv-candidates-v1 | All 637 [SV candidates](sv-candidates.md) |

## How the libraries build their test data

Osteosarc fetches, selects, packs and checks reads; each library keeps its own science:
which variants it tests, how it reads alleles, and what results it expects. Osteosarc
never imports Isovar, Topiary or Vaxrank.

| Library | Its test data | Built with |
| --- | --- | --- |
| Isovar | 311 fixtures of exact records, plus vaccine, fusion and SV cases | osteosarc.legacy_fixtures |
| Topiary | Variant, indel, fusion and pVACseq fixtures | osteosarc.regional_corpus |
| Vaxrank | 58 read cohorts | osteosarc.cohort_bundle |
| Varcode | Variants and SV records, no reads | a snapshot of the variant catalogue |

Isovar, Topiary and Vaxrank each accept the same recipe in their builders. From an
osteosarc checkout, this checks, with the network off, that all three select exactly
the same reads for the same reasons:

```sh
python -m scripts.check_fixture_consumers --isovar /path/to/isovar \
  --topiary /path/to/topiary --vaxrank /path/to/vaxrank
```

A change to a library's expected results needs its own review, even when its test data
rebuilds cleanly. The tests in tests/test_fixtures.py and tests/test_bundles.py build
recipes and BAMs offline, and cover each rule above.
