# Test data

Unit tests need small, real sets of reads that anyone can rebuild. Osteosarc makes
them three ways: a quick BAM of the reads around some variants; a bundle you make for
your library, with a balanced set of reads at each variant, every record pinned; and
openvax-v1, the bundle the OpenVax libraries share. A test reads a bundle's member with
one call, `osteosarc.bundle_file(bundle, member)`.

## A quick test BAM

```sh
osteosarc reads rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam --variant DYNC1H1-chr14-101980529 --padding 100 --to test-data
```

This streams just those reads from the public BAM and saves them, with an index, as
test-data/BG003082.Aligned.sortedByCoord.out.md.DYNC1H1-chr14-101980529.bam: the
BAM's name, then the variant. Give a sample ID instead of a file to
do the same for each of that sample's BAMs:

<!-- docs-check: skip (streams reads from several BAMs) -->
```sh
osteosarc reads T0_tumor --assay rna-seq --variant DYNC1H1-chr14-101980529 --padding 100 --to test-data
```

In Python, `data.extract_reads(file, variants=..., to="test-data")` does the same. See
[reads](reads.md) for regions, filters and mates.

## Make a bundle

Name the BAMs and the variants (or SVs):

<!-- docs-check: skip (makes the same bundle as the Python example below) -->
```sh
osteosarc test-data make dync1h1 rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam --variant DYNC1H1-chr14-101980529
```

This streams the reads around the variant, keeps a balanced set of them (see
[how reads are chosen](#how-reads-are-chosen)), and writes a verified bundle to the new
folder dync1h1, about 100 KB. Each member is one target in one BAM, named like the
files osteosarc reads --to writes:
BG003082.Aligned.sortedByCoord.out.md.DYNC1H1-chr14-101980529. A sample ID stands for
its indexed BAMs, narrowed with `--assay` or `--platform`; repeat `--variant`, and use
`--sv` for an SV candidate or fusion, as in
`osteosarc test-data make gabbr1 T2_tumor --assay wgs --sv GABBR1-SLC29A1`. Commit the
folder beside your tests, or make it in CI.

In Python:

```python
from osteosarc import Dataset

data = Dataset.open(offline=False)
data.make_bundle(
    "dync1h1",
    variants=["DYNC1H1-chr14-101980529"],
    files=["rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam"],
)
```

Variants can also be data.variants(...), and files can be samples, as in
data.samples["T0_tumor"]; caps change how many reads of each kind are kept.

## Use a bundle in tests

```python
import osteosarc

bam = osteosarc.bundle_file("dync1h1", "BG003082.Aligned.sortedByCoord.out.md.DYNC1H1-chr14-101980529")
```

That's an indexed BAM of exactly the member's records, exported once into the cache,
read-only, and reused, offline, by every later test run; `format="sam"` gives SAM
text. The bundle can be a folder or the name of a published bundle, such as
openvax-v1. An unknown member's error suggests close names, and
`osteosarc test-data list dync1h1` lists them all. To write members into a folder
instead, use `osteosarc test-data export`.

## Shared test data: openvax-v1

Isovar, Topiary, Varcode and Vaxrank can all test against one bundle of reads,
openvax-v1, chosen once, here. It holds:

- reads at every variant on the site, plus other alleles the libraries test, such as
  MAP2's allele before its correction and mitochondrial alleles on GRCh37;
- reads that join the breakends of all 7 RNA fusions and 8 of 14 DNA SVs (the other
  six have unresolved breakends or no joining reads in these BAMs, and stay listed,
  empty);
- every read in the libraries' current test files, record for record, so a library
  that moves onto openvax-v1 keeps every assertion it has.

<!-- docs-check: skip (downloads the bundle, 28 MB) -->
```python
bam = osteosarc.bundle_file("openvax-v1", "IPISRC044_tumor_T2_ucla.redux.DYNC1H1-chr14-101980529")
```

The first use downloads the bundle (28 MB), checks it against checksums that ship
with osteosarc, and keeps it in the cache; later uses work offline, and
`offline=True` makes sure they do. Members are named the BAM's name, then the target;
a library's own test files keep their names, such as isovar/chimeric/osteosarc-ont.sam.
From the command line, `osteosarc test-data list openvax-v1` lists the members and
`osteosarc test-data export openvax-v1 DIR --member NAME` writes them into a folder.

To check that your library's copies match, list them in a JSON file that maps
member names to your files, with paths relative to that JSON file: a BAM, SAM or SAM.gz path, or, for SAM lines kept inside
JSON, `{"json": path, "pointer": "/path/to/lines"}`. Every SAM line under the pointer
counts, whatever fields hold it, such as Isovar's sam and partner_sam; leave the
pointer out for the whole file. Then:

<!-- docs-check: skip (needs your own files) -->
```sh
osteosarc test-data check openvax-v1 fixtures.json
```

Records are compared as SAM text, repeats included, so formats and headers don't
matter. The command fails if any file differs.

## How reads are chosen

Bundles you make and openvax-v1 choose reads the same way. At each small variant,
osteosarc sorts every template (a
read with its mate) by what it shows across the allele: alt, ref, something else, or
nothing, when it doesn't span it. Indels are judged across any repeat they could
slide within. It keeps up to 20 alt, 10 ref, 5 other and 2 uncallable templates, in
an order fixed by a hash of each read's group and name, plus the two lowest-quality
alt templates. At a fusion or SV it keeps up to 50 templates with aligned bases within
1 kb of every breakend that join them: split and chimeric reads, pairs the aligner
didn't call proper, and reads spliced or deleted exactly from one breakend to another.
A read that simply runs across the breakends, a proper pair on either side, or a read
spliced between exons that merely lie near them doesn't count. A kept template keeps
all its records. Each record is pinned by checksum, so rebuilding gives the same
bytes, and any change upstream fails loudly. openvax-v1's reads come from the T2 tumor
RNA-seq and WGS BAMs for every target, from any BAM a library already uses there, and
from each RNA BAM the [SV candidates](sv-candidates.md) saw a junction in; a bundle you
make reads only the BAMs you name.

**Rebuilding openvax-v1.** As the libraries move their test reads here and delete their own
copies, each new version carries their members forward from the one before, pinned by
checksum. To make openvax-v2, copy osteosarc/data/bundles/openvax-v1.spec.json to
openvax-v2.spec.json, set its id to openvax-v2 (and its snapshot or targets, if they
change), then:

<!-- docs-check: skip (streams reads from about 40 BAMs; takes about an hour) -->
```sh
python scripts/shared_test_data/build.py openvax-v2 build --carry openvax-v1
```

The build finds those records in the public BAMs again, chooses the rest, builds the
bundle and packs it for a GitHub release. A library that needs different reads gives a
fresh list with `--required`, which replaces everything it had.

openvax-v1 itself came from the libraries' committed test files: the required_isovar.py,
required_topiary.py, required_vaxrank.py and required_varcode.py scripts in
scripts/shared_test_data list their SAM lines (for Varcode, the read names behind its
junctions) at the revision you give, and openvax-v1 as the revision means the commits
it was built from. Its published release record is never replaced by a rebuild.

## Recipes

For anything make doesn't cover (your own alleles or regions, named reads, other
selection rules), write a recipe: the JSON a bundle is built from. Every bundle keeps
its own as recipe.json, so one you made is a good start. Then:

<!-- docs-check: skip (needs your own recipe) -->
```sh
osteosarc test-data make bundle --recipe recipe.json
osteosarc --offline test-data make bundle --recipe recipe.json --source rna=archive.bam
```

The second uses a BAM you already have for the recipe's source rna. In Python,
`generate_bundle(recipe, "bundle", dataset=data)` does the same.

A recipe has targets (what each piece of test data is about), sources (the BAMs its
reads come from) and members (one target in one source, with a rule for picking
reads). Osteosarc only follows the rule: deciding which reads support an allele is up
to the library that uses them, and a member's reads are never an estimate of VAF.

```text
targets:   small_variant   one-based position, ref and alt, and where it came from
           sv              two or more breakends (zero-based, with orientation + / - / null)
           fixture         one of a library's test files, with a description
           unresolved      a reason; gets no reads
sources:   identity        which file: {"key": ...} or {"url": ...}
           assembly, sample, library, product   (null when unknown)
members:   target, source, regions (zero-based, half-open), and a policy:
           regional        reads overlapping the regions, optionally capped
           exact           specific records, by checksum and count, each with reasons
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
intact but not who made it. Exporting writes one sorted, indexed BAM per member (or SAM)
into any folder, named after the member; empty members become valid empty BAMs, and a
file already there with different contents is never replaced. A bundle's size limit is
64 MiB unless you set another. Installing osteosarc never downloads data.

## Target lists

`load_panel(name)` returns a shipped list of targets to build recipes from:

| Name | Holds |
| --- | --- |
| vaccine-loci-v1 | The vaccine target alleles |
| sv-regressions-v1 | Three RNA fusions and five candidate SVs used in regression tests |
| sv-candidates-v1 | All 637 [SV candidates](sv-candidates.md) |

## How the libraries use it

Osteosarc fetches, selects, packs and checks reads; each library keeps its own science:
which variants it tests, how it reads alleles, and what results it expects. Osteosarc
never imports Isovar, Topiary or Vaxrank. All four rely on openvax-v1 for their Sid
reads, and keep only their non-read test data (variants, expected results, recipes):

| Library | From openvax-v1 |
| --- | --- |
| Isovar (1.39.0) | Every Sid read file, exported under its original name |
| Topiary (5.75.0) | 14 read files; its rearrangement fixtures are checked against their members |
| Vaxrank (3.26.0) | Read cohorts, rebuilt byte for byte from members with its packaged recipe |
| Varcode (10.5.8) | No read files; its 3 junction fixtures are checked against their members |

A change to a library's expected results needs its own review, even when its test data
rebuilds cleanly. The tests in tests/test_fixtures.py and tests/test_bundles.py build
recipes and BAMs offline, and cover each rule above.
