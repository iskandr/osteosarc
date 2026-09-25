# OpenVax fixture adoption

Isovar, Topiary and Vaxrank build their test data with Osteosarc's fixture recipes,
bundles and checks, instead of each keeping its own code. Varcode uses only the
snapshot adapter and `Dataset`/`Cache`, without test BAMs. The shared code arrived
in Osteosarc 0.2.0 to 0.2.3 (PRs [#22](https://github.com/iskandr/osteosarc/pull/22)
to [#25](https://github.com/iskandr/osteosarc/pull/25)), and the libraries switched
over on 2026-09-23: [Isovar #351](https://github.com/openvax/isovar/pull/351),
[Topiary #380](https://github.com/openvax/topiary/pull/380),
[Vaxrank #503](https://github.com/openvax/vaxrank/pull/503) and
[Varcode #489](https://github.com/openvax/varcode/pull/489).

## Who owns what

Osteosarc owns fetching, selecting, packing and checking reads. Each library keeps
its own science: which variants to test, how to classify reads, and what results
to expect.

| Library and fixtures | Osteosarc provides | The library keeps |
| --- | --- | --- |
| Isovar: 311 exact-record fixtures | `legacy_fixtures` (keeps the original SAM order and counts) | The recipe and test assertions |
| Isovar: vaccine variants, mitochondrial and quality controls | `select_assigned_names`, `write_selected_names` | Read classification and chosen reads |
| Isovar: SVs, three fusions, empty controls | `select_window_segments` and indexed fetching | Breakends, exon windows, protein reconstruction |
| Topiary: variant, indel and RNA fixtures | `regional_corpus`, `fixture_assets` | Annotations and prediction expectations |
| Vaxrank: 58 read cohorts | `cohort_bundle` | The cohort recipe and ranking tests |
| Varcode: snapshot and VCF fixtures | The optional snapshot adapter | Alleles and protein tests |

Osteosarc never imports Isovar, Topiary or Vaxrank.

## One recipe, three builders

Each library's builder accepts the same recipe:

```text
--panel-recipe recipe.json --panel-source rna=archive.bam --output NEW_DIRECTORY --offline
```

- `python -m isovar.sid_data generate ...`
- `python -m scripts.generate_sid_fixtures ...` (Topiary)
- `python examples/osteosarc_test_data/build.py --cache CACHE ...` (Vaxrank)

To check that all three pick exactly the same reads for the same reasons, with
network access disabled, run this from an Osteosarc checkout:

```sh
python -m scripts.check_fixture_consumers --isovar /path/to/isovar \
  --topiary /path/to/topiary --vaxrank /path/to/vaxrank
```

Using a bundle doesn't need any of the other checkouts.

## The SV regression panel

`osteosarc/data/additional_sv_recipe.json` pins five candidate SVs: SV0461, SV0402,
SV0055 (near ATP8B5P), SV0175 (a GABBR1 duplication) and SV0499 (near IMMT). Their
251 reads come from four BAMs in the dataset: T1 and T2 Oxford Nanopore, T1 PacBio
and T1 short-read RNA-seq. The source distribution includes those reads and the
command that rebuilds them (`tests/data/additional_svs/README.md`).

The `sv-regressions-v1` panel adds three RNA fusions: PARD3B/CDKN2B, GABBR1/SLC29A1
and OTUD7A/FMN1.

## Changing versions

Each library pins its Osteosarc version and bumps its own version when the pin
changes. Before changing a pin, run the library's tests and check its built package
for the right fixtures and size. A change to expected scientific results needs its
own review, even when the fixtures rebuild cleanly.
