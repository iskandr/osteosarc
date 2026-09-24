# Changelog

Install or upgrade with `python -m pip install --upgrade osteosarc`, and check
your version with `python -m pip show osteosarc`. Pin both the package version and your
snapshot name for reproducible analyses. Full release notes are on
[GitHub](https://github.com/iskandr/osteosarc/releases).

## 0.2.6 (2026-09-24)

- Reorganized documentation with [key concepts](concepts.md), a
  [command-line guide](cli.md), a complete [API reference](api.md) and this
  changelog. The README is a shorter landing page that lists the main features
  ([#32](https://github.com/iskandr/osteosarc/issues/32)).
- The documentation example checker covers every published page and the strict
  build validates anchors ([#31](https://github.com/iskandr/osteosarc/issues/31)).

## 0.2.5 (2026-09-24)

- The [SV interest catalogue](sv-interest.md): 637 structural-variant nominations with
  original calls, breakend geometry, Ensembl 115 annotation, T2 expression and
  scoped RNA evidence. It loads offline with `load_sv_interest()` or
  `load_panel("sv-interest-v1")` ([#29](https://github.com/iskandr/osteosarc/pull/29)).
- The command line can now run the whole quickstart: `assets --sample`,
  `reads --variant ID --padding N` and `--version`
  ([#33](https://github.com/iskandr/osteosarc/pull/33)).

## 0.2.4 (2026-09-23)

- Partner recovery selects seed query names with `samtools view -N` before applying
  its record cap, so unrelated reads at dense loci can't exhaust it. Adds
  `ReadFilter(query_names=...)` ([#26](https://github.com/iskandr/osteosarc/issues/26)).

## 0.2.3 (2026-09-23)

- Shared historical fixture adapters for Isovar, Topiary and Vaxrank, and the
  additional SV research panel. See [OpenVax fixture adoption](fixture-migration.md)
  ([#15](https://github.com/iskandr/osteosarc/issues/15)).
- Python 3.9 support.

## 0.2.2 (2026-09-23)

- Portable fixture bundles: generate, pack, verify, list and export
  ([#14](https://github.com/iskandr/osteosarc/issues/14)).

## 0.2.1 (2026-09-23)

- Bounded mate and `SA` partner recovery with explicit receipts
  ([#13](https://github.com/iskandr/osteosarc/issues/13)).

## 0.2.0 (2026-09-23)

- Versioned [fixture recipes](fixtures.md) and one selection executor shared by the
  Python API, `Dataset` and CLI ([#12](https://github.com/iskandr/osteosarc/issues/12)).

## 0.1.4 (2026-09-22)

- Reconciled eight corrections after upstream source changes; historical snapshots
  keep their original IDs and corrections ([#10](https://github.com/iskandr/osteosarc/issues/10)).

## 0.1.3 (2026-09-21)

- `to_varcode()` converts `chrM` and `M` to Ensembl's `MT` and keeps the original
  contig name ([#8](https://github.com/iskandr/osteosarc/issues/8)).

## 0.1.2 (2026-09-21)

- Verified GRCh38 alleles for FAM157A and COL3A1, with explicit outcomes for MUC3A,
  OTUD4 and USH2A ([#5](https://github.com/iskandr/osteosarc/issues/5)).
- Clearer sample, timepoint and assay documentation.

## 0.1.1 (2026-09-20)

- Downloads use datacache; pysam is a standard dependency.
- Adds `describe_samples()`, `assets_for_sample()` and `extract_reads(variants=...)`.
- SAMtools capabilities are checked before acquisition
  ([#3](https://github.com/iskandr/osteosarc/issues/3)); malformed count rows no
  longer hide valid entries ([#4](https://github.com/iskandr/osteosarc/issues/4)).

## 0.1.0 (2026-09-19)

- First release: pinned metadata snapshots, file discovery, variants and vaccines,
  indexed read extraction, source corrections, the clinical timeline and the
  interactive explorer.
