# Changelog

Install or upgrade with `python -m pip install --upgrade osteosarc`, and check
your version with `osteosarc --version`. Pin both the package version and your
snapshot (by name or download date) for reproducible analyses. Full release notes are on
[GitHub](https://github.com/iskandr/osteosarc/releases).

## 0.7.0 (2026-09-24)

Every correction was re-checked against the reference genome and the site's data.
None only rewrites an equivalent form of a value, and these fixes follow from the
review:

- `tempus-timepoint` now says the Tempus tumor's variants match T0, so the T1 labels
  on the Tempus files are probably what's wrong.
- `provider-IPISRC044-T1-rna` also fixes the site's read-count rows, which switched
  to UCLA on 2026-09-21.
- `transcript-DCHS2` fixes only the typo, to the versionless NM_001142552 the site
  now uses; it no longer adds a version.
- New: `transcript-COL4A2` and `transcript-GTF3C5` fix two more accession typos.
- MAP2 keeps its published read counts, which already count reads carrying the real
  change.
- USH2A's retired entry is recorded as a confirmed duplicate, and the kept entry's
  sequence context is fixed. The moved Tempus alleles lose sequence context taken
  from the wrong position.
- `fam157a-withdrawn-protein` shows as fixed upstream where the site has the same
  note.
- The corrections page is rewritten in plain language.

## 0.6.1 (2026-09-24)

- The documentation is rewritten in plain language. [Key concepts](concepts.md) now
  says where every piece of data comes from: osteosarc.com, its S3 bucket, and the
  website's source repository on GitLab ([#41](https://github.com/iskandr/osteosarc/pull/41)).

## 0.6.0 (2026-09-24)

- `RecoveryPolicy(on_timeout="incomplete")` keeps the verified seed reads and
  partners already matched when a later partner query times out, instead of
  failing. The receipt says `status="incomplete"` and names the failed query; a
  later call retries it and resumes from verified extractions. Fixture recipes can
  request this through a source's `acquisition`, and bundles record the
  `incomplete` acquisition ([#28](https://github.com/iskandr/osteosarc/issues/28)).

## 0.5.0 (2026-09-24)

- The CLI is easier to explore with ([#38](https://github.com/iskandr/osteosarc/issues/38)):
  - `samples` and `specimens` show sequencing with the filter names, such as
    `rna-seq; wes; wgs; scrna-seq (ont, pacbio)`.
  - `samples --assay scrna-seq --platform ont` lists the samples with that sequencing;
    `describe_samples()` takes the same filters.
  - `assets` prints a table with complete file keys and a hint when rows are cut
    off; `--json` prints the full records.
  - Explorer tables never truncate file keys, and sizes read as KB, MB or GB.
- An assay, platform or tissue filter that no file uses raises an error listing the
  valid values, instead of an empty selection. A registry label such as
  `scRNA_ONT` names the filters to use instead.

## 0.4.0 (2026-09-24)

- Removed the CLI's positional snapshot argument (`osteosarc samples baseline`);
  choose a snapshot with `--snapshot`. Removed the empty `[reads]` install extra.
- The explorer's `assets` command accepts `sample=`, as in `assets sample=T1_tumor`.

## 0.3.0 (2026-09-24)

- Snapshots are indexed by download date, and the most recent is the default
  ([#35](https://github.com/iskandr/osteosarc/issues/35)). `Dataset.sync()` and
  `osteosarc sync` need no name: they save a snapshot named by the UTC date and
  reopen it for the rest of that day. `Dataset.open()` and every CLI command use the
  newest snapshot. `Dataset.snapshots()` and `osteosarc snapshots` list them.
  `Dataset.open(date="2026-09")` or `--snapshot 2026-09` picks the newest from a
  download date, month or year, and `Dataset.open(name)` or `--snapshot NAME` picks
  one exactly by name or ID prefix. Offline, `sync()` builds a snapshot from sources
  already in the cache.
- Named snapshots keep working. The CLI's leading snapshot argument
  (`osteosarc samples baseline`) is deprecated in favor of `--snapshot`.
- The documentation no longer uses a `baseline` snapshot name.

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
