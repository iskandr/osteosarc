# Changelog

Install or upgrade with `python -m pip install --upgrade osteosarc`, and check
your version with `osteosarc --version`. Pin both the package version and your
snapshot (by name or download date) for reproducible analyses. Full release notes are on
[GitHub](https://github.com/iskandr/osteosarc/releases).

## 0.11.2 (2026-09-26)

- **Rebuilding shared test data no longer needs the libraries' own copies**
  ([#63](https://github.com/iskandr/osteosarc/issues/63)). The build script's
  `--carry openvax-v1` keeps each library's members from the previous bundle, pinned by
  checksum, so the next version builds after a library deletes its test files; a
  `--required` list replaces all a library had. Fixture targets now record the regions
  they were planned from, so later carries plan them the same way.
- `osteosarc test-data` treats a published bundle's name as that bundle, even when a
  folder of the same name exists; write ./NAME for the folder.
- The required_*.py scripts take the library revision to read; openvax-v1 means the
  commits openvax-v1 was built from, and reproduces its inputs exactly. A rebuild never
  replaces a published release record.
- The cross-library builder check compares whichever libraries you give it.

## 0.11.1 (2026-09-26)

- `osteosarc test-data check` counts every SAM line under a JSON pointer, whatever
  fields hold it, so records kept as sam and partner_sam (Isovar's fusion corpus) no
  longer look half missing, and names or versions beside them aren't counted
  ([#61](https://github.com/iskandr/osteosarc/issues/61)). Pointers follow the JSON
  Pointer standard, a pointer to nothing says where it stopped, and a fixture that
  can't be read is reported without stopping the check of the others.

## 0.11.0 (2026-09-26)

The OpenVax libraries can take their test reads from one place: **openvax-v1**, chosen
once here ([#56](https://github.com/iskandr/osteosarc/issues/56)). See
[shared test data](test-data.md#shared-test-data-openvax-v1).

- **What it holds:** reads at every variant on the site plus other alleles the
  libraries test, reads joining the breakends of 7 RNA fusions and 8 DNA SVs, and
  every record in Isovar's, Topiary's, Vaxrank's and Varcode's current test files,
  exactly as they are.
- **Using it:** `osteosarc test-data list openvax-v1` downloads it the first time, and
  export and verify take its name too; check compares a library's own copies with
  it. In Python, `fetch_bundle("openvax-v1")` returns its folder.
- **How reads are chosen:** a new allele classifier sorts each read at a variant into
  alt, ref, other or uncallable, even for indels in repeats, and openvax-v1 keeps up
  to 20 alt, 10 ref, 5 other and 2 uncallable templates per variant and BAM, whole.
- **Export writes just the members.** `osteosarc test-data export BUNDLE DIR --member
  NAME` writes NAME.bam, with its index, into DIR, which may already exist; it no
  longer copies the whole bundle, and never replaces a different file. Manifests no
  longer record exports.
- **Recipes:** an exact member can give each record its own reason, and a target can
  be one of a library's test files.
- **Docs:** the example file names for reads saved with --to now match what it writes.

## 0.10.0 (2026-09-26)

A smaller, more consistent osteosarc, centered on exploring the data and making
test data ([#58](https://github.com/iskandr/osteosarc/issues/58)).

- **Test data in one command.** `osteosarc reads T1_tumor --assay rna-seq --variant ID
  --padding 100 --to tests/data` streams just those reads from each of a sample's
  BAMs, and saves each as a small indexed BAM named for its source and variant.
  `data.extract_reads(..., to="tests/data")` does the same in Python.
- **Command line.** `osteosarc timeline --around DATE` replaces `osteosarc on`, and
  `osteosarc test-data` (generate, list, verify, export) replaces `fixtures`; its
  select, pack and panel actions are gone (use the Python functions). `table` and
  `discover` are gone.
- **No silent overwrites.** Saving a download or reads into a folder never replaces
  a file with different contents; it stops and says so.
- **Python.**
  - Dataset loses `table` (use `parse`), `claims`, `timepoints`, `annotations`,
    `vaccine_names`, `pipeline_names`, `receipts`, `open_variants`, and the
    `generate_bundle` and `select_fixtures` shortcuts (use the functions).
  - The package's top-level names are the core API. Helpers such as
    `list_bucket`, `read_records` and `SNAPSHOT_SOURCES` now come from their modules.
- **Names.** The SV interest catalogue is now **SV candidates**: `load_sv_candidates()`
  and the `sv-candidates-v1` panel. FASTQ folder tables name assays and platforms
  as everywhere else, with a library column (gene expression, TCR, BCR, antibody
  tags, long reads).
- **Docs.**
  - Page addresses match their titles (command-line, python-api, snapshots,
    corrections, testing, map2, sv-candidates, test-data, openvax).
  - The three pages about moving code onto osteosarc are merged into two.
  - Every page has much less code formatting.

## 0.9.0 (2026-09-26)

The command line and the Python API now use the same two nouns, samples and files,
and every command prints something a person can read
([#46](https://github.com/iskandr/osteosarc/issues/46)–[#55](https://github.com/iskandr/osteosarc/issues/55)).

**Samples.** `data.samples` is the list of samples (tumor, organoid and blood), each
a `Sample` with its sequencing, BAMs and FASTQ folders. `data.samples["T1_tumor"]`
and `osteosarc samples T1_tumor` show one sample's files with their sizes and
whether they're downloaded, followed by the commands that fetch them;
`osteosarc samples --files` lists every sample's files. Sequencing now also comes
from the site's FASTQ table, so the four 2026 blood draws, which the registry leaves
blank, show single-cell and CITE-seq instead of "unknown".

**Files.** `data.assets`, `data.asset()`, `Asset` and `Assets` are now `data.files`,
`data.file()`, `File` and `Files`, and `osteosarc assets` is `osteosarc files`. With
no filters it summarizes the bucket by kind and folder; lists put BAMs first and name
each file's sample. `files.select(sample=...)` and `file.samples` link files to
samples, and scans and slides are a new `image` kind.

**Downloads.** `osteosarc downloads` and `data.downloads()` list what's on this
computer, with local paths, and `data.local_path(file)` finds one file's copy
without the network. `download(file, to=DIR)` and `osteosarc download FILE --to DIR`
put a file, and its index, in a folder under its own name.

**Timeline.** The chart has a row per treatment, named and grouped by class, over a
month axis, with rows for time points, samples, procedures, scans and each MRD
assay. Lab draws, DICOM studies and other frequent records are left out unless you
pass `--all` (`everything=True`), and the legend says how many.

**Command line.**

- Commands are grouped (browse, get data, snapshots, more) in `osteosarc` and
  `osteosarc --help`, and `osteosarc` on its own says which snapshot it's using.
- `variants`, `vaccines`, `corrections`, `sync`, `table` and `discover` print text;
  `--json` gives the old output. `variants ID` and `corrections ID` show one record.
- `curation` is now `corrections`. `specimens` and `timepoints` are gone: use
  `samples`.
- `download` and `reads` print plain paths (`reads --json` for the receipt).
- `osteosarc repl` replaces the `explore` mini-shell: it opens Python with the
  snapshot loaded as `data`.

**Python.** `repr(data)` lists what's there and how to get it; collections end with
how to narrow them. The old `data.samples` (what each source says about a file) is
now `data.claims`, with `file_ids`. `data.specimens`, `data.describe_samples()`,
`data.assets_for_sample()` and `data.explore()` are gone.

## 0.8.0 (2026-09-25)

A first look at the data is easier from both the command line and Python
([#44](https://github.com/iskandr/osteosarc/issues/44)):

- `osteosarc` on its own suggests where to start.
- Without a snapshot, commands say which cache they looked in and to run
  `osteosarc sync`. In a terminal, `osteosarc explore` offers to download one.
- In Python, `data`, variants, files, tables and the timeline show readable previews
  instead of `<object at 0x…>`. Text views such as `describe_samples()` display without
  quotes, and single variants, files and events no longer print their full source
  records.
- `data.summary()` shows what's in a snapshot and what to try next; `data.explore()`
  opens the interactive explorer.
- The README and docs home start with "Explore the data".

## 0.7.0 (2026-09-24)

Every correction was re-checked against the reference genome and the site's data.
None only rewrites an equivalent form of a value, and these fixes follow from the
review:

- The Tempus tumor's variants match T0, so the T1 labels on the Tempus files are
  probably what's wrong. `tempus-timepoint` now says so, and a new flag,
  `tempus-file-labels`, marks the files.
- `provider-IPISRC044-T1-rna` also fixes the site's read-count rows, which switched
  to UCLA on 2026-09-21.
- `transcript-DCHS2` fixes only the typo, to the versionless NM_001142552 the site
  now uses; it no longer adds a version.
- New: `transcript-COL4A2` and `transcript-GTF3C5` fix two more accession typos.
- MAP2 keeps its published read counts as an approximation: the site counts any
  large deletion there, which in the T0 tumor exome is always the real change.
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
  [command-line guide](command-line.md), a complete [API reference](python-api.md) and this
  changelog. The README is a shorter landing page that lists the main features
  ([#32](https://github.com/iskandr/osteosarc/issues/32)).
- The documentation example checker covers every published page and the strict
  build validates anchors ([#31](https://github.com/iskandr/osteosarc/issues/31)).

## 0.2.5 (2026-09-24)

- The [SV interest catalogue](sv-candidates.md): 637 structural-variant nominations with
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
  additional SV research panel. See [OpenVax fixture adoption](test-data.md#each-librarys-own-test-data)
  ([#15](https://github.com/iskandr/osteosarc/issues/15)).
- Python 3.9 support.

## 0.2.2 (2026-09-23)

- Portable fixture bundles: generate, pack, verify, list and export
  ([#14](https://github.com/iskandr/osteosarc/issues/14)).

## 0.2.1 (2026-09-23)

- Bounded mate and `SA` partner recovery with explicit receipts
  ([#13](https://github.com/iskandr/osteosarc/issues/13)).

## 0.2.0 (2026-09-23)

- Versioned [fixture recipes](test-data.md) and one selection executor shared by the
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
