# Changelog

Install or upgrade with `python -m pip install --upgrade osteosarc`, and check
your version with `osteosarc --version`. Pin both the package version and your
snapshot (by name or download date) for reproducible analyses. Full release notes are on
[GitHub](https://github.com/iskandr/osteosarc/releases).

## Unreleased

- `subset_bundle` and `test-data make --from` carve out members offline while
  preserving original source headers, acquisition receipts and selection results (#90).
- New `sv-regressions-v2` panel supplies the three Isovar fusions' retained sides;
  newly requested fusion bundles use it. Published v1 panels and bundles stay unchanged (#96).
- Historical retrieval cases have verified native coordinates, exact BAM contig names,
  mitochondrial reference identity and a guard against incompatible merges (#99).
- Reviewed Tempus evidence overlap is available on files and file selections, with
  explicit preference required to exclude one product. Both original labels and
  acquisition receipts remain intact; the true sample identity is still unresolved (#100).
- Python 3.14 support: the tests run on Python 3.9 to 3.14.
- Comprehensive vaccine/RNA test corpus: 51 vaccine-associated variants from the
  union of membership sources, a complete per-BAM support matrix, four
  DNA-supported SNVs with zero reported RNA alternate reads, a historical RNA
  disagreement case, frameshifts and SV evidence. Its
  openvax-v3 read specification is an unpublished candidate.
- Shared read builds can checkpoint completed source selections and set their
  supplementary-alignment depth and region limits explicitly in the spec.

## 0.14.4 (2026-09-28)

openvax-v2, the next bundle of reads the OpenVax libraries share
([#71](https://github.com/iskandr/osteosarc/issues/71),
[#88](https://github.com/iskandr/osteosarc/issues/88)). See
[shared test data](test-data.md#shared-test-data).

- It has openvax-v1's 993 members, the libraries' 429 record for record, and leaves
  out mates with no position of their own (116 records, none of them the libraries'),
  so it rebuilds without reading any BAM's unplaced reads.
- The reads it keeps bring their split alignments: 14 more records that their SA tags
  name. An alignment in a pileup, with more than 10,000 records at its first base, is
  left out.
- SV targets say which side of each breakend they keep, where that's known, and any
  sequence inserted at the junction. openvax-v1's recipe dropped the sides the SV panel
  and catalogue give; ATP5MG--KMT2A, FOXO3--STRADA-CCDC47 and TPST1--CRCP now give
  theirs, and DLG5's deletion is at the junction esvee resolved, with its 24-base
  insertion.
- Carrying a bundle's library records forward keeps a source's mates with no position
  only if one of those records lacks a position.
- Reading past `max_records`, or a recovery policy's, raises RecordLimitError (an
  IntegrityError). For regional reads the cache remembers it: asking again, even
  offline, gives the same answer without reading.

openvax-v1 stays published, and the libraries move to openvax-v2 at their own pace.

## 0.14.3 (2026-09-27)

Anyone can rebuild a bundle from its recipe
([#89](https://github.com/iskandr/osteosarc/issues/89)).

- A recipe's BAMs can be read through any snapshot that lists each one unchanged (the
  same URL, key, size and modification time as the recipe pins), not only through the
  snapshot it was made from, which only its builder had. `osteosarc test-data make
  --recipe` uses `--snapshot` if you give one, else the recipe's own snapshot if you
  have it, else your newest.

## 0.14.2 (2026-09-27)

A tidier command line, from running every command and reading what it printed, and one
place that decides where the OpenVax libraries' shared cache is.

- **The shared cache is datacache's.** osteosarc finds its cache with
  `datacache.get_cache_root("openvax", "OSTEOSARC_CACHE", "OPENVAX_DATA_CACHE")`, as the
  other OpenVax libraries can, so test data they share is downloaded once. It needs
  datacache 1.12.0, which no longer imports pandas.
- **Regions in a BAM's own genome build.** A Region's assembly can be None, meaning
  the build of the BAM it's read from, and `osteosarc reads FILE REGION` uses that
  when `--assembly` is left out (a sample, whose BAMs may differ, still needs it).
- `osteosarc test-data` on its own shows its actions and examples.
- **Tables.**
  - A cut cell ends with "…" instead of mid-word or mid-number.
  - The line under the headers sits under its columns and stops at the terminal's
    edge.
  - Command-line lists leave out columns empty in every row.
  - When even that is too wide, `osteosarc variants` leaves out its vaccines,
    found_by and corrections columns, and says so.
- **Wrapping.** A sample's notes and corrections wrap without splitting a word, and so
  do the one-week timeline's entries, which follow `--width`. `osteosarc files` shows
  each kind's two commonest formats, counts whole.
- **`osteosarc downloads`** lists each file's size and key, then its whole path. It no
  longer lists snapshots' metadata or pages of the bucket's listing.
- `--min-mapq` and `--exclude-flags` say what they do; "1 record" and "1 slide" are
  singular.

## 0.14.1 (2026-09-27)

Downloads survive a passing server error ([#86](https://github.com/iskandr/osteosarc/issues/86)).

- The checks of a remote object's identity, made before and after a download and when
  reading reads from a remote BAM, are tried again after a server error (500, 502,
  503, 504), too many requests (429) or a dropped connection: up to five tries, about
  1, 2, 4 and 8 seconds apart, or as long as the server asks, up to 30 s. A timeout,
  a TLS error or a refusal (other 4xx) still fails at once. The download itself was
  already tried again by datacache.
- A download pinned by checksum, such as a published bundle, no longer fails when
  those checks do: the checksum vouches for the bytes.

## 0.14.0 (2026-09-27)

Fast test data with no options ([#84](https://github.com/iskandr/osteosarc/issues/84)).

- **Mates with no position are left out by default.** `make_bundle`, `osteosarc
  test-data make` and `osteosarc reads --fetch-pairs` leave out mates that didn't
  align and have no position of their own, as STAR stores them. Finding those means
  reading every unplaced read of each BAM, most of the time an RNA-seq BAM takes:
  `test-data make` for one variant from T2_tumor's three RNA-seq BAMs took 15 s
  instead of 91 s, and for another 12 s instead of 19 minutes.
  `make_bundle(..., unplaced_mates=True)`, or `--unplaced-mates` on either command,
  keeps them (it replaces 0.13.2's `--placed-mates`). In Python,
  `extract_reads(fetch_pairs=True)` still fetches every mate, so existing recipes,
  such as openvax-v1's, rebuild the same. Bundles made with 0.13 read their BAMs
  again the first time they're made with 0.14 (their extractions are cached under
  the old setting), unless given `--unplaced-mates`.
- **Removed osteosarc.legacy_fixtures** ([#72](https://github.com/iskandr/osteosarc/issues/72)),
  which no OpenVax library imports any more. osteosarc.cohort_bundle stays.

## 0.13.3 (2026-09-26)

Faster bundle builds ([#82](https://github.com/iskandr/osteosarc/issues/82)). Rebuilding
openvax-v1 from cached reads takes 2.7 minutes and 1.4 GB instead of 31 minutes and
21 GB, and gives the same bundle, byte for byte.

- **Selection reads only what it can choose from.** For each BAM, osteosarc first reads
  every target's window from the extract's index, then keeps the templates found
  there, whole, and any a library names or pins. Before, it held all of the T2
  tumor RNA-seq extract's 13.5 million records.
- **Each source is counted once.** Making a bundle copied and counted every record of a
  source once for each of its members; it now counts them once, and keeps only the
  records the members pin (it still reads the whole extract to find them).
- Record checksums are about twice as fast.
- `select_fixtures(...).records` holds, for a source whose members all pin exact
  records, only those records.

## 0.13.2 (2026-09-26)

Faster read extraction ([#80](https://github.com/iskandr/osteosarc/issues/80)).

- **Leave out mates with no position.** `extract_reads(..., fetch_pairs=True,
  unplaced_mates=False)`, or `osteosarc reads --placed-mates`, gives what
  `fetch_pairs=True` does less the mates that didn't align and have no position of
  their own, as STAR stores them: finding those means reading every such read in the
  BAM. Six extractions from BG003082's RNA-seq BAM took 40 s instead of 68 s, the same
  reads less 28 unplaced mates. Where no read has such a mate it's a second or two
  slower. Bundles can leave them out too: `make_bundle(..., unplaced_mates=False)` or
  `osteosarc test-data make --placed-mates`. By default bundles keep every mate, as
  before, and openvax-v1 rebuilds the same.
- `osteosarc reads SAMPLE` reads a sample's BAMs four at a time, and the first error
  stops those not yet begun.
- An extraction's time limit covers all its samtools runs, however many.

## 0.13.1 (2026-09-26)

Faster commands ([#78](https://github.com/iskandr/osteosarc/issues/78)): three to eight
times faster once a snapshot has been used. On the current snapshot, `osteosarc reads`
on a cached result goes from 7.9 s to 2.8 s, `files` from 6.1 s to 1.6 s, `samples
T1_tumor` from 7.5 s to 1.8 s, `variants` from 1.8 s to 0.2 s and `--version` from 0.5 s
to 0.2 s.

- **A snapshot's catalogues are built once.** The file catalogue (some 400,000 bucket
  objects), its header and the variant catalogue are saved in the cache the first
  time, for that snapshot, osteosarc version and set of corrections. Loading one warns
  of the stale corrections to its sources, as building it does. Each user keeps their
  own copies, in a folder only they can write in, and a copy can hold only
  osteosarc's catalogue classes. Each kind keeps its three most recently used copies,
  and copies unused for 30 days are removed when another is saved.
- Looking files up by key, URL or sample, and looping over any collection, no longer
  goes through every file in Python.
- Corrections find their records through an index instead of scanning every record,
  and variants and the timeline list exactly the corrections that touched them.
- osteosarc imports pandas only to download a file. Checksums of local BAMs are
  remembered in records each user keeps for themselves; those from 0.13.0 (in the
  cache's digests folder) are no longer used and can be deleted. A read-only cache
  still works.

## 0.13.0 (2026-09-26)

A simpler site, and less code ([#75](https://github.com/iskandr/osteosarc/issues/75)).

- **The home page is the overview.** Concepts is folded into a short "How the data
  fits together" there, linking to each topic's page.
- **Removed, now that the OpenVax libraries take their reads from openvax-v1:**
  - osteosarc.regional_corpus;
  - the generators in legacy_fixtures (generate, pack, verify, export_fixture and
    their helpers) and cohort_bundle (generate_cohort_bundle, retrieval_regions);
    the functions Isovar and Vaxrank import stay;
  - osteosarc.bundles.generate_panel; `osteosarc test-data make --recipe`, or
    generate_bundle, does the same;
  - osteosarc.records.sam_digest, and SAM-text record identities in recipes;
  - osteosarc.reads.subset_templates, a random sample that make_bundle's balanced
    sets replace;
  - reading bundles made by osteosarc 0.2 or earlier, which now says to make them
    again.
- The frozen SV test recipe moves from the package into tests/data.

## 0.12.0 (2026-09-26)

Test data in a call or two ([#73](https://github.com/iskandr/osteosarc/issues/73)).

- **Make a bundle without a recipe.** `osteosarc test-data make DIR SAMPLE|FILE
  --variant ID --sv ID` and `data.make_bundle(dir, variants=..., files=...)` write a
  verified bundle with openvax-v1's selection: a balanced set of templates at each
  variant in each BAM, and templates joining each SV's breakends, every record pinned.
  Files can be BAMs or samples. They say which BAMs they skip and why, and stop,
  before streaming any reads, if a requested variant or SV can't be read from the BAMs
  given. `test-data make DIR --recipe FILE` replaces `test-data generate`, and fetches
  the recipe's sources from the snapshot it names, so a bundle's own recipe.json
  rebuilds it.
- **Use a member in one call.** `osteosarc.bundle_file(bundle, member)` returns it as a
  local indexed BAM (or SAM), exported once into the cache with the other members from
  its source BAM, read-only, and reused offline: the helper each OpenVax library had
  written for itself. An unknown member's error suggests close names.
- **One rule for bundles.** list_bundle, export_bundle, verify_bundle, check_fixtures
  and bundle_file all take a bundle folder or a published bundle's bare name, such as
  openvax-v1; a Path, or text with a slash, is always a folder. Offline, a bundle that
  isn't cached yet says how to get it.
- The package docstring, the `data` overview, `osteosarc --help` and the docs lead
  with these.

## 0.11.3 (2026-09-26)

Follow-ups from moving Isovar, Topiary, Vaxrank and Varcode onto openvax-v1, which
all four now use.

- **Receipts and bundles hold no local paths**
  ([#65](https://github.com/iskandr/osteosarc/issues/65)). Read receipts record a file
  by its place in the cache (objects/sha256/..., named by checksum) or else by its
  name, beside its checksum, and bundles leave out download times. So the same reads,
  extracted with the same tools from any cache, make the same bundle, and Isovar's and
  Vaxrank's recorded receipts lose their local paths too. Cache keys no longer depend
  on where the cache is, so reads extracted by earlier versions are fetched once more.
  openvax-v1 keeps its published bytes; the next version will be clean.
- `fetch_bundle(name, offline=True)` never downloads, and raises OfflineError if the
  bundle isn't cached ([#66](https://github.com/iskandr/osteosarc/issues/66)).
- `osteosarc test-data check --help` says fixture paths are relative to the fixtures
  file ([#67](https://github.com/iskandr/osteosarc/issues/67)).
- `check_fixtures` is importable from osteosarc and takes a published bundle's name,
  and the Python API page shows it and osteosarc.shared.published
  ([#68](https://github.com/iskandr/osteosarc/issues/68)).
- Export keeps the name of a member already named after its file: topiary/reads.bam
  is written as reads.bam, not reads.bam.bam; two members that would share a file
  are an error ([#69](https://github.com/iskandr/osteosarc/issues/69)).
- The cross-library builder check is gone, with the builders it compared.

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
[shared test data](test-data.md#shared-test-data).

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

- The documentation is rewritten in plain language. Key concepts now
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

- Reorganized documentation with key concepts, a
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
  additional SV research panel. See [OpenVax fixture adoption](test-data.md#how-the-libraries-use-it)
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
