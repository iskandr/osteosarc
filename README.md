# osteosarc

A shared Python API for the public [osteosarc.com](https://osteosarc.com/data/)
dataset: discover files, cache verified downloads, parse source tables, select
variants and vaccines, and extract indexed sequencing reads.

The initial implementation consolidates the acquisition patterns found in
Varcode, Isovar, Topiary, and Vaxrank. Those repositories have not yet been
migrated. See the [comparison and migration guide](docs/migration.md).

[Documentation and short examples](https://iskandr.github.io/osteosarc/) ·
[Consumer migration recipes](docs/consumers.md)

## Install

Python 3.10+, Linux or macOS, and `curl` are required. Indexed read extraction
also requires `samtools` on PATH and the `reads` extra.

```sh
python -m pip install -e '.[reads]'
```

## One dataset, one cache

```python
from osteosarc import Cache, Dataset

cache = Cache()  # the shared OpenVax cache (see below)
data = Dataset.sync("2026-09-18", cache=cache)  # explicit metadata download

# Later: verifies the saved snapshot and opens it without network access.
data = Dataset.open("2026-09-18", cache=cache)
```

`sync` acquires about 57 MB of metadata, including the full dated bucket index
and the timeline sources.
It does not download BAMs, FASTQs, or other large data files. Importing the
package also does no downloading.

Files are stored in the **shared OpenVax cache**, the same layout vaxrank's
downloader uses, so identical bytes are stored once across OpenVax tools:
`<root>/objects/sha256/<sha256><original suffixes>` (for example
`…9c2e.genes.results`). osteosarc keeps its receipts, snapshots, bindings and
extracted reads under `<root>/osteosarc/`. The root is `OPENVAX_DATA_CACHE`,
otherwise the platform cache directory for `openvax` (`~/Library/Caches/openvax`
on macOS, `$XDG_CACHE_HOME/openvax` or `~/.cache/openvax` on Linux).
`OSTEOSARC_CACHE` or `Cache(root)` selects an isolated cache instead.

Snapshots pin source URLs and SHA256 receipts. Downloads use atomic publication,
per-URL locks, bounded retries, and checksum verification on reuse. Interrupted
downloads restart; byte-range resume of whole objects is not implemented.
`Dataset.download` binds its first acquired receipt to that snapshot, so
refreshing a URL elsewhere cannot change previously acquired data. Use a new
name with `Dataset.sync(..., refresh=True)` to acquire fresh metadata. Existing
snapshot names cannot be overwritten. `Cache.fetch(url, refresh=True)` can
refresh standalone downloads while keeping old content available.

## Find samples and sequencing products

```python
data = Dataset.open("2026-09-18", offline=False)  # permit explicit acquisition

rna = data.assets.select(kind="alignment", assay="rna-seq", timepoint="T2")
ont = data.assets.select(kind="alignment", platform="ont")
vcfs = data.assets.select(kind="variants", format="vcf")
pvac = data.assets.select(prefix="neoantigen_prediction/pvactools/", format="tsv")
hla = data.assets.select(contains="hla", format="tsv")
raw_reads = data.assets.select(kind="reads", timepoint="T1")

sample_claims = data.samples          # table with source and asset IDs
timepoint_claims = data.timepoints    # published date precision is preserved

for asset in rna:
    print(asset.id, asset.key, asset.size, asset.conflicts)
```

An **asset is a processing product**, not an independent biological replicate.
An ID hashes the complete URL, so identical filenames in different directories
remain distinct. Sample claims come from the BAM viewer, consolidated public
metadata, VAF export, and data-page tables. Conflicting or absent values do not
match a metadata filter by default. Inspect `asset.claims`, `asset.conflicts`,
or `asset.values("timepoint")`; use `include_conflicts=True` to match any claim.
The catalog-wide viewer assembly label is never used to authorize read queries.

Assays are `rna-seq`, `scrna-seq`, `cite-seq`, `wgs`, and `wes`; known platforms
are `ont`, `pacbio`, and `illumina`. Platform is left unknown when it is not
explicit in the source. For paths with no published mapping, timepoint and a
few recognizable library IDs are retained as **inferences**, accessible with
`include_inferred=True`. Pooled blood libraries remain library-level records.

Every bucket object remains discoverable, including unclassified files.
The website listing has its own generation date and is not guaranteed current.
To discover newer files explicitly:

```python
from osteosarc import list_bucket

listing = list_bucket(cache, "ONT/", refresh=True)
# listing["files"]: complete key/size/mtime rows; listing["receipts"]: every page
```

Live listings are returned separately; they do not silently mutate a snapshot.
Missing or repeated continuation tokens and incomplete pagination raise errors.

## Variant sets, annotations, and vaccines

```python
site = data.variants()                  # every site entry, including unresolved
all_exported = data.variants("all")     # also includes count-export entries
vaccine_targets = data.variants("vaccine")  # site-reported vaccine count > 0
mrna = data.variants(vaccine="mRNA")     # vaccine-overlap JSON membership
ready = site.select(status="ready")
dynein = site.select(gene="DYNC1H1")     # may contain multiple alleles
detected = data.variants(pipeline="oncoanalyser")

print(data.vaccine_names, data.pipeline_names)
annotations = data.annotations         # original variant annotation records
peptides = data.vaccine_peptides("mRNA")
elispot = data.vaccines                 # original tested/status/response fields

# A different source's membership assertion, explicitly requested:
source_membership = data.variants(vaccine="JLF V3", vaccine_source="source_variants")
```

The variant index, vaccine-overlap JSON, and source variant JSON can disagree.
`Variant.vaccines` contains overlap-JSON membership joined only at an unambiguous
gene/locus; `annotations["source_vaccines"]` preserves source-JSON flags.
Disagreements and ambiguous joins remain in annotations. Counts from the site
table, peptide inclusion, pipeline detection, and ELISPOT measurements are kept
separately. A blank or untested measurement is not converted to a negative one.

Variants retain the published assembly and original VCF-style anchored allele.
`ready` means a unique literal allele with internally consistent coordinates;
it does **not** mean independently validated, somatic, normalized, or clinically
confirmed. Incomplete entries have a status such as `missing_literal_allele`,
`non_literal_allele`, or `ambiguous_literal_allele`. No allele is fabricated from
a gene symbol or protein label, and no coordinates are silently lifted over.

## Timelines, specimens, and the terminal explorer

```python
print(data.timeline.render(since="2024-05", until="2024-09", width=100))
print(data.timeline.around("2025-01-28", days=5).listing())
t2 = next(r for r in data.specimens if r["sample_id"] == "T2_tumor")
print(t2["date"], t2["site"], len(t2["assets"]), t2["corrections"])
mrd = data.measurements.select(source="mrd")
```

The timeline brings together every dated public source: treatments and doses,
procedures, imaging (events and DICOM studies), pathology, omics, time points,
specimens, MRD, flow-cytometry draws, and lab and cytometry dates. Run
`osteosarc explore baseline` to browse it interactively, or
`osteosarc timeline baseline --since 2024-05` for a chart. See
[timelines](docs/timeline.md) and the [guided tour](docs/tour.md).

## Corrections are central, optional, and drift-aware

All hand-written interpretation of the sources lives in `osteosarc/curation.py`.
Verified corrections are applied by default, for example the GRCh37 Tempus
counts, five relocated Tempus alleles, MAP2's observed allele, stale viewer
labels, and specimen dates and sites. Each touched object names its
corrections, and `Dataset.open(name, corrections=False)` gives the published
sources unchanged. Every load re-checks each correction against the sources, so
an upstream change makes it `stale` and unapplied rather than silently wrong.

```python
for row in data.corrections:
    print(row["status"], row["id"])
print(list(data.unrecognized))       # labels outside the vocabulary
```

`osteosarc curation baseline --strict` exits nonzero on drift. See
[corrections and drift](docs/curation.md).

## Download and parse

```python
counts = data.table("vafs")
counts_for_gene = counts.select(gene="SMC5")  # raw values, including "0" and "NA"
dna_fusions = data.table("dna_fusions")
rna_fusions = data.table("rna_fusions")
snvs = data.table("snv_top")

report = pvac[0]
path = data.download(report)     # explicit full-object download
table = data.parse(report)       # same cache; original column names/strings
frame = table.to_dataframe()     # optional pandas, no automatic numeric coercion

with data.open_variants(vcfs[0]) as calls:  # optional pysam
    for call in calls:
        print(call.contig, call.pos, call.ref, call.alts)
```

Built-in parsers cover TSV, CSV, JSON, and FASTA. `open_variants` preserves
VCF/BCF headers, INFO/FORMAT fields, sample genotypes, multiallelic records, and
symbolic alleles through pysam. Other objects (GTFs, matrices, RDS, imaging,
etc.) can be downloaded for their format-specific tools. Parsing a large text
table or FASTA currently materializes it in memory; use its downloaded path
with a streaming reader for very large files.

## Indexed read extraction

```python
from osteosarc import ReadFilter, Region

# Explicit source product; prefix/substring selectors can return several BAMs.
source = data.asset(
    "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam"
)
variant = data.variants()["DYNC1H1-chr14-101980529"]
subset = data.extract_reads(source, [variant.region(padding=100)])

# Region coordinates in Python are always zero-based and half-open.
region = Region("chr14", 101980528, 101980530, "GRCh38")
filtered = data.extract_reads(
    source, [region], filters=ReadFilter(min_mapq=20, exclude_flags=0x100 | 0x400)
)
with subset.open() as bam:
    for read in bam.fetch("chr14", 101980528, 101980530):
        print(read.query_name, read.cigarstring)
```

The extractor inspects reference lengths before querying, resolves contig
aliases, and refuses assembly conflicts, out-of-bounds intervals, ambiguous
contigs, and missing indexes. GRCh37 mitochondrial queries additionally require
the exact `reference_length` because hg19 and hs37d5 differ there. CRAM requires
a local indexed reference FASTA. Sparse/custom reference headers that cannot
establish the assembly are rejected.

SAMtools reads the indexed union of all intervals. Overlapping query intervals
do not cause repeated emission; genuine repeated source records are preserved.
Defaults retain secondary, supplementary, duplicate-marked, and low-MAPQ reads,
including original qualities, barcodes, UMIs, and other tags. Optional filters
are recorded in the request. `ReadFilter(barcodes=("...",))` selects the `CB`
tag (customizable via `barcode_tag`). `fetch_pairs=True` additionally retrieves
paired mates, but is not recovery of every supplementary record for a template.

The result contains an indexed BAM, source header, BED intervals, checksum
receipt, source/index identities, tool versions, and record count. It is cached
by the complete request. Remote identity is checked before and after new
extractions; cached results are verified offline without contacting the source.
No full-alignment scan is used as a fallback.

Local BAMs use the same implementation:

```python
from osteosarc import extract_reads, subset_templates

local = str(subset.path)  # any local indexed BAM; here, the subset from above
regional = extract_reads(local, [region], cache=cache)
fixture = subset_templates(regional, count=48, seed="regression-v1", cache=cache)
```

Template sampling is a separate operation, keyed by `(read group, query name)`.
It retains all available regional records for selected templates and marks the
receipt as a sampled fixture unsuitable for estimating VAF.

## Downstream libraries and CLI

```python
# Optional varcode / pyensembl adapter; the caller provides its reference.
from pyensembl import EnsemblRelease
variants = ready.to_varcode(genome=EnsemblRelease(95))
# Pass variants and subset.path to Isovar, then its products to Topiary/Vaxrank.
```

Osteosarc does not choose an Ensembl release, install a genome, calculate effects,
assemble RNA, score epitopes, infer somatic truth, or design a vaccine.

```sh
osteosarc --cache .cache/osteosarc sync baseline
osteosarc --cache .cache/osteosarc assets baseline --assay rna-seq --timepoint T2
osteosarc --cache .cache/osteosarc variants baseline --set vaccine --status ready
osteosarc --cache .cache/osteosarc variants baseline --vaccine 'JLF V3'
osteosarc --cache .cache/osteosarc table baseline dna_fusions
osteosarc --cache .cache/osteosarc reads baseline 'rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam' chr14:101980529-101980530 --assembly GRCh38
osteosarc --cache .cache/osteosarc discover 'neoantigen_prediction/pvactools/'
osteosarc --cache .cache/osteosarc timeline baseline --since 2024-05 --until 2024-09
osteosarc --cache .cache/osteosarc specimens baseline T1_tumor
osteosarc --cache .cache/osteosarc curation baseline --strict
osteosarc --cache .cache/osteosarc --no-corrections variants baseline --gene MAP2
```

CLI regions are **one-based inclusive**, matching SAMtools. Add the global
`--offline` flag to prohibit acquisition. `sync --source-revision <commit>` pins
every GitLab source-repository resource (variant JSON, BAM metadata, and the
timeline and specimen sources) to a full commit; site-served files are not
versioned by the site.

## Development and provenance

```sh
python -m pip install -e '.[test]'
ruff check osteosarc tests scripts
python -m pytest -q
python -m build --no-isolation
```

Tests are offline and use small public metadata excerpts plus synthetic indexed
BAMs. See [validation](docs/validation.md) for the bounded live check and
[API contracts](docs/design.md) for the design rationale.

Code is Apache-2.0. The public dataset is separately listed as CC0-1.0 in the
[AWS Open Data Registry](https://registry.opendata.aws/sid-osteosarc/). Cite the
dataset and the access date when using its data.
