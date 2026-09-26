# OpenVax libraries

Pass variants, reads and reports to the OpenVax libraries: Varcode for variant
effects, Isovar for RNA evidence and protein sequences, Topiary for epitope
predictions and Vaxrank for vaccine ranking. Install the ones you need. The
examples open your most recent snapshot (see [Get started](index.md#get-started)),
and the Varcode and Isovar ones need the Ensembl 95 human annotation:

```sh
pyensembl install --release 95 --species homo_sapiens
```

Start each example by opening your snapshot:

```python
from osteosarc import Dataset

data = Dataset.open(offline=False)
```

## Varcode

Turn variants into Varcode variants on the reference you choose:

```python
from pyensembl import EnsemblRelease

selected = data.variants("vaccine", status="ready")
native = selected.to_varcode(genome=EnsemblRelease(95))
for variant in native:
    print(variant, native.metadata[variant]["entries"][0]["id"])
```

This downloads no reference data. Every source entry is kept in the metadata, even
when Varcode merges two entries into one variant. Unusable alleles and genome
build mismatches raise an error.

By default, Varcode converts `chr1` to `1` and `chrM` (or `M`) to `MT`.
Use `variant.contig` for the annotation name and `variant.original_contig`
for the source name. If several entries become one variant, all their names
remain in `native.metadata[variant]["entries"]`.

For a reference with a custom name, say which build it is:

```python
from pyensembl import Genome

genome = Genome(reference_name="GRCh38-osteosarc-six-transcript-subset",
                annotation_name="fixture", gtf_path_or_url="subset.gtf")
custom = selected.to_varcode(genome=genome, assembly="GRCh38")
```

Annotating the variants later needs that reference's files, such as `subset.gtf`
above, to be installed and indexed.

If your GTF uses the same chromosome names as the source, such as `chrM`, turn
off both renaming options:

```python
custom = selected.to_varcode(
    genome=genome, assembly="GRCh38",
    convert_ucsc_contig_names=False, normalize_contig_names=False,
)
```

The names must still match `genome.contigs()`; PyEnsembl can change their case when
it indexes a GTF. Renaming never converts between genome builds, and names such as
`chrUn_KI270442v1` or `NC_012920.1` are left as they are, so your reference must
use them too.

## Isovar

Fetch the reads and hand the BAM to Isovar:

```python
from isovar import ReadCollector
from pyensembl import EnsemblRelease

source = data.file(
    "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam"
)
one = data.variants(ids=["DYNC1H1-chr14-101980529"], status="ready")
native = one.to_varcode(genome=EnsemblRelease(95))
subset = data.extract_reads(source, one.regions(padding=100))
with subset.open() as bam:
    evidence = ReadCollector().read_evidence_for_variant(native[0], bam)
    print(len(evidence.alt_reads), len(evidence.ref_reads))
```

To rebuild the mutant protein sequences:

```python
from isovar import run_isovar

with subset.open() as bam:
    results = list(run_isovar(native, bam))
```

Set transcript and read-collection options in Isovar as usual.

## Vaxrank

Check the BAM's genome build, then fetch reads with their mates:

```python
from osteosarc import Region

source = data.file(
    "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam"
)
info = data.inspect_alignment(source)
if info.assembly == "GRCh38":
    panel = [Region("chr14", 101980428, 101980630, "GRCh38")]
    corpus = data.extract_reads(source, panel, fetch_pairs=True)
    print(corpus.path, corpus.receipt["scope"])
```

Run the reads through Isovar, then rank with Vaxrank as usual. Fetching mates
needs a recent SAMtools; see [requirements](reads.md#requirements).

Compare your ranking with the peptides that were actually used:

```python
for peptide in data.vaccine_peptides("mRNA"):
    print(peptide["variant_id"], peptide["sequence"])
```

[Corrections](corrections.md) change some alleles, including a MAP2 vaccine target.
Use `corrections=False` to reproduce results from the published values.

## Topiary

Load a pVACseq report:

```python
from topiary import read_pvacseq

reports = data.files.select(prefix="neoantigen_prediction/pvactools/", format="tsv")
aggregated = reports.where(lambda f: f.key.endswith(".aggregated.tsv"))
if aggregated:
    predictions = read_pvacseq(data.download(aggregated[0]))
```

Or load RSEM expression:

```python
from topiary.rna.expression_loader import load_expression

rsem = data.files.where(lambda f: f.key.endswith(".genes.results"))
if rsem:
    expression = load_expression(data.download(rsem[0]))
```

Downloaded files keep their original extensions, so format detection works.
Header-only
reports remain empty reports.

For reproducible test data, see [test data](test-data.md).

## Move your own code to osteosarc

If your project downloads osteosarc.com data with its own code, switch the
downloading first, and change references, allele selection or analysis settings in
a separate step.

| Your code does | Use |
| --- | --- |
| Download the website's metadata | `Dataset.sync()` |
| Reopen it later | `Dataset.open()` |
| Parse the variants page and read counts | `data.variants()` |
| List a sample's BAMs and FASTQs | `data.samples["T1_tumor"].files` |
| Find files by assay, time point or folder | `data.files.select(...)` |
| Download a whole file | `data.download(file, to=DIR)` |
| Check a BAM's genome build | `data.inspect_alignment(file)` |
| Fetch the reads in some regions | `data.extract_reads(file, regions)` |
| Make Varcode variants | `variants.to_varcode(genome=...)` |

Keep your project's own reference releases, transcripts, read filters and scoring.
Before dropping an old read extractor, compare whole records, tags included, and
how often each appears: matching read counts aren't enough. Osteosarc fixes some
of the website's alleles by default; open a snapshot with corrections turned off to
compare against results from the published values (see [corrections](corrections.md)).

Files you already downloaded can go straight into the cache, so nothing is fetched
twice (see [snapshots and cache](snapshots.md#import-a-file-you-already-downloaded)).

### Names changed in 0.9

| Before 0.9 | Now |
| --- | --- |
| `data.assets`, `data.asset(key)`, `Asset` | `data.files`, `data.file(key)`, `File` |
| `data.specimens`, `data.describe_samples()` | `data.samples` |
| `data.assets_for_sample(...)` | `data.samples[ID].files` |
| the old `data.samples` | `data.claims` |
| `osteosarc explore` | `osteosarc repl` |
| `osteosarc assets`, `specimens`, `curation` | `osteosarc files`, `samples`, `corrections` |
