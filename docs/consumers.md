# Use other libraries

These examples use the `baseline` snapshot from [Get started](index.md).
Install the libraries you want to use. The Isovar examples also need an
indexed human Ensembl 95 reference:

```sh
pyensembl install --release 95 --species homo_sapiens
```

Start each example by opening your snapshot:

```python
from osteosarc import Dataset

data = Dataset.open("baseline", offline=False)
```

## Varcode

Convert selected alleles to native Varcode objects on your chosen reference:

```python
from pyensembl import EnsemblRelease

selected = data.variants("vaccine", status="ready")
native = selected.to_varcode(genome=EnsemblRelease(95))
for variant in native:
    print(variant, native.metadata[variant]["entries"][0]["id"])
```

The adapter downloads no reference data. It preserves snapshot provenance and
all source entries, including entries Varcode normalizes to the same allele.
Unresolved alleles and assembly mismatches raise errors.

For a custom reference name, declare the assembly:

```python
from pyensembl import Genome

genome = Genome(reference_name="GRCh38-osteosarc-six-transcript-subset",
                annotation_name="fixture", gtf_path_or_url="subset.gtf")
custom = selected.to_varcode(genome=genome, assembly="GRCh38")
```

Conversion preserves the genome's cache identity. Annotation later requires
the reference files, such as `subset.gtf` above, to be available and indexed.

## Isovar

Extract reads and pass the resulting alignment handle to Isovar:

```python
from isovar import ReadCollector
from pyensembl import EnsemblRelease

source = data.asset(
    "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.bam"
)
one = data.variants(ids=["DYNC1H1-chr14-101980529"], status="ready")
native = one.to_varcode(genome=EnsemblRelease(95))
subset = data.extract_reads(source, one.regions(padding=100))
with subset.open() as bam:
    evidence = ReadCollector().read_evidence_for_variant(native[0], bam)
    print(len(evidence.alt_reads), len(evidence.ref_reads))
```

To run protein reconstruction on those inputs:

```python
from isovar import run_isovar

with subset.open() as bam:
    results = list(run_isovar(native, bam))
```

Set transcript and read-collection options in Isovar as usual.

## Vaxrank

For a read corpus, inspect the reference and request paired mates when needed:

```python
from osteosarc import Region

source = data.asset(
    "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.bam"
)
info = data.inspect_alignment(source)
if info.assembly == "GRCh38":
    panel = [Region("chr14", 101980428, 101980630, "GRCh38")]
    corpus = data.extract_reads(source, panel, fetch_pairs=True)
    print(corpus.path, corpus.receipt["scope"])
```

Pass the resulting reads through Isovar, then use Vaxrank's existing predictor
and ranking configuration. See [read requirements](reads.md#requirements) for
paired-mate support.

Fetch published peptides to compare with your results:

```python
for peptide in data.vaccine_peptides("mRNA"):
    print(peptide["variant_id"], peptide["sequence"])
```

When updating fixtures, review [source corrections](curation.md), especially
MAP2's changed allele. Use `corrections=False` to reproduce the published inputs.

## Topiary

Load a pVAC report from the cache:

```python
from topiary import read_pvacseq

reports = data.assets.select(prefix="neoantigen_prediction/pvactools/", format="tsv")
aggregated = reports.where(lambda a: a.key.endswith(".aggregated.tsv"))
if aggregated:
    predictions = read_pvacseq(data.download(aggregated[0]))
```

Or load RSEM expression:

```python
from topiary.rna.expression_loader import load_expression

rsem = data.assets.where(lambda a: a.key.endswith(".genes.results"))
if rsem:
    expression = load_expression(data.download(rsem[0]))
```

Cached paths retain their file suffixes for format detection. Header-only
reports remain empty reports.

See [migration](migration.md) to replace existing download helpers and
[testing](validation.md) for integration coverage.
