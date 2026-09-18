# Recipes for the four consumers

These recipes replace acquisition and format conversion at the existing usage
sites. Keep reference releases, transcript selection, RNA interpretation,
peptide scoring, and expected scientific outputs in each consumer.

```python
from osteosarc import Dataset

data = Dataset.open("baseline", offline=False)
selected = data.variants("vaccine", status="ready")
print(selected.source["corrections"][:3])   # applied correction IDs travel with the selection
```

Verified corrections are on by default (see [corrections](curation.md)):

* The MAP2 vaccine target is the observed complex allele, not the published
  22-bp deletion.
* Five Tempus entries become `ready` at their literal alleles.
* The GRCh37 Tempus BAM's website counts are cleared.

Native Varcode metadata records the applied IDs. For comparisons with
published results, open a second Dataset with `corrections=False`.

## Varcode: selected alleles on your reference

```python
from pyensembl import EnsemblRelease

native = selected.to_varcode(genome=EnsemblRelease(95))
for variant in native:
    print(variant, native.metadata[variant]["entries"][0]["id"])
```

The caller selects the reference. The adapter itself downloads no annotation
data. Effect calculation may need your existing indexed reference installation.

For the custom-named partial genomes used in the fixture suites, keep their
unique names and state assembly explicitly:

```python
from pyensembl import Genome

# Stand-in for your fixture's custom subset genome; to_varcode reads no annotation files.
genome = Genome(reference_name="GRCh38-osteosarc-six-transcript-subset",
                annotation_name="fixture", gtf_path_or_url="subset.gtf")
native = selected.to_varcode(genome=genome, assembly="GRCh38")
print(native[0].ensembl.reference_name)
```

This preserves the genome's cache identity. A conflicting known reference name
is rejected. Native metadata retains the snapshot provenance and all original
entries if Varcode normalizes multiple entries to one allele. Unresolved
selected entries raise; they are never silently dropped by the adapter.

Migration sites: Varcode's osteosarc phasing `prepare.py` and `bounded_acquire.py`
can remove imports into Isovar's private test tree. Curated structural variants
and their expected effects remain Varcode fixtures.

## Isovar: native variants and an alignment handle

```python
from isovar import ReadCollector

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

The full pipeline accepts the same native inputs:

```python
from isovar import run_isovar

# Requires your selected reference's installed/indexed annotation resources.
with subset.open() as bam:
    results = list(run_isovar(native, bam))
```

Keep your current transcript whitelist, read collector, and protein-sequence
creator options. Migration sites: `tests/data/osteosarc/fetch_sources.py` and
`expansion/{inventory,discover,acquire}.py`. Keep liftover, fusion validation,
and the existing alternate-read stress-selection recipe local to Isovar.

## Topiary: original report paths

```python
from topiary import read_pvacseq

reports = data.assets.select(prefix="neoantigen_prediction/pvactools/", format="tsv")
aggregated = reports.where(lambda a: a.key.endswith(".aggregated.tsv"))
if aggregated:
    predictions = read_pvacseq(data.download(aggregated[0]))
```

```python
from topiary.rna.expression_loader import load_expression

rsem = data.assets.where(lambda a: a.key.endswith(".genes.results"))
if rsem:
    expression = load_expression(data.download(rsem[0]))
    raw = data.table(rsem[0])  # optional original strings/columns
```

The cache preserves filenames for Topiary's format detection. Header-only
reports remain empty reports. Topiary retains pVAC interpretation, numeric
coercion, expression overlays, and fixture row selection. Migration sites:
`scripts/osteosarc_variant_audit.py`, `scripts/osteosarc_rna_overlay.py`, and
`tests/data/pvacseq/osteosarc/regenerate.py`.

## Vaxrank: corpus acquisition and published comparator peptides

```python
from osteosarc import Region

# source is an explicitly selected RNA or DNA alignment product.
info = data.inspect_alignment(source)
if info.assembly == "GRCh38":
    panel = [Region("chr14", 101980428, 101980630, "GRCh38")]
    corpus = data.extract_reads(source, panel, fetch_pairs=True)
    print(corpus.receipt["scope"])
```

Survey headers before choosing native or explicitly lifted coordinates. Use
mate recovery for libraries where the corpus recipe requires it. For batches,
record per-source acquisition failures in the caller's manifest and continue
only according to that workflow's policy.

```python
published = data.vaccine_peptides("mRNA")
for peptide in published:
    print(peptide["variant_id"], peptide["sequence"])
```

Migration sites: `examples/osteosarc_read_corpus/build.py` and acquisition
portions of fixture support. Vaxrank's `run_vaxrank` continues to consume
Isovar results and its existing predictor/configuration objects; no additional
Osteosarc ranking wrapper is needed. Preserve reviewed extra calls, structural
events, source-selection rules, and biological expectations in Vaxrank.

## Tested handoffs

[`tests/test_consumer_usage.py`](https://github.com/iskandr/osteosarc/blob/main/tests/test_consumer_usage.py)
exercises native Varcode metadata/custom references, Isovar read evidence,
Topiary RSEM and pVAC loading (including empty reports), Vaxrank's header-survey
and paired-read acquisition pattern, and untested vaccine assay states.
The core suite also covers complete record preservation, indexed VCF access,
CRAM references, offline reuse, and conflicting source claims.

The dedicated consumer CI job installs Varcode 7.0.0, Isovar 1.17.0, and Topiary
5.55.1. It needs no genomes or prediction models. Full annotation, RNA protein
assembly, MHC prediction, and Vaxrank ranking belong to the consumer suites.
