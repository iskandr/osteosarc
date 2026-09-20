# osteosarc

Python tools for working with the public [osteosarc.com](https://osteosarc.com/data/)
dataset. Find sequencing files, look up variants and vaccine peptides, and fetch
reads around a variant without downloading an entire BAM.

[Documentation](https://iskandr.github.io/osteosarc/) ·
[Examples for Varcode, Isovar, Topiary, and Vaxrank](https://iskandr.github.io/osteosarc/consumers/)

## Install

```sh
python -m pip install osteosarc
```

Requires Python 3.10+ and Linux or macOS. Read extraction also needs
`samtools` on PATH; see the [read extraction guide](https://iskandr.github.io/osteosarc/reads/).

## Get started

Save a snapshot and see which samples and sequencing types are available:

```python
from osteosarc import Dataset

data = Dataset.sync("baseline")
print(data.describe_samples())

rna = data.assets_for_sample("T0_tumor", kind="alignment", assay="rna-seq")
for asset in rna:
    print(asset.key, asset.size)

for variant in data.variants("vaccine", status="ready"):
    print(variant.gene, variant.allele)
```

The first sync downloads about 57 MB of metadata. BAMs and FASTQs stay remote
until you request them. To use the saved snapshot later, without network access:

```python
data = Dataset.open("baseline")
```

Downloads use datacache and the shared OpenVax cache. Set `OSTEOSARC_CACHE` to choose a
separate directory. Snapshot names are fixed; use a new name with
`Dataset.sync("next-snapshot", refresh=True)` to fetch updated metadata.

## Fetch reads around a variant

```python
data = Dataset.open("baseline", offline=False)
targets = data.variants(ids=["DYNC1H1-chr14-101980529"], status="ready")
source = rna[
    "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam"
]
subset = data.extract_reads(
    source, variants=targets, padding=100,
)
print(subset.path)  # Local BAM, with an index
```

The library checks the alignment's assembly and caches the result for reuse.
See [read extraction](https://iskandr.github.io/osteosarc/reads/) for filters,
paired mates, and local BAMs.

## Use the command line

```sh
osteosarc sync baseline
osteosarc samples baseline
osteosarc variants baseline --gene MAP2
osteosarc timeline baseline --since 2024-05 --until 2024-09
osteosarc explore baseline
```

The explorer lets you browse specimens, files, variants, and the timeline.
Type `help` for commands and `quit` to leave.

## About the data

Documented source corrections are applied by default. Use
`Dataset.open("baseline", corrections=False)` to inspect the published values.
A variant marked `ready` has a usable literal allele; this is not independent
validation. See [corrections](https://iskandr.github.io/osteosarc/curation/) for
the changes and their evidence.

Code is Apache-2.0. The dataset is listed as CC0-1.0 in the
[AWS Open Data Registry](https://registry.opendata.aws/sid-osteosarc/).
Cite the dataset and access date when using it.

## Development

```sh
python -m pip install -e '.[test]'
ruff check osteosarc tests scripts
python -m pytest -q
```

See [testing](https://iskandr.github.io/osteosarc/validation/) for build and live-example checks.
