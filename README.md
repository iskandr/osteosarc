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

## Find samples and sequencing files

```python
from osteosarc import Dataset

data = Dataset.sync("baseline")
print(data.describe_samples())
```

`baseline` is a name you choose for the local metadata snapshot. The first sync
fetches about 57 MB; sequencing files stay remote until requested.

`T0_tumor` is the primary tumor specimen from the **T0 collection timepoint**
(2022-12-16). `T0_blood` is blood from that same timepoint. Sample type is the
`tissue` field; sequencing assay is a separate choice:

| Data available for `T0_tumor` | Assay filter |
| --- | --- |
| Bulk RNA sequencing | `rna-seq` |
| Bulk whole-exome DNA sequencing | `wes` |
| Bulk whole-genome DNA sequencing | `wgs` |

Find its bulk RNA alignments:

```python
rna = data.assets_for_sample("T0_tumor", kind="alignment", assay="rna-seq")
for asset in rna:
    print(asset.key)
```

`rna-seq` means **bulk RNA**; `scrna-seq` means **single-cell RNA**. A specimen
can have both, as `T1_tumor` does. `platform="ont"` or `"pacbio"` selects a
sequencing technology separately. See [sample IDs and assay names](https://iskandr.github.io/osteosarc/explore/)
for the full vocabulary and platform availability.

## Select variants

```python
targets = data.variants(gene="DYNC1H1", status="ready")
for variant in targets:
    print(variant.id, variant.allele)
```

Variant `status` describes whether its genomic allele is usable. `ready` means
one consistent chromosome, position, REF and ALT, with literal DNA bases.
Read support, somatic status and protein effects need separate analysis.
Omit the filter to include unresolved entries; see [all statuses](https://iskandr.github.io/osteosarc/variants/#variant-status).
Alleles are `(chromosome, one-based position, REF, ALT)`.

## Fetch reads around those variants

```python
source = rna[
    "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam"
]
reads = data.extract_reads(source, variants=targets, padding=100)
print(reads.path)  # Local indexed BAM
```

This retrieves overlapping reads from the chosen BAM and checks its assembly.
Use [Isovar](https://iskandr.github.io/osteosarc/consumers/#isovar) to classify
reference- and alternate-supporting reads.

Downloads use datacache and the shared OpenVax cache. The same extraction
request reuses its cached result. To reopen the snapshot offline:

```python
data = Dataset.open("baseline")
```

Use `offline=False` to acquire more data. Set `OSTEOSARC_CACHE` for a separate
cache directory. To refresh metadata, choose a new snapshot name with
`Dataset.sync("next-snapshot", refresh=True)`.

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
See [corrections](https://iskandr.github.io/osteosarc/curation/) for
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
