# Explore osteosarc data

Osteosarc gives Varcode, Isovar, Topiary, and Vaxrank one Python API for the
public [osteosarc.com](https://osteosarc.com/data/) dataset. Discover files,
keep verified local copies, select variants and vaccines, and extract indexed
reads. Your existing libraries handle annotation and analysis.

## Install

Python 3.10+, Linux or macOS, and `curl` are required. Install from the repository:

```sh
git clone https://github.com/iskandr/osteosarc.git
cd osteosarc
python -m pip install -e '.[reads]'
```

Read inspection/extraction also needs `samtools` on PATH. Pin a Git commit in
reproducible consumer environments; this project is not yet published to PyPI.

## Create a snapshot once

```python
from osteosarc import Dataset

data = Dataset.sync("baseline")
print(data.id)
print(len(data.assets), len(data.variants()))
```

This downloads about 57 MB of metadata, including the dated bucket inventory and the timeline sources.
Large sequencing files are acquired only when explicitly requested. Files go
into the shared OpenVax cache (`OPENVAX_DATA_CACHE`, else the platform `openvax`
cache directory) under `objects/sha256/`, the layout vaxrank also uses, so the
OpenVax projects share downloaded bytes. `OSTEOSARC_CACHE` selects an isolated
cache.

## Reopen offline

```python
data = Dataset.open("baseline")
for variant in data.variants(gene="SMC5"):
    print(variant.id, variant.status, variant.vaccines)
```

Opening verifies the pinned metadata and defaults to offline operation. To
allow new downloads or regional reads, reopen with `offline=False`.

```python
data = Dataset.open("baseline", offline=False)
counts = data.table("vafs").select(gene="SMC5")
print(counts.columns)
print(counts.rows[:2])
```

## Where to go next

| Task | Short examples |
| --- | --- |
| See the whole dataset through worked questions | [A guided tour](tour.md) |
| Browse treatments, imaging, MRD, and specimens over time | [Timelines and the explorer](timeline.md) |
| Find samples, timepoints, assays, and files | [Explore the catalog](explore.md) |
| Select alleles, annotations, vaccine targets, and peptides | [Variants and vaccines](variants.md) |
| Inspect assembly and extract an indexed region union | [Read extraction](reads.md) |
| Pass data to existing analysis libraries | [Consumer recipes](consumers.md) |
| Replace duplicated acquisition helpers | [Migration guide](migration.md) |
| Review, disable, or extend source corrections | [Corrections and drift](curation.md) |

## Command line

```sh
osteosarc sync baseline
osteosarc assets baseline --assay rna-seq --timepoint T2
osteosarc variants baseline --set vaccine --status ready
osteosarc --offline table baseline vafs
osteosarc timeline baseline --since 2024-05 --until 2024-09
osteosarc specimens baseline T2_tumor
osteosarc curation baseline
```

`osteosarc explore baseline` opens an interactive terminal explorer.

The CLI uses the same Dataset methods as Python. `osteosarc --help` lists the
commands. [API contracts](design.md) and [validation](validation.md) describe
source limitations and reproducibility checks.
