# Timelines, specimens, and the terminal explorer

A snapshot includes every public dated source. From them it builds one
chronological `Timeline`, a specimen registry linked to files, and a table of
measurements. All of it can be browsed as plain text in a terminal.

| Source | What it adds |
| --- | --- |
| `data/events.json` (+ the sheet export `timeline.csv`) | treatments (with ranges and dose numbers), procedures, imaging, pathology, omics, time points |
| `data/mrd.json` | Signatera, Northstar, and Personalis MRD values (`not detected` and below-LOQ kept distinct) |
| `scripts/data/samples-consolidated.tsv` | the specimen registry, linked to BAMs and FASTQ folders |
| `src/data/samples.json` | the homepage T0–T3 summary, used to cross-check the registry |
| `data/flow/manifest.json` | flow-cytometry blood-draw dates, with a `resolved` flag |
| `src/data/dicom-studies.json` | 136 DICOM imaging studies |
| `src/data/pathology-slides.json` | pathology specimens and slides |
| `data/lab_results.tsv`, `data/cytometry.tsv` | clinical labs and cytometry panels |

These add about 1.7 MB to `Dataset.sync`. Snapshots made before they existed
still open; their timeline raises a `SchemaError` that asks for a new snapshot.

## Timeline in Python

```python
from osteosarc import Dataset

data = Dataset.open("baseline")
timeline = data.timeline
print(len(timeline), timeline.lanes()[:6])
print(timeline.render(since="2024-05", until="2024-08", width=100))
```

The result looks like this, abbreviated:

```text
                              2024-05           2024-06         2024-07           2024-08
                              |                 |               |                 |
Time points                                       *T1
Specimens                                         3
Procedures: Biopsy                 *              *  *
Treatments: Radiation                                 ==============================
Treatments: Targeted therapy             2*==================*======================*===============
Treatments: Cancer vaccines    *   *                      *  *            *                *
Treatments: Steroids                   =====================
Imaging studies (DICOM)        *   *   *      2   *  *   * 4       *                     3      *2
MRD: Signatera                   o                   +         o        +          o        +
MRD: Personalis                                      +         +        +          +        +
Lab draws                      * * *     **  *    *                                *    *   *   *22

2024-05-01 .. 2024-08-31; one column = 1.8 days.  * event  2-9/# several  = range  > ongoing  MRD: + detected  o not detected  ~ below LOQ
```

Filters compose, and events keep their sources:

```python
vaccines = timeline.select(lane="Cancer vaccines")
print(vaccines.listing())                     # date, lane, label (with dose number), source
t1 = timeline.select(timepoint="T1")          # only where a source states the timepoint
week = timeline.around("2024-06-11", days=3)  # everything within three days
mrd = timeline.select(category="MRD", since="2025")
print(mrd.render(width=90))
event = vaccines[0]
print(event.date, event.end, event.precision, event.source, event.details)
```

Dates keep their published precision, and `event.precision` reports it.
Open-ended treatments (the website ends them at its build date) are marked
`open_end` and drawn with `>`. A timepoint is attached only when a source
states one. Procedure titles such as "T4/5" are vertebral levels, not
timepoints. `event.corrections` names any correction that flagged or edited
the underlying record. Source rows whose dates cannot be read are left off the
chart and listed in `timeline.source["undated"]` rather than failing the
timeline. Dates given to `select`, `around`, and `render` must be `YYYY`,
`YYYY-MM`, or `YYYY-MM-DD`; anything else raises `ValueError`.

## Specimens

```python
for row in data.specimens:
    print(row["sample_id"], row["timepoint"], row["date"], row["site"], len(row["assets"]))

t2 = next(r for r in data.specimens if r["sample_id"] == "T2_tumor")
print(t2["specimen"], t2["assays"], t2["corrections"])
print(t2["assets"][:3])          # BAM keys, via the consolidated metadata display names
print(t2["fastq_folders"][:3])
print(t2["disagreements"])       # other sources that give a different date or site
```

Registry dates and sites are cross-checked against the homepage summary and
the timeline's time points. With corrections off, T1–T3 show several
disagreements. With them on (see [corrections](curation.md)), those resolve to
the evidence-backed values, and each row names the correction that changed it.
The four PBMC capture-date specimens stay flagged rather than edited.

Measurements keep raw strings, units, and an explicit kind (`numeric`,
`not_detected`, `below_loq`, `text`, `missing`). A below-LOQ MRD value is the
limit itself, for example `"25.0"` with kind `below_loq`:

```python
labs = data.measurements.select(source="labs", measurement="WBC")
print(labs.rows[:3])
print(data.measurements.select(source="mrd").rows[:3])
```

## Command line

```sh
osteosarc timeline baseline --width 110                          # the whole history
osteosarc timeline baseline --since 2024-05 --until 2024-09      # zoom
osteosarc timeline baseline --lane MRD --since 2025              # some lanes
osteosarc timeline baseline --contains vaccine --list            # one line per event
osteosarc on baseline 2025-01-28 --days 5                        # what happened around a date
osteosarc specimens baseline                                     # the registry
osteosarc specimens baseline T1_tumor                            # one specimen: files, notes, nearby events
osteosarc timeline baseline --json --lane Biopsy                 # machine-readable events
```

## Interactive explorer

```sh
osteosarc explore baseline
```

```text
snapshot baseline (4580057f4269)
corrections: 29 applied
timeline: 787 events, 2022-10-06 .. 2026-09-16, 31 lanes
specimens: 21
variants: 182 on site, 177 ready, 3 missing_literal_allele, 2 non_literal_allele
Type help (or ?) for commands, quit to exit.
osteosarc> zoom 2024-05 2024-09
osteosarc> only MRD
osteosarc> events vaccine
osteosarc> on 2024-06-11 3
osteosarc> specimen T1_tumor
osteosarc> assets kind=alignment timepoint=T3
osteosarc> variants MAP2
osteosarc> corrections tempus-grch37-counts
osteosarc> quit
```

| Command | Does |
| --- | --- |
| `timeline [SINCE [UNTIL]] [lane=TEXT]`, `zoom`, `only TEXT`, `reset` | Draw, zoom, and filter the chart. The zoom and lane filter persist until `reset`. |
| `lanes` | Lane names with event counts |
| `events [TEXT]` | Events in the current zoom (all lanes) whose label contains `TEXT` |
| `on DATE [DAYS]` | Everything near a date |
| `specimens`, `specimen ID` | The registry; one specimen's files, corrections, and neighbouring events |
| `assets key=value ...` | `Assets.select` filters, e.g. `assay=rna-seq timepoint=T2` |
| `variants [GENE] [status=ready] [vaccine=mRNA]` | Catalogue entries with alleles and corrections |
| `corrections [ID]`, `summary` | Correction status and evidence; snapshot overview |

The same views are functions in `osteosarc.explore`, such as
`specimen_view(data, "T1_tumor")`. Scripts and notebooks can print exactly
what the shell shows.
