# Timeline

Treatments, procedures, scans, sample collection, MRD and lab results on one
timeline. The examples open your most recent snapshot; see
[Get started](index.md#get-started).

## Draw it

```sh
osteosarc timeline --since 2024-05 --until 2024-08
```

```python
from osteosarc import Dataset

data = Dataset.open()
print(data.timeline.render(since="2024-05", until="2024-08", width=90))
```

Each treatment gets a row, under its group, with rows for each MRD assay below:

```text
                                      2024
                                      May          Jun          Jul          Aug
Time points                                          T1
Samples collected                                    *
Imaging                               *     **    *        **                     *     *
Procedures
  Biopsy                                 *           * *
Radiation
  Proton therapy                                       =======================
Targeted therapy
  DeltaRex-G                                 ============================================
  ...
Cancer vaccines
  Peptide neoantigen vaccine (JLFv2)  *  *                *           *
MRD
  Signatera                             o              +      o      +       o      +
```

A star is an event or dose, = a treatment that continues, and > one still going.
For MRD, + is detected, o not detected, and ~ below the limit of quantification.
The whole timeline takes a month per column; a shorter window gets more.

Lab draws, DICOM studies, cytometry, flow draws, pathology slides and sequencing
runs happen too often to chart, so they're left out; the legend says how many.
Add `--all` (or `everything=True`) to include them.

## What happened around a date

```sh
osteosarc timeline --around 2025-01-28 --days 5
```

```python
print(data.timeline.around("2025-01-28", days=5).listing())
```

This lists every event within five days, one per line, with any
[corrections](corrections.md).

## Filter events

```python
mrd = data.timeline.select(category="MRD", since="2025")
print(mrd.render(width=90))
print(data.timeline.select(lane="Cancer vaccines").listing())
```

Dates can be a year, a month or a day, as in `since="2025-01"`; a date the site
gives as a month stays a month. On the command line, `--lane MRD` and
`--contains` filter the same way, and `--list` prints one line per event.

## Measurements

```python
signatera = data.measurements.select(source="mrd", measurement="Signatera")
print([(r["date"], r["kind"], r["value"]) for r in signatera.rows[-5:]])
```

Lab, cytometry and MRD values keep their published text and units, with a kind:
numeric, not detected, below the limit of quantification (the value is then the
lab's limit, not a measurement), text or missing.

## Where the dates come from

The site's events and timeline sheet, sample registry, T0–T3 summary, MRD
results, flow-cytometry list, imaging and pathology indexes, and lab and
cytometry tables. The Samples collected row marks when each [sample](samples.md) was
collected.
