# Timeline

See treatments, procedures, scans, sample collection, MRD and lab results on one
timeline. The examples open your most recent snapshot; see
[Get started](index.md#get-started) to save one.

## Draw the timeline

```python
from osteosarc import Dataset

data = Dataset.open()
print(data.timeline.render(since="2024-05", until="2024-08", width=90))
```

Each treatment gets its own row, under its group (chemotherapy, immunotherapy,
cancer vaccines and so on), above rows for each MRD assay:

```text
                                      2024
                                      May          Jun          Jul          Aug
Time points                                          T1
Samples collected                                    *
Imaging                               *     **    *        **                     *     *
Pathology                                *           * *
Procedures
  Biopsy                                 *           * *
  Other                                    *
Radiation
  Proton therapy                                       =======================
Chemotherapy
  Trabectedin                                 *
Targeted therapy
  DeltaRex-G                                 ============================================
  ...
Cancer vaccines
  Peptide neoantigen vaccine (JLFv2)  *  *                *           *
  Peptide neoantigen vaccine (CeGaT)                        *                      *
MRD
  Signatera                             o              +      o      +       o      +
  Northstar                                                   +      +       +      +
```

`*` is an event or dose, `=` a treatment that continues over time, and `>` one
still going. MRD uses `+` for detected, `o` for not detected, and `~` for below the
limit of quantification. The whole timeline fits a month per column; a shorter
window gets several columns per month.

Lab draws, DICOM studies, cytometry panels, flow draws, pathology slides and
sequencing runs happen too often to chart usefully, so they're left out; the
legend says how many. `render(everything=True)` includes them, as does asking for
one with `select`.

## View events around a date

```python
print(data.timeline.around("2025-01-28", days=5).listing())
```

`listing()` prints one line per event, with its date or date range, lane and any
[corrections](corrections.md). Each event keeps the record it came from.

## Filter events

```python
vaccines = data.timeline.select(lane="Cancer vaccines")
print(vaccines.listing())

mrd = data.timeline.select(category="MRD", since="2025")
print(mrd.render(width=90))

print(data.timeline.select(timepoint="T1").listing())
event = vaccines[0]
print(event.date, event.precision, event.source, event.details)
```

Dates can be `YYYY`, `YYYY-MM` or `YYYY-MM-DD`. A date published as just a month
stays a month. Events have a time point only where the site gives one. Rows with a
date that can't be read are in `data.timeline.source["undated"]`.

## Read measurements

```python
labs = data.measurements.select(source="labs", measurement="WBC")
print(labs.rows[:3])
signatera = data.measurements.select(source="mrd", measurement="Signatera")
print([(r["date"], r["kind"], r["value"]) for r in signatera.rows[-5:]])
```

Values keep their published text and units. `kind` is `numeric`, `not_detected`,
`below_loq`, `text` or `missing`; for `below_loq`, the value is the lab's reporting
limit, not a measurement.

## Samples on the timeline

The "Samples collected" row marks each [sample](samples.md)'s collection date.
`data.samples["T2_tumor"]` shows a sample's files, and
`data.timeline.around("2025-01-28").listing()` what happened around it. A sample's
`disagreements` lists dates or sites that the site's other pages give differently.

## Use the terminal

```sh
osteosarc timeline
osteosarc timeline --since 2024-05 --until 2024-09
osteosarc timeline --lane MRD --since 2025 --list
osteosarc timeline --all
osteosarc on 2025-01-28 --days 5
```

`--list` prints one line per event and `--json` the event records. `--lane` and
`--contains` also chart lanes that are otherwise left out, such as
`--lane "Lab draws"`.

## Sources

Dates come from the site's events and timeline sheet, sample registry, T0–T3
summary, MRD results, flow-cytometry list, imaging and pathology indexes, and lab
and cytometry tables. A snapshot saved before the timeline existed can't show it;
run `osteosarc sync --refresh` for a new one.
