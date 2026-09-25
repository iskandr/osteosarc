# Browse the timeline and specimens

See treatments, procedures, imaging, specimens, MRD and lab results on one
timeline, and look up each specimen's files. The examples open your most recent snapshot;
see [Get started](index.md#get-started) to save one.

## View events around a date

```python
from osteosarc import Dataset

data = Dataset.open()
print(data.timeline.around("2025-01-28", days=5).listing())
```

The timeline puts treatments, procedures, imaging, pathology, sample collection,
MRD and lab results on one axis. Each event keeps the record it came from.

## Draw a timeline

```python
print(data.timeline.render(since="2024-05", until="2024-08", width=100))
```

The text chart uses `*` for events, `=` for ranges, and `>` for ongoing treatment.
MRD uses `+` for detected, `o` for not detected, and `~` for below the limit of
quantification.

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
stays a month. Events have a timepoint only where the site gives one. Rows with a
date that can't be read are in `data.timeline.source["undated"]`.

## Find a specimen's files

```python
for row in data.specimens:
    print(row["sample_id"], row["date"], row["site"])

t2 = next(r for r in data.specimens if r["sample_id"] == "T2_tumor")
print(t2["assets"][:3])
print(t2["fastq_folders"][:3])
print(t2["corrections"], t2["disagreements"])
```

Each specimen's date and site are checked against the site's other pages, and
the row lists any disagreements and [corrections](curation.md).

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

## Use the terminal

```sh
osteosarc timeline --since 2024-05 --until 2024-09
osteosarc timeline --lane MRD --since 2025
osteosarc on 2025-01-28 --days 5
osteosarc specimens T2_tumor
osteosarc explore
```

Inside the explorer:

```text
osteosarc> zoom 2024-05 2024-09
osteosarc> only MRD
osteosarc> events vaccine
osteosarc> specimen T1_tumor
osteosarc> assets kind=alignment timepoint=T3
osteosarc> variants MAP2
osteosarc> reset
osteosarc> quit
```

`zoom` and `only` persist until `reset`. Type `help` for all commands.
Add `--json` to the `timeline` CLI command for machine-readable events. The
[command-line guide](cli.md#interactive-explorer) lists every explorer command.

## Sources

Dates come from the site's events and timeline sheet, specimen registry, T0–T3
summary, MRD results, flow-cytometry list, imaging and pathology indexes, and lab
and cytometry tables. A snapshot saved before the timeline existed can't show it;
run `osteosarc sync --refresh` for a new one.
