# Browse the timeline and specimens

See treatments, procedures, imaging, specimens, MRD and lab results on one
timeline, and look up each specimen's files. These examples use the `baseline`
snapshot from [Get started](index.md).

## View events around a date

```python
from osteosarc import Dataset

data = Dataset.open("baseline")
print(data.timeline.around("2025-01-28", days=5).listing())
```

The timeline combines treatments, procedures, imaging, pathology, sampling,
MRD, and lab dates. Each event keeps its source record.

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

Filters accept dates as `YYYY`, `YYYY-MM`, or `YYYY-MM-DD`. Published date
precision is preserved. Timepoints appear only where a source states them.
Rows without readable dates are listed in `data.timeline.source["undated"]`.

## Find a specimen's files

```python
for row in data.specimens:
    print(row["sample_id"], row["date"], row["site"])

t2 = next(r for r in data.specimens if r["sample_id"] == "T2_tumor")
print(t2["assets"][:3])
print(t2["fastq_folders"][:3])
print(t2["corrections"], t2["disagreements"])
```

Specimen dates and sites are checked against other public sources. The rows
record disagreements and any [corrections](curation.md) applied.

## Read measurements

```python
labs = data.measurements.select(source="labs", measurement="WBC")
print(labs.rows[:3])
signatera = data.measurements.select(source="mrd", measurement="Signatera")
print([(r["date"], r["kind"], r["value"]) for r in signatera.rows[-5:]])
```

Values keep their raw strings and units. `kind` distinguishes `numeric`,
`not_detected`, `below_loq`, `text`, and `missing`. For `below_loq`, the value
is the reported limit, not a measured concentration.

## Use the terminal

```sh
osteosarc timeline baseline --since 2024-05 --until 2024-09
osteosarc timeline baseline --lane MRD --since 2025
osteosarc on baseline 2025-01-28 --days 5
osteosarc specimens baseline T2_tumor
osteosarc explore baseline
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

Dates come from the site's event JSON and timeline sheet, specimen registry,
T0–T3 summary, MRD records, flow-cytometry manifest, DICOM and pathology indexes,
and lab and cytometry tables. Snapshots predating these sources need to be
recreated under a new name to use the timeline.
