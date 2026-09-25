"""Dated events from every public timeline source, plus plain-text views.

Events keep their published date precision ("2024-06" stays a month). A
timepoint (T0-T3) is attached only where a source states it; titles such as
"T4/5" are vertebral levels, never timepoints. Measurements (MRD, labs,
cytometry) keep their raw values; per-date summaries appear on the timeline.
"""

from __future__ import annotations

import calendar
import json
import re
import shutil
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import date, timedelta

from .cache import stable_id
from .display import Text
from .models import Collection

#: Display order of lanes; unlisted lanes follow alphabetically within their category.
LANE_ORDER = (
    "Time points", "Specimens", "Procedures: Surgery", "Procedures: Biopsy", "Treatments: Surgery",
    "Treatments: Radiation", "Treatments: Chemotherapy", "Treatments: Targeted therapy",
    "Treatments: Immunotherapy", "Treatments: Cancer vaccines", "Treatments: Steroids",
    "Treatments: Antibiotics", "Imaging", "Imaging studies (DICOM)", "Pathology", "Pathology slides",
    "Omics: Genomics", "Omics: Transcriptomics", "Omics: Functional testing",
    "MRD: Signatera", "MRD: Northstar", "MRD: Personalis", "Flow cytometry draws",
    "Cytometry panels", "Lab draws", "Symptoms")
CATEGORY_ORDER = ("Time points", "Specimens", "Procedures", "Treatments", "Imaging", "Pathology",
                  "Omics", "MRD", "Blood draws", "Symptoms")


def normalize_date(text):
    """ISO date at its published precision; accepts YYYYMMDD and MM/DD/YYYY."""
    text = (text or "").strip()
    if re.fullmatch(r"\d{8}", text):
        return f"{text[:4]}-{text[4:6]}-{text[6:]}"
    match = re.fullmatch(r"(\d{1,2})/(\d{1,2})/(\d{4})", text)
    if match:
        return f"{match[3]}-{int(match[1]):02d}-{int(match[2]):02d}"
    if re.fullmatch(r"\d{4}(-\d{2}(-\d{2})?)?", text):
        return text
    return None


def to_day(text, *, end=False):
    """First (or last, with end=True) day covered by a YYYY[-MM[-DD]] string."""
    parts = [int(p) for p in text.split("-")]
    year, month = parts[0], parts[1] if len(parts) > 1 else (12 if end else 1)
    day = parts[2] if len(parts) > 2 else (calendar.monthrange(year, month)[1] if end else 1)
    return date(year, month, day)


@dataclass(frozen=True)
class Event:
    """One dated source record. ``end`` is set for ranges (e.g. a treatment)."""

    id: str
    date: str
    lane: str
    category: str
    label: str
    source: str
    end: str | None = None
    kind: str = "event"
    subcategory: str = ""
    track: str = ""
    timepoint: str | None = None
    value: str | None = None
    open_end: bool = False
    links: tuple[str, ...] = ()
    corrections: tuple[str, ...] = ()
    details: dict = field(default_factory=dict, compare=False, repr=False)

    @property
    def first_day(self):
        return to_day(self.date)

    @property
    def last_day(self):
        return to_day(self.end or self.date, end=True)

    @property
    def precision(self):
        return ("year", "month", "day")[self.date.count("-")]


def parse_day(text, *, end=False):
    """First (or last) day of a YYYY, YYYY-MM or YYYY-MM-DD date; ValueError otherwise."""
    day = valid_date(str(text))
    if day is None:
        raise ValueError(f"Expected a date as YYYY, YYYY-MM or YYYY-MM-DD, not {text!r}")
    return to_day(day, end=end)


def _window(since, until):
    return (parse_day(since) if since else None, parse_day(until, end=True) if until else None)


class Timeline(Collection):
    """Chronological events with composable filters and text rendering."""

    def overview(self):
        """How many events, their date range, and how many lanes."""
        if not len(self):
            return "no events"
        last = max(e.last_day for e in self)
        return f"{len(self)} events, {self[0].date} .. {last}, {len(self.lanes())} lanes"

    def __repr__(self):
        from .display import preview
        if not len(self):
            return "Timeline: no events"
        return preview(f"Timeline: {self.overview()}. .render() draws it; .listing() lists every event.",
                       self, ("date", "lane", "event"),
                       lambda e: dict(date=e.date + (f"..{e.end}" if e.end else ""), lane=e.lane,
                                      event=e.label + (f" = {e.value}" if e.kind == "event" and e.value else "")))

    def select(self, *, lane=None, category=None, kind=None, source=None, track=None,
               timepoint=None, contains=None, since=None, until=None):
        """Filter events. lane/contains match case-insensitive substrings;
        since/until accept YYYY, YYYY-MM or YYYY-MM-DD and keep overlapping ranges."""
        lo, hi = _window(since, until)
        text = contains.lower() if contains else None
        return self.where(lambda e:
                          (lane is None or lane.lower() in e.lane.lower())
                          and (category is None or e.category == category)
                          and (kind is None or e.kind == kind)
                          and (source is None or e.source == source)
                          and (track is None or e.track == track)
                          and (timepoint is None or e.timepoint == timepoint)
                          and (text is None or text in f"{e.label} {e.lane} {e.subcategory}".lower())
                          and (lo is None or e.last_day >= lo) and (hi is None or e.first_day <= hi))

    def around(self, day, days=7):
        center = parse_day(day)
        lo, hi = center - timedelta(days=days), center + timedelta(days=days)
        return self.where(lambda e: e.last_day >= lo and e.first_day <= hi)

    def lanes(self):
        order = {name: i for i, name in enumerate(LANE_ORDER)}
        categories = {name: i for i, name in enumerate(CATEGORY_ORDER)}
        seen = {e.lane: e.category for e in self}
        return sorted(seen, key=lambda lane: (categories.get(seen[lane], len(categories)),
                                              order.get(lane, len(order)), lane))

    def listing(self):
        """One line per event: dates, lane, label, source, and any corrections."""
        lines = []
        for e in self:
            when = e.date + (f"..{e.end}" + ("+" if e.open_end else "") if e.end else "")
            label = e.label + (f" = {e.value}" if e.kind == "event" and e.value else "")
            notes = f" (corrections: {', '.join(e.corrections)})" if e.corrections else ""
            lines.append(f"{when:<23} {e.lane[:30]:<30} {label}  [{e.source}]{notes}")
        return Text("\n".join(lines))

    def render(self, *, width=None, since=None, until=None, legend=True):
        """ASCII lanes over a shared date axis.

        *  one event  2-9 events in one column  #  ten or more  = range  > ongoing
        +  MRD detected  o  not detected  ~  below limit of quantification
        """
        events = list(self.select(since=since, until=until)) if (since or until) else list(self)
        if not events:
            return "(no events)"
        lo, hi = _window(since, until)
        lo = lo or min(e.first_day for e in events)
        hi = hi or max(e.last_day for e in events)
        lanes = Timeline(events).lanes()
        label_width = min(30, max(len(lane) for lane in lanes))
        width = width or shutil.get_terminal_size((110, 24)).columns
        plot = max(20, width - label_width - 2)
        span = max(1, (hi - lo).days)

        def column(day):
            return min(plot - 1, max(0, round((day - lo).days * (plot - 1) / span)))

        rows = [_axis(lo, hi, plot, column, label_width)]
        for lane in lanes:
            cells, counts, marks = [" "] * plot, Counter(), {}
            for e in (e for e in events if e.lane == lane):
                start, stop = column(max(e.first_day, lo)), column(min(e.last_day, hi))
                if e.end:
                    for c in range(start, stop + 1):
                        cells[c] = "="
                    if e.open_end and e.last_day >= hi:
                        cells[stop] = ">"
                elif e.kind == "measurement" and e.category == "MRD":
                    marks[start] = max(marks.get(start, "?"), _MRD_MARK[e.details.get("value_kind")],
                                       key=_MRD_RANK.index)
                else:
                    counts[start] += 1
            for c, n in counts.items():
                cells[c] = "*" if n == 1 else str(n) if n < 10 else "#"
            for c, mark in marks.items():
                cells[c] = mark
            if lane == "Time points":
                for e in (e for e in events if e.lane == lane and e.timepoint):
                    c = column(e.first_day)
                    if c + 1 + len(e.timepoint) <= plot and all(x == " " for x in cells[c + 1:c + 1 + len(e.timepoint)]):
                        cells[c + 1:c + 1 + len(e.timepoint)] = e.timepoint
            rows.append(f"{lane[:label_width]:<{label_width}}  {''.join(cells).rstrip()}")
        if legend:
            days = span / max(1, plot - 1)
            rows.append("")
            rows.append(f"{lo} .. {hi}; one column = {days:.1f} days.  * event  2-9/# several  "
                        "= range  > ongoing  MRD: + detected  o not detected  ~ below LOQ")
        return Text("\n".join(rows))


_MRD_MARK = {"numeric": "+", "not_detected": "o", "below_loq": "~", "missing": "?"}
_MRD_RANK = ["?", "o", "~", "+"]


def _axis(lo, hi, plot, column, label_width):
    """Year (or month, for short windows) labels above tick marks."""
    labels, ticks = [" "] * plot, [" "] * plot
    months = (hi - lo).days <= 540
    day = date(lo.year, lo.month, 1) if months else date(lo.year, 1, 1)
    while day <= hi:
        if day >= lo:
            c = column(day)
            text = day.strftime("%Y-%m" if months else "%Y")
            if all(x == " " for x in labels[max(0, c - 1):c + len(text) + 1]) and c + len(text) <= plot:
                labels[c:c + len(text)] = text
                ticks[c] = "|"
        day = date(day.year + (day.month == 12), day.month % 12 + 1, 1) if months else date(day.year + 1, 1, 1)
    pad = " " * (label_width + 2)
    return pad + "".join(labels).rstrip() + "\n" + pad + "".join(ticks).rstrip()


# ---------------------------------------------------------------------------
# Building events from snapshot sources
# ---------------------------------------------------------------------------

def _event(source, raw, occurrence, **fields):
    key = json.dumps(raw, sort_keys=True, default=str)
    return Event(id=stable_id([source, key, occurrence])[:16], source=source, details=dict(raw), **fields)


def valid_date(text):
    """A normalized date that also exists on the calendar, else None."""
    day = normalize_date(text) if isinstance(text, str) else None
    try:
        return day if day and to_day(day) else None
    except ValueError:
        return None


def _dated(source, rows, field, marks=None, undated=None):
    """Yield (row, occurrence, date, correction IDs); rows without a usable date go to undated."""
    seen = Counter()
    for i, row in enumerate(rows):
        key = json.dumps(row, sort_keys=True, default=str)
        seen[key] += 1
        day = field(row) if callable(field) else valid_date(row.get(field))
        if day is None:
            if undated is not None:
                undated.append(dict(source=source, record=row))
            continue
        yield row, seen[key] - 1, day, tuple((marks or {}).get(i, ()))


def events_from_site(document, sheet_rows=(), *, marks=None, sheet_marks=None, undated=None):
    """events.json rows; dose and dose number joined from the sheet export when unique."""
    doses = defaultdict(list)
    for i, row in enumerate(sheet_rows):
        start = normalize_date(row.get("Start date"))
        doses[(row.get("Title", "").strip(), start, row.get("Category", "").strip())].append(
            (row, tuple((sheet_marks or {}).get(i, ()))))
    build_end = document.get("date_range", {}).get("end")
    names = {t["id"]: t["name"] for t in document.get("timelines", ())}
    result = []
    for row, n, day, corrections in _dated("events", document["events"], "date", marks, undated):
        category, subcategory = row.get("category", ""), row.get("subcategory", "")
        group, title = row.get("group", ""), row.get("title", "")
        timepoint = None
        if category == "Omics" and group == "Time point" and re.fullmatch(r"T\d+", title):
            lane, timepoint = "Time points", title
            category = "Time points"
        elif category == "Treatments" and group == "Cancer vaccines":
            lane = "Treatments: Cancer vaccines"
        elif category in ("Treatments", "Procedures", "Omics") and subcategory:
            lane = f"{category}: {subcategory}"
        else:
            lane = category or "Other"
        matches = doses.get((title.strip(), day, row.get("category", "")), [])
        dose, sheet_corrections = matches[0] if len(matches) == 1 else ({}, ())
        corrections = tuple(dict.fromkeys((*corrections, *sheet_corrections)))
        end = valid_date(row.get("end_date"))
        result.append(_event(
            "events", row, n, date=day, end=end, lane=lane, category=category or "Other",
            subcategory=subcategory, label=title + (f" ({group})" if group and group not in (
                subcategory, title, "Time point") else ""),
            track=names.get(row.get("timeline"), row.get("timeline", "")), timepoint=timepoint,
            open_end=bool(end) and end == build_end, corrections=corrections,
            value=" ".join(x for x in (dose.get("Dose", ""), f"#{dose['Dose #']}" if dose.get("Dose #") else "")
                           if x) or None))
    return result


def events_from_mrd(document, *, marks=None, undated=None):
    """MRD values: the raw value as a string, with its assay unit and an explicit kind."""
    assays = document.get("assays", {})
    result = []
    for row, n, day, corrections in _dated("mrd", document["measurements"], "date", marks, undated):
        assay = assays.get(row.get("assay"), {})
        name, unit = assay.get("name", row.get("assay", "")), assay.get("unit", "")
        value = row.get("value")
        if isinstance(value, dict) and "below_loq" in value:
            kind, raw, text = "below_loq", str(value["below_loq"]), f"<LOQ ({value['below_loq']} {unit})"
        elif value == "not_detected":
            kind, raw, text = "not_detected", "not_detected", "not detected"
        elif value in (None, ""):
            kind, raw, text = "missing", "", "missing"
        else:
            kind, raw, text = "numeric", str(value), f"{value} {unit}".strip()
        result.append(_event("mrd", dict(row, value_kind=kind, unit=unit, assay_name=name), n, date=day,
                             lane=f"MRD: {name}", category="MRD", kind="measurement",
                             label=f"{name} {text}", value=raw, corrections=corrections))
    return result


def events_from_measurements(rows, marks, source, lane, *, undated=None):
    """One event per date summarizing measurements (labs, cytometry)."""
    by_date = defaultdict(list)
    for row, _, day, corrections in _dated(source, rows, "date", marks, undated):
        by_date[day].append((row, corrections))
    result = []
    for day, items in sorted(by_date.items()):
        flagged = sum(r.get("out_of_range", "").upper() == "TRUE" for r, _ in items)
        categories = sorted({r.get("category", "") for r, _ in items} - {""})
        label = f"{len(items)} measurements" + (f", {flagged} out of range" if flagged else "")
        corrections = tuple(dict.fromkeys(c for _, ids in items for c in ids))
        result.append(_event(source, dict(date=day, measurements=len(items), out_of_range=flagged,
                                          categories=categories), 0,
                             date=day, lane=lane, category="Blood draws", kind="measurement", label=label,
                             corrections=corrections))
    return result


def events_from_flow(manifest, *, marks=None, undated=None):
    return [_event("flow", row, n, date=day, lane="Flow cytometry draws", category="Blood draws",
                   label=f"Flow draw {row.get('label', day)}"
                         + ("" if row.get("resolved", True) else " (date provisional)"),
                   links=(row["id"],) if row.get("id") else (), corrections=corrections)
            for row, n, day, corrections in _dated("flow", manifest["samples"], "draw_date", marks, undated)]


def events_from_imaging(document, *, marks=None, undated=None):
    result = []
    for row, n, day, corrections in _dated("imaging", document["studies"], "studyDate", marks, undated):
        summary = {k: row.get(k) for k in ("slug", "studyInstanceUID", "title", "modality",
                                           "seriesCount", "instanceCount", "viewerUrl", "impression")}
        result.append(_event("imaging", summary, n, date=day, lane="Imaging studies (DICOM)",
                             category="Imaging", subcategory=row.get("modality", ""),
                             label=f"{row.get('modality', '')} {row.get('title', '')}".strip(),
                             corrections=corrections))
    return result


def _specimen_date(group):
    match = re.match(r"(\d{4}-\d{2}-\d{2})_", group.get("specimen", ""))
    return valid_date(match[1]) if match else None


def events_from_pathology(document, *, marks=None, undated=None):
    result = []
    for group, n, day, corrections in _dated("pathology", document["groups"], _specimen_date, marks, undated):
        kind = group["specimen"].split("_", 1)[1]
        result.append(_event("pathology", {k: group.get(k) for k in ("id", "title", "specimen", "ids")}, n,
                             date=day, lane="Pathology slides", category="Pathology", subcategory=kind,
                             label=f"{group.get('title', group['specimen'])} ({len(group.get('ids', ()))} slides)",
                             links=tuple(group.get("ids", ())), corrections=corrections))
    return result


def events_from_specimens(rows, *, marks=None, undated=None):
    return [_event("specimens", row, n, date=day, lane="Specimens", category="Specimens", kind="specimen",
                   timepoint=row.get("timepoint") or None,
                   label=f"{row['sample_id']}: {row.get('tissue_source', '')} ({row.get('collection_site', '')})",
                   links=(row["sample_id"],), corrections=corrections)
            for row, n, day, corrections in _dated("specimens", rows, "collection_date", marks, undated)]
