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
import textwrap
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from datetime import date, timedelta

from .cache import stable_id
from .display import Text
from .models import Collection

#: Display order of lanes; unlisted lanes follow alphabetically within their category.
LANE_ORDER = (
    "Time points", "Samples", "Procedures: Surgery", "Procedures: Biopsy", "Treatments: Surgery",
    "Treatments: Radiation", "Treatments: Chemotherapy", "Treatments: Targeted therapy",
    "Treatments: Immunotherapy", "Treatments: Cancer vaccines", "Treatments: Steroids",
    "Treatments: Antibiotics", "Imaging", "Imaging studies (DICOM)", "Pathology", "Pathology slides",
    "Omics: Genomics", "Omics: Transcriptomics", "Omics: Functional testing",
    "MRD: Signatera", "MRD: Northstar", "MRD: Personalis", "Flow cytometry draws",
    "Cytometry panels", "Lab draws", "Symptoms")
CATEGORY_ORDER = ("Time points", "Samples", "Procedures", "Treatments", "Imaging", "Pathology",
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
        return preview(f"Timeline: {self.overview()}. .render() charts it, .listing() lists every event, "
                       ".select(since=, until=, lane=, contains=) narrows it.",
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
        """One line per event: dates, lane, label, and any corrections."""
        width = max((len(e.lane) for e in self), default=0)
        lines = []
        for e in self:
            when = e.date + (f"..{e.end}" + ("+" if e.open_end else "") if e.end else "")
            label = e.label + (f" = {e.value}" if e.kind == "event" and e.value else "")
            notes = f"  (corrections: {', '.join(e.corrections)})" if e.corrections else ""
            lines.append(f"{when:<23} {e.lane:<{width}}  {label}{notes}".rstrip())
        return Text("\n".join(lines))

    def render(self, *, width=None, since=None, until=None, everything=False, legend=True):
        """A chart with one row per treatment and rows for time points, samples,
        procedures, scans and each MRD assay, over a month axis.

        Frequent records (lab draws, DICOM studies, cytometry, flow draws,
        slides and sequencing runs) are left out unless everything=True; the
        legend says how many. Every event is in listing().
        """
        events = list(self.select(since=since, until=until)) if (since or until) else list(self)
        if not events:
            return Text("(no events)")
        rows, hidden = _chart_rows(events, everything)
        if not any(label for _, label, _ in rows):
            rows, hidden = _chart_rows(events, True)
        lo, hi = _window(since, until)
        lo = date((lo or min(e.first_day for e in events)).year, (lo or min(e.first_day for e in events)).month, 1)
        hi = hi or max(e.last_day for e in events)
        labels = [("  " if heading is None and _indented(rows, i) else "") + label
                  for i, (heading, label, _) in enumerate(rows)]
        width = width or shutil.get_terminal_size((110, 24)).columns
        months = (hi.year - lo.year) * 12 + hi.month - lo.month + 1
        # Labels get what a column per month leaves, between 20 and 42 characters.
        label_width = min(max(len(label) for label in labels), 42, max(20, width - 2 - months))
        axis = _MonthAxis(lo, hi, max(24, width - label_width - 2))
        out = [" " * (label_width + 2) + line for line in axis.lines()]
        for (heading, _, items), label in zip(rows, labels):
            if heading is not None:
                out.append(heading)
                continue
            if len(label) > label_width:
                label = label[:label_width - 2] + ".."
            out.append(f"{label:<{label_width}}  {axis.draw(items, lo, hi)}".rstrip())
        if legend:
            out += ["", *textwrap.wrap("* event or dose   = continuing   > still ongoing   "
                                       "MRD: + detected  o not detected  ~ below LOQ", width,
                                       break_on_hyphens=False),
                    axis.scale()]
            if hidden:
                out += textwrap.wrap("Not shown: " + ", ".join(f"{n} {lane}" for lane, n in hidden)
                                     + ". Add --all (everything=True in Python) to include them.", width,
                                     break_on_hyphens=False)
        return Text("\n".join(out))


#: Lanes of frequent records, charted only with everything=True (with every
#: Omics lane: the sequencing runs, which follow the samples).
RECORD_LANES = ("Imaging studies (DICOM)", "Pathology slides", "Flow cytometry draws", "Cytometry panels",
                "Lab draws")
#: Chart names for the rows above the headed sections, in display order.
OVERVIEW_ROWS = {"Time points": "Time points", "Samples": "Samples collected", "Imaging": "Imaging",
                 "Pathology": "Pathology", "Symptoms": "Symptoms"}
#: Names of treatment groups that would otherwise read like a procedure row.
TREATMENT_HEADINGS = {"Surgery": "Surgical treatments", "Other": "Other treatments"}


def _chart_rows(events, everything):
    """Chart rows as (heading, label, events): a heading row has no label or events.

    Returns the rows and the (lane, count) pairs left out.
    """
    by_lane = defaultdict(list)
    for e in events:
        by_lane[e.lane].append(e)
    overview, procedures, treatments, mrd, records, hidden = [], [], defaultdict(dict), [], [], []
    for lane in Timeline(events).lanes():
        items = by_lane[lane]
        category = items[0].category
        if category == "Treatments":
            group = lane.split(": ", 1)[1] if ": " in lane else "Other"
            for e in items:
                treatments[group].setdefault(e.details.get("title") or e.label, []).append(e)
        elif category == "Procedures":
            procedures.append((None, lane.split(": ", 1)[1] if ": " in lane else "Other", items))
        elif category == "MRD":
            mrd.append((None, lane.split(": ", 1)[-1], items))
        elif lane in RECORD_LANES or category == "Omics":
            if everything:
                records.append((None, lane, items))
            else:
                hidden.append((lane, len(items)))
        else:
            overview.append((None, OVERVIEW_ROWS.get(lane, lane), items))
    order = list(OVERVIEW_ROWS.values())
    overview.sort(key=lambda row: order.index(row[1]) if row[1] in order else len(order))
    rows = list(overview)
    if procedures:
        rows += [("Procedures", "", ()), *sorted(procedures, key=lambda row: row[1] == "Other")]
    groups = [lane.split(": ", 1)[1] for lane in LANE_ORDER if lane.startswith("Treatments: ")]
    for group in sorted(treatments, key=lambda g: (g == "Other", groups.index(g) if g in groups else len(groups), g)):
        rows.append((TREATMENT_HEADINGS.get(group, group), "", ()))
        rows += [(None, name, items) for name, items in
                 sorted(treatments[group].items(), key=lambda item: min(e.first_day for e in item[1]))]
    if mrd:
        rows += [("MRD", "", ()), *mrd]
    if records:
        rows += [("Records", "", ()), *records]
    return rows, hidden


def _indented(rows, i):
    """Whether row i sits under a heading."""
    return any(heading is not None for heading, _, _ in rows[:i])


_MRD_MARK = {"numeric": "+", "not_detected": "o", "below_loq": "~", "missing": "?"}
_MRD_RANK = ["?", "o", "~", "+"]


class _MonthAxis:
    """Whole months across the plot: several columns per month, or several months per column."""

    def __init__(self, lo, hi, plot):
        self.lo = lo
        self.months = (hi.year - lo.year) * 12 + hi.month - lo.month + 1
        self.per_month = max(1, plot // self.months)
        self.per_column = max(1, -(-self.months // plot))
        self.width = (self.months * self.per_month if self.per_column == 1
                      else -(-self.months // self.per_column))

    def column(self, day):
        month = (day.year - self.lo.year) * 12 + day.month - self.lo.month
        if self.per_column > 1:
            return month // self.per_column
        days = calendar.monthrange(day.year, day.month)[1]
        return month * self.per_month + min(self.per_month - 1, (day.day - 1) * self.per_month // days)

    def _starts(self):
        for i in range(self.months):
            year, month = divmod(self.lo.month - 1 + i, 12)
            yield date(self.lo.year + year, month + 1, 1)

    def lines(self):
        years, months = [" "] * self.width, [" "] * self.width

        def put(cells, column, text):
            if column + len(text) <= self.width and all(c == " " for c in cells[max(0, column - 1):column + len(text)]):
                cells[column:column + len(text)] = text
                return True
            return False

        for i, start in enumerate(self._starts()):
            column = self.column(start)
            # The first month is labeled only when its year's label ends before January's.
            if start.month == 1 or (i == 0 and (start.month > 1 and self.column(
                    date(start.year + 1, 1, 1)) > column + 4 or self.months - i <= 12 - start.month)):
                put(years, column, str(start.year))
            if self.per_column > 1:
                if start.month == 1:
                    months[column] = "|"
            elif self.per_month >= 4:
                put(months, column, calendar.month_abbr[start.month])
            else:
                months[column] = calendar.month_abbr[start.month][0]
        return ["".join(years).rstrip(), "".join(months).rstrip()]

    def scale(self):
        if self.per_column > 1:
            return f"Each column is {self.per_column} months."
        if self.per_month == 1:
            return "Each column is a month."
        return f"Each month is {self.per_month} columns."

    def draw(self, events, lo, hi):
        cells = [" "] * self.width
        marks = {}
        for e in events:
            start, stop = self.column(max(e.first_day, lo)), self.column(min(e.last_day, hi))
            if e.end:
                for c in range(start, stop + 1):
                    if cells[c] == " ":
                        cells[c] = "="
                if e.open_end and e.last_day >= hi:
                    cells[stop] = ">"
        for e in events:
            start = self.column(max(e.first_day, lo))
            if e.kind == "measurement" and e.category == "MRD":
                marks[start] = max(marks.get(start, "?"), _MRD_MARK[e.details.get("value_kind")],
                                   key=_MRD_RANK.index)
            elif not e.end:
                cells[start] = "*"
        for c, mark in marks.items():
            cells[c] = mark
        for e in events:
            if e.lane == "Time points" and e.timepoint:
                c = self.column(max(e.first_day, lo))
                text = e.timepoint
                if c + len(text) <= self.width and all(x in " *" for x in cells[c:c + len(text)]) and (
                        c == 0 or cells[c - 1] == " "):
                    cells[c:c + len(text)] = text
        return "".join(cells)


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
    return [_event("specimens", row, n, date=day, lane="Samples", category="Samples", kind="sample",
                   timepoint=row.get("timepoint") or None,
                   label=f"{row['sample_id']}: {row.get('tissue_source', '')} ({row.get('collection_site', '')})",
                   links=(row["sample_id"],), corrections=corrections)
            for row, n, day, corrections in _dated("specimens", rows, "collection_date", marks, undated)]
