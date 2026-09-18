"""Plain-text views of a Dataset and an interactive terminal explorer.

    osteosarc explore baseline

Every view is also a function returning a string, so the same output is
available from Python, the CLI subcommands, and the interactive shell.
"""

from __future__ import annotations

import cmd
import shlex
import shutil
from collections import Counter, defaultdict

from .errors import OsteosarcError
from .timeline import _window


def table(rows, columns, *, width=None, limit=None):
    """Fixed-width text table; long cells are truncated to fit the terminal."""
    if limit is not None and (not isinstance(limit, int) or isinstance(limit, bool) or limit < 0):
        raise ValueError("limit must be a nonnegative integer")
    rows = list(rows)
    shown = rows if limit is None else rows[:limit]
    cells = [[_text(row.get(c)) for c in columns] for row in shown]
    widths = [max([len(c)] + [len(r[i]) for r in cells]) for i, c in enumerate(columns)]
    width = width or shutil.get_terminal_size((120, 24)).columns
    # Shrink the widest columns to fit, but never the first (identifying) column.
    while sum(widths) + 2 * (len(widths) - 1) > width and max(widths[1:], default=0) > 12:
        widest = max(range(1, len(widths)), key=widths.__getitem__)
        widths[widest] -= 1
    line = "  ".join
    out = [line(c[:w].ljust(w) for c, w in zip(columns, widths)).rstrip(),
           line("-" * w for w in widths)]
    out += [line(cell[:w].ljust(w) for cell, w in zip(row, widths)).rstrip() for row in cells]
    if limit is not None and len(rows) > limit:
        out.append(f"... {len(rows) - limit} more")
    return "\n".join(out)


def _text(value):
    if value is None:
        return ""
    if isinstance(value, (tuple, list)):
        return "; ".join(_text(v) for v in value)
    return str(value)


def specimens_view(data, *, width=None):
    rows = [dict(row, bams=len(row["assets"]), fastq_dirs=len(row["fastq_folders"]),
                 notes=("disagrees: " + ", ".join(sorted({d["field"] for d in row["disagreements"]}))
                        if row["disagreements"] else "") + (" corrected: " + ", ".join(row["corrections"])
                                                            if row["corrections"] else ""))
            for row in data.specimens]
    return table(rows, ("sample_id", "timepoint", "date", "tissue", "site", "specimen", "assays",
                        "bams", "fastq_dirs", "notes"), width=width)


def specimen_view(data, sample_id, *, days=7, width=None):
    matches = [r for r in data.specimens if r["sample_id"] == sample_id]
    if not matches:
        known = ", ".join(r["sample_id"] for r in data.specimens)
        raise KeyError(f"No specimen {sample_id!r}; known: {known}")
    row = matches[0]
    lines = [f"{row['sample_id']}  {row['timepoint'] or '-'}  {row['date']}  {row['tissue']}  "
             f"{row['site']}  {row['specimen']}",
             f"vendors: {_text(row['vendors']) or '-'}   assays: {_text(row['assays']) or '-'}"]
    for item in row["disagreements"]:
        lines.append(f"  ! {item['source']} gives {item['field']} {item['value']!r}")
    summaries = {r["id"]: r["summary"] for r in data.corrections} if row["corrections"] else {}
    for correction in row["corrections"]:
        lines.append(f"  corrected by {correction}: {summaries[correction]}")
    by_assay = defaultdict(list)
    for key in row["assets"]:
        asset = data.asset(key)
        by_assay[_text(asset.values("assay")) or "unknown assay"].append(asset)
    lines.append(f"\nBAMs ({len(row['assets'])}):")
    for assay, assets in sorted(by_assay.items()):
        for asset in assets:
            size = f"{asset.size / 1e9:.1f} GB" if asset.size else "size ?"
            lines.append(f"  {assay:<10} {size:>9}  {asset.key}")
    if row["unmatched_bams"]:
        lines.append(f"  (named in the registry but not found: {_text(row['unmatched_bams'])})")
    if row["fastq_folders"]:
        lines.append(f"\nFASTQ folders ({len(row['fastq_folders'])}):")
        lines += [f"  {folder}" for folder in row["fastq_folders"]]
    if row["date"]:
        nearby = data.timeline.around(row["date"], days)
        lines.append(f"\nTimeline within {days} days of {row['date']}:")
        lines.append(nearby.listing() or "  (none)")
    return "\n".join(lines)


def assets_view(data, *, limit=40, width=None, **filters):
    selected = data.assets.select(**filters)
    rows = [dict(key=a.key, kind=a.kind, size=f"{a.size / 1e9:.2f} GB" if a.size else "",
                 timepoint=_text(a.values("timepoint")), assay=_text(a.values("assay")),
                 tissue=_text(a.values("tissue")), provider=_text(a.values("provider")),
                 notes=_text(a.metadata.get("corrections")))
            for a in selected]
    return f"{len(selected)} assets\n" + table(rows, ("timepoint", "assay", "tissue", "provider", "size",
                                                      "key", "notes"), width=width, limit=limit)


def variants_view(data, *, width=None, limit=None, **filters):
    rows = []
    for v in data.variants("all" if filters.pop("all", False) else "site", **filters):
        allele = f"{v.alleles[0][0]}:{v.alleles[0][1]} {v.alleles[0][2][:12]}>{v.alleles[0][3][:12]}" \
            if len(v.alleles) == 1 else ""
        rows.append(dict(id=v.id, gene=v.gene, status=v.status, allele=allele,
                         vaccines=_text(v.vaccines), pipelines=_text(v.pipelines),
                         corrections=_text(v.annotations.get("corrections"))))
    return table(rows, ("gene", "id", "status", "allele", "vaccines", "pipelines", "corrections"),
                 width=width, limit=limit)


def corrections_view(data, *, width=None):
    rows = [dict(r, touches=", ".join(sorted({c["source"] for c in r["changes"]})),
                 records=sum(c["records"] for c in r["changes"])) for r in data.corrections]
    return table(rows, ("status", "action", "id", "touches", "records", "summary"), width=width)


def summary_view(data):
    lines = [f"snapshot {data.manifest['name']} ({data.id[:12]})"]
    statuses = Counter(r["status"] for r in data.corrections)
    lines.append("corrections: " + ", ".join(f"{n} {s}" for s, n in sorted(statuses.items())))
    try:
        timeline = data.timeline
        first, last = timeline[0].date, max(e.last_day for e in timeline)
        lines.append(f"timeline: {len(timeline)} events, {first} .. {last}, {len(timeline.lanes())} lanes")
        lines.append(f"specimens: {len(data.specimens)}")
    except OsteosarcError as error:
        lines.append(f"timeline: unavailable ({error})")
    variants = data.variants()
    lines.append(f"variants: {len(variants)} on site, "
                 + ", ".join(f"{n} {s}" for s, n in Counter(v.status for v in variants).most_common()))
    kinds = Counter(a.kind for a in data.assets)
    lines.append("assets: " + ", ".join(f"{n:,} {k}" for k, n in kinds.most_common()))
    return "\n".join(lines)


BOOLEAN_OPTIONS = {"include_conflicts", "include_inferred", "all", "on_site", "vaccinated"}
INTEGER_OPTIONS = {"limit"}


def options(words):
    """key=value words as keyword arguments, with booleans and integers converted."""
    result = {}
    for word in words:
        if "=" not in word:
            raise ValueError(f"Expected key=value, not {word!r}")
        key, value = word.split("=", 1)
        if key in BOOLEAN_OPTIONS:
            if value.lower() not in ("true", "false", "yes", "no", "1", "0"):
                raise ValueError(f"{key} must be true or false, not {value!r}")
            result[key] = value.lower() in ("true", "yes", "1")
        elif key in INTEGER_OPTIONS:
            result[key] = int(value)
        else:
            result[key] = value
    return result


class Explorer(cmd.Cmd):
    """Interactive browsing of one snapshot. Type help <command> for details."""

    prompt = "osteosarc> "

    def __init__(self, data, *, width=None, stdin=None, stdout=None):
        super().__init__(stdin=stdin, stdout=stdout)
        if stdin is not None:
            self.use_rawinput = False
        self.data, self.width = data, width
        self.since = self.until = self.lane = None
        self.intro = summary_view(data) + "\nType help (or ?) for commands, quit to exit."

    def _print(self, text):
        self.stdout.write(text + "\n")

    def onecmd(self, line):
        try:
            return super().onecmd(line)
        except (OsteosarcError, KeyError, ValueError, TypeError, IndexError) as error:
            self._print(f"error: {error}")
            return False

    def emptyline(self):
        return False

    def do_summary(self, arg):
        """summary -- snapshot, correction, timeline, variant and asset counts."""
        self._print(summary_view(self.data))

    def do_timeline(self, arg):
        """timeline [SINCE [UNTIL]] [lane=TEXT] -- ASCII chart; dates are YYYY, YYYY-MM or YYYY-MM-DD.
        With no arguments, redraws the current zoom (see zoom, only, reset)."""
        words = shlex.split(arg)
        options = dict(w.split("=", 1) for w in words if "=" in w)
        dates = [w for w in words if "=" not in w]
        if dates:
            since, until = dates[0], dates[1] if len(dates) > 1 else None
            _window(since, until)  # validate before changing the zoom
            self.since, self.until = since, until
        if "lane" in options:
            self.lane = options["lane"] or None
        events = self.data.timeline.select(lane=self.lane) if self.lane else self.data.timeline
        self._print(events.render(since=self.since, until=self.until, width=self.width))

    def do_zoom(self, arg):
        """zoom SINCE [UNTIL] -- set the date window and redraw."""
        if not arg.strip():
            raise ValueError("zoom needs a start date, e.g. zoom 2024-05 2024-09")
        self.do_timeline(arg)

    def do_only(self, arg):
        """only TEXT -- show lanes whose name contains TEXT (e.g. only MRD); 'only' alone clears."""
        self.lane = arg.strip() or None
        self.do_timeline("")

    def do_reset(self, arg):
        """reset -- show the whole timeline and all lanes."""
        self.since = self.until = self.lane = None
        self.do_timeline("")

    def do_lanes(self, arg):
        """lanes -- lane names with event counts."""
        counts = Counter(e.lane for e in self.data.timeline)
        for lane in self.data.timeline.lanes():
            self._print(f"{counts[lane]:>5}  {lane}")

    def do_events(self, arg):
        """events [TEXT] -- list events in the current zoom (all lanes) whose label or lane contains TEXT."""
        events = self.data.timeline.select(contains=arg.strip() or None, since=self.since, until=self.until)
        self._print(events.listing() or "(no events)")

    def do_on(self, arg):
        """on DATE [DAYS] -- everything within DAYS (default 7) of DATE."""
        words = arg.split()
        if not words:
            raise ValueError("on needs a date, e.g. on 2024-06-11")
        days = int(words[1]) if len(words) > 1 else 7
        self._print(self.data.timeline.around(words[0], days).listing() or "(no events)")

    def do_specimens(self, arg):
        """specimens -- the specimen registry with linked BAM and FASTQ counts."""
        self._print(specimens_view(self.data, width=self.width))

    def do_specimen(self, arg):
        """specimen SAMPLE_ID -- details, linked files, disagreements and nearby events."""
        self._print(specimen_view(self.data, arg.strip(), width=self.width))

    def do_assets(self, arg):
        """assets [key=value ...] -- e.g. assets kind=alignment timepoint=T1 assay=rna-seq."""
        self._print(assets_view(self.data, width=self.width, **options(shlex.split(arg))))

    def do_variants(self, arg):
        """variants [GENE] [status=ready] [vaccine=mRNA] [pipeline=...] -- catalogue entries."""
        words = shlex.split(arg)
        filters = options(w for w in words if "=" in w)
        genes = [w for w in words if "=" not in w]
        if genes:
            filters["gene"] = genes[0]
        self._print(variants_view(self.data, width=self.width, **filters))

    def do_corrections(self, arg):
        """corrections [ID] -- status of every correction, or one correction's details."""
        if arg.strip():
            row = next((r for r in self.data.corrections if r["id"] == arg.strip()), None)
            if row is None:
                raise KeyError(f"No correction {arg.strip()!r}")
            self._print(f"{row['id']} [{row['status']}] {row['summary']}")
            for change in row["changes"]:
                self._print(f"  {change['state']:>15}  {change['source']}: {change['match']} "
                            f"({change['records']} records)")
            for item in row["evidence"]:
                self._print(f"  evidence: {item}")
            return
        self._print(corrections_view(self.data, width=self.width))

    def do_quit(self, arg):
        """quit -- leave the explorer."""
        return True

    do_exit = do_quit
    do_EOF = do_quit
