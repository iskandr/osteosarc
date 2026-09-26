"""Plain-text views of a Dataset, shared by the CLI and Python reprs.

Every view is a function returning a string, so the command line and a
Python session show the same thing.
"""

from __future__ import annotations

import shutil
import textwrap
from collections import Counter, defaultdict
from pathlib import PurePosixPath

from .curation import ASSAY_NAMES, ASSAYS
from .errors import OsteosarcError
from .urls import BUCKET, S3_BUCKET

#: Display order of file kinds.
KINDS = ("alignment", "reads", "variants", "expression", "annotation", "table", "reference", "index",
         "image", "other")
KIND_NAMES = {"alignment": "BAM/CRAM aligned reads", "reads": "FASTQ raw reads", "variants": "VCF/BCF calls",
              "expression": "expression matrices and counts", "annotation": "GTF/GFF/BED annotation",
              "table": "CSV/TSV/JSON tables", "reference": "reference sequences",
              "index": "indexes for BAMs, VCFs and FASTAs", "image": "scans, microscopy and slides",
              "other": "logs, raw signal and the rest"}


def table(rows, columns, *, width=None, limit=None, wrap=False, fixed=(), total=None):
    """Fixed-width text table; long cells are truncated to fit the terminal.

    Columns in fixed, such as file keys people copy, are never shortened. total
    is the full row count when rows holds only the first few.
    """
    if limit is not None and (not isinstance(limit, int) or isinstance(limit, bool) or limit < 0):
        raise ValueError("limit must be a nonnegative integer")
    rows = list(rows)
    shown = rows if limit is None else rows[:limit]
    cells = [[_text(row.get(c)) for c in columns] for row in shown]
    widths = [max([len(c)] + [len(r[i]) for r in cells]) for i, c in enumerate(columns)]
    width = width or shutil.get_terminal_size((120, 24)).columns
    # Shrink the widest columns to fit, but never the first (identifying) column.
    shrinkable = [i for i in range(1, len(widths)) if columns[i] not in fixed]
    while sum(widths) + 2 * (len(widths) - 1) > width and max((widths[i] for i in shrinkable), default=0) > 12:
        widest = max(shrinkable, key=widths.__getitem__)
        widths[widest] -= 1
    line = "  ".join
    out = [line(c[:w].ljust(w) for c, w in zip(columns, widths)).rstrip(),
           line("-" * w for w in widths)]
    for row in cells:
        parts = [textwrap.wrap(cell, width=w, break_on_hyphens=False) or [""] if wrap else [cell[:w]]
                 for cell, w in zip(row, widths)]
        for i in range(max(map(len, parts), default=0)):
            out.append(line((part[i] if i < len(part) else "").ljust(w)
                            for part, w in zip(parts, widths)).rstrip())
    total = len(rows) if total is None else total
    if limit is not None and total > limit:
        out.append(f"... {total - limit:,} more")
    return "\n".join(out)


def _text(value):
    if value is None:
        return ""
    if isinstance(value, (tuple, list)):
        return ", ".join(_text(v) for v in value)
    return str(value)


def size_text(n):
    if not n:
        return ""
    for unit, scale in (("TB", 1e12), ("GB", 1e9), ("MB", 1e6), ("KB", 1e3)):
        if n >= scale:
            return f"{n / scale:.1f} {unit}"
    return f"{n} B"


def sequencing_text(pairs):
    """(assay, platform) pairs as 'rna-seq, wgs, scrna-seq (ont, pacbio)'."""
    platforms = defaultdict(set)
    for assay, platform in pairs:
        platforms[assay].update([platform] if platform else [])
    order = {name: i for i, name in enumerate(ASSAY_NAMES)}
    return ", ".join(assay + (f" ({', '.join(sorted(platforms[assay]))})" if platforms[assay] else "")
                     for assay in sorted(platforms, key=lambda a: (order.get(a, len(order)), a)))


def allele_text(variant):
    """chrom:pos REF>ALT, with long alleles shortened; blank without a single allele."""
    if len(variant.alleles) != 1:
        return ""
    chrom, pos, ref, alt = variant.alleles[0]
    return f"{chrom}:{pos} {_short(ref)}>{_short(alt)}"


def _short(allele, n=12):
    return allele if len(allele) <= n else allele[:n - 2] + ".."


def snapshot_line(data):
    when = (data.downloaded or "")[:16].replace("T", " ")
    return (f"snapshot {data.name} ({data.id[:12]}), "
            + (f"downloaded {when} UTC" if when else "download time unknown"))


def regions_text(regions):
    """Zero-based region records as one-based samtools text: 'chr2:10-20', or '3 regions from chr2:10-20'."""
    if not regions:
        return ""
    first = f"{regions[0]['contig']}:{regions[0]['start'] + 1}-{regions[0]['end']}"
    return first if len(regions) == 1 else f"{len(regions)} regions from {first}"


def plural(n, word):
    return f"{n:,} {word}" + ("" if n == 1 else "s")


def hints(lines):
    """Aligned 'command   what it does' lines."""
    width = max(len(command) for command, _ in lines)
    return "\n".join(f"  {command:<{width}}  {what}" for command, what in lines)


# ---------------------------------------------------------------------------
# Samples
# ---------------------------------------------------------------------------

def samples_view(samples, *, width=None, footer=None):
    if not len(samples):
        return "(no matching samples)"
    tissues = Counter(s.tissue or "unknown" for s in samples).most_common()
    title = (f"{plural(len(samples), tissues[0][0] + ' sample')}" if len(tissues) == 1 else
             f"{len(samples)} samples: " + ", ".join(f"{n} {t}" for t, n in tissues))
    rows = [dict(sample=s.id, timepoint=s.timepoint or "", date=s.date, where=s.site,
                 sequencing=sequencing_text(s.sequencing) or "none recorded", BAMs=len(s.bams),
                 **{"FASTQ folders": len(s.fastq_folders)}) for s in samples]
    shown = table(rows, ("sample", "timepoint", "date", "where", "sequencing", "BAMs", "FASTQ folders"),
                  width=width, wrap=True)
    footer = footer if footer is not None else hints([
        ('data.samples["T1_tumor"]', "one sample's files, with how to get them"),
        ('.select(tissue="blood", assay="cite-seq")', "filter by timepoint, tissue, assay or platform")])
    return f"{title}\n\n{shown}\n\n{footer}"


def _fastq_folders(sample, files):
    """One row per FASTQ folder: its library label, provider, and files and bytes in the bucket."""
    labels = {}
    for row in sample.details.get("fastqs", ()):
        labels.setdefault(row["folder"], (row.get("assay") or "", row.get("provider") or ""))
    counts, sizes = Counter(), Counter()
    if files is not None:
        folders = set(sample.fastq_folders)
        for file in files:
            parts = file.key.split("/")
            for depth in range(len(parts) - 1, 0, -1):
                folder = "/".join(parts[:depth])
                if folder in folders:
                    counts[folder] += 1
                    sizes[folder] += file.size or 0
                    break
    rows = []
    for folder in sample.fastq_folders:
        label, provider = labels.get(folder, ("", ""))
        rows.append(dict(library=label, provider=provider, files=counts[folder] if files is not None else "",
                         size=size_text(sizes[folder]), folder=folder + "/"))
    order = {name: i for i, name in enumerate(ASSAY_NAMES)}
    return sorted(rows, key=lambda r: (order.get(ASSAYS.get(r["library"], (r["library"],))[0], len(order)),
                                       r["library"], r["folder"]))


def _bam_rows(sample, data, local):
    rows = []
    for key in sample.bams:
        try:
            file = data.file(key) if data is not None else None
        except KeyError:
            file = None
        rows.append(dict(
            assay=_text(file.values("assay")) if file else "", platform=_text(file.values("platform")) if file else "",
            provider=_text(file.values("provider")) if file else "", size=size_text(file.size) if file else "",
            local="yes" if file and file.url in local else "", key=key, indexed=bool(file and file.index_urls)))
    order = {name: i for i, name in enumerate(ASSAY_NAMES)}
    return sorted(rows, key=lambda r: (order.get(r["assay"].split(", ")[0], len(order)), r["key"]))


def sample_files_view(sample, data, *, width=None, files=None, local=None):
    """A sample's BAMs and FASTQ folders as two tables.

    files (the sample's files) and local (URLs with a local copy) save looking
    them up again when showing many samples.
    """
    local = data.local_urls() if local is None else local
    bams = _bam_rows(sample, data, local)
    lines = []
    if bams:
        downloaded = sum(r["local"] == "yes" for r in bams)
        lines.append(f"BAMs, aligned reads ({len(bams)}"
                     + (f", {downloaded} downloaded" if downloaded else "") + "):")
        lines.append(table(bams, ("assay", "platform", "provider", "size", "local", "key"),
                           width=width, fixed=("key",)))
    else:
        lines.append("BAMs: none")
    if sample.missing_bams:
        lines.append("The site names BAMs the bucket doesn't have: " + "; ".join(sample.missing_bams))
    folders = _fastq_folders(sample, sample.files if files is None else files)
    lines.append("")
    if folders:
        lines.append(f"FASTQ folders, raw reads ({len(folders)}):")
        lines.append(table(folders, ("library", "provider", "files", "size", "folder"), width=width,
                           fixed=("folder",)))
    else:
        lines.append("FASTQ folders: none")
    return "\n".join(lines), bams, folders


def sample_view(sample, data=None, *, width=None, python=False):
    """One sample: where it came from, its files, and commands that fetch them
    (Python calls instead of shell commands with python=True)."""
    what = ", ".join(x for x in (sample.description or sample.tissue, sample.site, sample.date) if x)
    lines = [f"{sample.id}: {what}" + (f" (time point {sample.timepoint})" if sample.timepoint else ""),
             f"Sequencing: {sequencing_text(sample.sequencing) or 'none recorded'}"]
    if sample.providers:
        lines.append(f"Providers: {', '.join(sample.providers)}")
    if sample.notes:
        lines.append(f"Note from the site: {sample.notes}")
    for item in sample.disagreements:
        lines.append(f"The {item['source']} source gives {item['field']} {item['value']}.")
    if sample.corrections and data is not None:
        lines += corrected_lines(data, sample.corrections)
    if data is None:
        lines += ["", f"BAMs: {len(sample.bams)}, FASTQ folders: {len(sample.fastq_folders)}"]
        return "\n".join(lines)
    files, bams, folders = sample_files_view(sample, data, width=width)
    lines += ["", files, "", "Get the data:", get_data_hints(data, bams, folders, python=python), ""]
    if python:
        lines.append(f'All of its files: data.samples["{sample.id}"].files')
        if sample.date:
            lines.append(f'What happened around then: data.timeline.around("{sample.date}").listing()')
    else:
        lines.append(f'In Python: data.samples["{sample.id}"].files')
        if sample.date:
            lines.append(f"What happened around then: osteosarc on {sample.date}")
    return "\n".join(lines)


def get_data_hints(data, bams, folders, *, python=False):
    """Commands, written out, that fetch these BAMs and FASTQ folders."""
    lines = []
    if bams:
        bam = next((b for b in bams if b["indexed"]), bams[0])
        variant = example_variant(data)
        if variant and bam["indexed"]:
            lines.append((f'data.extract_reads("{bam["key"]}", variants=data.variants(ids="{variant}"), padding=100)'
                          if python else f"osteosarc reads {bam['key']} --variant {variant} --padding 100",
                          "reads around a variant, streamed into a small local BAM"))
        lines.append((f'data.download("{bam["key"]}", to=".")' if python else
                      f"osteosarc download {bam['key']} --to .",
                      f"the whole BAM ({bam['size'] or 'size unknown'})"
                      + (" and its index" if bam["indexed"] else "") + ", into this folder"))
    if folders:
        folder = folders[0]["folder"]
        lines.append((f'data.files.select(prefix="{folder}")' if python else f"osteosarc files --prefix {folder}",
                      "the files in a FASTQ folder"))
        if data._download_header.get("download_base", BUCKET) == BUCKET:  # not for a mirror
            lines.append((f"aws s3 cp --recursive --no-sign-request {S3_BUCKET}{folder} "
                          f"{PurePosixPath(folder).name}/", "a whole folder, with the AWS CLI (in a shell)"))
    if not lines:
        return "  (no files)"
    return "\n".join(f"  {command}\n      {what}" for command, what in lines)


def example_variant(data):
    """A catalogue variant with a ready allele for example commands (MAP2's, when present)."""
    try:
        ready = data.variants(status="ready")
    except OsteosarcError:
        return None
    return next((v.id for v in ready if v.gene == "MAP2"), ready[0].id if len(ready) else None)


def corrected_lines(data, ids):
    """Each correction's summary, wrapped, under its ID."""
    summaries = {r["id"]: r["summary"] for r in data.corrections}
    return [textwrap.fill(f"Corrected by {i}: {summaries.get(i, '')}", 100, subsequent_indent="  ")
            for i in ids]


def _first_sentence(text):
    end = text.find(". ")
    return text if end < 0 else text[:end + 1]


def all_sample_files_view(samples, data, *, width=None):
    by_sample, local = defaultdict(list), data.local_urls()
    for file in data.files:  # one pass, rather than one per sample
        for sample in file.samples:
            by_sample[sample].append(file)
    parts = []
    for sample in samples:
        what = ", ".join(x for x in (sample.description or sample.tissue, sample.site, sample.date) if x)
        files, _, _ = sample_files_view(sample, data, width=width, files=by_sample[sample.id], local=local)
        parts.append(f"== {sample.id}: {what}\n{files}")
    parts.append("Get a file: osteosarc download KEY --to .   Reads in a region: osteosarc reads KEY REGION"
                 "\nOne sample, with commands written out: osteosarc samples SAMPLE")
    return "\n\n".join(parts)


# ---------------------------------------------------------------------------
# Files
# ---------------------------------------------------------------------------

def files_overview(data, *, width=None):
    files = [f for f in data.files if not f.key.startswith("site/")]
    kinds, sizes, formats = Counter(), Counter(), defaultdict(Counter)
    tops, top_sizes = Counter(), Counter()
    for file in files:
        kinds[file.kind] += 1
        sizes[file.kind] += file.size or 0
        formats[file.kind][file.format or "(none)"] += 1
        top = file.key.split("/")[0] + ("/" if "/" in file.key else "")
        tops[top] += 1
        top_sizes[top] += file.size or 0
    by_kind = [dict(kind=kind, files=f"{kinds[kind]:,}", size=size_text(sizes[kind]),
                    formats=", ".join(f"{fmt} {n:,}" for fmt, n in formats[kind].most_common(3)),
                    what=KIND_NAMES.get(kind, "")) for kind in KINDS if kinds[kind]]
    folders = [dict(folder=top, files=f"{tops[top]:,}", size=size_text(top_sizes[top]))
               for top, _ in top_sizes.most_common(12)]
    tables = len(data.files) - len(files)
    return "\n".join([
        f"{len(files):,} files in {S3_BUCKET} ({size_text(sum(sizes.values()))}), "
        f"plus {tables} of the site's own tables",
        "",
        table(by_kind, ("kind", "files", "size", "what", "formats"), width=width),
        "",
        "Largest top-level folders:",
        table(folders, ("folder", "files", "size"), width=width),
        "",
        "Go further:",
        hints([("osteosarc files --kind alignment", "every BAM"),
               ("osteosarc files --sample T1_tumor", "one sample's BAMs and FASTQs (see osteosarc samples)"),
               ("osteosarc files --prefix hudson_lab/", "one folder"),
               ("osteosarc files --kind table --prefix site/", "the site's own tables"),
               ("osteosarc files --downloaded", "files already on this computer")])])


def files_view(files, data, *, limit=50, width=None, more="Use --limit N to show more, or --json for records."):
    local = data.local_urls()
    order = {kind: i for i, kind in enumerate(KINDS)}
    files = sorted(files, key=lambda f: (order.get(f.kind, len(order)), f.key))
    rows = [dict(kind=f.kind, sample=_text(f.samples), assay=_text(f.values("assay")),
                 provider=_text(f.values("provider")), size=size_text(f.size),
                 local="yes" if f.url in local else "", notes=_text(f.metadata.get("corrections")), key=f.key)
            for f in (files if limit is None else files[:limit])]
    columns = [c for c in ("kind", "sample", "assay", "provider") if any(r[c] for r in rows) or c == "kind"]
    if len({f.kind for f in files}) == 1:
        columns.remove("kind")
    columns += ["size", "local", *(("notes",) if any(r["notes"] for r in rows) else ()), "key"]
    total = sum(f.size or 0 for f in files)
    shown = table(rows, columns, width=width, limit=limit, total=len(files), fixed=("key",))
    return (f"{plural(len(files), 'file')} ({size_text(total) or '0 B'})\n{shown}"
            + (f"\n{more}" if limit is not None and len(files) > limit else ""))


def downloads_view(rows, cache_root, *, width=None):
    files = [r for r in rows if r["kind"] == "file"]
    reads = [r for r in rows if r["kind"] == "reads"]
    out = []
    if files:
        out.append(f"Downloaded files ({len(files)}, {size_text(sum(r['size'] or 0 for r in files))}):")
        out.append(table([dict(r, size=size_text(r["size"])) for r in files], ("key", "size", "path"),
                         width=width, fixed=("key", "path")))
        out.append("Cached files are named by content. Put one in a folder under its own name with\n"
                   "osteosarc download KEY --to DIR (a read-only hard link, so no second copy).")
    else:
        out.append("Downloaded files: none yet (osteosarc download FILE)")
    out.append("")
    if reads:
        out.append(f"Extracted reads ({len(reads)}):")
        out.append(table([dict(r, size=size_text(r["size"])) for r in reads], ("key", "regions", "size", "path"),
                         width=width, fixed=("path", "regions")))
        out.append("The same request reuses its extract; osteosarc reads prints the path.")
    else:
        out.append("Extracted reads: none yet (osteosarc reads FILE REGION)")
    out.append(f"\nCache: {cache_root} (set OSTEOSARC_CACHE to use another).")
    return "\n".join(out)


# ---------------------------------------------------------------------------
# Variants, vaccines and corrections
# ---------------------------------------------------------------------------

def variants_view(variants, *, width=None, limit=None):
    if not len(variants):
        return "(no matching variants)"
    rows = [dict(id=v.id, status=v.status, allele=allele_text(v),
                 protein=v.annotations.get("protein_change") or "", vaccines=_text(v.vaccines),
                 found_by=_text(v.pipelines), corrections=_text(v.annotations.get("corrections")))
            for v in variants]
    ready = sum(v.status == "ready" for v in variants)
    columns = ("id", *(("status",) if ready < len(variants) else ()), "allele", "protein", "vaccines",
               "found_by", *(("corrections",) if any(r["corrections"] for r in rows) else ()))
    return (f"{len(variants)} variants, {ready} with a ready allele\n"
            + table(rows, columns, width=width, limit=limit, fixed=("id", "allele")))


def variant_view(variant, data, *, width=None):
    lines = [f"{variant.id}: {variant.gene}, {variant.status}"]
    if len(variant.alleles) == 1:
        chrom, pos, ref, alt = variant.alleles[0]
        lines.append(f"Allele: {chrom}:{pos} {ref}>{alt} ({variant.assembly})")
    elif variant.alleles:
        lines.append(f"Candidate alleles ({variant.assembly}): "
                     + "; ".join(f"{c}:{p} {r}>{a}" for c, p, r, a in variant.alleles))
    effect = ", ".join(x for x in (variant.annotations.get("protein_change"), variant.annotations.get("consequence"),
                                   variant.annotations.get("variant_type")) if x)
    if effect:
        lines.append(f"Effect: {effect}")
    lines.append(f"Vaccines: {_text(variant.vaccines) or 'none'}")
    lines.append(f"Found by: {_text(variant.pipelines) or 'no pipeline recorded'}")
    if variant.annotations.get("corrections"):
        lines += corrected_lines(data, variant.annotations["corrections"])
    counts = [r for r in data.vafs if r.get("variant_id") == variant.id]
    if counts:
        rows = sorted(({"sample": r["sample_label"], "alt/total": f"{r['alt_reads']}/{r['total_reads']}",
                        "VAF": r["vaf"], "BAM": r["bam_file"]} for r in counts
                       if r.get("total_reads") not in (None, "")), key=lambda r: r["sample"])
        lines += ["", f"Read counts published by the site ({len(rows)}):",
                  table(rows, ("sample", "alt/total", "VAF", "BAM"), width=width)]
    if variant.status == "ready":
        lines += ["", "Reads around it: osteosarc reads BAM_KEY --variant "
                      f"{variant.id} --padding 100   (BAM keys: osteosarc files --kind alignment)"]
    return "\n".join(lines)


def vaccines_view(data, *, width=None):
    names = data.vaccine_names
    rows = []
    for row in data.vaccines:
        included = [n for n in names if (row.get("vaccines") or {}).get(n)]
        rows.append(dict(gene=row.get("gene"), mutation=row.get("mutation"), vaccines=", ".join(included),
                         ELISPOT=(row.get("elispot_status") or "").replace("_", " ")))
    return (f"{len(rows)} vaccine targets across {len(names)} vaccines ({', '.join(names)})\n"
            + table(rows, ("gene", "mutation", "vaccines", "ELISPOT"), width=width))


def corrections_view(data, *, width=None):
    rows = [dict(r, summary=_first_sentence(r["summary"])) for r in data.corrections]
    statuses = Counter(r["status"] for r in rows)
    return (f"{len(rows)} corrections: " + ", ".join(f"{n} {s}" for s, n in statuses.most_common()) + "\n"
            + table(rows, ("id", "status", "action", "summary"), width=width, fixed=("id",))
            + "\n\nOne correction with its evidence: osteosarc corrections ID")


def correction_view(row):
    lines = [f"{row['id']} ({row['status']}, {row['action']})", "", textwrap.fill(row["summary"], 88), ""]
    lines.append("Records it checks or changes:")
    for change in row["changes"]:
        lines.append(f"  {change['source']}: {change['match']} ({change['records']} records, {change['state']})")
    lines += ["", "Evidence:", *(f"  - {item}" for item in row["evidence"])]
    return "\n".join(lines)


def summary_view(data):
    lines = [snapshot_line(data)]
    try:
        samples = data.samples
        lines.append(f"samples: {len(samples)} (" + ", ".join(
            f"{n} {t}" for t, n in Counter(s.tissue or "unknown" for s in samples).most_common()) + ")")
        lines.append(f"timeline: {data.timeline.overview()}")
    except OsteosarcError as error:
        lines.append(f"samples and timeline: unavailable ({error})")
    kinds = Counter(f.kind for f in data.files)
    lines.append(f"files: {len(data.files):,} (" + ", ".join(
        f"{kinds[k]:,} {k}" for k in KINDS if kinds[k]) + ")")
    variants = data.variants()
    lines.append(f"variants: {len(variants)} on the site, "
                 + ", ".join(f"{n} {s}" for s, n in Counter(v.status for v in variants).most_common()))
    statuses = Counter(r["status"] for r in data.corrections)
    lines.append("corrections: " + ", ".join(f"{n} {s}" for s, n in statuses.most_common()))
    return "\n".join(lines)
