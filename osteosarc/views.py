"""Plain-text views of a Dataset, shared by the CLI and Python reprs.

Every view is a function returning a string, so the command line and a
Python session show the same thing.
"""

from __future__ import annotations

import shlex
import shutil
import textwrap
from collections import Counter, defaultdict
from pathlib import PurePosixPath

from .catalog import fastq_prefix, matches_fastq
from .curation import ASSAY_NAMES, ASSAYS, normalize_provider
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


def table(rows, columns, *, width=None, limit=None, wrap=False, fixed=(), optional=(), drop_empty=False,
          total=None):
    """Fixed-width text table; long cells are cut, with an ellipsis, to fit the terminal.

    Columns in fixed, such as file keys people copy, are never shortened, and the
    line under the headers stops at the terminal's edge. With drop_empty, columns
    empty in every row shown are left out. Columns in optional are left out, the
    last first, while the table is too wide even with every column at its
    narrowest, with a line saying so. total is the full row count when rows holds
    only the first few.
    """
    if limit is not None and (not isinstance(limit, int) or isinstance(limit, bool) or limit < 0):
        raise ValueError("limit must be a nonnegative integer")
    rows = list(rows)
    shown = rows if limit is None else rows[:limit]
    width = width or shutil.get_terminal_size((120, 24)).columns
    text = {c: [_text(row.get(c)) for row in shown] for c in columns}
    natural = {c: max([len(c)] + [len(t) for t in text[c]]) for c in columns}
    columns = [c for i, c in enumerate(columns) if not drop_empty or i == 0 or any(text[c])]

    def narrowest(columns):
        return sum(natural[c] if i == 0 or c in fixed else min(12, natural[c]) for i, c in enumerate(columns)) \
            + 2 * (len(columns) - 1)
    dropped = []
    for c in reversed(optional):
        if c in columns and narrowest(columns) > width:
            columns = [other for other in columns if other != c]
            dropped.insert(0, c)
    widths = [natural[c] for c in columns]
    # Shrink the widest columns to fit, but never the first (identifying) column.
    shrinkable = [i for i in range(1, len(widths)) if columns[i] not in fixed]
    while sum(widths) + 2 * (len(widths) - 1) > width and max((widths[i] for i in shrinkable), default=0) > 12:
        widest = max(shrinkable, key=widths.__getitem__)
        widths[widest] -= 1
    rule, used = [], 0
    for w in widths:  # under each column, as far as the terminal's edge
        if used >= width:
            break
        rule.append("-" * min(w, width - used))
        used += w + 2
    line = "  ".join
    out = [line(_cut(c, w).ljust(w) for c, w in zip(columns, widths)).rstrip(), line(rule)]
    for i in range(len(shown)):
        cells = [text[c][i] for c in columns]
        parts = [textwrap.wrap(cell, width=w, break_on_hyphens=False) or [""] if wrap else [_cut(cell, w)]
                 for cell, w in zip(cells, widths)]
        for j in range(max(map(len, parts), default=0)):
            out.append(line((part[j] if j < len(part) else "").ljust(w) for part, w in zip(parts, widths)).rstrip())
    total = len(rows) if total is None else total
    if limit is not None and total > limit:
        out.append(f"... {total - limit:,} more")
    if dropped:
        out.append(f"(Too narrow for {' and '.join(dropped)}: widen the terminal, or see --json.)")
    return "\n".join(out)


def prose(text, width=None):
    """A paragraph wrapped to the terminal (at most 100 columns), later lines indented."""
    width = min(width or shutil.get_terminal_size((100, 24)).columns, 100)
    return textwrap.fill(text, width, subsequent_indent="  ", break_on_hyphens=False, break_long_words=False)


def _cut(text, width):
    return text if len(text) <= width else text[:width - 1] + "…"


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


def plural(n, word, words=None):
    """"1 BAM", "3 BAMs"; words is the plural when it isn't word + "s"."""
    return f"{n:,} {word if n == 1 else words or word + 's'}"


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
                 **{"FASTQ locations": len(s.fastq_folders)}) for s in samples]
    shown = table(rows, ("sample", "timepoint", "date", "where", "sequencing", "BAMs", "FASTQ locations"),
                  width=width, wrap=True)
    footer = footer if footer is not None else hints([
        ('data.samples["T1_tumor"]', "one sample's files, with how to get them"),
        ('.select(tissue="blood", assay="cite-seq")', "filter by timepoint, tissue, assay or platform")])
    return f"{title}\n\n{shown}\n\n{footer}"


#: What kind of single-cell library a site label names; bulk libraries need no name.
LIBRARIES = {"scRNA_GEX": "gene expression", "scRNA": "gene expression", "Tumor scRNA": "gene expression",
             "Blood scRNA": "gene expression", "scRNA_TCR": "αβ TCR", "scRNA_TCRgd": "γδ TCR", "scRNA_BCR": "BCR",
             "CITE": "antibody tags", "scRNA_ONT": "long reads", "scRNA ONT": "long reads", "PacBio": "long reads"}


def _fastq_folders(sample, files):
    """One row per published FASTQ location, counting only matching raw reads."""
    labels = {}
    for row in sample.details.get("fastqs", ()):
        labels.setdefault(row["folder"], (row.get("assay") or "", row.get("provider") or ""))
    reads = [f for f in files or () if f.kind == "reads"]
    rows = []
    for folder in sample.fastq_folders:
        label, provider = labels.get(folder, ("", ""))
        assay, platform = ASSAYS.get(label, (label, None))
        matched = [f for f in reads if matches_fastq(f.key, folder)]
        prefix = fastq_prefix(folder)
        directory = not matched or all(f.key.startswith(prefix + "/") for f in matched)
        rows.append(dict(assay=assay, platform=platform or "", library=LIBRARIES.get(label, ""),
                         provider=normalize_provider(provider) if provider else "",
                         files=len(matched) if files is not None else "",
                         size=size_text(sum(f.size or 0 for f in matched)),
                         folder=prefix + ("/" if directory else ""),
                         published_location=folder, directory=directory,
                         keys=tuple(f.key for f in matched)))
    order = {name: i for i, name in enumerate(ASSAY_NAMES)}
    kinds = {name: i for i, name in enumerate(dict.fromkeys(LIBRARIES.values()))}
    return sorted(rows, key=lambda r: (order.get(r["assay"], len(order)), r["assay"],
                                       kinds.get(r["library"], -1), r["folder"]))


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


def sample_files_view(sample, data, *, width=None, local=None):
    """A sample's BAMs and FASTQ folders as two tables.

    local (URLs with a local copy) saves looking them up again when showing many
    samples.
    """
    local = data.local_urls() if local is None else local
    bams = _bam_rows(sample, data, local)
    lines = []
    if bams:
        downloaded = sum(r["local"] == "yes" for r in bams)
        lines.append(f"BAMs, aligned reads ({len(bams)}"
                     + (f", {downloaded} downloaded" if downloaded else "") + "):")
        lines.append(table(bams, ("assay", "platform", "provider", "size", "local", "key"),
                           width=width, fixed=("key",), drop_empty=True))
    else:
        lines.append("BAMs: none")
    if sample.missing_bams:
        lines.append(prose("The site names BAMs the bucket doesn't have: " + "; ".join(sample.missing_bams), width))
    folders = _fastq_folders(sample, data.files.select(sample=sample.id))
    lines.append("")
    if folders:
        lines.append(f"FASTQ locations, raw reads ({len(folders)}; folders or filename prefixes):")
        lines.append(table(folders, ("assay", "platform", "library", "provider", "files", "size", "folder"),
                           width=width, fixed=("folder", "library"), drop_empty=True))
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
        lines.append(prose(f"Note from the site: {sample.notes}", width))
    for item in sample.disagreements:
        lines.append(prose(f"The {item['source']} source gives {item['field']} {item['value']}.", width))
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
        location = folders[0]
        folder = location["folder"]
        lines.append((f'data.files.select(prefix="{folder}")' if python else f"osteosarc files --prefix {folder}",
                      "files under this FASTQ location"))
        if data._download_header.get("download_base", BUCKET) == BUCKET:  # not for a mirror
            if location.get("directory", folder.endswith("/")):
                lines.append((f"aws s3 cp --recursive --no-sign-request {shlex.quote(S3_BUCKET + folder)} "
                              f"{shlex.quote(PurePosixPath(folder).name + '/')}",
                              "a whole folder, with the AWS CLI (in a shell)"))
            elif location["keys"]:
                parent = str(PurePosixPath(folder).parent) + "/"
                includes = " ".join("--include " + shlex.quote(key[len(parent):]) for key in location["keys"])
                lines.append((f"aws s3 cp --recursive --no-sign-request {shlex.quote(S3_BUCKET + parent)} "
                              f"./ --exclude '*' {includes}",
                              "only this location's FASTQs, with the AWS CLI (in a shell)"))
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
    summaries = {c.id: c.summary for c in data.curation.corrections}  # no need to evaluate them
    return [prose(f"Corrected by {i}: {summaries.get(i, '')}")
            for i in ids]


def _first_sentence(text):
    end = text.find(". ")
    return text if end < 0 else text[:end + 1]


def all_sample_files_view(samples, data, *, width=None):
    local = data.local_urls()
    parts = []
    for sample in samples:
        what = ", ".join(x for x in (sample.description or sample.tissue, sample.site, sample.date) if x)
        files, _, _ = sample_files_view(sample, data, width=width, local=local)
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
                    # The two commonest, whole: a count cut short would mislead.
                    formats=", ".join([f"{fmt} {n:,}" for fmt, n in formats[kind].most_common(2)]
                                      + (["…"] if len(formats[kind]) > 2 else [])),
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
    columns = ["kind", "sample", "assay", "provider", "size", "local", "notes", "key"]
    if len({f.kind for f in files}) == 1:
        columns.remove("kind")
    total = sum(f.size or 0 for f in files)
    shown = table(rows, columns, width=width, limit=limit, total=len(files), fixed=("key",), drop_empty=True)
    return (f"{plural(len(files), 'file')} ({size_text(total) or '0 B'})\n{shown}"
            + (f"\n{more}" if limit is not None and len(files) > limit else ""))


def downloads_view(rows, cache_root):
    """Downloaded files and extracted reads: each one's size and key (and regions),
    and on the next line where it is, whole, to copy."""
    files = [r for r in rows if r["kind"] == "file"]
    reads = [r for r in rows if r["kind"] == "reads"]

    def listing(rows, *labels):
        sizes = [size_text(r["size"]) for r in rows]
        pad = max(map(len, sizes), default=0)
        return "\n".join(f"{size.rjust(pad)}  {'  '.join(r[label] for label in labels)}\n{' ' * (pad + 2)}{r['path']}"
                         for size, r in zip(sizes, rows))
    out = []
    if files:
        out.append(f"Downloaded files ({len(files)}, {size_text(sum(r['size'] or 0 for r in files))}):")
        out.append(listing(files, "key"))
        out.append("Cached files are named by content. Put one in a folder under its own name with\n"
                   "osteosarc download KEY --to DIR (a read-only hard link, so no second copy).")
    else:
        out.append("Downloaded files: none yet (osteosarc download FILE)")
    out.append("")
    if reads:
        out.append(f"Extracted reads ({len(reads)}):")
        out.append(listing(reads, "key", "regions"))
        out.append("The same request reuses its extract; osteosarc reads prints its path.")
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
               "found_by", "corrections")
    return (f"{len(variants)} variants, {ready} with a ready allele\n"
            + table(rows, columns, width=width, limit=limit, fixed=("id", "allele"), drop_empty=True,
                    optional=("vaccines", "found_by", "corrections")))


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
        lines += ["", f"Reads around it: osteosarc reads BAM_KEY --variant {variant.id} --padding 100",
                  "  (BAM keys: osteosarc files --kind alignment, or a sample ID for all its BAMs)"]
    return "\n".join(lines)


def vaccines_view(data, *, width=None):
    vaccines = list(data.vaccines)
    names = data._json("vaccine_overlap").get("vaccine_names") or list(dict.fromkeys(
        name for row in vaccines for name in (row.get("vaccines") or {})))  # the site's own list and order
    rows = []
    for row in vaccines:
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
        lines.append(f"  {change['source']}: {change['match']} ({plural(change['records'], 'record')}, {change['state']})")
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
