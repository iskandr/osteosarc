"""Command-line access using exactly the same Dataset API as Python callers."""

import argparse
import json
import os
import subprocess
import sys
from dataclasses import asdict

from . import __version__
from .cache import Cache
from .dataset import DATE_SELECTOR, Dataset
from .discovery import list_bucket
from .errors import NoSnapshotsError, OsteosarcError
from .models import Region
from .reads import ReadFilter

#: Every command, grouped as `osteosarc` and `osteosarc --help` show them.
COMMANDS = (
    ("Browse", (
        ("samples [SAMPLE]", "Tumor, organoid and blood samples; one sample's files and how to get them"),
        ("files", "Files in the S3 bucket: an overview, or a list with --kind, --sample or --prefix"),
        ("variants [ID]", "The variant catalogue, with alleles, vaccines and read counts"),
        ("vaccines", "Vaccine targets and ELISPOT results"),
        ("timeline", "Treatments, procedures, scans and MRD over time"),
        ("on DATE", "Everything within a week of a date"),
        ("corrections [ID]", "Known problems in the website's data, and the fixes applied"))),
    ("Get data", (
        ("download FILE", "Download a whole file (--to DIR puts it, and its index, in DIR)"),
        ("reads FILE REGION", "Stream the reads in a region, or around a variant, into a small local BAM"),
        ("downloads", "What's already on this computer, and where"),
        ("table NAME", "Print one of the site's tables, or a bucket CSV/TSV, as TSV"))),
    ("Snapshots of the website's metadata", (
        ("sync", "Download the current metadata (about 57 MB); commands use the newest"),
        ("snapshots", "Saved snapshots, by download date"))),
    ("More", (
        ("repl", "Python with the newest snapshot loaded as `data`"),
        ("discover PREFIX", "List a bucket folder live, without a snapshot"),
        ("fixtures ...", "Build and check read fixtures for tests"))),
)


def command_guide():
    width = max(len(name) for _, commands in COMMANDS for name, _ in commands)
    lines = []
    for group, commands in COMMANDS:
        lines.append(f"{group}:")
        lines += [f"  {name:<{width}}  {what}" for name, what in commands]
        lines.append("")
    lines.append("Run osteosarc COMMAND --help for its options. Docs: https://iskandr.github.io/osteosarc/")
    return "\n".join(lines)


def parser():
    root = argparse.ArgumentParser(
        prog="osteosarc", formatter_class=argparse.RawDescriptionHelpFormatter,
        description="Browse and download the public osteosarc.com dataset.\n\n" + command_guide())
    root.add_argument("--version", action="version", version=f"osteosarc {__version__}")
    root.add_argument("--cache", help="Cache root (default: OSTEOSARC_CACHE, else the shared OpenVax "
                                      "cache: OPENVAX_DATA_CACHE or the platform 'openvax' cache)")
    root.add_argument("--offline", action="store_true", help="Forbid any network access")
    root.add_argument("--no-corrections", action="store_true",
                      help="Use the website's data unchanged (see osteosarc corrections)")
    commands = root.add_subparsers(dest="command", required=True, metavar="COMMAND")
    # Commands that read a snapshot take --snapshot; the newest is the default.
    snapshot = argparse.ArgumentParser(add_help=False)
    snapshot.add_argument("--snapshot", help="A UTC download date, month or year for the newest snapshot "
                                             "downloaded then (2026-09-24, 2026-09), or a snapshot name or ID "
                                             "prefix (default: the most recent; see osteosarc snapshots)")
    described = {name.split()[0]: what for _, group in COMMANDS for name, what in group}

    def command(name, *, data=True, epilog=None):
        return commands.add_parser(name, parents=[snapshot] if data else [], description=described[name],
                                   epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)

    samples = command("samples", epilog="examples:\n  osteosarc samples\n  osteosarc samples T1_tumor\n"
                                        "  osteosarc samples --tissue blood --assay cite-seq\n"
                                        "  osteosarc samples --files")
    samples.add_argument("sample", nargs="?", help="One sample, such as T1_tumor: its files, and commands to get them")
    samples.add_argument("--timepoint", help="T0, T1, T2 or T3")
    samples.add_argument("--tissue", help="tumor, blood or organoid")
    samples.add_argument("--assay", help="rna-seq, wes, wgs, scrna-seq or cite-seq")
    samples.add_argument("--platform", help="illumina, ont or pacbio")
    samples.add_argument("--files", action="store_true", help="List every sample's BAMs and FASTQ folders")
    samples.add_argument("--json", action="store_true", help="Sample records as JSON")

    files = command("files", epilog="With no filters, a summary of the bucket by kind and folder.\n\n"
                                    "examples:\n  osteosarc files --kind alignment\n"
                                    "  osteosarc files --sample T1_tumor\n  osteosarc files --prefix hudson_lab/\n"
                                    "  osteosarc files --downloaded")
    files.add_argument("--sample", help="A sample's BAMs and FASTQ files, such as T0_tumor (see osteosarc samples)")
    files.add_argument("--kind", help="alignment, reads, variants, expression, annotation, table, reference, "
                                      "index, image or other")
    files.add_argument("--format", help="File format, such as bam, fastq or vcf")
    files.add_argument("--prefix", help="Keys starting with this, such as vendor/tempus/")
    files.add_argument("--contains", help="Keys containing this text")
    files.add_argument("--timepoint", help="T0, T1, T2 or T3, as the sources state it")
    files.add_argument("--assay", help="rna-seq, wes, wgs, scrna-seq or cite-seq")
    files.add_argument("--platform", help="illumina, ont or pacbio")
    files.add_argument("--tissue", help="tumor, blood or organoid")
    files.add_argument("--provider", help="Who made the file, such as BostonGene or Tempus")
    files.add_argument("--library", help="A library ID, such as BG009368")
    files.add_argument("--include-conflicts", action="store_true",
                       help="Also match files whose sources disagree about a filtered field")
    files.add_argument("--include-inferred", action="store_true",
                       help="Also match fields inferred from file paths")
    files.add_argument("--downloaded", action="store_true", help="Only files already on this computer")
    files.add_argument("--limit", type=int, default=50, help="Rows to show (default 50)")
    files.add_argument("--json", action="store_true", help="File records as JSON (up to --limit)")

    variants = command("variants", epilog="examples:\n  osteosarc variants --gene MAP2\n"
                                          "  osteosarc variants MAP2-chr2-209694768\n"
                                          "  osteosarc variants --vaccine mRNA --status ready")
    variants.add_argument("id", nargs="?", help="One variant, with its read counts")
    variants.add_argument("--set", choices=("site", "all", "vaccine"), default="site",
                          help="site: the variants page (default); all: also count-only entries; "
                               "vaccine: vaccine targets")
    variants.add_argument("--gene")
    variants.add_argument("--vaccine", help="In this vaccine, such as mRNA or 'JLF V2'")
    variants.add_argument("--pipeline", help="Found by this pipeline, such as 'Mutect2 2024'")
    variants.add_argument("--status", help="ready (has one allele), or another status")
    variants.add_argument("--vaccine-source", choices=("overlap", "source_variants"), default="overlap",
                          help="Which source says a variant is in a vaccine (default: overlap)")
    variants.add_argument("--json", action="store_true", help="Variant records as JSON")

    vaccines = command("vaccines")
    vaccines.add_argument("--json", action="store_true", help="Vaccine-overlap records as JSON")

    timeline = command("timeline", epilog="examples:\n  osteosarc timeline\n"
                                          "  osteosarc timeline --since 2024-05 --until 2024-09\n"
                                          "  osteosarc timeline --lane MRD --list")
    timeline.add_argument("--since", help="YYYY, YYYY-MM or YYYY-MM-DD")
    timeline.add_argument("--until", help="YYYY, YYYY-MM or YYYY-MM-DD")
    timeline.add_argument("--lane", help="Only lanes whose name contains this text, such as MRD or Lab")
    timeline.add_argument("--contains", help="Only events whose name contains this text")
    timeline.add_argument("--all", action="store_true", help="Also chart lab draws, DICOM studies, cytometry, "
                                                             "flow draws, slides and sequencing runs")
    timeline.add_argument("--width", type=int, help="Chart width (default: the terminal's)")
    timeline.add_argument("--list", action="store_true", help="One line per event instead of a chart")
    timeline.add_argument("--json", action="store_true", help="Event records as JSON")
    around = command("on")
    around.add_argument("date", help="YYYY-MM-DD")
    around.add_argument("--days", type=int, default=7, help="Days either side (default 7)")
    corrections = command("corrections")
    corrections.add_argument("id", nargs="?", help="One correction, with its evidence")
    corrections.add_argument("--strict", action="store_true",
                             help="Exit 1 if a correction is stale or the data uses unknown labels")
    corrections.add_argument("--json", action="store_true", help="Correction reports as JSON")

    download = command("download", epilog="examples:\n  osteosarc download "
                       "rna-seq/reprocessed/BG003082/BG003082.Aligned.sortedByCoord.out.md.bam --to .\n"
                       "  osteosarc download vafs")
    download.add_argument("file", help="A file's key (see osteosarc files), URL or ID")
    download.add_argument("--to", metavar="DIR", help="Also put the file, and its index, in DIR under its own name")
    download.add_argument("--refresh", action="store_true",
                          help="Download again, for a file this snapshot hasn't downloaded yet")
    reads = command("reads", epilog="examples:\n"
                    "  osteosarc reads KEY --variant MAP2-chr2-209694768 --padding 100\n"
                    "  osteosarc reads KEY chr17:7661779-7687538 --assembly GRCh38\n\n"
                    "Prints the path of the BAM it writes; the same request reuses it.")
    reads.add_argument("file", help="An indexed BAM/CRAM: its key (see osteosarc files), URL or ID")
    reads.add_argument("regions", nargs="*", help="contig:start-end, one-based inclusive (or use --variant)")
    reads.add_argument("--variant", action="append", default=[], metavar="ID",
                       help="A catalogue variant with a ready allele; repeat for several")
    reads.add_argument("--padding", type=int, default=0, help="Bases added on each side of each --variant")
    reads.add_argument("--assembly", help="Assembly of explicit regions, such as GRCh38")
    reads.add_argument("--reference", help="Local indexed FASTA (required for CRAM)")
    reads.add_argument("--index", help="Explicit index path/URL when the bucket has none")
    reads.add_argument("--reference-length", type=int, help="Expected contig length (mitochondrial queries)")
    reads.add_argument("--min-mapq", type=int, default=0)
    reads.add_argument("--exclude-flags", type=lambda s: int(s, 0), default=0)
    reads.add_argument("--fetch-pairs", action="store_true", help="Also retrieve paired mates outside the regions")
    reads.add_argument("--recover-linked", action="store_true", help="Bounded mate and SA-linked recovery")
    reads.add_argument("--json", action="store_true", help="The extract's paths and receipt as JSON")
    downloads = command("downloads")
    downloads.add_argument("--json", action="store_true")
    table = command("table", epilog="examples:\n  osteosarc table vafs\n  osteosarc table dna_fusions --json")
    table.add_argument("file", help="A site table's name (see osteosarc files --prefix site/) or a CSV/TSV key")
    table.add_argument("--json", action="store_true", help="Rows as JSON")

    sync = command("sync", data=False)
    sync.add_argument("name", nargs="?", help="Optional name (default: today's UTC date)")
    sync.add_argument("--refresh", action="store_true", help="Download again even if today's snapshot exists")
    sync.add_argument("--source-revision", help="Pin every public site-repository source to this commit")
    sync.add_argument("--json", action="store_true")
    snapshots = command("snapshots", data=False)
    snapshots.add_argument("--json", action="store_true")
    command("repl")
    discover = command("discover", data=False)
    discover.add_argument("prefix")
    discover.add_argument("--refresh", action="store_true")
    discover.add_argument("--json", action="store_true")
    fixtures = command("fixtures", data=False)
    actions = fixtures.add_subparsers(dest="fixture_command", required=True)
    select = actions.add_parser("select", help="Return record membership and inclusion reasons")
    select.add_argument("recipe")
    select.add_argument("--source", action="append", default=[], metavar="ID=LOCAL_BAM")
    panel = actions.add_parser("panel", help="Print a shipped named target panel")
    panel.add_argument("name")
    for name in ("generate", "pack"):
        action = actions.add_parser(name, help="Select and publish a portable bundle")
        action.add_argument("recipe")
        action.add_argument("output")
        action.add_argument("--source", action="append", default=[], metavar="ID=LOCAL_BAM")
        action.add_argument("--header-policy", choices=("full", "compact"), default="full")
        action.add_argument("--size-budget", type=int, default=64 * 1024 * 1024)
    for name in ("verify", "list", "export"):
        action = actions.add_parser(name, help="Work with a bundle entirely offline")
        action.add_argument("bundle")
        if name == "export":
            action.add_argument("output")
            action.add_argument("--member", action="append")
            action.add_argument("--format", choices=("bam", "sam", "sam.gz"), default="bam")
        if name == "verify":
            action.add_argument("--sha256", help="Pinned manifest digest")
    return root


def pinned_sources(revision):
    """Every site-repository source URL, rewritten from main to one commit."""
    from .catalog import SNAPSHOT_SOURCES, SOURCE_REPO, TIMELINE_SOURCES
    if any(c not in "0123456789abcdefABCDEF" for c in revision) or len(revision) != 40:
        raise ValueError("--source-revision must be a full 40-character commit SHA")
    return {key: url.replace(SOURCE_REPO, SOURCE_REPO.replace("/main/", f"/{revision}/"))
            for key, url in {**SNAPSHOT_SOURCES, **TIMELINE_SOURCES}.items() if url.startswith(SOURCE_REPO)}


def read_targets(dataset, args):
    """Parse the reads command's regions or catalogue variants; Dataset.extract_reads checks the rest."""
    if args.variant and (args.regions or args.assembly or args.reference_length):
        raise ValueError("Supply either regions (with --assembly) or --variant, not both; "
                         "variants carry their own assembly")
    if args.variant:
        variants = dataset.variants("all", ids=args.variant)
        missing = sorted(set(args.variant) - {v.id for v in variants})
        if missing:
            raise ValueError(f"Unknown variant ID(s): {', '.join(missing)}; "
                             "list them with `osteosarc variants --set all"
                             + (f" --snapshot {args.snapshot}" if args.snapshot else "") + "`")
        return dict(variants=variants, padding=args.padding)
    if args.regions and args.assembly is None:
        raise ValueError("--assembly is required with explicit regions")
    return dict(regions=[Region.from_samtools(r, assembly=args.assembly, reference_length=args.reference_length)
                         for r in args.regions], padding=args.padding)


def snapshot_selector(args):
    """Dataset.open's arguments for --snapshot: a download date, or a name or ID prefix."""
    if args.snapshot is not None and DATE_SELECTOR.fullmatch(args.snapshot):
        return dict(date=args.snapshot)
    return dict(name=args.snapshot)


def snapshots_view(rows):
    if not rows:
        return "No saved snapshots; run `osteosarc sync`."
    from .views import table
    shown = [{"name": r["name"], "downloaded (UTC)": (r["downloaded"] or "")[:16].replace("T", " "),
              "id": r["id"][:12]} for r in rows]
    return (table(shown, ["name", "downloaded (UTC)", "id"])
            + f"\n\nCommands use {rows[0]['name']} unless you pass --snapshot.")


def start(cache):
    """What `osteosarc` alone prints: the snapshot in use, and every command."""
    rows = list(Dataset.snapshots(cache=cache))
    if rows:
        when = (rows[0]["downloaded"] or "")[:16].replace("T", " ")
        status = f"Using snapshot {rows[0]['name']}" + (f", downloaded {when} UTC." if when else ".")
    else:
        status = "No data yet: start with `osteosarc sync`, which downloads the website's metadata (about 57 MB)."
    return (f"osteosarc {__version__}: browse and download the public osteosarc.com dataset.\n{status}\n\n"
            + command_guide())


def open_snapshot(args, cache, online):
    """Open the chosen snapshot; without any, say what to run, or offer to sync for the REPL."""
    try:
        return Dataset.open(**snapshot_selector(args), cache=cache, offline=not online,
                            corrections=not args.no_corrections)
    except NoSnapshotsError as error:
        where = f"No snapshot yet in {error.root}."
        if args.command == "repl" and args.snapshot is None and not args.offline and sys.stdin.isatty():
            try:
                answer = input(f"{where} Download the website's metadata now (about 57 MB)? [Y/n] ")
            except (EOFError, KeyboardInterrupt):
                print()
                answer = "n"
            if answer.strip().lower() in ("", "y", "yes"):
                dataset = Dataset.sync(cache=cache, corrections=not args.no_corrections)
                print(f"Saved snapshot {dataset.name}.")
                return dataset
        raise FileNotFoundError(f"{where} Run `osteosarc sync` to download the website's metadata "
                                "(about 57 MB).") from None


def repl(dataset):
    """Python with the snapshot as `data`: IPython when it's installed, else the standard console."""
    import osteosarc
    namespace = {"data": dataset, "osteosarc": osteosarc, "Dataset": Dataset}
    banner = f"{dataset!r}\n\n`data` is this snapshot, and help(data) documents it. Leave with Ctrl-D."
    try:
        from IPython import start_ipython
    except ImportError:
        import code
        try:
            import readline
            import rlcompleter
            readline.set_completer(rlcompleter.Completer(namespace).complete)
            readline.parse_and_bind("tab: complete")
        except ImportError:
            pass
        code.interact(banner=banner, local=namespace, exitmsg="")
        return
    print(banner)
    start_ipython(argv=["--no-banner"], user_ns=namespace)


def print_json(value):
    print(json.dumps(value, indent=2, default=lambda x: asdict(x)))


def browse(args, dataset):
    """The Browse commands: text by default, JSON with --json. Returns the exit code."""
    from . import views
    if args.command == "samples":
        if args.sample:
            if args.timepoint or args.tissue or args.assay or args.platform or args.files:
                raise ValueError("Filters and --files are for several samples; leave them out to see one")
            sample = dataset.samples[args.sample]
            print_json(asdict(sample)) if args.json else print(views.sample_view(sample, dataset))
            return 0
        samples = dataset.samples.select(timepoint=args.timepoint, tissue=args.tissue, assay=args.assay,
                                         platform=args.platform)
        if args.json:
            print_json(samples.to_records())
        elif args.files:
            print(views.all_sample_files_view(samples, dataset))
        else:
            print(views.samples_view(samples, footer=views.hints([
                ("osteosarc samples T1_tumor", "one sample's files, with commands to get them"),
                ("osteosarc samples --files", "every sample's files"),
                ("osteosarc samples --tissue blood --assay cite-seq",
                 "filter by timepoint, tissue, assay or platform")])))
        return 0
    if args.command == "files":
        if args.limit < 0:
            raise ValueError("--limit must be nonnegative")
        filters = {name: getattr(args, name) for name in
                   ("sample", "kind", "format", "prefix", "contains", "timepoint", "assay", "platform", "tissue",
                    "provider", "library", "include_conflicts", "include_inferred")}
        if args.sample is not None:
            dataset.samples[args.sample]  # an unknown sample is an error, not an empty list
        if not any(filters.values()) and not args.downloaded and not args.json:
            print(views.files_overview(dataset))
            return 0
        selection = dataset.files.select(**filters)
        if args.downloaded:
            local = dataset.local_urls()
            selection = selection.where(lambda f: f.url in local)
        if args.json:
            print_json(dict(total=len(selection), files=selection[:args.limit].to_records()))
        else:
            print(views.files_view(selection, dataset, limit=args.limit))
        return 0
    if args.command == "variants":
        if args.id:
            variant = dataset.variants("all").select(ids=args.id)
            if not variant:
                raise KeyError(f"No variant {args.id!r}; list them with `osteosarc variants --set all`")
            print_json(asdict(variant[0])) if args.json else print(views.variant_view(variant[0], dataset))
            return 0
        variants = dataset.variants(args.set, **{name: getattr(args, name) for name in
                                                 ("gene", "vaccine", "pipeline", "status", "vaccine_source")})
        if args.json:
            print_json(variants.to_records())
        else:
            print(views.variants_view(variants) + "\n\nOne variant, with its read counts: osteosarc variants ID")
        return 0
    if args.command == "vaccines":
        print_json(list(dataset.vaccines)) if args.json else print(views.vaccines_view(dataset))
        return 0
    if args.command == "timeline":
        events = dataset.timeline.select(lane=args.lane, contains=args.contains, since=args.since, until=args.until)
        if args.json:
            print_json(events.to_records())
        elif args.list:
            print(events.listing() or "(no events)")
        else:
            # Asking for a lane or text shows it even when it's one the chart leaves out.
            print(events.render(width=args.width, since=args.since, until=args.until,
                                everything=args.all or bool(args.lane or args.contains)))
        return 0
    if args.command == "on":
        print(dataset.timeline.around(args.date, args.days).listing() or "(no events)")
        return 0
    rows = list(dataset.corrections)
    if args.id:
        row = next((r for r in rows if r["id"] == args.id), None)
        if row is None:
            raise KeyError(f"No correction {args.id!r}; list them with `osteosarc corrections`")
        print_json(row) if args.json else print(views.correction_view(row))
        return 0
    unrecognized = list(dataset.unrecognized)
    if args.json:
        print_json(dict(corrections=rows, unrecognized=unrecognized))
    else:
        print(views.corrections_view(dataset))
        if unrecognized:
            print(f"\n{len(unrecognized)} labels in the data aren't in osteosarc's vocabulary: "
                  + ", ".join(f"{r['field']}={r['value']!r}" for r in unrecognized[:5])
                  + (" ..." if len(unrecognized) > 5 else ""))
    # Drift is judged on the evaluation, even when this run disables corrections.
    stale = [r["id"] for r in rows if r["evaluation"] == "stale"]
    if args.strict and (stale or unrecognized):
        print(f"osteosarc: stale corrections {stale}; {len(unrecognized)} unrecognized source labels",
              file=sys.stderr)
        return 1
    return 0


def get_data(args, dataset):
    """download, reads, downloads and table. Returns the exit code."""
    from . import views
    if args.command == "download":
        file = dataset.file(args.file)
        print(dataset.download(file, to=args.to, refresh=args.refresh))
        if args.to and file.index_urls:
            print(dataset.download(file.index_urls[0], to=args.to))  # already placed; this names it
    elif args.command == "reads":
        subset = dataset.extract_reads(args.file, **read_targets(dataset, args),
                                       reference=args.reference, index=args.index,
                                       filters=ReadFilter(args.min_mapq, args.exclude_flags),
                                       fetch_pairs=args.fetch_pairs,
                                       recovery={} if args.recover_linked else None)
        if args.json:
            print_json(dict(path=str(subset.path), index=str(subset.index_path), receipt=subset.receipt))
        else:
            print(subset.path)
    elif args.command == "downloads":
        rows = list(dataset.downloads())
        print_json(rows) if args.json else print(views.downloads_view(rows, dataset.cache.root))
    else:
        parsed = dataset.parse(args.file)
        if args.json or not hasattr(parsed, "rows"):
            print_json(list(parsed) if hasattr(parsed, "rows") else parsed)
        else:
            print("\t".join(parsed.columns))
            for row in parsed:
                print("\t".join("" if row.get(c) is None else str(row.get(c)) for c in parsed.columns))
    return 0


def main(argv=None):
    root = parser()
    arguments = sys.argv[1:] if argv is None else argv
    if not arguments:
        print(start(Cache(offline=True)))
        return 0
    args, extra = root.parse_known_args(arguments)
    # Before Python 3.13, argparse binds the optional regions positional before any
    # option, so regions written after an option arrive here as extra arguments.
    if extra and (args.command != "reads" or any(item.startswith("-") for item in extra)):
        root.error(f"unrecognized arguments: {' '.join(extra)}")
    if extra:
        args.regions = [*args.regions, *extra]
    cache = Cache(args.cache, offline=args.offline)
    try:
        if args.command == "fixtures":
            return fixtures(args, cache)
        if args.command == "sync":
            sources = pinned_sources(args.source_revision) if args.source_revision else None
            dataset = Dataset.sync(args.name, cache=cache, refresh=args.refresh, sources=sources,
                                   corrections=not args.no_corrections)
            if args.json:
                print_json(dict(snapshot=dataset.name, id=dataset.id, downloaded=dataset.downloaded))
            else:
                print(f"Saved snapshot {dataset.name} ({dataset.id[:12]}). Try `osteosarc samples`.")
            return 0
        if args.command == "snapshots":
            rows = list(Dataset.snapshots(cache=cache))
            print_json(rows) if args.json else print(snapshots_view(rows))
            return 0
        if args.command == "discover":
            listing = list_bucket(cache, args.prefix, refresh=args.refresh)
            if args.json:
                print_json(listing)
            else:
                from .views import size_text, table
                rows = [dict(size=size_text(size), modified=modified[:10], key=key)
                        for key, size, modified in listing["files"]]
                print(f"{len(rows):,} files under {args.prefix} (live, not from a snapshot)\n"
                      + table(rows, ("size", "modified", "key"), fixed=("key",)) if rows
                      else f"Nothing under {args.prefix}")
            return 0
        # Browsing stays offline; commands that fetch bytes may use the network unless --offline.
        online = args.command in ("download", "table", "reads", "repl") and not args.offline
        dataset = open_snapshot(args, cache, online)
        if args.command == "repl":
            repl(dataset)
            return 0
        if args.command in ("download", "reads", "downloads", "table"):
            return get_data(args, dataset)
        return browse(args, dataset)
    except BrokenPipeError:
        # Output piped into head and the like: stop quietly, as other command-line tools do.
        os.dup2(os.open(os.devnull, os.O_WRONLY), sys.stdout.fileno())
        return 1
    except (OsteosarcError, ValueError, KeyError, OSError, subprocess.SubprocessError) as error:
        # A KeyError's str() is its repr; show its message as written.
        message = error.args[0] if isinstance(error, KeyError) and error.args else error
        print(f"osteosarc: {message}", file=sys.stderr)
        if isinstance(error, subprocess.CalledProcessError) and error.stderr:
            print(error.stderr.decode(errors="replace") if isinstance(error.stderr, bytes) else error.stderr,
                  file=sys.stderr)
        return 1


def fixtures(args, cache):
    from .fixtures import load_panel, select_fixtures
    if args.fixture_command == "panel":
        value = load_panel(args.name)
    elif args.fixture_command in ("select", "generate", "pack"):
        from pathlib import Path

        from .bundles import generate_bundle, pack_bundle
        recipe = json.loads(Path(args.recipe).read_text())
        sources = dict(item.split("=", 1) for item in args.source)
        if args.fixture_command == "select":
            value = select_fixtures(recipe, sources).manifest
        elif args.fixture_command == "generate":
            value = generate_bundle(recipe, args.output, sources=sources, cache=cache,
                                    header_policy=args.header_policy, size_budget=args.size_budget)
        else:
            value = pack_bundle(select_fixtures(recipe, sources), args.output,
                                header_policy=args.header_policy, size_budget=args.size_budget)
    else:
        from .bundles import export_bundle, list_bundle, verify_bundle
        if args.fixture_command == "verify":
            value = verify_bundle(args.bundle, sha256=args.sha256)
        elif args.fixture_command == "list":
            value = list_bundle(args.bundle)
        else:
            value = export_bundle(args.bundle, args.output, members=args.member, format=args.format)
    print_json(value)
    return 0
