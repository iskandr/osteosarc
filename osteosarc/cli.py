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
        ("timeline", "Treatments, procedures, scans and MRD over time (--around DATE for one week)"),
        ("corrections [ID]", "Known problems in the website's data, and the fixes applied"))),
    ("Get data", (
        ("reads FILE|SAMPLE", "The reads around variants or in regions, as a small BAM: test data in seconds"),
        ("download FILE", "A whole file (--to DIR puts it, and its index, in DIR)"),
        ("downloads", "What's already on this computer, and where"))),
    ("Test data for libraries", (
        ("test-data ...", "Shared test reads (openvax-v1): list, export, and check your copies"),)),
    ("Snapshots of the website's metadata", (
        ("sync", "Download the current metadata (about 57 MB); commands use the newest"),
        ("snapshots", "Saved snapshots, by download date"),
        ("repl", "Python with the newest snapshot loaded as `data`"))),
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
    described = {name.split()[0].split("|")[0]: what for _, group in COMMANDS for name, what in group}

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
                                          "  osteosarc timeline --around 2025-01-28\n"
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
    timeline.add_argument("--around", metavar="DATE", help="List everything within a week of this date")
    timeline.add_argument("--days", type=int, default=7, help="Days either side of --around (default 7)")
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
                    "  osteosarc reads T1_tumor --assay rna-seq --variant MAP2-chr2-209694768 --padding 100 --to tests/data\n"
                    "  osteosarc reads KEY chr17:7661779-7687538 --assembly GRCh38\n\n"
                    "Streams only those reads (never the whole BAM) and prints the path of each small BAM it\n"
                    "writes; asking again reuses it.")
    reads.add_argument("file", metavar="FILE|SAMPLE",
                       help="An indexed BAM: its key (see osteosarc files), URL or ID; or a sample ID, for each of "
                            "its BAMs")
    reads.add_argument("--assay", help="With a sample: only its BAMs of this assay, such as rna-seq")
    reads.add_argument("--platform", help="With a sample: only its BAMs from this platform, such as ont")
    reads.add_argument("--to", metavar="DIR", help="Also put each BAM and its index in DIR, named for the source "
                                                   "file and what was asked for")
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
    sync = command("sync", data=False)
    sync.add_argument("name", nargs="?", help="Optional name (default: today's UTC date)")
    sync.add_argument("--refresh", action="store_true", help="Download again even if today's snapshot exists")
    sync.add_argument("--source-revision", help="Pin every public site-repository source to this commit")
    sync.add_argument("--json", action="store_true")
    snapshots = command("snapshots", data=False)
    snapshots.add_argument("--json", action="store_true")
    command("repl")
    test_data = command("test-data", data=False, epilog="examples:\n"
                        "  osteosarc test-data list openvax-v1\n"
                        "  osteosarc test-data export openvax-v1 tests/data --member DYNC1H1-rna\n"
                        "  osteosarc test-data check openvax-v1 fixtures.json\n"
                        "  osteosarc test-data generate recipe.json bundle")
    actions = test_data.add_subparsers(dest="test_data_command", required=True, metavar="ACTION")
    generate = actions.add_parser("generate", help="Fetch a recipe's reads and write a bundle")
    generate.add_argument("recipe")
    generate.add_argument("output", help="A new directory")
    generate.add_argument("--source", action="append", default=[], metavar="ID=LOCAL_BAM",
                          help="Use a BAM you already have for one of the recipe's sources")
    generate.add_argument("--size-budget", type=int, default=64 * 1024 * 1024)
    generate.add_argument("--header-policy", choices=("full", "compact"), default="full",
                          help="compact keeps only the header lines the records need")
    for name, what in (("list", "Each member of a bundle, with its records and why"),
                       ("verify", "Check a bundle's files, records and indexes"),
                       ("export", "Write a bundle's members as indexed BAMs (or SAM)")):
        action = actions.add_parser(name, help=what)
        action.add_argument("bundle", help="A bundle directory, or a published bundle such as openvax-v1")
        if name == "export":
            action.add_argument("output", help="A new directory")
            action.add_argument("--member", action="append", help="Only these members; repeat for several")
            action.add_argument("--format", choices=("bam", "sam", "sam.gz"), default="bam")
        if name == "verify":
            action.add_argument("--sha256", help="The manifest checksum you expect")
    check = actions.add_parser("check", help="Compare your library's test files with a bundle's members")
    check.add_argument("bundle", help="A bundle directory, or a published bundle such as openvax-v1")
    check.add_argument("fixtures", help='JSON mapping member names to your files: a BAM, SAM or SAM.gz path, '
                                        'or {"json": path, "pointer": "/path/to/lines"}')
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
        found = dataset.variants("all", ids=args.variant)
        missing = sorted(set(args.variant) - {v.id for v in found})
        if missing:
            raise ValueError(f"Unknown variant ID(s): {', '.join(missing)}; "
                             "list them with `osteosarc variants --set all"
                             + (f" --snapshot {args.snapshot}" if args.snapshot else "") + "`")
        order = {v: i for i, v in enumerate(dict.fromkeys(args.variant))}
        variants = type(found)(sorted(found, key=lambda v: order[v.id]), source=found.source)  # as given
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
        if args.days != 7 and not args.around:
            raise ValueError("--days sets the window around --around DATE; give a date too")
        if args.around:
            events = events.around(args.around, args.days)
            print_json(events.to_records()) if args.json else print(events.listing() or "(no events)")
            return 0
        if args.json:
            print_json(events.to_records())
        elif args.list:
            print(events.listing() or "(no events)")
        else:
            # Asking for a lane or text shows it even when it's one the chart leaves out.
            print(events.render(width=args.width, since=args.since, until=args.until,
                                everything=args.all or bool(args.lane or args.contains)))
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


def is_sample(dataset, name):
    """Whether a reads command's FILE|SAMPLE names a sample (file keys and URLs have slashes)."""
    return "/" not in name and name in {sample.id for sample in dataset.samples}


def read_sources(dataset, args):
    """The BAMs a reads command reads from: one file, or a sample's BAMs."""
    if not is_sample(dataset, args.file):
        if args.assay or args.platform:
            raise ValueError("--assay and --platform choose among a sample's BAMs; give a sample ID")
        return [dataset.file(args.file)]
    if args.index:
        raise ValueError("--index belongs to one BAM; give a file's key, not a sample, to use it")
    bams = dataset.samples[args.file].files.select(kind="alignment", assay=args.assay, platform=args.platform)
    indexed = [f for f in bams if f.index_urls]
    for file in bams:
        if not file.index_urls:
            print(f"skipping {file.key}: it has no index in the bucket", file=sys.stderr)
    if not indexed:
        raise ValueError(f"{args.file} has no indexed BAMs" + (" of that kind" if args.assay or args.platform else ""))
    return indexed


def get_data(args, dataset):
    """reads, download and downloads. Returns the exit code."""
    from . import views
    from .errors import CoordinateError
    if args.command == "download":
        file = dataset.file(args.file)
        print(dataset.download(file, to=args.to, refresh=args.refresh))
        if args.to and file.index_urls:
            print(dataset.download(file.index_urls[0], to=args.to))  # already placed; this names it
    elif args.command == "reads":
        targets, sources = read_targets(dataset, args), read_sources(dataset, args)
        results = []
        for file in sources:
            try:
                subset = dataset.extract_reads(file, **targets, reference=args.reference, index=args.index,
                                               filters=ReadFilter(args.min_mapq, args.exclude_flags),
                                               fetch_pairs=args.fetch_pairs,
                                               recovery={} if args.recover_linked else None, to=args.to)
            except CoordinateError as error:
                if len(sources) == 1:
                    raise
                print(f"skipping {file.key}: {error}", file=sys.stderr)  # e.g. a GRCh37 BAM
                continue
            results.append(dict(file=file.key, path=str(subset.path), index=str(subset.index_path),
                                receipt=subset.receipt))
            if not args.json:
                print(subset.path)
        if not results:
            raise ValueError(f"No reads extracted: every BAM of {args.file} was skipped")
        if args.json:
            # A sample always gives a list, however many of its BAMs it has.
            print_json(results if is_sample(dataset, args.file) else results[0])
    else:
        rows = list(dataset.downloads())
        print_json(rows) if args.json else print(views.downloads_view(rows, dataset.cache.root))
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
        if args.command == "test-data":
            return test_data(args, cache)
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
        # Browsing stays offline; commands that fetch bytes may use the network unless --offline.
        online = args.command in ("download", "reads", "repl") and not args.offline
        dataset = open_snapshot(args, cache, online)
        if args.command == "repl":
            repl(dataset)
            return 0
        if args.command in ("download", "reads", "downloads"):
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


def test_data(args, cache):
    """test-data generate, list, verify, export and check."""
    from pathlib import Path

    from .bundles import export_bundle, generate_bundle, list_bundle, verify_bundle
    from .shared import check_fixtures, fetch_bundle, read_json
    action = args.test_data_command
    if action == "generate":
        recipe = read_json(args.recipe)
        sources = dict(item.split("=", 1) for item in args.source)
        print_json(generate_bundle(recipe, args.output, sources=sources, cache=cache,
                                   size_budget=args.size_budget, header_policy=args.header_policy))
        return 0
    if not Path(args.bundle).exists():
        args.bundle = fetch_bundle(args.bundle, cache=cache)
    if action == "check":
        manifest = Path(args.fixtures)
        problems = check_fixtures(args.bundle, read_json(manifest), root=manifest.parent)
        if problems:
            print_json(problems)
            print(f"osteosarc: {len(problems)} of your files differ from the bundle", file=sys.stderr)
            return 1
        print(f"Every file in {manifest} matches the bundle.")
        return 0
    if action == "verify":
        value = verify_bundle(args.bundle, sha256=args.sha256)
    elif action == "list":
        value = list_bundle(args.bundle)
    else:
        value = export_bundle(args.bundle, args.output, members=args.member, format=args.format)
    print_json(value)
    return 0
