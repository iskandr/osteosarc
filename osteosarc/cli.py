# PYTHON_ARGCOMPLETE_OK
"""Command-line access using exactly the same Dataset API as Python callers."""

import argparse
import json
import os
import subprocess
import sys
from dataclasses import asdict

from datacache.cli_help import HelpArgumentParser, HelpFormatter, enable_argcomplete

from ._version import __version__

_DEFAULT_DERIVED_UNUSED_DAYS = 30

#: Every command, grouped as `osteosarc` and `osteosarc --help` show them.
COMMANDS = (
    ("Browse", (
        ("samples [SAMPLE]", "Tumor, organoid and blood samples; one sample's files and how to get them"),
        ("files", "Files in the S3 bucket: an overview, or a list by --kind, --sample or --prefix"),
        ("variants [ID]", "The variant catalogue, with alleles, vaccines and read counts"),
        ("vaccines", "Vaccine targets and ELISPOT results"),
        ("timeline", "Treatments, procedures, scans and MRD over time (--around DATE for one week)"),
        ("corrections [ID]", "Known problems in the website's data, and the fixes applied"))),
    ("Get data", (
        ("reads FILE|SAMPLE", "The reads around variants or in regions, as a small BAM: test data in seconds"),
        ("download FILE", "A whole file (--to DIR puts it, and its index, in DIR)"),
        ("downloads", "What's already on this computer, and where"),
        ("cache ...", "Inspect cache directories or prune derived/ read subsets and leftovers"))),
    ("Test data for libraries", (
        ("test-data ...", "Bundles of test reads: make your own, or use the shared openvax-v2"),)),
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
    root = HelpArgumentParser(
        prog="osteosarc",
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
                                   epilog=epilog, formatter_class=HelpFormatter)

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
    reads.add_argument("--assembly", help="Assembly of explicit regions, such as GRCh38 (default: a BAM's own)")
    reads.add_argument("--reference", help="Local indexed FASTA (required for CRAM)")
    reads.add_argument("--index", help="Explicit index path/URL when the bucket has none")
    reads.add_argument("--reference-length", type=int, help="Expected contig length (mitochondrial queries)")
    reads.add_argument("--min-mapq", type=int, default=0, help="Only reads with at least this mapping quality")
    reads.add_argument("--exclude-flags", type=lambda s: int(s, 0), default=0,
                       help="Leave out reads with any of these SAM flags, such as 0x400 (duplicates)")
    mates = reads.add_mutually_exclusive_group()
    mates.add_argument("--fetch-pairs", action="store_true",
                       help="Also retrieve mates outside the regions, except those with no position")
    mates.add_argument("--unplaced-mates", action="store_true",
                       help="Also retrieve mates outside the regions, even those with no position (slower)")
    mates.add_argument("--recover-linked", action="store_true", help="Bounded mate and SA-linked recovery")
    reads.add_argument("--json", action="store_true", help="The extract's paths and receipt as JSON")
    downloads = command("downloads")
    downloads.add_argument("--json", action="store_true")
    cache_cmd = command("cache", data=False, epilog="examples:\n"
                        "  osteosarc cache dirs\n"
                        "  osteosarc cache prune --days 30\n"
                        "  osteosarc cache prune --max-gb 2 --yes")
    cache_actions = cache_cmd.add_subparsers(dest="cache_command", metavar="ACTION")
    cache_dirs = cache_actions.add_parser(
        "dirs", help="Show cache directories, disk usage, and offline state (read-only)")
    cache_dirs.add_argument("--json", action="store_true", help="Machine-readable JSON")
    cache_prune = cache_actions.add_parser(
        "prune", help="Prune unused derived/ read subsets and interrupted staging leftovers (dry-run by default)")
    cache_prune.add_argument(
        "--days", type=float, default=None,
        help=f"Remove derived entries older than this many days (default: {_DEFAULT_DERIVED_UNUSED_DAYS} unless --max-gb/--all)")
    cache_prune.add_argument(
        "--max-gb", type=float, default=None,
        help="Cap remaining derived/ cache size at this many GB (evicts oldest first)")
    cache_prune.add_argument(
        "--all", action="store_true",
        help="Prune all derived/ read subsets and leftovers regardless of age")
    cache_prune.add_argument(
        "--yes", action="store_true",
        help="Actually delete candidates (default is dry-run)")
    cache_prune.add_argument("--json", action="store_true", help="Machine-readable JSON")
    sync = command("sync", data=False)
    sync.add_argument("name", nargs="?", help="Optional name (default: today's UTC date)")
    sync.add_argument("--refresh", action="store_true", help="Download again even if today's snapshot exists")
    sync.add_argument("--source-revision", help="Pin every public site-repository source to this commit")
    sync.add_argument("--json", action="store_true")
    snapshots = command("snapshots", data=False)
    snapshots.add_argument("--json", action="store_true")
    command("repl")
    test_data = command("test-data", data=False, epilog="examples:\n"
                        "  osteosarc test-data make tests/data/dync1h1 T0_tumor --assay rna-seq "
                        "--variant DYNC1H1-chr14-101980529\n"
                        "  osteosarc test-data list openvax-v2\n"
                        "  osteosarc test-data export openvax-v2 tests/data --member "
                        "IPISRC044_tumor_T2_ucla.redux.DYNC1H1-chr14-101980529\n"
                        "  osteosarc test-data check openvax-v2 fixtures.json")
    actions = test_data.add_subparsers(dest="test_data_command", metavar="ACTION")
    make = actions.add_parser(
        "make", parents=[snapshot], formatter_class=HelpFormatter,
        help="Make a bundle of test reads for variants and SVs in some BAMs (or from a recipe)",
        description="Make a bundle of test reads in a new folder. For each variant, each BAM gives a\n"
                    "balanced set of templates (up to 20 alt, 10 ref, 5 other and 2 not spanning it); for\n"
                    "each SV or fusion, up to 50 templates that join its breakends. Every record is pinned\n"
                    "by checksum, so anyone can rebuild the bundle and check it offline.",
        epilog="examples:\n"
               "  osteosarc test-data make tests/data/dync1h1 T0_tumor --assay rna-seq "
               "--variant DYNC1H1-chr14-101980529\n"
               "  osteosarc test-data make tests/data/gabbr1 T2_tumor --assay wgs --sv GABBR1-SLC29A1\n"
               "  osteosarc test-data make bundle --recipe recipe.json")
    make.add_argument("output", help="A new folder for the bundle")
    make.add_argument("sources", nargs="*", metavar="FILE|SAMPLE",
                      help="BAMs to read: a file's key, URL or ID, or a sample ID for each of its indexed BAMs")
    make.add_argument("--variant", action="append", default=[], metavar="ID",
                      help="A catalogue variant with a ready allele; repeat for several")
    make.add_argument("--sv", action="append", default=[], metavar="ID",
                      help="An SV candidate (such as SV0461) or SV regression target; repeat for several")
    make.add_argument("--assay", help="With a sample: only its BAMs of this assay, such as rna-seq")
    make.add_argument("--platform", help="With a sample: only its BAMs from this platform, such as ont")
    make.add_argument("--recipe", metavar="FILE",
                      help="Make the bundle from a recipe instead, such as another bundle's recipe.json. Its "
                           "BAMs are read through --snapshot, else the recipe's own snapshot if you have it, else "
                           "your newest, which must list each BAM unchanged; or give them with --source")
    make.add_argument("--source", action="append", default=[], metavar="ID=LOCAL_BAM",
                      help="With --recipe: use a BAM you already have for one of its sources")
    make.add_argument("--from", dest="from_bundle", metavar="BUNDLE",
                      help="Subset an existing bundle, keeping its original headers and acquisition receipts")
    make.add_argument("--member", action="append", metavar="PREFIX",
                      help="With --from: keep members whose names start with this prefix; repeat for several")
    make.add_argument("--header-policy", choices=("full", "compact"), default="full",
                      help="compact keeps only the header lines the records need")
    make.add_argument("--unplaced-mates", action="store_true",
                      help="Also keep mates with no position (slower: reads each BAM's unplaced reads)")
    make.add_argument("--size-budget", type=int, default=64 * 1024 * 1024, help="Largest bundle, in bytes")
    make.add_argument("--json", action="store_true", help="The bundle's manifest as JSON")
    for name, what in (("list", "Each member of a bundle, with its records and why"),
                       ("verify", "Check a bundle's files, records and indexes"),
                       ("export", "Write members into a folder as indexed BAMs (or SAM), named after them")):
        action = actions.add_parser(name, help=what)
        action.add_argument("bundle", help="A published bundle such as openvax-v2, or a bundle folder's path")
        if name == "export":
            action.add_argument("output", help="Folder to write the members into; it may already exist")
            action.add_argument("--member", action="append", help="Only these members; repeat for several")
            action.add_argument("--format", choices=("bam", "sam", "sam.gz"), default="bam")
        if name == "verify":
            action.add_argument("--sha256", help="The manifest checksum you expect")
        if name in ("list", "verify"):
            action.add_argument("--json", action="store_true",
                                help="Every member's records and reasons as JSON" if name == "list"
                                else "The bundle's manifest as JSON")
    check = actions.add_parser("check", help="Compare your library's test files with a bundle's members")
    check.add_argument("bundle", help="A published bundle such as openvax-v2, or a bundle folder's path")
    check.add_argument("fixtures", help='JSON mapping member names to your files: a BAM, SAM or SAM.gz path, '
                                        'or {"json": path, "pointer": "/path/to/lines"}; paths are relative '
                                        "to this JSON file's folder")
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
        raise ValueError("Supply either regions or --variant, not both; variants carry their own assembly")
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
    if args.regions and args.assembly is None and is_sample(dataset, args.file):
        raise ValueError("Give --assembly with regions and a sample: its BAMs may use different genome builds")
    # Without --assembly, a BAM's regions are in its own genome build (Region's None).
    from .models import Region
    return dict(regions=[Region.from_samtools(r, assembly=args.assembly, reference_length=args.reference_length)
                         for r in args.regions], padding=args.padding)


def snapshot_selector(args):
    """Dataset.open's arguments for --snapshot: a download date, or a name or ID prefix."""
    from .dataset import DATE_SELECTOR
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
    from .dataset import Dataset
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
    from .dataset import Dataset
    from .errors import NoSnapshotsError
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

    from .dataset import Dataset
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
            print_json(events.to_records()) if args.json else print(events.listing(args.width) or "(no events)")
            return 0
        if args.json:
            print_json(events.to_records())
        elif args.list:
            print(events.listing(args.width) or "(no events)")
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


def bams(dataset, name, *, assay=None, platform=None):
    """A FILE|SAMPLE argument's BAMs: that file, or a sample's indexed BAMs (of this assay and platform)."""
    from .shared import alignment_files
    if not is_sample(dataset, name):
        if assay or platform:
            raise ValueError("--assay and --platform choose among a sample's BAMs; give a sample ID")
        return [dataset.file(name)]
    return alignment_files(dataset, name, assay=assay, platform=platform,
                           warn=lambda text: print(text, file=sys.stderr))


def read_sources(dataset, args):
    """The BAMs a reads command reads from: one file, or a sample's BAMs."""
    if is_sample(dataset, args.file) and args.index:
        raise ValueError("--index belongs to one BAM; give a file's key, not a sample, to use it")
    return bams(dataset, args.file, assay=args.assay, platform=args.platform)


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
        from concurrent.futures import ThreadPoolExecutor, as_completed

        from .reads import ReadFilter
        sources = read_sources(dataset, args)
        targets = read_targets(dataset, args)

        def extract(file):
            return dataset.extract_reads(file, **targets, reference=args.reference, index=args.index,
                                         filters=ReadFilter(args.min_mapq, args.exclude_flags),
                                         fetch_pairs=args.fetch_pairs or args.unplaced_mates,
                                         unplaced_mates=args.unplaced_mates,
                                         recovery={} if args.recover_linked else None, to=args.to)
        results, finished, shown = [], {}, 0
        # A sample's BAMs are read four at a time and reported in order; the first error
        # stops those not yet begun.
        pool = ThreadPoolExecutor(max_workers=max(1, min(4, len(sources))))
        try:
            order = {pool.submit(extract, file): i for i, file in enumerate(sources)}
            for extracted in as_completed(order):
                try:
                    finished[order[extracted]] = extracted.result()
                except CoordinateError as error:
                    if len(sources) == 1:
                        raise
                    finished[order[extracted]] = error
                while shown in finished:
                    file, subset = sources[shown], finished.pop(shown)
                    shown += 1
                    if isinstance(subset, CoordinateError):
                        print(f"skipping {file.key}: {subset}", file=sys.stderr)  # e.g. a GRCh37 BAM
                        continue
                    results.append(dict(file=file.key, path=str(subset.path), index=str(subset.index_path),
                                        receipt=subset.receipt))
                    if not args.json:
                        print(subset.path, flush=True)
        except BaseException:
            for i in sorted(finished):  # finished after one that failed: still say where they are
                if not isinstance(finished[i], CoordinateError) and not args.json:
                    print(finished[i].path, flush=True)
            raise
        finally:
            # BAMs not yet begun are cancelled; any being read finish (and are cached)
            # before the process exits.
            pool.shutdown(wait=False, cancel_futures=True)
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
    enable_argcomplete(root)
    arguments = sys.argv[1:] if argv is None else argv
    if not arguments:
        from .cache import Cache
        print(start(Cache(offline=True)))
        return 0
    args, extra = root.parse_known_args(arguments)
    if args.command == "test-data" and args.test_data_command is None:
        root.parse_args([*arguments, "--help"])  # its actions, and examples (exits)
    if args.command == "cache" and args.cache_command is None:
        root.parse_args([*arguments, "--help"])
    # Before Python 3.13, argparse binds an optional list positional (reads' regions,
    # test-data make's BAMs) before any option, so items written after an option
    # arrive here as extra arguments.
    trailing = ("regions" if args.command == "reads" else
                "sources" if args.command == "test-data" and args.test_data_command == "make" else None)
    if extra and (trailing is None or any(item.startswith("-") for item in extra)):
        root.error(f"unrecognized arguments: {' '.join(extra)}")
    if extra:
        setattr(args, trailing, [*getattr(args, trailing), *extra])
    from .cache import Cache, is_offline_env
    from .errors import OsteosarcError
    offline_enabled = bool(args.offline or is_offline_env())
    cache = Cache(args.cache, offline=offline_enabled)
    try:
        if args.command == "cache":
            return cache_action(args, cache)
        if args.command == "test-data":
            return test_data(args, cache)
        if args.command == "sync":
            from .dataset import Dataset
            sources = pinned_sources(args.source_revision) if args.source_revision else None
            dataset = Dataset.sync(args.name, cache=cache, refresh=args.refresh, sources=sources,
                                   corrections=not args.no_corrections)
            if args.json:
                print_json(dict(snapshot=dataset.name, id=dataset.id, downloaded=dataset.downloaded))
            else:
                print(f"Saved snapshot {dataset.name} ({dataset.id[:12]}). Try `osteosarc samples`.")
            return 0
        if args.command == "snapshots":
            from .dataset import Dataset
            rows = list(Dataset.snapshots(cache=cache))
            print_json(rows) if args.json else print(snapshots_view(rows))
            return 0
        # Browsing stays offline; commands that fetch bytes may use the network unless --offline.
        online = args.command in ("download", "reads", "repl") and not offline_enabled
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


def cache_action(args, cache):
    """cache dirs and cache prune."""
    from .cache import (
        DERIVED_UNUSED_DAYS,
        _entry_size_bytes,
        offline_env_source,
        prune_derived,
    )
    from .views import size_text

    if args.cache_command == "dirs":
        derived = cache.workspace / "derived"
        snapshots_dir = cache.workspace / "snapshots"
        derived_bytes = _entry_size_bytes(derived)
        objects_bytes = _entry_size_bytes(cache.objects)
        workspace_bytes = _entry_size_bytes(cache.workspace)
        offline_src = offline_env_source()
        payload = {
            "package": "osteosarc",
            "version": __version__,
            "root": str(cache.root),
            "objects_dir": str(cache.objects),
            "workspace_dir": str(cache.workspace),
            "derived_dir": str(derived),
            "snapshots_dir": str(snapshots_dir),
            "exists": cache.root.exists(),
            "objects_bytes": objects_bytes,
            "workspace_bytes": workspace_bytes,
            "derived_bytes": derived_bytes,
            "offline": cache.offline,
            "offline_source": offline_src or ("--offline" if args.offline else None),
        }
        if args.json:
            print_json(payload)
            return 0
        print(f"osteosarc cache locations (v{__version__})")
        print(f"  root      : {cache.root}")
        print(f"  objects   : {cache.objects} [{size_text(objects_bytes) if cache.objects.exists() else 'not created'}]")
        print(f"  workspace : {cache.workspace} [{size_text(workspace_bytes) if cache.workspace.exists() else 'not created'}]")
        print(f"  derived   : {derived} [{size_text(derived_bytes) if derived.exists() else 'not created'}]")
        if cache.offline:
            detail = f" ({payload['offline_source']})" if payload["offline_source"] else ""
            print(f"  offline   : enabled{detail}")
        else:
            print("  offline   : disabled")
        return 0

    if args.all:
        days = 0
        max_bytes = None
        leftover_seconds = 0
    else:
        if args.max_gb is not None and args.max_gb < 0:
            raise ValueError("--max-gb must be nonnegative")
        max_bytes = None if args.max_gb is None else int(args.max_gb * 1_000_000_000)
        days = DERIVED_UNUSED_DAYS if (args.days is None and max_bytes is None) else args.days
        leftover_seconds = 24 * 3600

    report = prune_derived(
        cache.workspace,
        days=days,
        max_bytes=max_bytes,
        leftover_seconds=leftover_seconds,
        dry_run=not args.yes,
    )
    if args.json:
        print_json(report)
        return 0

    print(f"osteosarc derived cache: {report['derived_dir']}")
    if not report["candidates"]:
        print(f"Nothing to prune ({report['kept_count']} kept, {size_text(report['kept_bytes'])}).")
        return 0
    for item in report["candidates"]:
        tag = "LEFTOVER" if item["kind"] == "leftover" else "DERIVED "
        print(f"  DELETE  {tag}  {size_text(item['size_bytes']):>10}  {item['name']}")
    if report["dry_run"]:
        print(
            f"\nDry run — would remove {len(report['candidates'])} item(s) "
            f"({size_text(report['reclaimable_bytes'])}), keeping {report['kept_count']}. "
            "Re-run with --yes to delete."
        )
    else:
        print(
            f"\nDeleted {len(report['deleted'])} item(s), freed {size_text(report['freed_bytes'])}; "
            f"{report['kept_count']} kept."
        )
    return 0


def test_data(args, cache):
    """test-data make, list, verify, export and check."""
    from pathlib import Path

    from .bundles import (
        export_bundle,
        generate_bundle,
        list_bundle,
        sources_to_read,
        subset_bundle,
        verify_bundle,
    )
    from .fixtures import validate_recipe
    from .shared import bundle_folder, check_fixtures, read_json
    from .views import plural, table
    action = args.test_data_command
    if action == "make":
        if args.member and not args.from_bundle:
            raise ValueError("--member goes with --from")
        if args.from_bundle:
            if (args.recipe or args.source or args.sources or args.variant or args.sv or args.assay
                    or args.platform or args.unplaced_mates or args.header_policy != "full"):
                raise ValueError("--from preserves the existing bundle's selections and header policy; "
                                 "use only --member and --size-budget")
            if not args.member:
                raise ValueError("--from requires at least one --member prefix")
            root = bundle_folder(args.from_bundle, cache=cache)
            available = list_bundle(root)
            for prefix in args.member:
                if not any(name.startswith(prefix) for name in available):
                    raise ValueError(f"No bundle member starts with {prefix!r}")
            names = [name for name in available if any(name.startswith(p) for p in args.member)]
            manifest = subset_bundle(root, names, args.output, size_budget=args.size_budget)
        elif args.recipe:
            if args.sources or args.variant or args.sv or args.assay or args.platform or args.unplaced_mates:
                raise ValueError("--recipe makes the bundle from the recipe alone, which says how each BAM is read; "
                                 "drop the BAMs, --variant, --sv and --unplaced-mates")
            recipe = validate_recipe(read_json(args.recipe))
            given = dict(item.split("=", 1) for item in args.source)
            dataset = None
            if set(sources_to_read(recipe)) - set(given):
                from .dataset import Dataset
                # Its BAMs are found through a snapshot: --snapshot, else the recipe's own if
                # it's here, else the newest. Any will do in which each BAM is unchanged.
                named = recipe["snapshot"] if isinstance(recipe.get("snapshot"), dict) else {}
                own = [row["id"] for row in Dataset.snapshots(cache=cache) if row["id"] == named.get("id")
                       or (named.get("id") is None and row["name"] == named.get("name"))]
                if args.snapshot is None and own:
                    args.snapshot = own[0]  # by ID: a name like 2026-09-25 would be read as a date
                dataset = open_snapshot(args, cache, not args.offline)
            manifest = generate_bundle(recipe, args.output, sources=given, cache=cache, dataset=dataset,
                                       size_budget=args.size_budget, header_policy=args.header_policy)
        else:
            if args.source:
                raise ValueError("--source goes with --recipe")
            if not args.sources:
                raise ValueError("Name the BAMs to read: a file's key, or a sample ID (see osteosarc samples)")
            dataset = open_snapshot(args, cache, not args.offline)
            files = [f for name in args.sources for f in bams(dataset, name, assay=args.assay, platform=args.platform)]
            print(f"Reading {plural(len(files), 'BAM')}; the first time, this streams the reads it needs "
                  "from each.", file=sys.stderr)
            from .shared import make_bundle
            manifest = make_bundle(dataset, args.output, variants=args.variant, svs=args.sv, files=files,
                                   unplaced_mates=args.unplaced_mates,
                                   size_budget=args.size_budget, header_policy=args.header_policy,
                                   log=lambda text: print(text, file=sys.stderr))
        if args.json:
            print_json(manifest)
            return 0
        rows = [dict(member=name, status=m["status"], records=m["record_count"])
                for name, m in manifest["members"].items()]
        print(bundle_summary(args.output, manifest) + "\n\n" + table(rows, ["status", "records", "member"],
                                                                      fixed=("member",)))
        print(f'\nIn a test: osteosarc.bundle_file("{args.output}", MEMBER) gives a member as an indexed BAM.')
        return 0
    args.bundle = bundle_folder(args.bundle, cache=cache)
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
        manifest = verify_bundle(args.bundle, sha256=args.sha256)
        if args.json:
            print_json(manifest)
        else:
            print(bundle_summary(args.bundle, manifest) + "\nEvery file, record and index checks out.")
    elif action == "list":
        members = list_bundle(args.bundle)
        if args.json:
            print_json(members)
        else:
            rows = [dict(member=name, status=m["status"], records=m["record_count"]) for name, m in members.items()]
            print(table(rows, ["status", "records", "member"], fixed=("member",)))
    else:
        for path in export_bundle(args.bundle, args.output, members=args.member, format=args.format).values():
            print(path)
    return 0


def bundle_summary(where, manifest):
    from .views import plural, size_text
    members = manifest["members"].values()
    with_reads = sum(1 for m in members if m["record_count"])
    records = sum(s["record_count"] for s in manifest["sources"].values())
    return (f"{where}: {plural(len(members), 'member')} ({with_reads} with reads), "
            f"{plural(records, 'record')} from {plural(len(manifest['sources']), 'BAM')}, "
            f"{size_text(manifest['total_size_bytes'])}")
