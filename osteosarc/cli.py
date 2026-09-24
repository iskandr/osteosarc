"""Command-line access using exactly the same Dataset API as Python callers."""

import argparse
import json
import subprocess
import sys
from dataclasses import asdict

from . import __version__
from .cache import Cache
from .dataset import DATE_SELECTOR, Dataset
from .discovery import list_bucket
from .errors import OsteosarcError
from .models import Region
from .reads import ReadFilter


def parser():
    root = argparse.ArgumentParser(prog="osteosarc", description=__doc__)
    root.add_argument("--version", action="version", version=f"osteosarc {__version__}")
    root.add_argument("--cache", help="Cache root (default: OSTEOSARC_CACHE, else the shared OpenVax "
                                      "cache: OPENVAX_DATA_CACHE or the platform 'openvax' cache)")
    root.add_argument("--offline", action="store_true", help="Forbid any network acquisition")
    root.add_argument("--no-corrections", action="store_true",
                      help="Use the published sources unchanged (see osteosarc.curation)")
    commands = root.add_subparsers(dest="command", required=True)
    # Commands that read a snapshot take --snapshot; the newest is the default.
    snapshot = argparse.ArgumentParser(add_help=False)
    snapshot.add_argument("--snapshot", help="A UTC download date, month or year for the newest snapshot "
                                             "downloaded then (2026-09-24, 2026-09), or a snapshot name or ID "
                                             "prefix (default: the most recent; see snapshots)")

    def data_command(name, **kwargs):
        return commands.add_parser(name, parents=[snapshot], **kwargs)

    sync = commands.add_parser("sync", help="Save the website's current metadata as a dated snapshot")
    sync.add_argument("name", nargs="?", help="Optional name (default: today's UTC date)")
    sync.add_argument("--refresh", action="store_true", help="Download again even if today's snapshot exists")
    sync.add_argument("--source-revision", help="Pin every public site-repository source to this commit")
    snapshots = commands.add_parser("snapshots", help="Saved snapshots by download date; the newest is the default")
    snapshots.add_argument("--json", action="store_true")
    assets = data_command("assets", help="List objects without downloading their data")
    assets.add_argument("--sample", help="Only files linked to this specimen, such as T0_tumor (see samples)")
    for field in ("kind", "format", "prefix", "contains", "timepoint", "assay", "platform", "tissue", "provider", "library"):
        assets.add_argument("--" + field)
    assets.add_argument("--include-conflicts", action="store_true")
    assets.add_argument("--include-inferred", action="store_true")
    assets.add_argument("--limit", type=int, default=50)
    variants = data_command("variants", help="List source-reported variant entries")
    variants.add_argument("--set", choices=("site", "all", "vaccine"), default="site")
    for field in ("gene", "vaccine", "pipeline", "status"):
        variants.add_argument("--" + field)
    variants.add_argument("--vaccine-source", choices=("overlap", "source_variants"), default="overlap")
    curation = data_command("curation", help="Report corrections and unrecognized source labels")
    curation.add_argument("--strict", action="store_true",
                          help="Exit 1 if any correction is stale or any source label is unrecognized")
    timeline = data_command("timeline", help="ASCII timeline of every dated source (or --list)")
    timeline.add_argument("--since", help="YYYY, YYYY-MM or YYYY-MM-DD")
    timeline.add_argument("--until")
    timeline.add_argument("--lane", help="Only lanes whose name contains this text")
    timeline.add_argument("--contains", help="Only events whose label contains this text")
    timeline.add_argument("--width", type=int)
    timeline.add_argument("--list", action="store_true", help="One line per event instead of a chart")
    timeline.add_argument("--json", action="store_true", help="Event records as JSON")
    around = data_command("on", help="Events within some days of a date")
    around.add_argument("date")
    around.add_argument("--days", type=int, default=7)
    specimens = data_command("specimens", help="Specimen registry, or one specimen's details")
    specimens.add_argument("sample_id", nargs="?")
    specimens.add_argument("--json", action="store_true")
    data_command("explore", help="Interactive terminal explorer")
    samples = data_command("samples", help="Readable specimen and sequencing overview")
    samples.add_argument("--timepoint")
    samples.add_argument("--tissue")
    samples.add_argument("--json", action="store_true", help="Original source-attributed sample claims")
    for name, description in (("timepoints", "Published timepoint and date pairs"),
                              ("vaccines", "Vaccine-overlap rows with ELISPOT results")):
        data_command(name, help=description)
    table = data_command("table", help="Parse a named site table or bucket table")
    table.add_argument("asset")
    download = data_command("download", help="Explicitly fetch one full data object")
    download.add_argument("asset")
    download.add_argument("--refresh", action="store_true")
    reads = data_command("reads", help="Extract reads around variants or regions to a cached BAM")
    reads.add_argument("asset", help="File key, URL or asset ID of an indexed BAM/CRAM in the dataset")
    reads.add_argument("regions", nargs="*", help="contig:start-end, one-based inclusive (or use --variant)")
    reads.add_argument("--variant", action="append", default=[], metavar="ID",
                       help="Catalogue variant ID with a ready allele; repeat for several")
    reads.add_argument("--padding", type=int, default=0, help="Bases added on each side of each --variant")
    reads.add_argument("--assembly", help="Assembly of explicit regions, such as GRCh38")
    reads.add_argument("--reference", help="Local indexed FASTA (required for CRAM)")
    reads.add_argument("--index", help="Explicit index path/URL when absent from catalog")
    reads.add_argument("--reference-length", type=int, help="Expected contig length (mitochondrial queries)")
    reads.add_argument("--min-mapq", type=int, default=0)
    reads.add_argument("--exclude-flags", type=lambda s: int(s, 0), default=0)
    reads.add_argument("--fetch-pairs", action="store_true", help="Also retrieve paired mates outside the regions")
    reads.add_argument("--recover-linked", action="store_true", help="Bounded mate and SA-linked recovery")
    discover = commands.add_parser("discover", help="Explicitly list a live S3 prefix")
    discover.add_argument("prefix")
    discover.add_argument("--refresh", action="store_true")
    fixtures = commands.add_parser("fixtures", help="Execute pinned read-fixture recipes")
    actions = fixtures.add_subparsers(dest="fixture_command", required=True)
    select = actions.add_parser("select", help="Return record membership and inclusion reasons")
    select.add_argument("recipe")
    select.add_argument("--source", action="append", default=[], metavar="ID=LOCAL_BAM")
    panel = actions.add_parser("panel", help="Print a shipped named target panel")
    panel.add_argument("name")
    for name in ("generate", "pack"):
        command = actions.add_parser(name, help="Select and publish a portable bundle")
        command.add_argument("recipe")
        command.add_argument("output")
        command.add_argument("--source", action="append", default=[], metavar="ID=LOCAL_BAM")
        command.add_argument("--header-policy", choices=("full", "compact"), default="full")
        command.add_argument("--size-budget", type=int, default=64 * 1024 * 1024)
    for name in ("verify", "list", "export"):
        command = actions.add_parser(name, help="Work with a bundle entirely offline")
        command.add_argument("bundle")
        if name == "export":
            command.add_argument("output")
            command.add_argument("--member", action="append")
            command.add_argument("--format", choices=("bam", "sam", "sam.gz"), default="bam")
        if name == "verify":
            command.add_argument("--sha256", help="Pinned manifest digest")
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
    from .explore import table
    shown = [{"name": r["name"], "downloaded (UTC)": (r["downloaded"] or "")[:16].replace("T", " "),
              "id": r["id"][:12]} for r in rows]
    return (table(shown, ["name", "downloaded (UTC)", "id"])
            + f"\n\nCommands use {rows[0]['name']} unless you pass --snapshot.")


def main(argv=None):
    root = parser()
    args, extra = root.parse_known_args(argv)
    # Before Python 3.13, argparse binds the optional regions positional before any
    # option, so regions written after an option arrive here as extra arguments.
    if extra and (args.command != "reads" or any(item.startswith("-") for item in extra)):
        root.error(f"unrecognized arguments: {' '.join(extra)}")
    if extra:
        args.regions = [*args.regions, *extra]
    cache = Cache(args.cache, offline=args.offline)
    try:
        if args.command == "fixtures":
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
        elif args.command == "sync":
            sources = pinned_sources(args.source_revision) if args.source_revision else None
            dataset = Dataset.sync(args.name, cache=cache, refresh=args.refresh, sources=sources,
                                   corrections=not args.no_corrections)
            value = dict(snapshot=dataset.name, id=dataset.id, downloaded=dataset.downloaded,
                         assets=len(dataset.assets), variants=len(dataset.variants()))
        elif args.command == "snapshots":
            rows = list(Dataset.snapshots(cache=cache))
            if not args.json:
                print(snapshots_view(rows))
                return 0
            value = rows
        elif args.command == "discover":
            value = list_bucket(cache, args.prefix, refresh=args.refresh)
        else:
            # Listing and parsing pinned metadata remain offline automatically;
            # commands that acquire new bytes opt in unless --offline is set.
            online = args.command in ("download", "table", "reads") and not args.offline
            dataset = Dataset.open(**snapshot_selector(args), cache=cache, offline=not online,
                                   corrections=not args.no_corrections)
            if args.command == "assets":
                filters = {name: getattr(args, name) for name in
                    ("kind", "format", "prefix", "contains", "timepoint", "assay", "platform", "tissue", "provider", "library", "include_conflicts", "include_inferred")}
                selection = (dataset.assets_for_sample(args.sample, **filters) if args.sample
                             else dataset.assets.select(**filters))
                if args.limit < 0:
                    raise ValueError("--limit must be nonnegative")
                value = dict(total=len(selection), assets=selection[:args.limit].to_records())
            elif args.command == "variants":
                value = dataset.variants(**{name: getattr(args, name) for name in
                    ("set", "gene", "vaccine", "pipeline", "status", "vaccine_source")}).to_records()
            elif args.command in ("timeline", "on", "specimens", "explore"):
                from . import explore
                if args.command == "explore":
                    explore.Explorer(dataset).cmdloop()
                    return 0
                if args.command == "on":
                    print(dataset.timeline.around(args.date, args.days).listing() or "(no events)")
                    return 0
                if args.command == "specimens":
                    if args.json:
                        value = list(dataset.specimens)
                    else:
                        print(explore.specimen_view(dataset, args.sample_id) if args.sample_id
                              else explore.specimens_view(dataset))
                        return 0
                else:
                    events = dataset.timeline.select(lane=args.lane, contains=args.contains,
                                                     since=args.since, until=args.until)
                    if args.json:
                        value = events.to_records()
                    else:
                        print(events.listing() if args.list else
                              events.render(width=args.width, since=args.since, until=args.until))
                        return 0
            elif args.command == "curation":
                value = dict(corrections=list(dataset.corrections), unrecognized=list(dataset.unrecognized))
                # Drift is judged on the evaluation, even when this run disables corrections.
                drift = [r["id"] for r in value["corrections"] if r["evaluation"] == "stale"]
                if args.strict and (drift or value["unrecognized"]):
                    print(json.dumps(value, indent=2))
                    print(f"osteosarc: stale corrections {drift}; "
                          f"{len(value['unrecognized'])} unrecognized source labels", file=sys.stderr)
                    return 1
            elif args.command == "samples":
                if args.json:
                    if args.timepoint or args.tissue:
                        raise ValueError("--json returns original sample claims; filters apply to the specimen overview")
                    value = list(dataset.samples)
                else:
                    print(dataset.describe_samples(timepoint=args.timepoint, tissue=args.tissue))
                    return 0
            elif args.command in ("timepoints", "vaccines"):
                value = list(getattr(dataset, args.command))
            elif args.command == "download":
                value = str(dataset.download(args.asset, refresh=args.refresh))
            elif args.command == "table":
                parsed = dataset.parse(args.asset)
                value = list(parsed) if hasattr(parsed, "rows") else parsed
            else:
                subset = dataset.extract_reads(args.asset, **read_targets(dataset, args),
                                               reference=args.reference, index=args.index,
                                               filters=ReadFilter(args.min_mapq, args.exclude_flags),
                                               fetch_pairs=args.fetch_pairs,
                                               recovery={} if args.recover_linked else None)
                value = dict(path=str(subset.path), index=str(subset.index_path), receipt=subset.receipt)
        print(json.dumps(value, indent=2, default=lambda x: asdict(x)))
        return 0
    except (OsteosarcError, ValueError, KeyError, OSError, subprocess.SubprocessError) as error:
        print(f"osteosarc: {error}", file=sys.stderr)
        if isinstance(error, subprocess.CalledProcessError) and error.stderr:
            print(error.stderr.decode(errors="replace") if isinstance(error.stderr, bytes) else error.stderr,
                  file=sys.stderr)
        return 1
