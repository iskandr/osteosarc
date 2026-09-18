"""Command-line access using exactly the same Dataset API as Python callers."""

import argparse
import json
import subprocess
import sys
from dataclasses import asdict

from .cache import Cache
from .dataset import Dataset
from .discovery import list_bucket
from .errors import OsteosarcError
from .models import Region
from .reads import ReadFilter


def parser():
    root = argparse.ArgumentParser(description=__doc__)
    root.add_argument("--cache", help="Shared cache directory (or OSTEOSARC_CACHE)")
    root.add_argument("--offline", action="store_true", help="Forbid any network acquisition")
    commands = root.add_subparsers(dest="command", required=True)
    sync = commands.add_parser("sync", help="Acquire metadata into a named snapshot")
    sync.add_argument("snapshot")
    sync.add_argument("--refresh", action="store_true")
    sync.add_argument("--source-revision", help="Public site repository commit instead of main")
    assets = commands.add_parser("assets", help="List objects without downloading their data")
    assets.add_argument("snapshot")
    for field in ("kind", "format", "prefix", "contains", "timepoint", "assay", "platform", "tissue", "provider", "library"):
        assets.add_argument("--" + field)
    assets.add_argument("--include-conflicts", action="store_true")
    assets.add_argument("--include-inferred", action="store_true")
    assets.add_argument("--limit", type=int, default=50)
    variants = commands.add_parser("variants", help="List source-reported variant entries")
    variants.add_argument("snapshot")
    variants.add_argument("--set", choices=("site", "all", "vaccine"), default="site")
    for field in ("gene", "vaccine", "pipeline", "status"):
        variants.add_argument("--" + field)
    variants.add_argument("--vaccine-source", choices=("overlap", "source_variants"), default="overlap")
    for name in ("samples", "timepoints", "vaccines"):
        command = commands.add_parser(name)
        command.add_argument("snapshot")
    table = commands.add_parser("table", help="Parse a named site table or bucket table")
    table.add_argument("snapshot")
    table.add_argument("asset")
    download = commands.add_parser("download", help="Explicitly fetch one full data object")
    download.add_argument("snapshot")
    download.add_argument("asset")
    download.add_argument("--refresh", action="store_true")
    reads = commands.add_parser("reads", help="Extract indexed regions to a cached BAM")
    reads.add_argument("snapshot")
    reads.add_argument("asset")
    reads.add_argument("regions", nargs="+", help="contig:start-end, one-based inclusive")
    reads.add_argument("--assembly", required=True)
    reads.add_argument("--reference", help="Local indexed FASTA (required for CRAM)")
    reads.add_argument("--index", help="Explicit index path/URL when absent from catalog")
    reads.add_argument("--reference-length", type=int, help="Expected contig length (mitochondrial queries)")
    reads.add_argument("--min-mapq", type=int, default=0)
    reads.add_argument("--exclude-flags", type=lambda s: int(s, 0), default=0)
    discover = commands.add_parser("discover", help="Explicitly list a live S3 prefix")
    discover.add_argument("prefix")
    discover.add_argument("--refresh", action="store_true")
    return root


def main(argv=None):
    args = parser().parse_args(argv)
    cache = Cache(args.cache, offline=args.offline)
    try:
        if args.command == "sync":
            sources = None
            if args.source_revision:
                from .catalog import SNAPSHOT_SOURCES, SOURCE_REPO
                revision = args.source_revision
                if any(c not in "0123456789abcdefABCDEF" for c in revision) or len(revision) != 40:
                    raise ValueError("--source-revision must be a full 40-character commit SHA")
                sources = {key: url.replace(SOURCE_REPO, SOURCE_REPO.replace("/main/", f"/{revision}/"))
                           for key, url in SNAPSHOT_SOURCES.items() if url.startswith(SOURCE_REPO)}
            dataset = Dataset.sync(args.snapshot, cache=cache, refresh=args.refresh, sources=sources)
            value = dict(snapshot=dataset.manifest["name"], id=dataset.id,
                         assets=len(dataset.assets), variants=len(dataset.variants()))
        elif args.command == "discover":
            value = list_bucket(cache, args.prefix, refresh=args.refresh)
        else:
            # Listing and parsing pinned metadata remain offline automatically;
            # commands that acquire new bytes opt in unless --offline is set.
            online = args.command in ("download", "table", "reads") and not args.offline
            dataset = Dataset.open(args.snapshot, cache=cache, offline=not online)
            if args.command == "assets":
                selection = dataset.assets.select(**{name: getattr(args, name) for name in
                    ("kind", "format", "prefix", "contains", "timepoint", "assay", "platform", "tissue", "provider", "library", "include_conflicts", "include_inferred")})
                if args.limit < 0:
                    raise ValueError("--limit must be nonnegative")
                value = dict(total=len(selection), assets=selection[:args.limit].to_records())
            elif args.command == "variants":
                value = dataset.variants(**{name: getattr(args, name) for name in
                    ("set", "gene", "vaccine", "pipeline", "status", "vaccine_source")}).to_records()
            elif args.command in ("samples", "timepoints", "vaccines"):
                value = list(getattr(dataset, args.command))
            elif args.command == "download":
                value = str(dataset.download(args.asset, refresh=args.refresh))
            elif args.command == "table":
                parsed = dataset.parse(args.asset)
                value = list(parsed) if hasattr(parsed, "rows") else parsed
            else:
                regions = [Region.from_samtools(r, assembly=args.assembly,
                                               reference_length=args.reference_length) for r in args.regions]
                subset = dataset.extract_reads(args.asset, regions, reference=args.reference, index=args.index,
                                               filters=ReadFilter(args.min_mapq, args.exclude_flags))
                value = dict(path=str(subset.path), index=str(subset.index_path), receipt=subset.receipt)
        print(json.dumps(value, indent=2, default=lambda x: asdict(x)))
        return 0
    except (OsteosarcError, ValueError, KeyError, OSError, subprocess.SubprocessError) as error:
        print(f"osteosarc: {error}", file=sys.stderr)
        if isinstance(error, subprocess.CalledProcessError) and error.stderr:
            print(error.stderr.decode(errors="replace") if isinstance(error.stderr, bytes) else error.stderr,
                  file=sys.stderr)
        return 1
