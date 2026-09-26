"""Build shared test data from its spec: the frozen recipe, the bundle, and its release record.

python scripts/shared_test_data/build.py openvax-v2 OUT_DIR --carry openvax-v1 [--required FILE ...]

1. Selects every member's records from the snapshot that
   osteosarc/data/bundles/NAME.spec.json pins (osteosarc.shared.build_shared_recipe),
   keeping each library's required records exactly: those it had in the bundles
   given with --carry, except where a --required file (made by a required_*.py
   script here) gives that library's records afresh. The first run streams reads
   from the public BAMs, which takes about an hour for openvax-v1; reruns reuse
   the cache.
2. Builds the bundle, with that frozen recipe in it, into OUT_DIR/NAME. The
   bundle takes the same records from the public BAMs, or fails if any changed
   upstream.
3. Packs the bundle as OUT_DIR/NAME.tar.gz and writes
   osteosarc/data/bundles/NAME.release.json, pointing at the GitHub release NAME.

Then publish the archive:

    gh release create NAME OUT_DIR/NAME.tar.gz --title "Test data: NAME"
"""

import argparse
import json
import sys
from pathlib import Path

# This checkout's osteosarc, whatever is installed: the spec and release record live in it.
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from osteosarc import Dataset, generate_bundle  # noqa: E402
from osteosarc.cache import stable_id  # noqa: E402
from osteosarc.shared import (  # noqa: E402
    BUNDLES,
    build_shared_recipe,
    bundle_fixtures,
    fetch_bundle,
    load_required,
    pack_release,
    read_json,
)

REPOSITORY = "https://github.com/iskandr/osteosarc"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("name", help="A spec in osteosarc/data/bundles, such as openvax-v1")
    parser.add_argument("output", type=Path, help="Folder for the bundle and its archive")
    parser.add_argument("--carry", action="append", default=[], metavar="BUNDLE",
                        help="Keep the library members of this published bundle (or folder); repeat for several")
    parser.add_argument("--required", action="append", default=[], metavar="FILE",
                        help="A library's required records, replacing what --carry gives for that library")
    parser.add_argument("--size-budget", type=int, default=1024 * 1024 * 1024)
    parser.add_argument("--cache", help="Cache directory (default: the usual one)")
    args = parser.parse_args()

    spec = read_json(BUNDLES / f"{args.name}.spec.json")
    dataset = Dataset.open(spec["snapshot"]["name"], cache=args.cache, offline=False)
    carried = {}
    for bundle in args.carry:
        carried.update(bundle_fixtures(bundle if Path(bundle).is_dir() else fetch_bundle(bundle, cache=args.cache)))
    required = load_required(args.required)
    fresh = {subset["consumer"] for subset in required.values()}
    kept = {name: s for name, s in carried.items() if s["consumer"] not in fresh}
    if clash := sorted(set(kept) & set(required)):
        parser.error(f"carried and required both name {', '.join(clash)}")
    required.update(kept)
    recipe = build_shared_recipe(spec, dataset, required=required, log=lambda text: print(text, file=sys.stderr))
    args.output.mkdir(parents=True, exist_ok=True)
    bundle = args.output / args.name
    generate_bundle(recipe, bundle, dataset=dataset, size_budget=args.size_budget)
    archive = args.output / f"{args.name}.tar.gz"
    release = pack_release(bundle, archive)
    release = dict(url=f"{REPOSITORY}/releases/download/{args.name}/{archive.name}", spec_sha256=stable_id(spec),
                   **release)
    (BUNDLES / f"{args.name}.release.json").write_text(json.dumps(release, indent=1) + "\n")
    print(json.dumps(dict(bundle=str(bundle), archive=str(archive), **release), indent=1))


if __name__ == "__main__":
    main()
