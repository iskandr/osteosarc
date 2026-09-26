"""Build shared test data from its spec: the frozen recipe, the bundle, and its release record.

python scripts/shared_test_data/build.py openvax-v1 OUT_DIR --required required-isovar.json.gz ...

1. Selects every member's records from the snapshot that
   osteosarc/data/bundles/NAME.spec.json pins (osteosarc.shared.build_shared_recipe),
   keeping each library's required records (made by the required_*.py scripts
   here) exactly. The first run streams reads from the public BAMs, which takes
   about an hour for openvax-v1; reruns reuse the cache.
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
from osteosarc.shared import BUNDLES, build_shared_recipe, pack_release, read_json  # noqa: E402

REPOSITORY = "https://github.com/iskandr/osteosarc"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("name", help="A spec in osteosarc/data/bundles, such as openvax-v1")
    parser.add_argument("output", type=Path, help="Folder for the bundle and its archive")
    parser.add_argument("--required", action="append", default=[], metavar="FILE",
                        help="A library's required records; repeat for several")
    parser.add_argument("--size-budget", type=int, default=1024 * 1024 * 1024)
    parser.add_argument("--cache", help="Cache directory (default: the usual one)")
    args = parser.parse_args()

    spec = read_json(BUNDLES / f"{args.name}.spec.json")
    dataset = Dataset.open(spec["snapshot"]["name"], cache=args.cache, offline=False)
    recipe = build_shared_recipe(spec, dataset, required=args.required, log=lambda text: print(text, file=sys.stderr))
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
