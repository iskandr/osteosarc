"""Build shared test data from its spec: the frozen recipe, the bundle, and its release record.

python scripts/shared_test_data/build.py openvax-v2 OUT_DIR --carry openvax-v1 [--required FILE ...]

1. Selects every member's records from the snapshot that
   osteosarc/data/bundles/NAME.spec.json pins (osteosarc.shared.build_shared_recipe),
   keeping each library's required records exactly: those it had in the bundle
   given with --carry, pinned by checksum, except for a library a --required
   file (made by a required_*.py script here) lists afresh. The first run
   streams reads from the public BAMs, which takes about an hour for openvax-v1;
   reruns reuse the cache.
2. Builds the bundle, with that frozen recipe in it, into OUT_DIR/NAME. The
   bundle takes the same records from the public BAMs, or fails if any changed
   upstream.
3. Packs the bundle as OUT_DIR/NAME.tar.gz and writes its release record,
   pointing at the GitHub release NAME, to OUT_DIR/NAME.release.json, and to
   osteosarc/data/bundles/NAME.release.json unless NAME is already published
   (a published record is never replaced: users' downloads are checked against it).

Then publish the archive:

    gh release create NAME OUT_DIR/NAME.tar.gz --title "Test data: NAME" --latest=false
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
    bundle_folder,
    merge_required,
    pack_release,
    read_json,
)

REPOSITORY = "https://github.com/iskandr/osteosarc"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("name", help="A spec in osteosarc/data/bundles, such as openvax-v2")
    parser.add_argument("output", type=Path, help="Folder for the bundle and its archive")
    parser.add_argument("--carry", metavar="BUNDLE",
                        help="The previous bundle, whose library members to keep: a published name, or a folder's path")
    parser.add_argument("--required", action="append", default=[], metavar="FILE",
                        help="A library's required records, replacing all that --carry gives for that library")
    parser.add_argument("--size-budget", type=int, default=1024 * 1024 * 1024)
    parser.add_argument("--cache", help="Cache directory (default: the usual one)")
    args = parser.parse_args()

    spec = read_json(BUNDLES / f"{args.name}.spec.json")
    if spec["id"] != args.name:
        parser.error(f"{args.name}.spec.json has id {spec['id']!r}; a new version needs its own id")
    dataset = Dataset.open(spec["snapshot"]["name"], cache=args.cache, offline=False)
    carried = bundle_fixtures(bundle_folder(args.carry, cache=args.cache)) if args.carry else {}
    required = merge_required(carried, args.required)
    recipe = build_shared_recipe(spec, dataset, required=required, log=lambda text: print(text, file=sys.stderr))
    args.output.mkdir(parents=True, exist_ok=True)
    bundle = args.output / args.name
    generate_bundle(recipe, bundle, dataset=dataset, size_budget=args.size_budget)
    archive = args.output / f"{args.name}.tar.gz"
    release = pack_release(bundle, archive)
    release = dict(url=f"{REPOSITORY}/releases/download/{args.name}/{archive.name}", spec_sha256=stable_id(spec),
                   **release)
    text = json.dumps(release, indent=1) + "\n"
    (args.output / f"{args.name}.release.json").write_text(text)
    published = BUNDLES / f"{args.name}.release.json"
    if published.exists():
        print(f"{published} already pins the published {args.name}; left as it is", file=sys.stderr)
    else:
        published.write_text(text)
    print(json.dumps(dict(bundle=str(bundle), archive=str(archive), **release), indent=1))


if __name__ == "__main__":
    main()
