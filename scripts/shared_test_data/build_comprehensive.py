"""Build the unpublished comprehensive corpus without changing release records.

python scripts/shared_test_data/build_comprehensive.py build/comprehensive

The recipe is checkpointed before packing. Interrupted read acquisition is cached;
rerunning reuses it. --offline proves that a completed build is reproducible locally.
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from osteosarc import Dataset, generate_bundle, verify_bundle  # noqa: E402
from osteosarc.cache import digest, stable_id  # noqa: E402
from osteosarc.provenance import evidence_overlaps  # noqa: E402
from osteosarc.shared import (  # noqa: E402
    build_shared_recipe,
    bundle_fixtures,
    bundle_folder,
    pack_release,
    read_json,
)

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "tests/data/comprehensive"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--cache")
    parser.add_argument("--offline", action="store_true")
    parser.add_argument("--workers", type=int, default=6)
    args = parser.parse_args()
    spec = read_json(INPUT / "spec.json")
    data = Dataset.open(spec["snapshot"]["name"], cache=args.cache, offline=args.offline)
    args.output.mkdir(parents=True, exist_ok=True)
    recipe_path = args.output / "recipe.json"
    stamp_path = args.output / "inputs.json"
    stamp = dict(spec_sha256=stable_id(spec), snapshot_id=data.id, carry="openvax-v2",
                 carry_regions="record_starts",
                 evidence_overlaps_sha256=digest(ROOT / "osteosarc/data/evidence-overlaps.json"))
    carried = bundle_fixtures(bundle_folder("openvax-v2", cache=data.cache), compact_regions=True)
    if recipe_path.exists():
        if read_json(stamp_path) != stamp:
            raise ValueError("Checkpoint inputs changed; choose a new output folder")
        recipe = read_json(recipe_path)
    else:
        recipe = build_shared_recipe(spec, data, required=carried, workers=args.workers,
                                     checkpoint=args.output / "selections",
                                     log=lambda message: print(message, flush=True))
        for source in recipe["sources"].values():
            if overlaps := evidence_overlaps(source["identity"]["url"]):
                source["evidence_overlaps"] = overlaps
        stamp_path.write_text(json.dumps(stamp, indent=2) + "\n")
        recipe_path.write_text(json.dumps(recipe, sort_keys=True, indent=1) + "\n")
    bundle = args.output / spec["id"]
    if not bundle.exists():
        generate_bundle(recipe, bundle, dataset=data, size_budget=1024**3)
    manifest = verify_bundle(bundle)
    if manifest["recipe_sha256"] != stable_id(recipe):
        raise ValueError("Existing bundle differs from the checkpoint")
    for name, subset in carried.items():
        if manifest["members"][name]["records"] != subset["records"]:
            raise ValueError(f"Carried fixture changed: {name}")
    archive = args.output / (spec["id"] + ".tar.gz")
    record = dict(pack_release(bundle, archive), **stamp, published=False,
                  library_fixtures_preserved=len(carried))
    (args.output / "build.json").write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(dict(bundle=str(bundle), archive=str(archive), **record), indent=2))


if __name__ == "__main__":
    main()
