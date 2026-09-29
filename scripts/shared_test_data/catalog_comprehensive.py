"""Freeze/rebuild the comprehensive evidence catalogue, independently of read sampling.

python scripts/shared_test_data/catalog_comprehensive.py --freeze
python scripts/shared_test_data/catalog_comprehensive.py  # deterministic, offline
"""

import argparse
import gzip
import json
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from osteosarc import Dataset  # noqa: E402
from osteosarc.cache import digest  # noqa: E402
from osteosarc.corpus import evidence_catalog  # noqa: E402
from osteosarc.shared import read_json  # noqa: E402
from osteosarc.sv_candidates import load_sv_candidates  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
FOLDER = ROOT / "tests/data/comprehensive"


def write_gzip(path, value):
    raw = (json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False) + "\n").encode()
    path.write_bytes(gzip.compress(raw, mtime=0))


def freeze(folder):
    spec = read_json(folder / "spec.json")
    data = Dataset.open(spec["snapshot"]["name"])
    if data.id != spec["snapshot"]["id"]:
        raise ValueError("Snapshot ID differs from the pinned spec")
    aliases = read_json(folder / "count-source-aliases.json")
    bams = [f for f in data.files if f.format == "bam"]
    by_name = defaultdict(list)
    for file in bams:
        by_name[file.key.rsplit("/", 1)[-1]].append(file)
    count_sources = {}
    for name in sorted({r["bam_file"] for r in data.vafs}):
        files = [data.file(aliases[name])] if name in aliases else by_name[name]
        count_sources[name] = dict(
            identity_status="sample_catalog_alias" if name in aliases else
                            "unique_basename" if len(files) == 1 else "ambiguous" if files else "unresolved",
            files=[dict(id=f.id, key=f.key, url=f.url, samples=list(f.samples),
                        conflicts=f.conflicts) for f in files])
    samples = [{key: row[key] for key in ("id", "timepoint", "date", "tissue", "site", "description", "sequencing")}
               for row in data.samples.to_records()]
    source_files = ("variant_index", "source_variants", "vaccine_overlap", "vafs", "vaf_columns", "bams",
                    "bam_metadata", "specimens")
    inputs = dict(snapshot=spec["snapshot"], variants=data.variants("all").to_records(),
                  vaccine_rows=list(data.vaccines), counts=list(data.vafs), count_sources=count_sources,
                  samples=samples, structural_candidates=load_sv_candidates(),
                  historical_rna=read_json(folder / "historical-rna.json"),
                  provenance=dict(sources={k: {f: data.manifest["sources"][k].get(f) for f in
                                               ("url", "sha256", "size", "retrieved_at")}
                                           for k in source_files},
                                  corrections="osteosarc curated snapshot; corrections retained per row/variant",
                                  aliases=aliases, aliases_basis="Explicit T1/T2/T3/CD45neg CellRanger sample paths in bams.json",
                                  curation_sha256=digest(ROOT / "osteosarc/curation.py"),
                                  allele_resolutions_sha256=digest(ROOT / "osteosarc/data/allele_resolutions.json")))
    write_gzip(folder / "inputs.json.gz", inputs)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--freeze", action="store_true")
    parser.add_argument("--folder", type=Path, default=FOLDER)
    args = parser.parse_args()
    if args.freeze:
        freeze(args.folder)
    catalog = evidence_catalog(read_json(args.folder / "inputs.json.gz"))
    write_gzip(args.folder / "catalog.json.gz", catalog)
    checksums = {name: digest(args.folder / name) for name in
                 ("inputs.json.gz", "catalog.json.gz", "spec.json", "count-source-aliases.json", "historical-rna.json")}
    (args.folder / "checksums.json").write_text(json.dumps(checksums, indent=2, sort_keys=True) + "\n")
    (args.folder / "summary.json").write_text(json.dumps(catalog["summary"], indent=2, sort_keys=True) + "\n")
    print(json.dumps(catalog["summary"], indent=2))


if __name__ == "__main__":
    main()
