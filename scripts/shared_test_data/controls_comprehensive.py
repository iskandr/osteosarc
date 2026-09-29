"""Recount and preserve all reads at the negative and discordant control SNVs.

python scripts/shared_test_data/controls_comprehensive.py build/comprehensive

These unsampled one-base extracts make the independent count audit reproducible
from the delivered data, without the much larger cached regional extracts.
"""

import argparse
import json
import shutil
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from osteosarc import Dataset, Region  # noqa: E402
from osteosarc.cache import digest  # noqa: E402
from osteosarc.corpus import AUDIT_SNVS, snv_read_support  # noqa: E402
from osteosarc.reads import resolve_regions  # noqa: E402
from osteosarc.shared import read_json  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "tests/data/comprehensive"


def recount(output, cache=None, offline=False):
    spec = read_json(INPUT / "spec.json")
    catalog = read_json(INPUT / "catalog.json.gz")
    data = Dataset.open(spec["snapshot"]["name"], cache=cache, offline=offline)
    if data.id != spec["snapshot"]["id"]:
        raise ValueError("Snapshot differs from the spec")
    destination = output / "controls"
    destination.mkdir(parents=True, exist_ok=True)
    regions = [Region(catalog["variants"][vid]["allele"][0], catalog["variants"][vid]["allele"][1] - 1,
                      catalog["variants"][vid]["allele"][1], "GRCh38") for vid in AUDIT_SNVS]

    def one(url):
        file = data.file(url)
        label = data.short_name(file)
        subset = data.extract_reads(file, regions, fetch_pairs=False, timeout=600)
        bam_path = destination / (label + ".bam")
        shutil.copyfile(subset.path, bam_path)
        shutil.copyfile(str(subset.path) + ".bai", str(bam_path) + ".bai")
        sha256 = digest(bam_path)
        rows = []
        with subset.open() as bam:
            for vid in AUDIT_SNVS:
                chrom, pos, ref, alt = catalog["variants"][vid]["allele"]
                region = resolve_regions([Region(chrom, pos - 1, pos, "GRCh38")], bam.header.to_dict())[0]
                evidence = snv_read_support(bam.fetch(region.contig, region.start, region.end), pos, ref, alt)
                rows.append(dict(variant_id=vid, source=label, source_url=file.url, samples=list(file.samples),
                                 assay=file.resolved("assay"), source_extract_sha256=sha256,
                                 extract_file=bam_path.relative_to(output).as_posix(), **evidence))
        return label, rows, subset.receipt

    rows, receipts = [], {}
    with ThreadPoolExecutor(max_workers=4) as pool:
        for future in as_completed([pool.submit(one, url) for url in spec["sources"]["all_targets"]]):
            label, found, receipt = future.result()
            rows.extend(found)
            receipts[label] = receipt
            print(f"Recounted {label}: " + ", ".join(f"{r['variant_id'].split('-chr')[0]} "
                  f"alt={r['templates']['alt']} ref={r['templates']['ref']}" for r in found), flush=True)
    result = dict(snapshot=spec["snapshot"], catalog_sha256=digest(INPUT / "catalog.json.gz"),
                  audits=sorted(rows, key=lambda r: (r["source"], r["variant_id"])),
                  acquisition=receipts, files={p.name: digest(p) for p in sorted(destination.glob("*.bam*"))},
                  interpretation="Unsampled regional reads; quality-filtered template counts, not transcript abundance.")
    (destination / "audit.json").write_text(json.dumps(result, sort_keys=True, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--cache")
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    recount(args.output, args.cache, args.offline)
