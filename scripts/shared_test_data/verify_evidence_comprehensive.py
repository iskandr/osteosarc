"""Verify the complete evidence catalogue and unsampled controls offline.

python scripts/shared_test_data/verify_evidence_comprehensive.py build/comprehensive

This check needs neither the source BAMs nor the larger balanced read bundle.
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import pysam  # noqa: E402

from osteosarc import Region  # noqa: E402
from osteosarc.cache import digest  # noqa: E402
from osteosarc.corpus import (  # noqa: E402
    AUDIT_SNVS,
    NEGATIVE_CONTROLS,
    RNA_ASSAYS,
    evidence_catalog,
    snv_read_support,
)
from osteosarc.reads import resolve_regions  # noqa: E402
from osteosarc.shared import read_json  # noqa: E402

INPUT = Path(__file__).resolve().parents[2] / "tests/data/comprehensive"


def verify(output):
    for name, checksum in read_json(INPUT / "checksums.json").items():
        if digest(INPUT / name) != checksum:
            raise ValueError(f"Evidence input changed: {name}")
    catalog = read_json(INPUT / "catalog.json.gz")
    if evidence_catalog(read_json(INPUT / "inputs.json.gz")) != catalog:
        raise ValueError("Catalogue is not reproducible from the frozen inputs")
    spec = read_json(INPUT / "spec.json")
    control = read_json(output / "controls/audit.json")
    if control["snapshot"] != spec["snapshot"] or control["catalog_sha256"] != digest(INPUT / "catalog.json.gz"):
        raise ValueError("Controls refer to different evidence inputs")
    for name, checksum in control["files"].items():
        if digest(output / "controls" / name) != checksum:
            raise ValueError(f"Control file changed: {name}")
    expected = {(url, vid) for url in spec["sources"]["all_targets"] for vid in AUDIT_SNVS}
    rows = control["audits"]
    if len(rows) != len(expected) or {(r["source_url"], r["variant_id"]) for r in rows} != expected:
        raise ValueError("Control audit is missing or duplicates variant/source pairs")
    for row in rows:
        chrom, pos, ref, alt = catalog["variants"][row["variant_id"]]["allele"]
        path = output / row["extract_file"]
        if digest(path) != row["source_extract_sha256"]:
            raise ValueError(f"Wrong source extract: {path}")
        with pysam.AlignmentFile(str(path)) as bam:
            region = resolve_regions([Region(chrom, pos - 1, pos, "GRCh38")], bam.header.to_dict())[0]
            actual = snv_read_support(bam.fetch(region.contig, region.start, region.end), pos, ref, alt)
        if any(row[k] != value for k, value in actual.items()):
            raise ValueError(f"Control recount differs: {row['source']} / {row['variant_id']}")
    rna = {f["url"] for s in catalog["count_sources"].values() if s["assay_type"] in RNA_ASSAYS for f in s["files"]}
    negative = [r for r in rows if r["source_url"] in rna and r["variant_id"] in NEGATIVE_CONTROLS]
    if len(negative) != len(rna) * len(NEGATIVE_CONTROLS):
        raise ValueError("Not all reported RNA sources were independently recounted")
    discrepancies = [r for r in negative if r["templates"]["alt"] > 0]
    report = dict(schema_version=1, component="evidence-and-unsampled-controls", snapshot=spec["snapshot"],
                  catalog_sha256=digest(INPUT / "catalog.json.gz"), controls_audit_sha256=digest(output / "controls/audit.json"),
                  catalogue_rebuilt_exactly=True, control_recounts_verified=len(rows),
                  tumor_rna_negative_recounts=len(negative), published_zero_rna_discrepancies=discrepancies,
                  negative_controls=catalog["negative_controls"], known_rna_disagreements=catalog["known_rna_disagreements"],
                  summary=catalog["summary"],
                  broader_read_bundle="Verified separately by audit_comprehensive.py; not certified by this report.")
    (output / "evidence-audit.json").write_text(json.dumps(report, sort_keys=True, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ("catalogue_rebuilt_exactly", "control_recounts_verified",
                                            "tumor_rna_negative_recounts", "published_zero_rna_discrepancies")}, indent=2))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    verify(parser.parse_args().output)
