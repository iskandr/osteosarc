"""Verify the built corpus and recount negative-control SNVs before sampling.

python scripts/shared_test_data/audit_comprehensive.py build/comprehensive

Recounts the unsampled control extracts, not the balanced fixture BAMs. No
network is used. Run controls_comprehensive.py first to preserve those extracts.
"""

import argparse
import copy
import gzip
import hashlib
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from osteosarc import Dataset, Region, verify_bundle  # noqa: E402
from osteosarc.cache import digest, stable_id  # noqa: E402
from osteosarc.corpus import (  # noqa: E402
    AUDIT_SNVS,
    NEGATIVE_CONTROLS,
    RNA_ASSAYS,
    snv_read_support,
    verify_control_extracts,
)
from osteosarc.reads import resolve_regions  # noqa: E402
from osteosarc.reference import reference_sequence  # noqa: E402
from osteosarc.shared import bundle_fixtures, bundle_folder, read_json  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "tests/data/comprehensive"


def audit(output, cache=None):
    spec = read_json(INPUT / "spec.json")
    catalog = read_json(INPUT / "catalog.json.gz")
    control_data = read_json(output / "controls/audit.json")
    if control_data["snapshot"] != spec["snapshot"] or control_data["catalog_sha256"] != digest(INPUT / "catalog.json.gz"):
        raise ValueError("Control extracts refer to other inputs")
    verify_control_extracts(output, control_data)
    bundle = output / spec["id"]
    manifest = verify_bundle(bundle)
    recipe = read_json(bundle / "recipe.json")
    if recipe["snapshot"] != spec["snapshot"] or stable_id(recipe) != manifest["recipe_sha256"]:
        raise ValueError("Bundle snapshot or recipe changed")
    data = Dataset.open(spec["snapshot"]["name"], cache=cache, offline=True)
    if data.id != recipe["snapshot"]["id"]:
        raise ValueError("Audit snapshot differs")
    previous_bundle = bundle_folder("openvax-v2", cache=data.cache)
    previous_recipe = read_json(previous_bundle / "recipe.json")
    carried = bundle_fixtures(previous_bundle)
    for name, required in carried.items():
        if manifest["members"][name]["records"] != required["records"]:
            raise ValueError(f"Historical library records changed: {name}")
    targets = recipe["targets"]
    # Exact-record rebuilding needs no website snapshot when the public index
    # URLs accompany the frozen source identities. The original recipe stays intact.
    portable = copy.deepcopy(recipe)
    for source in portable["sources"].values():
        file = data.file(source["identity"]["key"])
        if not file.index_urls:
            raise ValueError(f"No public index for {file.key}")
        source["identity"].update(index_urls=list(file.index_urls), format=file.format)
    (output / "portable-recipe.json").write_text(json.dumps(portable, sort_keys=True, indent=1) + "\n")
    contexts = {}
    for name, target in targets.items():
        context = target.get("window", {}).get("context")
        if context is None:
            continue
        sequence = reference_sequence(target["contig"], context["start"], context["end"], target["assembly"],
                                      cache=data.cache, reference_length=target.get("reference_length"))
        if hashlib.sha256(sequence.encode()).hexdigest() != context["sha256"]:
            raise ValueError(f"Reference context changed: {name}")
        offset = target["position"] - 1 - context["start"]
        if sequence[offset:offset + len(target["ref"])] != target["ref"]:
            raise ValueError(f"REF mismatch: {name}")
        contexts[name] = dict(contig=target["contig"], assembly=target["assembly"], **context, sequence=sequence)
    (output / "reference-contexts.json.gz").write_bytes(gzip.compress(
        (json.dumps(contexts, sort_keys=True, separators=(",", ":")) + "\n").encode(), mtime=0))
    vaccinated = {vid for vid, v in catalog["variants"].items() if v["vaccine"]["included"]}
    by_pair = {(m["source"], m["target"]): (name, m) for name, m in recipe["members"].items()}
    expected_sources = set(spec["sources"]["all_targets"])
    actual_sources = {s["identity"]["url"] for s in recipe["sources"].values()}
    if expected_sources - actual_sources:
        raise ValueError("An all-target source is absent")
    for sid, source in recipe["sources"].items():
        if source["identity"]["url"] not in expected_sources:
            continue
        for vid in vaccinated:
            if (sid, vid) not in by_pair:
                raise ValueError(f"Vaccine allele not covered in source plan: {sid} / {vid}")
            if targets[vid]["kind"] != "small_variant" or "unavailable" in targets[vid]["window"]:
                raise ValueError(f"Vaccine allele is not reference validated: {vid}")
    members = []
    for name, member in recipe["members"].items():
        target = targets[member["target"]]
        members.append(dict(name=name, target=member["target"], source=member["source"],
                            target_kind=target["kind"], observed_before_sampling=member.get("observed"),
                            selected_record_count=manifest["members"][name]["record_count"],
                            selected_record_sha256=stable_id(manifest["members"][name]["records"]),
                            status=manifest["members"][name]["status"],
                            interpretation="Geometric joining templates, not exact DNA-allele confirmation."
                            if target["kind"] == "sv" else "Balanced fixtures are not abundance estimates."))
    audits = []
    import pysam
    if {a["source_url"] for a in control_data["audits"]} != expected_sources:
        raise ValueError("Control extracts do not cover every all-target source")
    for sid, source in recipe["sources"].items():
        if source["identity"]["url"] not in expected_sources:
            continue
        saved = [a for a in control_data["audits"] if a["source_url"] == source["identity"]["url"]]
        if len(saved) != len(AUDIT_SNVS) or {a["variant_id"] for a in saved} != set(AUDIT_SNVS):
            raise ValueError(f"Control audit has missing/duplicate rows: {sid}")
        with pysam.AlignmentFile(str(output / saved[0]["extract_file"])) as bam:
            for vid in AUDIT_SNVS:
                chrom, pos, ref, alt = catalog["variants"][vid]["allele"]
                region = resolve_regions([Region(chrom, pos - 1, pos, "GRCh38")], bam.header.to_dict())[0]
                evidence = snv_read_support(bam.fetch(region.contig, region.start, region.end), pos, ref, alt)
                original = next(a for a in saved if a["variant_id"] == vid)
                if any(original[k] != v for k, v in evidence.items()):
                    raise ValueError(f"Independent control counts changed: {sid} / {vid}")
                audits.append(original)
        print(f"Audited {sid}", flush=True)
    # Prove that every reported tumor RNA product was independently checked.
    expected_rna = {f["url"] for s in catalog["count_sources"].values() if s["assay_type"] in RNA_ASSAYS
                    for f in s["files"]}
    # Do not certify controls across two labels for the reviewed Tempus overlap.
    from osteosarc.models import Files
    Files(data.file(url) for url in expected_rna).require_no_known_overlaps()
    for vid in NEGATIVE_CONTROLS:
        if expected_rna - {a["source_url"] for a in audits if a["variant_id"] == vid}:
            raise ValueError(f"Negative control RNA audit incomplete: {vid}")
    discrepancies = [dict(variant_id=a["variant_id"], source=a["source"], alt_templates=a["templates"]["alt"])
                     for a in audits if a["source_url"] in expected_rna and a["templates"]["alt"] > 0]
    historical = []
    for row in catalog["known_rna_disagreements"]:
        original = previous_recipe["members"][row["member"]]["policy"]["records"]
        current = [name for name, member in recipe["members"].items()
                   if recipe["sources"][member["source"]]["identity"]["url"] == row["source_url"]
                   and all(manifest["members"][name]["records"].get(key, 0) >= count
                           for key, count in original.items())]
        if not current:
            raise ValueError(f"Historical RNA disagreement records were lost: {row['member']}")
        historical.append(dict(row, current_read_members=sorted(current),
                               historical_selected_records_preserved=sum(original.values())))
    report = dict(schema_version=1, snapshot=spec["snapshot"], manifest_sha256=digest(bundle / "manifest.json"),
                  catalog_sha256=digest(INPUT / "catalog.json.gz"), suitable_for_abundance=False,
                  controls_audit_sha256=digest(output / "controls/audit.json"),
                  portable_recipe_sha256=digest(output / "portable-recipe.json"),
                  reference_contexts_sha256=digest(output / "reference-contexts.json.gz"),
                  library_fixtures_preserved=len(carried), vaccine_variants=len(vaccinated),
                  verified_reference_contexts=len(contexts),
                  all_target_sources=len(expected_sources), source_count=len(recipe["sources"]),
                  member_count=len(members), statuses=dict(Counter(m["status"] for m in members)),
                  unresolved_targets={k: t for k, t in targets.items() if t["kind"] == "unresolved"},
                  unavailable_allele_windows={k: t["window"] for k, t in targets.items() if "unavailable" in t.get("window", {})},
                  sources=recipe["sources"], members=members, snv_control_audit=audits,
                  published_zero_rna_discrepancies=discrepancies,
                  known_historical_rna_disagreements=historical,
                  supplementary_policy=spec["selection"],
                  supplementary_limitations="MAPQ-zero, malformed, absent-contig SA entries and regions exceeding split_depth are excluded; see acquisition log.")
    (output / "read-audit.json.gz").write_bytes(gzip.compress(
        (json.dumps(report, sort_keys=True, separators=(",", ":")) + "\n").encode(), mtime=0))
    summary = {k: report[k] for k in ("vaccine_variants", "all_target_sources", "source_count", "member_count", "statuses",
                                     "library_fixtures_preserved", "unresolved_targets", "unavailable_allele_windows",
                                     "published_zero_rna_discrepancies", "known_historical_rna_disagreements")}
    (output / "audit-summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--cache")
    args = parser.parse_args()
    audit(args.output, args.cache)
