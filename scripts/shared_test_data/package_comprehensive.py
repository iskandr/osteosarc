"""Pack the evidence, verified reads, audit and reproduction code in one archive."""

import argparse
import gzip
import io
import json
import subprocess
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

from osteosarc import verify_bundle  # noqa: E402
from osteosarc.cache import digest  # noqa: E402
from osteosarc.shared import read_json  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
INPUT = ROOT / "tests/data/comprehensive"


def package(output, *, evidence_only=False):
    from scripts.shared_test_data.verify_evidence_comprehensive import verify
    evidence_audit = verify(output)
    spec = read_json(INPUT / "spec.json")
    bundle = output / spec["id"]
    if not evidence_only:
        verify_bundle(bundle)
        audit = read_json(output / "read-audit.json.gz")
        if audit["manifest_sha256"] != digest(bundle / "manifest.json"):
            raise ValueError("Audit is for another read bundle")
        if audit["catalog_sha256"] != digest(INPUT / "catalog.json.gz"):
            raise ValueError("Audit is for another evidence catalogue")
        for name, field in (("controls/audit.json", "controls_audit_sha256"),
                            ("portable-recipe.json", "portable_recipe_sha256"),
                            ("reference-contexts.json.gz", "reference_contexts_sha256")):
            if audit[field] != digest(output / name):
                raise ValueError(f"Audited artifact changed: {name}")
    for name, checksum in read_json(output / "controls/audit.json")["files"].items():
        if digest(output / "controls" / name) != checksum:
            raise ValueError(f"Audited control extract changed: {name}")
    for name, checksum in read_json(INPUT / "checksums.json").items():
        if digest(INPUT / name) != checksum:
            raise ValueError(f"Evidence input changed: {name}")
    files = {"evidence/" + p.name: p for p in INPUT.iterdir() if p.is_file()}
    discrepancies = evidence_audit["published_zero_rna_discrepancies"]
    if discrepancies:
        rna_result = (
            f"Independent recounts found high-quality ALT templates in {len(discrepancies)} "
            "control/source comparisons despite published zero RNA ALT counts. "
            "These controls are not uniformly RNA-negative; see published_zero_rna_discrepancies "
            "in audit/evidence-audit.json for the affected variants, sources and counts. ")
    else:
        rna_result = (
            f"The {len(evidence_audit['negative_controls'])} controls have zero high-quality ALT templates "
            f"in all {evidence_audit['tumor_rna_products_checked']} tested tumor RNA products. ")
    readme_path = output / ("README-evidence.md" if evidence_only else "README-dataset.md")
    readme_path.write_text(
        "# Osteosarc comprehensive test dataset\n\n"
        + ("This archive contains the complete evidence catalogue and independently verified, unsampled control BAMs. "
           "The larger balanced openvax-v3 read bundle is not included. Its build specification and code are included.\n\n"
           if evidence_only else "This archive contains the evidence catalogue, unsampled controls, and verified balanced read bundle.\n\n")
        + "See evidence/README.md for provenance, counting units, known exceptions and reproduction details. "
          "manifest.json pins every delivered file by SHA-256. Dataset: CC0-1.0; code: Apache-2.0.\n\n"
          "From this extracted directory, with Python and the dependencies in reproduction/pyproject.toml installed:\n\n"
          "```sh\npython reproduction/scripts/shared_test_data/verify_evidence_comprehensive.py .\n```\n\n"
          f"This reproduces the catalogue and all {evidence_audit['control_recounts_verified']} control recounts offline. "
          "It needs no source BAM downloads or website snapshot. "
        + rna_result
        + "Coverage limitations, BTN3A3's blood-RNA summary exception, and NR2F2's historical RNA evidence are preserved explicitly.\n")
    files["README.md"] = readme_path
    if not evidence_only:
        files.update({"reads/" + str(p.relative_to(bundle)): p for p in bundle.rglob("*") if p.is_file()})
    files.update({"controls/" + p.name: p for p in (output / "controls").iterdir() if p.is_file()})
    files["audit/evidence-audit.json"] = output / "evidence-audit.json"
    if not evidence_only:
        files.update({"audit/" + name: output / name for name in
                  ("read-audit.json.gz", "audit-summary.json", "build.json", "portable-recipe.json",
                   "reference-contexts.json.gz")})
    files.update({"reproduction/" + str(p.relative_to(ROOT)): p for p in
                  [*sorted((ROOT / "osteosarc").glob("*.py")),
                   *sorted((ROOT / "scripts/shared_test_data").glob("*comprehensive.py")),
                   ROOT / "pyproject.toml", ROOT / "README.md", ROOT / "LICENSE", ROOT / "MANIFEST.in"]})
    for folder in (ROOT / "osteosarc/data", ROOT / "tests"):
        files.update({"reproduction/" + str(p.relative_to(ROOT)): p for p in folder.rglob("*")
                      if p.is_file() and "__pycache__" not in p.parts})
    suffix = "-evidence" if evidence_only else ""
    manifest = dict(schema_version=1, id="osteosarc-comprehensive-2026-09-28" + suffix, published=False,
                    component="evidence-and-unsampled-controls" if evidence_only else "evidence-controls-and-balanced-reads",
                    snapshot=spec["snapshot"], redistribution=spec["redistribution"],
                    git_base=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                    files={name: dict(sha256=digest(path), size_bytes=path.stat().st_size)
                           for name, path in sorted(files.items())})
    raw_manifest = (json.dumps(manifest, sort_keys=True, indent=2) + "\n").encode()
    manifest_path = output / ("dataset" + suffix + "-manifest.json")
    manifest_path.write_bytes(raw_manifest)
    archive = output / (manifest["id"] + ".tar.gz")
    with archive.open("wb") as raw, gzip.GzipFile(filename="", fileobj=raw, mode="wb", mtime=0) as compressed, \
            tarfile.open(fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT) as tar:
        info = tarfile.TarInfo("manifest.json")
        info.size, info.mode = len(raw_manifest), 0o644
        tar.addfile(info, io.BytesIO(raw_manifest))
        for name, path in sorted(files.items()):
            info = tarfile.TarInfo(name)
            info.size, info.mode = path.stat().st_size, 0o644
            with path.open("rb") as handle:
                tar.addfile(info, handle)
    record = dict(archive=archive.name, sha256=digest(archive), size_bytes=archive.stat().st_size,
                  manifest_sha256=digest(manifest_path), published=False)
    (output / ("dataset" + suffix + "-archive.json")).write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--evidence-only", action="store_true", help="Package the completed evidence and unsampled controls independently of the larger read build")
    args = parser.parse_args()
    package(args.output, evidence_only=args.evidence_only)
