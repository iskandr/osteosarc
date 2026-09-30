"""Receipt attribution and truthful packaging, using small indexed control BAMs."""

import copy
import gzip
import json
import tarfile
from pathlib import Path

import pysam
import pytest

from osteosarc import IntegrityError
from osteosarc.cache import digest
from osteosarc.corpus import snv_read_support, verify_control_extracts
from scripts.shared_test_data import (
    audit_comprehensive,
    package_comprehensive,
    verify_evidence_comprehensive,
)

INPUT = Path(__file__).parent / "data/comprehensive"


@pytest.fixture
def controls(bam, tmp_path):
    output = tmp_path / "output"
    directory = output / "controls"
    directory.mkdir(parents=True)
    control = dict(snapshot=json.loads((INPUT / "spec.json").read_text())["snapshot"],
                   catalog_sha256=digest(INPUT / "catalog.json.gz"), files={}, audits=[], acquisition={})
    for label in ("rna-a", "rna-b"):
        path = directory / (label + ".bam")
        with pysam.AlignmentFile(bam) as source, pysam.AlignmentFile(path, "wb", template=source) as out:
            for read in source:
                read.query_name = label + ":" + read.query_name
                out.write(read)
        pysam.index(str(path))
        url = "https://example.test/" + path.name
        sha256, index_sha256 = digest(path), digest(str(path) + ".bai")
        control["files"].update({path.name: sha256, path.name + ".bai": index_sha256})
        control["acquisition"][label] = dict(
            request=dict(source=url, snapshot_id=control["snapshot"]["id"]),
            files={"reads.bam": sha256, "reads.bam.bai": index_sha256})
        control["audits"].append(dict(source=label, source_url=url, variant_id="control",
            extract_file="controls/" + path.name, source_extract_sha256=sha256))
    return output, control


def test_controls_bind_the_delivered_files_to_receipts(controls):
    output, control = controls
    verify_control_extracts(output, control)


@pytest.mark.parametrize("change, message", [
    ("urls", "source URL"), ("bam", "Control BAM"), ("index", "Control index"),
    ("missing_receipt", "acquisition receipts"), ("snapshot", "snapshot"),
    ("missing_inventory", "Control BAM"), ("conflicting_rows", "Conflicting control"),
])
def test_recountable_files_cannot_override_the_original_receipts(controls, change, message):
    output, control = controls
    a, b = control["audits"]
    if change == "urls":
        a["source_url"], b["source_url"] = b["source_url"], a["source_url"]
    elif change == "bam":
        # The substituted extract still agrees with its inventory and row hash.
        a["extract_file"], a["source_extract_sha256"] = b["extract_file"], b["source_extract_sha256"]
    elif change == "index":
        path = output / (a["extract_file"] + ".bai")
        path.write_bytes(path.read_bytes() + b"changed")
        control["files"][path.name] = digest(path)
    elif change == "missing_receipt":
        del control["acquisition"][a["source"]]
    elif change == "snapshot":
        control["acquisition"][a["source"]]["request"]["snapshot_id"] = "another-snapshot"
    elif change == "missing_inventory":
        del control["files"][Path(a["extract_file"]).name]
    elif change == "conflicting_rows":
        control["audits"].append(dict(a, extract_file=b["extract_file"]))
    with pytest.raises(IntegrityError, match=message):
        verify_control_extracts(output, control)


@pytest.mark.parametrize("audit", [verify_evidence_comprehensive.verify, audit_comprehensive.audit])
def test_both_audits_reject_swapped_source_urls_before_recounting(controls, audit):
    output, control = controls
    a, b = control["audits"]
    a["source_url"], b["source_url"] = b["source_url"], a["source_url"]
    (output / "controls/audit.json").write_text(json.dumps(control))
    with pytest.raises(IntegrityError, match="source URL differs"):
        audit(output)


@pytest.mark.parametrize("has_alt", [False, True])
def test_packaged_readme_matches_actual_control_recounts(controls, tmp_path, monkeypatch, has_alt):
    output, control = controls
    root = tmp_path / "reproduction"
    inputs = root / "tests/data/comprehensive"
    inputs.mkdir(parents=True)
    # Only the catalogue fixture is stubbed; receipt checks, indexed recounts,
    # evidence report, archive and its README are all produced by the real code.
    ref, alt = ("C", "A") if has_alt else ("A", "C")
    catalog = dict(variants={"control": dict(allele=["chr1", 101, ref, alt])},
        count_sources={row["source"]: dict(assay_type="RNA", files=[dict(url=row["source_url"])])
                       for row in control["audits"]},
        negative_controls=[dict(variant_id="control")], known_rna_disagreements=[], summary={})
    (inputs / "catalog.json.gz").write_bytes(gzip.compress(json.dumps(catalog).encode()))
    (inputs / "inputs.json.gz").write_bytes(gzip.compress(b"{}"))
    spec = dict(id="test-reads", snapshot=control["snapshot"], redistribution=dict(license="CC0-1.0"),
                sources=dict(all_targets=[row["source_url"] for row in control["audits"]]))
    (inputs / "spec.json").write_text(json.dumps(spec))
    (inputs / "checksums.json").write_text(json.dumps({p.name: digest(p) for p in inputs.iterdir()}))
    for name in ("pyproject.toml", "README.md", "LICENSE", "MANIFEST.in"):
        (root / name).write_text("test fixture\n")
    control["catalog_sha256"] = digest(inputs / "catalog.json.gz")
    for row in control["audits"]:
        with pysam.AlignmentFile(output / row["extract_file"]) as source:
            row.update(snv_read_support(source.fetch("chr1", 100, 101), 101, ref, alt))
    (output / "controls/audit.json").write_text(json.dumps(control))
    monkeypatch.setattr(verify_evidence_comprehensive, "INPUT", inputs)
    monkeypatch.setattr(verify_evidence_comprehensive, "AUDIT_SNVS", ("control",))
    monkeypatch.setattr(verify_evidence_comprehensive, "NEGATIVE_CONTROLS", ("control",))
    monkeypatch.setattr(verify_evidence_comprehensive, "evidence_catalog", lambda inputs: copy.deepcopy(catalog))
    monkeypatch.setattr(package_comprehensive, "INPUT", inputs)
    monkeypatch.setattr(package_comprehensive, "ROOT", root)
    monkeypatch.setattr(package_comprehensive.subprocess, "check_output", lambda *a, **kw: "review-test\n")
    package_comprehensive.package(output, evidence_only=True)
    archive, = output.glob("*.tar.gz")
    with tarfile.open(archive) as tar:
        readme = tar.extractfile("README.md").read().decode()
        report = json.load(tar.extractfile("audit/evidence-audit.json"))
    assert "all 2 control recounts" in readme
    assert report["tumor_rna_products_checked"] == 2
    if has_alt:
        assert len(report["published_zero_rna_discrepancies"]) == 2
        assert all(row["templates"]["alt"] == 2 for row in report["published_zero_rna_discrepancies"])
        assert "not uniformly RNA-negative" in readme
        assert "zero high-quality ALT templates" not in readme
    else:
        assert report["published_zero_rna_discrepancies"] == []
        assert "zero high-quality ALT templates in all 2 tested tumor RNA products" in readme
