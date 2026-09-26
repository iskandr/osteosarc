"""The shared interest set must retain unproductive and unassessed SVs too."""

import json

from osteosarc import load_panel, load_sv_candidates, validate_recipe
from osteosarc.cli import main


def test_interest_set_keeps_requested_candidates_and_all_nomination_outcomes():
    catalogue = load_sv_candidates()
    targets = catalogue["targets"]
    requested = {"SV0203", "SV0089", "SV0085", "SV0078", "SV0368", "SV0381",
                 "SV0281", "SV0111", "SV0499", "SV0172", "SV0377"}
    assert requested <= targets.keys()
    assert len(targets) == 637
    assert {"SV0634", "SV0635", "SV0636", "SV0637"} <= targets.keys()
    assert any(t["kind"] == "unresolved" for t in targets.values())
    assert any(not t["dna_samples"] for t in targets.values())
    validate_recipe(dict(schema_version=1, id="all-SVs", targets=targets, sources={}, members={}))
    assert targets["SV0203"]["dna_samples"] == ["T1", "T2"]
    assert targets["SV0089"]["dna_samples"] == ["T0"]
    assert targets["SV0368"]["dna_samples"] == ["T1-organoid"]


def test_api_panel_and_cli_preserve_identical_metadata_and_return_copies(capsys):
    catalogue = load_sv_candidates()
    expected = catalogue["targets"]
    assert load_panel("sv-candidates-v1") == expected
    assert main(["--offline", "fixtures", "panel", "sv-candidates-v1"]) == 0
    assert json.loads(capsys.readouterr().out) == expected
    catalogue["targets"].clear()
    assert len(load_sv_candidates()["targets"]) == 637


def test_expression_and_missing_rna_remain_sample_and_source_scoped():
    targets = load_sv_candidates()["targets"]
    t = targets["SV0203"]
    expression, = [r for r in t["gene_expression"] if r["gene"] == "MACROD2"]
    assert expression["sample_id"] == "T2"
    assert expression["tpm"] == 1.40
    assert not expression["allele_specific"]
    t1, = [r for r in t["rna_evidence"] if r["source_id"] == "T1-ONT-tagged"]
    t3, = [r for r in t["rna_evidence"] if r["source_id"] == "T3-ONT-tagged"]
    assert t1["split_path_templates"] == 0
    assert t3["split_path_templates"] is None and t3["status"] == "not_complete"
    # Same adjacency reported by separate callers is not two independent events.
    assert targets["SV0179"]["adjacency_group_id"] == targets["SV0546"]["adjacency_group_id"]
