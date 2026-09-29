import copy
import json
from pathlib import Path

import pytest

from osteosarc import File, IntegrityError
from osteosarc.cache import digest
from osteosarc.models import Files
from osteosarc.provenance import audit_alignment_overlap, evidence_overlaps

DATA = Path(__file__).parent / "data/evidence-overlap"
REGISTRY = Path(__file__).parents[1] / "osteosarc/data/evidence-overlaps.json"


def test_real_overlap_is_reproducible_including_qualities_and_multiplicity():
    group, = json.loads(REGISTRY.read_text())["groups"]
    paths = [DATA / source["regional_bam"] for source in group["sources"]]
    receipts = json.loads((DATA / "receipts.json").read_text())
    for source, path in zip(group["sources"], paths):
        assert digest(path) == source["regional_bam_sha256"]
        assert receipts[source["id"]]["request"]["source"] == source["url"]
        assert receipts[source["id"]]["files"]["reads.bam"] == digest(path)
    observed = audit_alignment_overlap(*paths)
    assert observed == group["audit"]
    assert observed["shared_records_excluding_tags"] == 40684
    assert observed["different_a"] == observed["different_b"] == 12
    assert observed["shared_names"] == observed["union_names"] == 18971


def test_selection_preserves_labels_requires_explicit_preference_and_reports_exclusion():
    group, = json.loads(REGISTRY.read_text())["groups"]
    sources = group["sources"]
    files = Files([File(s["id"], s["url"].split(".com/", 1)[1], s["url"], "alignment", "bam",
                        metadata=dict(original_label=s["id"])) for s in sources], source={"snapshot_id": "pinned"})
    before = copy.deepcopy([f.metadata for f in files])
    assert len(files.overlap_groups) == 1
    with pytest.raises(IntegrityError, match="known overlapping evidence"):
        files.require_no_known_overlaps()
    with pytest.raises(ValueError, match="exactly one"):
        files.without_known_overlaps()
    with pytest.raises(ValueError, match="exactly one"):
        files.without_known_overlaps(prefer=[f.url for f in files])
    selected = files.without_known_overlaps(prefer=[files[1].url])
    assert list(selected) == [files[1]]
    assert selected.require_no_known_overlaps() is selected
    assert selected.source["snapshot_id"] == "pinned"
    assert selected.source["excluded_evidence_overlaps"][0]["url"] == files[0].url
    assert [f.metadata for f in files] == before
    # The same sample name in a different processing product is not enough.
    assert evidence_overlaps(files[0].url.replace("out.md.bam", "out.bam")) == []
    found = files[0].evidence_overlaps
    found[0]["status"] = "modified"
    assert files[0].evidence_overlaps[0]["status"] == "biological_identity_unresolved"
