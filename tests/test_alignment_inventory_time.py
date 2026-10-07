from collections import Counter
from dataclasses import replace
from email.utils import formatdate
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from threading import Thread

import pytest

from osteosarc import (
    Cache,
    File,
    IntegrityError,
    RecoveryPolicy,
    Region,
    count_reads,
    extract_reads,
)
from osteosarc.cache import _check_inventory_time, stable_id, write_json
from osteosarc.read_receipts import write_read_receipt
from osteosarc.reads import inspect_alignment

EPOCH = 1784586023
REGIONS = [Region("chr1", 100, 160, "GRCh38")]


@pytest.fixture
def remote_bam(bam):
    """Serve actual BAM bytes with a controllable HTTP modification date."""
    state = dict(last_modified=formatdate(EPOCH, usegmt=True), requests=[])

    class Handler(SimpleHTTPRequestHandler):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, directory=str(bam.parent), **kwargs)

        def send_header(self, keyword, value):
            if keyword.lower() == "last-modified":
                value = state["last_modified"]
                if value is None:
                    return
            super().send_header(keyword, value)

        def do_HEAD(self):
            state["requests"].append("HEAD")
            super().do_HEAD()

        def do_GET(self):
            state["requests"].append("GET")
            super().do_GET()

        def log_message(self, *args):
            pass

    server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    thread = Thread(target=server.serve_forever, daemon=True)
    thread.start()
    source = File("remote", "input.bam", f"http://127.0.0.1:{server.server_port}/input.bam",
                  "alignment", "bam", size=bam.stat().st_size, modified=EPOCH)
    try:
        yield source, bam, state
    finally:
        server.shutdown()
        server.server_close()
        thread.join()


def acquire(operation, source, bam, cache):
    if operation == "inspect":
        return inspect_alignment(source, cache=cache, snapshot_id="frozen")
    kwargs = dict(cache=cache, index=str(bam) + ".bai", snapshot_id="frozen", max_records=6)
    if operation == "count":
        return count_reads(source, REGIONS, **kwargs)
    if operation == "recover":
        kwargs.pop("max_records")
        kwargs["recovery"] = RecoveryPolicy(max_rounds=1, max_records=100)
    return extract_reads(source, REGIONS, **kwargs)


@pytest.mark.parametrize("operation", ["inspect", "extract", "recover", "count"])
@pytest.mark.parametrize("offset", [-86400, 86400])
def test_remote_inventory_date_mismatch_fails_before_reading(remote_bam, tmp_path, operation, offset):
    source, bam, state = remote_bam
    cache = Cache(tmp_path / "cache")
    with pytest.raises(IntegrityError, match="inventory lists.*new snapshot"):
        acquire(operation, replace(source, modified=EPOCH + offset), bam, cache)
    assert "HEAD" in state["requests"]
    assert "GET" not in state["requests"]
    assert not list(cache.workspace.rglob("receipt.json"))


@pytest.mark.parametrize("offset", [0, 0.5])
def test_matching_inventory_date_preserves_exact_regional_reads(remote_bam, tmp_path, offset):
    import pysam

    source, bam, state = remote_bam
    subset = acquire("extract", replace(source, modified=EPOCH + offset), bam, Cache(tmp_path / "cache"))
    assert subset.receipt["inventory_modification"] == "matched"
    assert subset.receipt["header_receipt"]["inventory_modification"] == "matched"
    assert subset.receipt["request"]["source_modified"] == EPOCH + offset
    with pysam.AlignmentFile(bam) as original, subset.open() as selected:
        expected = Counter(r.to_string() for r in original if r.reference_start < 160)
        assert Counter(r.to_string() for r in selected) == expected
    assert subset.receipt["records"] == 6
    assert "GET" in state["requests"]


@pytest.mark.parametrize("modified,last_modified,status", [
    (None, formatdate(EPOCH, usegmt=True), "inventory_timestamp_unavailable"),
    ("unknown", formatdate(EPOCH, usegmt=True), "inventory_timestamp_unavailable"),
    (EPOCH, None, "last_modified_unavailable"),
])
def test_missing_timestamp_evidence_is_explicit(remote_bam, tmp_path, modified, last_modified, status):
    source, bam, state = remote_bam
    state["last_modified"] = last_modified
    subset = acquire("extract", replace(source, modified=modified), bam, Cache(tmp_path / "cache"))
    assert subset.receipt["inventory_modification"] == status
    assert subset.receipt["header_receipt"]["inventory_modification"] == status
    assert subset.receipt["remote_identity"]["last-modified"] == last_modified


def test_invalid_server_date_is_an_integrity_error(remote_bam, tmp_path):
    source, bam, state = remote_bam
    state["last_modified"] = "not a date"
    cache = Cache(tmp_path / "cache")
    with pytest.raises(IntegrityError, match="Invalid Last-Modified"):
        acquire("extract", source, bam, cache)
    assert "GET" not in state["requests"]
    assert not list(cache.workspace.rglob("receipt.json"))


@pytest.mark.parametrize("operation", ["inspect", "extract", "count", "legacy_count"])
@pytest.mark.parametrize("legacy", [False, True])
def test_offline_cached_date_mismatch_is_rejected(remote_bam, tmp_path, monkeypatch, operation, legacy):
    import osteosarc.reads as reads

    source, bam, state = remote_bam
    cache = Cache(tmp_path / "cache")
    result = acquire("extract" if operation == "legacy_count" else operation, source, bam, cache)
    receipt = result.receipt
    directory = (cache.workspace / "derived" / stable_id(receipt["request"])
                 if operation == "count" else result.path.parent)
    # Model a result accepted by the old implementation under a mismatched inventory.
    receipt["request"]["source_modified"] += 86400
    if legacy:
        receipt.pop("inventory_modification")
    if operation == "inspect":
        write_json(directory / "receipt.json", receipt)
    else:
        write_read_receipt(directory / "receipt.json", receipt, cache.workspace)
    directory.rename(directory.parent / stable_id(receipt["request"]))
    monkeypatch.setattr(reads, "_remote_identity", lambda *a, **k: pytest.fail("offline HTTP access"))
    monkeypatch.setattr(reads, "_run", lambda *a, **k: pytest.fail("offline SAMtools execution"))
    with pytest.raises(IntegrityError, match="inventory lists.*new snapshot"):
        acquire("count" if operation == "legacy_count" else operation,
                replace(source, modified=EPOCH + 86400), bam, Cache(cache.root, offline=True))


def test_matching_inventory_count_is_reusable_offline_without_bams(remote_bam, tmp_path, monkeypatch):
    import osteosarc.reads as reads

    source, bam, state = remote_bam
    cache = Cache(tmp_path / "cache")
    count = acquire("count", source, bam, cache)
    assert count.records == 6
    assert count.receipt["inventory_modification"] == "matched"
    assert not list((cache.workspace / "derived").rglob("*.bam*"))
    monkeypatch.setattr(reads, "_remote_identity", lambda *a, **k: pytest.fail("offline HTTP access"))
    assert acquire("count", source, bam, Cache(cache.root, offline=True)) == count


def test_size_check_still_precedes_reading(remote_bam, tmp_path):
    source, bam, state = remote_bam
    with pytest.raises(IntegrityError, match="size differs"):
        acquire("inspect", replace(source, size=source.size + 1), bam, Cache(tmp_path / "cache"))
    assert "GET" not in state["requests"]


@pytest.mark.parametrize("modified", [float("nan"), float("inf"), 10 ** 1000])
def test_invalid_numeric_inventory_times_are_rejected(modified):
    with pytest.raises(IntegrityError, match="Invalid inventory modification timestamp"):
        _check_inventory_time(modified, formatdate(EPOCH, usegmt=True), source="source.bam")


def test_obsolete_http_date_is_interpreted_as_utc():
    assert _check_inventory_time(EPOCH, "Mon Jul 20 22:20:23 2026", source="source.bam") == "matched"
