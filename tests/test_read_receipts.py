import json

import pytest

from osteosarc import Cache, IntegrityError, ReadFilter, Region, extract_reads
from osteosarc.cache import digest, stable_id, write_json
from osteosarc.read_receipts import (
    compact_read_cache,
    query_names_asset,
    read_read_receipt,
    read_receipt_files,
    write_read_receipt,
)


def test_nested_receipts_share_one_asset_without_changing_requests(tmp_path):
    names = [f"read-{i}" for i in range(1000)]
    request = dict(filters=dict(query_names=names), regions=["chr1:1-20"])
    receipt = dict(request=request, acquisition=[dict(request=request) for _ in range(12)])
    before = stable_id(receipt)
    path = tmp_path / "receipt.json"
    write_read_receipt(path, receipt, tmp_path)
    assert stable_id(receipt) == before
    assert read_read_receipt(path, tmp_path) == receipt
    assert len(list((tmp_path / "query-names").glob("*.txt"))) == 1
    assert path.stat().st_size < len(json.dumps(receipt)) / 10
    files = read_receipt_files(path, tmp_path)
    assert len(files) == 2 and files[path] == digest(path)
    other = tmp_path / "second.json"
    write_read_receipt(other, dict(request=request), tmp_path)
    assert set(read_receipt_files(other, tmp_path)) - {other} == set(files) - {path}
    # Expanded lists remain independent mutable values, as in old JSON receipts.
    expanded = read_read_receipt(path, tmp_path)
    expanded["acquisition"][0]["request"]["filters"]["query_names"].append("new")
    assert expanded["request"]["filters"]["query_names"] == names


@pytest.mark.parametrize("damage", ["missing", "changed", "count", "unresolved", "traversal"])
def test_asset_corruption_is_rejected_before_reusing_receipt(tmp_path, damage):
    path = tmp_path / "receipt.json"
    write_read_receipt(path, dict(filters=dict(query_names=["a", "b"])), tmp_path)
    saved = json.loads(path.read_text())
    checksum = next(iter(saved["query_name_assets"]))
    asset = tmp_path / "query-names" / (checksum + ".txt")
    if damage == "missing":
        asset.unlink()
    elif damage == "changed":
        asset.write_text("c\nb\n")
    else:
        if damage == "count":
            saved["query_name_assets"][checksum] = 3
        elif damage == "unresolved":
            saved["receipt"]["filters"]["query_names"]["sha256"] = "0" * 64
        else:
            saved["query_name_assets"] = {"../outside": 2}
        write_json(path, saved)
    with pytest.raises(IntegrityError):
        read_read_receipt(path, tmp_path)
    with pytest.raises(IntegrityError):
        read_receipt_files(path, tmp_path)


def test_two_extractions_share_asset_and_legacy_cache_compacts_offline(bam, tmp_path, monkeypatch):
    import osteosarc.reads as reads
    cache = Cache(tmp_path / "cache", offline=True)
    filters = ReadFilter(query_names=["primary", "absent"])
    first = extract_reads(bam, [Region("chr1", 100, 160, "GRCh38")], cache=cache, filters=filters)
    second = extract_reads(bam, [Region("chr1", 100, 170, "GRCh38")], cache=cache, filters=filters)
    assert first.path != second.path
    asset, = (cache.workspace / "query-names").glob("*.txt")
    for subset in [first, second]:
        assert not (subset.path.parent / "query-names.txt").exists()
        assert subset.receipt["request"]["filters"]["query_names"] == list(filters.query_names)
        assert subset.receipt["command"][subset.receipt["command"].index("-N") + 1].endswith(asset.name)
    # Simulate old expanded disk receipts and duplicate -N files.
    pins = {p: digest(p) for sub in [first, second] for p in [sub.path, sub.index_path]}
    for subset in [first, second]:
        write_json(subset.receipt_path, subset.receipt)
        (subset.path.parent / "query-names.txt").write_bytes(asset.read_bytes())
    monkeypatch.setattr(reads, "_run", lambda *a: pytest.fail("offline reuse must not run samtools"))
    assert extract_reads(bam, [Region("chr1", 100, 160, "GRCh38")], cache=cache, filters=filters).receipt == first.receipt
    assert compact_read_cache(cache) == 2
    for subset in [first, second]:
        assert read_read_receipt(subset.receipt_path, cache.workspace) == subset.receipt
        assert (subset.path.parent / "query-names.txt").stat().st_ino == asset.stat().st_ino
    assert {p: digest(p) for p in pins} == pins
    assert extract_reads(bam, [Region("chr1", 100, 160, "GRCh38")], cache=cache, filters=filters).receipt == first.receipt
    asset.write_text("corruption\n")
    with pytest.raises(IntegrityError):
        extract_reads(bam, [Region("chr1", 100, 160, "GRCh38")], cache=cache, filters=filters)


def test_failed_asset_publication_leaves_no_partial_asset(tmp_path, monkeypatch):
    import osteosarc.read_receipts as module

    def interrupted(*args):
        raise OSError("interrupted")

    monkeypatch.setattr(module.os, "replace", interrupted)
    with pytest.raises(OSError, match="interrupted"):
        query_names_asset(["a"], tmp_path)
    assert not list((tmp_path / "query-names").iterdir())
