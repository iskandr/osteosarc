import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pytest

from osteosarc import Cache, IntegrityError, OfflineError
from osteosarc.cache import stable_id


@pytest.fixture
def download_transport(monkeypatch):
    """Exercise curl's file/receipt boundary without network access in CI."""
    state = dict(body=b"original bytes", calls=0, fail=False)

    def run(command, **kwargs):
        state["calls"] += 1
        Path(command[command.index("--output") + 1]).write_bytes(state["body"])
        Path(command[command.index("--dump-header") + 1]).write_text(
            'HTTP/1.1 302 Found\r\nETag: "redirect"\r\n\r\n'
            'HTTP/2 200\r\nETag: "opaque-multipart-2"\r\nLast-Modified: yesterday\r\n\r\n')
        if state["fail"]:
            raise subprocess.CalledProcessError(22, command, stderr=b"failed")
        return subprocess.CompletedProcess(command, 0)

    monkeypatch.setattr("osteosarc.cache.subprocess.run", run)
    return state


def test_refresh_preserves_old_content_and_offline_verifies(tmp_path, download_transport):
    cache = Cache(tmp_path)
    url = "https://example.test/a/data.tsv"
    first = cache.fetch(url)
    assert first.etag == '"opaque-multipart-2"'
    assert first.md5 is None
    assert cache.fetch(url) == first
    assert download_transport["calls"] == 1
    download_transport["body"] = b"new bytes"
    second = cache.fetch(url, refresh=True)
    assert first.sha256 != second.sha256
    assert cache.path(first).read_bytes() == b"original bytes"
    offline = Cache(tmp_path, offline=True)
    assert offline.fetch(url) == second
    with pytest.raises(OfflineError):
        offline.fetch(url, refresh=True)
    with pytest.raises(OfflineError):
        offline.fetch("https://example.test/missing")
    cache.path(second).write_bytes(b"modified")
    with pytest.raises(IntegrityError):
        offline.fetch(url)


def test_failed_and_mismatched_downloads_never_publish(tmp_path, download_transport):
    cache = Cache(tmp_path)
    url = "https://example.test/file"
    for kwargs in ({"sha256": "0" * 64}, {"md5": "0" * 32}, {"size": 9}, {"max_bytes": 1}):
        with pytest.raises(IntegrityError):
            cache.fetch(url, **kwargs)
        assert not (cache.workspace / "urls" / (stable_id(url) + ".json")).exists()
    download_transport["fail"] = True
    with pytest.raises(subprocess.CalledProcessError):
        cache.fetch(url)
    assert not list((cache.workspace / "staging").iterdir())


def test_same_basename_and_concurrent_writers(tmp_path, download_transport):
    cache = Cache(tmp_path)
    first_url = "https://example.test/a/data.tsv"
    with ThreadPoolExecutor(max_workers=4) as executor:
        receipts = list(executor.map(cache.fetch, [first_url] * 4))
    assert len({r.sha256 for r in receipts}) == 1
    assert download_transport["calls"] == 1
    download_transport["body"] = b"a different file"
    second = cache.fetch("https://example.test/b/data.tsv")
    assert cache.path(receipts[0]) != cache.path(second)
    pointer = cache.workspace / "urls" / (stable_id(first_url) + ".json")
    record = json.loads(pointer.read_text())
    record["filename"] = "../escape"
    pointer.write_text(json.dumps(record))
    with pytest.raises(IntegrityError):
        cache.fetch(first_url)


def test_pruned_object_is_fetched_again_online_but_not_offline(tmp_path, download_transport):
    cache = Cache(tmp_path)
    url = "https://example.test/a/data.tsv"
    first = cache.fetch(url)
    path = cache.path(first)
    path.unlink()
    with pytest.raises(OfflineError):
        Cache(tmp_path, offline=True).fetch(url)
    restored = cache.fetch(url, sha256=first.sha256)
    assert restored.sha256 == first.sha256 and cache.path(restored) == path
    assert path.read_bytes() == b"original bytes" and download_transport["calls"] == 2


def test_files_are_shareable_and_local_digests_are_remembered(tmp_path, download_transport, monkeypatch):
    import os
    import stat

    from osteosarc.cache import write_json
    old = os.umask(0o022)
    try:
        write_json(tmp_path / "x.json", {})
        receipt = Cache(tmp_path).fetch("https://example.test/shared.tsv")
        assert stat.S_IMODE((tmp_path / "x.json").stat().st_mode) == 0o644
        assert stat.S_IMODE(Cache(tmp_path).path(receipt).stat().st_mode) == 0o644
    finally:
        os.umask(old)
    local = tmp_path / "local.bam"
    local.write_bytes(b"alignment")
    cache, calls = Cache(tmp_path / "cache"), []
    import osteosarc.cache as module
    original = module.digest
    monkeypatch.setattr(module, "digest", lambda *a, **k: calls.append(a) or original(*a, **k))
    first = cache.file_digest(local)
    assert Cache(tmp_path / "cache").file_digest(local) == first and len(calls) == 1
    local.write_bytes(b"rewritten alignment")
    assert cache.file_digest(local) != first and len(calls) == 2


def test_objects_use_the_shared_openvax_layout(tmp_path, download_transport, monkeypatch):
    from osteosarc.cache import default_root
    cache = Cache(tmp_path)
    receipt = cache.fetch("https://example.test/dir/sample.genes.results")
    path = cache.path(receipt)
    # The same path vaxrank's downloader uses: objects/sha256/<sha><suffixes>.
    assert path == tmp_path / "objects" / "sha256" / (receipt.sha256 + ".genes.results")
    assert receipt.filename == "sample.genes.results"
    assert not (tmp_path / "urls").exists() and (tmp_path / "osteosarc" / "urls").is_dir()
    monkeypatch.delenv("OSTEOSARC_CACHE", raising=False)
    monkeypatch.setenv("OPENVAX_DATA_CACHE", str(tmp_path / "openvax"))
    assert Cache().root == (tmp_path / "openvax").resolve()
    monkeypatch.setenv("OSTEOSARC_CACHE", str(tmp_path / "isolated"))
    assert Cache().root == (tmp_path / "isolated").resolve()
    monkeypatch.delenv("OPENVAX_DATA_CACHE")
    assert default_root().name == "openvax"


def test_objects_written_by_another_openvax_tool_are_reused(tmp_path, download_transport):
    import hashlib
    body = b"shared bytes"
    sha = hashlib.sha256(body).hexdigest()
    shared = tmp_path / "objects" / "sha256" / (sha + ".bam.bai")
    shared.parent.mkdir(parents=True)
    shared.write_bytes(body)  # e.g. vaxrank's datacache published this object
    source = tmp_path / "elsewhere.bam.bai"
    source.write_bytes(body)
    inode = shared.stat().st_ino
    receipt = Cache(tmp_path).import_file(source, "https://example.test/x.bam.bai", sha256=sha)
    assert Cache(tmp_path).path(receipt) == shared and shared.read_bytes() == body
    assert shared.stat().st_ino == inode  # reused, not rewritten
