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
        assert not (tmp_path / "urls" / (stable_id(url) + ".json")).exists()
    download_transport["fail"] = True
    with pytest.raises(subprocess.CalledProcessError):
        cache.fetch(url)
    assert not list((tmp_path / "staging").iterdir())


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
    pointer = tmp_path / "urls" / (stable_id(first_url) + ".json")
    record = json.loads(pointer.read_text())
    record["filename"] = "../escape"
    pointer.write_text(json.dumps(record))
    with pytest.raises(IntegrityError):
        cache.fetch(first_url)
