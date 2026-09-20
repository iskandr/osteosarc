import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import datacache
import pytest
import requests

from osteosarc import Cache, IntegrityError, OfflineError, OsteosarcError
from osteosarc.cache import stable_id


@pytest.fixture
def download_transport(monkeypatch):
    """Run the real datacache downloader over an in-memory HTTP response."""
    state = dict(body=b"original bytes", calls=0, fail=False, etag='"opaque-multipart-2"')

    class Response:
        def __init__(self, body=b"", fail=False):
            self.body, self.fail = body, fail
            self.headers = requests.structures.CaseInsensitiveDict({
                "ETag": state["etag"], "Last-Modified": "yesterday",
                "Content-Length": str(len(state["body"]))})

        def raise_for_status(self):
            if self.fail:
                response = requests.Response()
                response.status_code = 404
                raise requests.HTTPError("failed", response=response)

        def iter_content(self, chunk_size):
            for i in range(0, len(self.body), chunk_size):
                yield self.body[i:i + chunk_size]

        def close(self):
            pass

        def __enter__(self):
            return self

        def __exit__(self, *args):
            self.close()

    def get(url, **kwargs):
        state["calls"] += 1
        return Response(state["body"], state["fail"])

    monkeypatch.setattr("requests.get", get)
    monkeypatch.setattr("requests.head", lambda *a, **k: Response())
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
    with pytest.raises(OsteosarcError, match="Cannot download"):
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


def test_datacache_objects_are_reused_in_both_directions(tmp_path, download_transport):
    import hashlib

    cache = Cache(tmp_path)
    backend = datacache.Cache(cache_root=cache.objects)
    url = "https://example.test/shared.bam.bai"
    sha = hashlib.sha256(download_transport["body"]).hexdigest()
    path = Path(backend.fetch(url, filename=sha + ".bam.bai", expected_sha256=sha))
    inode = path.stat().st_ino
    receipt = Cache(tmp_path, offline=True).fetch(url, sha256=sha)
    assert cache.path(receipt) == path
    assert download_transport["calls"] == 1
    assert path.stat().st_ino == inode

    other_url = "https://example.test/other.tsv"
    receipt = cache.fetch(other_url)
    shared = cache.path(receipt)
    inode = shared.stat().st_ino
    assert Path(backend.fetch(other_url, filename=shared.name,
                              expected_sha256=receipt.sha256)) == shared
    assert shared.stat().st_ino == inode
    assert download_transport["calls"] == 2


def test_compressed_downloads_keep_the_original_bytes(tmp_path, download_transport):
    import gzip

    download_transport["body"] = gzip.compress(b"##fileformat=VCFv4.2\n")
    cache = Cache(tmp_path)
    receipt = cache.fetch("https://example.test/calls.vcf.gz?revision=1")
    assert cache.path(receipt).read_bytes() == download_transport["body"]
    assert cache.path(receipt).name.endswith(".vcf.gz")


def test_changed_http_identity_does_not_replace_existing_receipt(tmp_path, download_transport, monkeypatch):
    cache = Cache(tmp_path)
    url = "https://example.test/data.tsv"
    old = cache.fetch(url)
    original_get = requests.get

    def changed_get(*args, **kwargs):
        response = original_get(*args, **kwargs)
        download_transport["etag"] = '"changed-during-transfer"'
        return response

    monkeypatch.setattr(requests, "get", changed_get)
    with pytest.raises(IntegrityError, match="changed during download"):
        cache.fetch(url, refresh=True)
    assert Cache(tmp_path, offline=True).fetch(url) == old
    assert not list((cache.workspace / "staging").iterdir())


def test_oversized_stream_aborts_and_cleans_datacache_staging(tmp_path, download_transport):
    download_transport["body"] = b"x" * (3 * 1024 * 1024)
    cache = Cache(tmp_path)
    with pytest.raises(IntegrityError, match="exceeds"):
        cache.fetch("https://example.test/large.bin", max_bytes=10)
    assert not list((cache.workspace / "staging").iterdir())
    assert not list(cache.objects.glob("*"))


def test_download_works_when_server_does_not_support_head(tmp_path, download_transport, monkeypatch):
    response = requests.Response()
    response.status_code = 405
    response._content_consumed = True
    monkeypatch.setattr(requests, "head", lambda *a, **k: response)
    cache = Cache(tmp_path)
    receipt = cache.fetch("https://example.test/get-only.tsv")
    assert receipt.etag is None and receipt.last_modified is None
    assert cache.path(receipt).read_bytes() == download_transport["body"]


def test_transient_failure_is_retried_by_datacache(tmp_path, download_transport, monkeypatch):
    original_get = requests.get
    attempts = []

    def get(*args, **kwargs):
        attempts.append(args)
        if len(attempts) == 1:
            response = requests.Response()
            response.status_code = 503
            raise requests.HTTPError("temporarily unavailable", response=response)
        return original_get(*args, **kwargs)

    monkeypatch.setattr(requests, "get", get)
    monkeypatch.setattr("datacache.download.time.sleep", lambda *a: None)
    cache = Cache(tmp_path)
    receipt = cache.fetch("https://example.test/retried.tsv")
    assert len(attempts) == 2
    assert cache.path(receipt).read_bytes() == download_transport["body"]
