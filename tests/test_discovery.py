from urllib.parse import parse_qs, urlsplit

import pytest

from osteosarc import Cache, SchemaError
from osteosarc.discovery import list_bucket


def page(prefix, key, *, token=None, truncated=False):
    return ('<ListBucketResult xmlns="http://s3.amazonaws.com/doc/2006-03-01/">'
            f'<Prefix>{prefix}</Prefix><IsTruncated>{str(truncated).lower()}</IsTruncated>'
            f'<Contents><Key>{key}</Key><Size>42</Size><LastModified>date</LastModified></Contents>'
            + (f'<NextContinuationToken>{token}</NextContinuationToken>' if token else '')
            + '</ListBucketResult>')


def test_pagination_and_tokens_are_complete(tmp_path, monkeypatch):
    cache = Cache(tmp_path / "cache", offline=True)
    calls = []

    def fetch(url, **kwargs):
        query = parse_qs(urlsplit(url).query)
        calls.append(query)
        first = "continuation-token" not in query
        path = tmp_path / "page.xml"
        path.write_text(page("a/", "a/one" if first else "a/two", token="t+/=" if first else None, truncated=first))
        return cache.import_file(path, url)

    monkeypatch.setattr(cache, "fetch", fetch)
    listing = list_bucket(cache, "a/")
    assert [r[0] for r in listing["files"]] == ["a/one", "a/two"]
    assert calls[1]["continuation-token"] == ["t+/="]
    assert len(listing["receipts"]) == 2


@pytest.mark.parametrize("payload", [page("wrong/", "wrong/file"), page("a/", "a/file", truncated=True),
                                      page("a/", "outside/file")])
def test_invalid_listings_do_not_return_partial_success(tmp_path, monkeypatch, payload):
    cache = Cache(tmp_path / "cache")
    path = tmp_path / "page.xml"
    path.write_text(payload)
    monkeypatch.setattr(cache, "fetch", lambda url, **kw: cache.import_file(path, url))
    with pytest.raises(SchemaError):
        list_bucket(cache, "a/")
