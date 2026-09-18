"""Explicit, paginated discovery of the public S3 bucket (no credentials)."""

from urllib.parse import urlencode
from xml.etree import ElementTree

from .catalog import BUCKET
from .errors import SchemaError

NS = {"s": "http://s3.amazonaws.com/doc/2006-03-01/"}


def list_bucket(cache, prefix, *, refresh=False, max_pages=1000):
    """Return a complete prefix listing and receipts, or raise on truncation.

    Website metadata is dated. This separately records live S3 object claims;
    it never downloads the objects themselves. Existing pages are reused unless
    refresh=True. Use a new dataset snapshot when incorporating newer listings.
    """
    files, receipts, seen, token = [], [], set(), None
    for _ in range(max_pages):
        query = {"list-type": "2", "prefix": prefix, "max-keys": "1000"}
        if token:
            query["continuation-token"] = token
        receipt = cache.fetch(BUCKET + "?" + urlencode(query), refresh=refresh, max_bytes=8_000_000)
        receipts.append(receipt.to_dict())
        root = ElementTree.fromstring(cache.path(receipt).read_bytes())
        if root.findtext("s:Prefix", default="", namespaces=NS) != prefix:
            raise SchemaError("S3 returned a different prefix")
        for row in root.findall("s:Contents", NS):
            key = row.findtext("s:Key", namespaces=NS)
            if key is None or not key.startswith(prefix):
                raise SchemaError("S3 returned a key outside the requested prefix")
            files.append([key, int(row.findtext("s:Size", namespaces=NS)),
                          row.findtext("s:LastModified", namespaces=NS)])
        truncated = root.findtext("s:IsTruncated", namespaces=NS)
        if truncated == "false":
            if len({r[0] for r in files}) != len(files):
                raise SchemaError("S3 listing contains duplicate keys")
            return dict(download_base=BUCKET, files=files, prefix=prefix, receipts=receipts)
        token = root.findtext("s:NextContinuationToken", namespaces=NS)
        if truncated != "true" or not token or token in seen:
            raise SchemaError("Missing or repeated S3 continuation token")
        seen.add(token)
    raise SchemaError(f"S3 listing exceeded {max_pages} pages; no complete inventory returned")
