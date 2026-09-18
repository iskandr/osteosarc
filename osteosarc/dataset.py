"""One consistent entry point for a reproducible snapshot of osteosarc data."""

from __future__ import annotations

import json
from dataclasses import asdict
from datetime import datetime, timezone
from functools import cached_property
from pathlib import Path

from .cache import Cache, Receipt, file_lock, stable_id, write_json
from .catalog import (
    SNAPSHOT_SOURCES,
    build_assets,
    parse_data_paths,
)
from .errors import IntegrityError, SchemaError
from .models import Asset
from .parsing import PARSE_FORMATS, Table, parse_file, parse_table, parse_variants, read_text


class Dataset:
    """A pinned set of metadata receipts plus lazily parsed public resources.

    Use Dataset.sync(name) explicitly to acquire metadata, then Dataset.open
    to reopen it offline. Full data objects are downloaded only on request.
    """

    def __init__(self, cache, manifest):
        self.cache = cache
        self.manifest = manifest
        if manifest.get("schema_version") != 1:
            raise SchemaError("Unsupported snapshot schema")
        if stable_id(manifest["sources"]) != manifest.get("id"):
            raise IntegrityError("Snapshot source receipts were modified")
        missing = set(SNAPSHOT_SOURCES) - set(manifest["sources"])
        if missing:
            raise SchemaError(f"Missing snapshot sources: {sorted(missing)}")
        for receipt in manifest["sources"].values():
            cache.path(receipt)

    @staticmethod
    def _snapshot_path(cache, name):
        if not name or Path(name).name != name or name in (".", ".."):
            raise ValueError("Snapshot name must be one path component")
        return cache.root / "snapshots" / (name + ".json")

    @classmethod
    def sync(cls, name, *, cache=None, refresh=False, sources=None):
        """Create a named snapshot, acquiring metadata only (~55 MB currently).

        Existing names cannot be overwritten. refresh=True fetches new source
        bytes for a NEW snapshot. sources may override endpoint URLs (for a
        pinned public source-repository commit or an internal mirror).
        """
        cache = cache if isinstance(cache, Cache) else Cache(cache)
        path = cls._snapshot_path(cache, name)
        with file_lock(path.with_suffix(".lock")):
            if path.exists():
                if refresh or sources:
                    raise FileExistsError("Choose a new snapshot name to change its sources")
                return cls.open(name, cache=cache, offline=cache.offline)
            urls = dict(SNAPSHOT_SOURCES)
            if sources:
                unknown = set(sources) - set(urls)
                if unknown:
                    raise KeyError(f"Unknown metadata sources: {sorted(unknown)}")
                urls.update(sources)
            receipts = {key: cache.fetch(url, refresh=refresh, max_bytes=256_000_000).to_dict()
                        for key, url in urls.items()}
            manifest = dict(schema_version=1, name=name, id=stable_id(receipts), sources=receipts,
                            created_at=datetime.now(timezone.utc).isoformat())
            dataset = cls(cache, manifest)
            # Validate identity joins before making a snapshot discoverable.
            dataset.variants()
            dataset.assets
            write_json(path, manifest)
            return dataset

    @classmethod
    def open(cls, name, *, cache=None, offline=True):
        """Open a saved snapshot; offline by default, including later downloads."""
        if isinstance(cache, Cache):
            cache = Cache(cache.root, offline=offline, timeout=cache.timeout)
        else:
            cache = Cache(cache, offline=offline)
        path = cls._snapshot_path(cache, name)
        return cls(cache, json.loads(path.read_text()))

    @property
    def id(self):
        return self.manifest["id"]

    def source_path(self, name):
        """Return verified bytes for a metadata source in this snapshot."""
        return self.cache.path(self.manifest["sources"][name])

    def _json(self, name):
        return json.loads(self.source_path(name).read_text())

    @cached_property
    def vafs(self):
        """Original count/annotation rows; no inference from missing counts."""
        return parse_table(read_text(self.source_path("vafs")),
                           source=self.manifest["sources"]["vafs"])

    @cached_property
    def assets(self):
        """All listed files plus metadata-only catalog objects and site tables."""
        return build_assets(self._json("bucket"), self._json("bams"),
                            parse_table(read_text(self.source_path("bam_metadata"))), self.vafs,
                            parse_data_paths(read_text(self.source_path("data_page"))))

    @cached_property
    def _variants(self):
        return parse_variants(read_text(self.source_path("variant_index")), self.vafs,
                              source_variants=self._json("source_variants"),
                              vaccine_overlap=self._json("vaccine_overlap"),
                              source={"snapshot_id": self.id, "receipts": self.manifest["sources"]})

    def variants(self, set="site", **filters):
        """Select 'site', 'all' (includes count-export entries), or 'vaccine'.

        Named vaccine membership defaults to the vaccine-overlap JSON; use
        vaccine_source='source_variants' for the separate source JSON flags.
        Pipeline detection comes from source JSON. The site table's vaccine
        count is kept independently. Unresolved
        alleles are retained unless status='ready' is requested explicitly.
        """
        if set not in ("site", "all", "vaccine"):
            raise ValueError("Variant set must be 'site', 'all', or 'vaccine'")
        variants = self._variants
        if set != "all":
            variants = variants.select(on_site=True)
        if set == "vaccine":
            variants = variants.select(vaccinated=True)
        return variants.select(**filters)

    @property
    def vaccines(self):
        """Published vaccine-overlap rows including unmodified ELISPOT states."""
        return Table(self._json("vaccine_overlap")["mutations"],
                     source=self.manifest["sources"]["vaccine_overlap"])

    @property
    def vaccine_names(self):
        return tuple(self._json("vaccine_overlap")["vaccine_names"])

    def vaccine_peptides(self, vaccine=None):
        """Published peptide sequences and experiments, keyed by exact variant ID."""
        rows = []
        for variant in self._variants:
            for peptide in variant.annotations["source_record"].get("vaccine_peptides", []):
                if vaccine is None or vaccine in peptide.get("in_vaccines", []):
                    rows.append(dict(peptide, variant_id=variant.id, gene=variant.gene))
        return Table(rows, source=self.manifest["sources"]["source_variants"])

    @property
    def annotations(self):
        """Original variant annotation records, including peptides and validation."""
        return Table(self._json("source_variants"), source=self.manifest["sources"]["source_variants"])

    @property
    def pipeline_names(self):
        return tuple(sorted({p for v in self._variants for p in
                             v.annotations["source_record"].get("detection", {})}))

    @property
    def samples(self):
        """Source-attributed sample/library claims with asset IDs.

        These rows intentionally do not turn processing products or pooled
        libraries into independent specimens. Unresolved identities stay visible.
        """
        groups = {}
        for asset in self.assets:
            for claim in asset.claims:
                key = stable_id(asdict(claim))
                groups.setdefault(key, dict(asdict(claim), id=key, asset_ids=[]))["asset_ids"].append(asset.id)
        return Table(groups.values(), source={"snapshot_id": self.id})

    @property
    def timepoints(self):
        """Published timepoint/date pairs; dates retain their source precision."""
        rows = {(c.timepoint, c.date, c.source) for a in self.assets for c in a.claims
                if c.basis == "published" and (c.timepoint or c.date)}
        return Table((dict(timepoint=t, date=d, source=s) for t, d, s in
                      sorted(rows, key=lambda r: tuple(x or "" for x in r))))

    def asset(self, key_or_id):
        """Resolve an exact bucket key, URL, resource name, or asset ID."""
        matches = [a for a in self.assets if key_or_id in
                   (a.id, a.key, a.url, a.metadata.get("resource"))]
        if len(matches) != 1:
            raise KeyError(f"Expected one asset for {key_or_id!r}, found {len(matches)}")
        return matches[0]

    def download(self, asset, *, refresh=False, verify_size=True):
        """Fetch one full object and bind its receipt to this metadata snapshot.

        After the first acquisition, URL refreshes elsewhere cannot change this
        snapshot's bytes. Use a new snapshot to acquire a newer object. The
        binding is populated lazily, so metadata-only use downloads no data.
        """
        asset = asset if isinstance(asset, Asset) else self.asset(asset)
        if asset.id != stable_id(asset.url):
            raise ValueError("Asset ID must be the stable hash of its URL")
        # Metadata tables already pinned by this snapshot never drift to latest.
        for receipt in self.manifest["sources"].values():
            if receipt["url"] == asset.url:
                if refresh:
                    raise ValueError("Create a new snapshot to refresh pinned metadata")
                return self.cache.path(receipt)
        binding = self.cache.root / "bindings" / self.id / (asset.id + ".json")
        with file_lock(binding.with_suffix(".lock")):
            if binding.exists():
                if refresh:
                    raise ValueError("Create a new snapshot to refresh an acquired object")
                receipt = Receipt(**json.loads(binding.read_text()))
                if receipt.url != asset.url:
                    raise IntegrityError("Snapshot object binding has the wrong URL")
                return self.cache.path(receipt)
            md5s = {r.get("md5sum") for r in asset.metadata.get("metadata_rows", []) if r.get("md5sum")}
            if len(md5s) > 1:
                raise IntegrityError("Conflicting published MD5 claims")
            receipt = self.cache.fetch(asset.url, refresh=refresh,
                                       size=asset.size if verify_size else None,
                                       md5=next(iter(md5s), None))
            write_json(binding, receipt.to_dict())
            return self.cache.path(receipt)

    def parse(self, asset):
        """Download explicitly selected data and parse it with original columns."""
        asset = asset if isinstance(asset, Asset) else self.asset(asset)
        if asset.format not in PARSE_FORMATS:
            raise ValueError(f"No built-in parser for {asset.format!r}; download the original asset")
        return parse_file(self.download(asset), format=asset.format,
                          source=dict(url=asset.url, snapshot_id=self.id))

    def table(self, asset):
        """Read any CSV/TSV asset, including named site tables and RSEM output."""
        asset = asset if isinstance(asset, Asset) else self.asset(asset)
        if asset.format not in ("csv", "tsv"):
            raise ValueError("table requires a CSV or TSV asset")
        return self.parse(asset)

    def inspect_alignment(self, asset, **kwargs):
        """Cache the alignment header to inspect assembly before choosing regions."""
        from .reads import inspect_alignment
        asset = asset if isinstance(asset, Asset) else self.asset(asset)
        return inspect_alignment(asset, cache=self.cache, snapshot_id=self.id, **kwargs)

    def extract_reads(self, asset, regions, **kwargs):
        """Extract an indexed region union; see osteosarc.extract_reads."""
        from .reads import extract_reads
        asset = asset if isinstance(asset, Asset) else self.asset(asset)
        return extract_reads(asset, regions, cache=self.cache, snapshot_id=self.id, **kwargs)

    def open_variants(self, asset):
        """Download one VCF/BCF and its listed index, returning pysam.VariantFile.

        Use as a context manager. Header, caller annotations, multiallelic
        records, symbolic alleles and per-sample FORMAT fields are preserved.
        Indexed fetch is available when an index is published; otherwise iterate.
        """
        import pysam
        asset = asset if isinstance(asset, Asset) else self.asset(asset)
        if asset.format not in ("vcf", "bcf"):
            raise ValueError("open_variants requires a VCF or BCF")
        path = self.download(asset)
        index = None
        if asset.index_urls:
            index = self.download(asset.index_urls[0])
        return pysam.VariantFile(str(path), index_filename=str(index) if index else None)

    def receipts(self):
        return {name: Receipt(**value) for name, value in self.manifest["sources"].items()}
