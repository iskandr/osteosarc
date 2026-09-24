"""One consistent entry point for a reproducible snapshot of osteosarc data."""

from __future__ import annotations

import copy
import json
import re
import warnings
from collections import Counter, defaultdict
from dataclasses import asdict
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from functools import cached_property
from pathlib import Path
from urllib.parse import urlsplit

from .cache import Cache, Receipt, file_lock, stable_id, write_json
from .catalog import (
    BUCKET,
    SNAPSHOT_SOURCES,
    TABLE_SOURCES,
    TIMELINE_SOURCES,
    build_assets,
    data_page_rows,
    object_key,
    parse_data_paths,
)
from .curation import CORRECTIONS, Curation, normalize_tissue, unrecognized_values
from .errors import CoordinateError, IntegrityError, OfflineError, SchemaError
from .models import Asset, Region
from .parsing import (
    PARSE_FORMATS,
    Table,
    parse_file,
    parse_table,
    parse_variant_index,
    parse_variants,
    read_text,
)


def _check_inventory_time(asset, receipt):
    """Refuse bytes newer (or older) than the object this snapshot's inventory listed."""
    if not isinstance(asset.modified, (int, float)) or not receipt.last_modified:
        return
    served = parsedate_to_datetime(receipt.last_modified).timestamp()
    if abs(served - asset.modified) > 1:
        raise IntegrityError(
            f"{asset.key} was modified at {receipt.last_modified}, not at the time this snapshot's "
            f"inventory lists ({datetime.fromtimestamp(asset.modified, timezone.utc).isoformat()}); "
            "create a new snapshot to use the current object")


DATE_SELECTOR = re.compile(r"\d{4}(-\d{2}(-\d{2})?)?")
DATED_NAME = re.compile(r"\d{4}-\d{2}-\d{2}(\.\d+)?")


def downloaded_at(manifest):
    """UTC time of a snapshot's latest source download (its creation time if it has none)."""
    return max((r["retrieved_at"] for r in manifest["sources"].values()), default=manifest.get("created_at"))


def choose_snapshot(rows, name=None, *, date=None):
    """Name of the snapshot to open, given rows sorted newest first.

    With neither argument, the newest. name is an exact name or, if no snapshot
    has that name, a unique ID prefix of six or more hexadecimal characters.
    date is a UTC year, month or day and selects the newest snapshot
    downloaded then.
    """
    rows = list(rows)
    if name is not None and date is not None:
        raise ValueError("Choose a snapshot by name or by date, not both")
    if date is not None and not DATE_SELECTOR.fullmatch(date):
        raise ValueError("date must be YYYY, YYYY-MM or YYYY-MM-DD")
    if not rows:
        raise FileNotFoundError("No saved snapshots; run Dataset.sync() or `osteosarc sync` first")
    if date is not None:
        matches, wanted = [r for r in rows if (r["downloaded"] or "").startswith(date)], f"downloaded in {date} (UTC)"
    elif name is not None:
        matches, wanted = [r for r in rows if r["name"] == name], f"named {name!r}"
        if not matches and re.fullmatch(r"[0-9a-fA-F]{6,64}", name):
            matches = [r for r in rows if r["id"].startswith(name.lower())]
            if len({r["id"] for r in matches}) > 1:
                raise FileNotFoundError(f"Snapshot ID prefix {name!r} is ambiguous")
    else:
        return rows[0]["name"]
    if not matches:
        raise FileNotFoundError(f"No snapshot {wanted}; saved: " + ", ".join(r["name"] for r in rows))
    return matches[0]["name"]


class Dataset:
    """A pinned set of metadata receipts plus lazily parsed public resources.

    Use Dataset.sync() explicitly to acquire metadata, then Dataset.open()
    to reopen the most recent snapshot offline. Full data objects are
    downloaded only on request.

    corrections=True applies the verified corrections in osteosarc.curation;
    False uses the published sources unchanged; a sequence supplies your own.
    Either way, data.corrections reports every correction's status.
    """

    def __init__(self, cache, manifest, *, corrections=True):
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
        chosen = CORRECTIONS if corrections in (True, False) else tuple(corrections)
        self.curation = Curation(chosen, self._published, enabled=corrections is not False)
        self._columns, self._table_diagnostics, self._bucket_header = {}, {}, None

    @staticmethod
    def _snapshot_path(cache, name):
        if not name or Path(name).name != name or name in (".", ".."):
            raise ValueError("Snapshot name must be one path component")
        return cache.workspace / "snapshots" / (name + ".json")

    @classmethod
    def snapshots(cls, *, cache=None):
        """Saved snapshots, most recently downloaded first.

        ``downloaded`` is the UTC time of the snapshot's latest source download.
        It can precede ``created`` when a snapshot reused cached sources.
        Unreadable manifests are skipped with a warning.
        """
        cache = cache if isinstance(cache, Cache) else Cache(cache, offline=True)
        rows = []
        for path in sorted((cache.workspace / "snapshots").glob("*.json")):
            try:
                manifest = json.loads(path.read_text())
                rows.append(dict(name=path.stem, downloaded=downloaded_at(manifest),
                                 created=manifest.get("created_at"), id=manifest["id"]))
            except (OSError, ValueError, KeyError, TypeError, AttributeError) as error:
                warnings.warn(f"Skipping unreadable snapshot {path.name}: {error!r}", stacklevel=2)
        rows.sort(key=lambda r: (r["downloaded"] or "", r["created"] or "", r["name"]), reverse=True)
        return Table(rows, columns=["name", "downloaded", "created", "id"])

    @classmethod
    def sync(cls, name=None, *, cache=None, refresh=False, sources=None, corrections=True):
        """Save the website's current metadata as a snapshot (~57 MB currently).

        Without a name, the snapshot is named by its UTC download date, such as
        2026-09-24, and calling sync() again that day reopens it; refresh=True
        downloads another (2026-09-24.2). Offline, the snapshot is built from
        source bytes already in the cache and dated by their download.

        A named snapshot is created once: repeating the name reopens it, and
        refresh=True fetches new source bytes for a new name. sources overrides
        endpoint URLs (a pinned source-repository commit or a mirror) and needs
        a name, so it never becomes the dated default.
        """
        cache = cache if isinstance(cache, Cache) else Cache(cache)
        if name is not None:
            return cls._create(name, cache, refresh=refresh, sources=sources, corrections=corrections)
        if sources:
            raise ValueError("Name a snapshot with custom sources, as in Dataset.sync('rev-abc123', sources=...)")
        if refresh and cache.offline:
            raise OfflineError("Cannot download new metadata in offline mode")
        with file_lock(cache.workspace / "snapshots" / ".dated.lock"):
            dated = [r for r in cls.snapshots(cache=cache) if DATED_NAME.fullmatch(r["name"])]
            today = datetime.now(timezone.utc).date().isoformat()
            current = [r for r in dated if (r["downloaded"] or "").startswith(today)]
            if current and not refresh:
                return cls.open(current[0]["name"], cache=cache, offline=cache.offline, corrections=corrections)
            receipts = cls._receipts(cache, refresh=not cache.offline)
            same = [r for r in dated if r["id"] == stable_id(receipts)]
            if same and not refresh:
                return cls.open(same[0]["name"], cache=cache, offline=cache.offline, corrections=corrections)
            # Name by the actual download date, which can differ from today's if
            # the download crossed midnight or came from the cache offline.
            day = max(r["retrieved_at"] for r in receipts.values())[:10]
            taken = {path.stem for path in (cache.workspace / "snapshots").glob("*.json")}
            name = day if day not in taken else next(
                f"{day}.{n}" for n in range(2, len(taken) + 3) if f"{day}.{n}" not in taken)
            path = cls._snapshot_path(cache, name)
            with file_lock(path.with_suffix(".lock")):
                return cls._save(cache, path, name, receipts, corrections)

    @classmethod
    def _receipts(cls, cache, *, refresh, sources=None):
        urls = dict(SNAPSHOT_SOURCES, **TIMELINE_SOURCES)
        if sources:
            unknown = set(sources) - set(urls)
            if unknown:
                raise KeyError(f"Unknown metadata sources: {sorted(unknown)}")
            urls.update(sources)
        return {key: cache.fetch(url, refresh=refresh, max_bytes=256_000_000).to_dict()
                for key, url in urls.items()}

    @classmethod
    def _save(cls, cache, path, name, receipts, corrections):
        manifest = dict(schema_version=1, name=name, id=stable_id(receipts), sources=receipts,
                        created_at=datetime.now(timezone.utc).isoformat())
        dataset = cls(cache, manifest, corrections=corrections)
        # Validate identity joins before making a snapshot discoverable.
        dataset.variants()
        dataset.assets
        write_json(path, manifest)
        return dataset

    @classmethod
    def _create(cls, name, cache, *, refresh, sources, corrections):
        path = cls._snapshot_path(cache, name)
        with file_lock(path.with_suffix(".lock")):
            if path.exists():
                if refresh or sources:
                    raise FileExistsError("Choose a new snapshot name to change its sources")
                return cls.open(name, cache=cache, offline=cache.offline, corrections=corrections)
            receipts = cls._receipts(cache, refresh=refresh, sources=sources)
            return cls._save(cache, path, name, receipts, corrections)

    @classmethod
    def open(cls, name=None, *, date=None, cache=None, offline=True, corrections=True):
        """Open a saved snapshot, the most recently downloaded by default.

        name is an exact snapshot name or a snapshot ID prefix of six or more
        characters. date (2026, 2026-09 or 2026-09-24) opens the newest snapshot
        downloaded then, in UTC. Offline by default, including later downloads.
        """
        if isinstance(cache, Cache):
            cache = Cache(cache.root, offline=offline, timeout=cache.timeout)
        else:
            cache = Cache(cache, offline=offline)
        if name is None or date is not None or not cls._snapshot_path(cache, name).exists():
            name = choose_snapshot(cls.snapshots(cache=cache), name, date=date)
        path = cls._snapshot_path(cache, name)
        return cls(cache, json.loads(path.read_text()), corrections=corrections)

    @property
    def name(self):
        return self.manifest["name"]

    @property
    def downloaded(self):
        """UTC time of this snapshot's latest source download."""
        return downloaded_at(self.manifest)

    @property
    def id(self):
        return self.manifest["id"]

    def source_path(self, name):
        """Return verified bytes for a metadata source in this snapshot."""
        return self.cache.path(self.manifest["sources"][name])

    def _json(self, name):
        return json.loads(self.source_path(name).read_text())

    def _published(self, name):
        """Published records of one source as a list of dicts (for curation)."""
        if name in ("vafs", "bam_metadata"):
            table = parse_table(read_text(self.source_path(name)), strict=name != "vafs")
            self._columns[name] = table.columns
            self._table_diagnostics[name] = table.diagnostics
            return list(table.rows)
        if name == "variant_index":
            return parse_variant_index(read_text(self.source_path(name)))
        if name == "vaccine_overlap":
            return self._json(name)["mutations"]
        if name == "bams":
            return [row for category in self._json(name)["categories"] for row in category["bams"]]
        if name == "bucket":
            listing = self._json(name)
            self._bucket_header = {k: v for k, v in listing.items() if k != "files"}
            records = []
            for row in listing["files"]:
                try:
                    records.append(dict(key=row[0], size=int(row[1]), modified=row[2]))
                except (IndexError, TypeError, ValueError) as error:
                    raise SchemaError(f"Malformed bucket object: {row!r}") from error
            return records
        if name == "source_variants":
            return self._json(name)
        if name not in TIMELINE_SOURCES:
            raise KeyError(f"No correctable source named {name!r}")
        if name not in self.manifest["sources"]:
            return None  # an older snapshot: corrections for this source are unavailable
        if name in ("specimens", "fastqs", "labs", "cytometry"):
            return list(parse_table(read_text(self.source_path(name))).rows)
        if name == "events_sheet":
            return list(parse_table(read_text(self.source_path(name)), delimiter=",").rows)
        field = {"events": "events", "mrd": "measurements", "timepoint_summary": "rows",
                 "flow": "samples", "imaging": "studies", "pathology": "groups"}[name]
        return self._json(name)[field]

    def _require_timeline(self):
        missing = sorted(set(TIMELINE_SOURCES) - set(self.manifest["sources"]))
        if missing:
            raise SchemaError(f"Snapshot {self.manifest['name']!r} predates the timeline sources "
                              f"({', '.join(missing)}); create a new one with Dataset.sync(new_name)")

    @cached_property
    def timeline(self):
        """Every dated event: treatments, procedures, imaging, pathology, omics,
        specimens, MRD, flow draws, and lab/cytometry draw dates."""
        from .timeline import (
            Timeline,
            events_from_flow,
            events_from_imaging,
            events_from_measurements,
            events_from_mrd,
            events_from_pathology,
            events_from_site,
            events_from_specimens,
        )
        self._require_timeline()
        records, undated = self.curation.records, []
        events, event_marks = records("events")
        mrd, mrd_marks = records("mrd")
        sheet, sheet_marks = records("events_sheet")
        items = [*events_from_site(dict(self._json("events"), events=events), sheet,
                                   marks=event_marks, sheet_marks=sheet_marks, undated=undated),
                 *events_from_mrd(dict(self._json("mrd"), measurements=mrd), marks=mrd_marks, undated=undated),
                 *events_from_flow(dict(samples=records("flow")[0]), marks=records("flow")[1], undated=undated),
                 *events_from_imaging(dict(studies=records("imaging")[0]), marks=records("imaging")[1],
                                      undated=undated),
                 *events_from_pathology(dict(groups=records("pathology")[0]), marks=records("pathology")[1],
                                        undated=undated),
                 *events_from_specimens(records("specimens")[0], marks=records("specimens")[1], undated=undated),
                 *events_from_measurements(*records("labs"), "labs", "Lab draws", undated=undated),
                 *events_from_measurements(*records("cytometry"), "cytometry", "Cytometry panels",
                                           undated=undated)]
        # Rows whose dates cannot be read are kept out of the chart but reported here.
        return Timeline(sorted(items, key=lambda e: (e.first_day, e.lane, e.label, e.id)),
                        source={"snapshot_id": self.id, "corrections": self.curation.applied(),
                                "undated": undated})

    @cached_property
    def specimens(self):
        """The sample registry: one row per biological specimen, linked to its
        BAMs (via consolidated metadata) and FASTQ folders, with any dates or
        sites that other sources state differently."""
        self._require_timeline()
        rows, touched = self.curation.records("specimens")
        metadata = {r["display_name"]: r for r in self.curation.records("bam_metadata")[0]}
        fastqs, related = defaultdict(list), defaultdict(list)
        fastq_rows, fastq_touched = self.curation.records("fastqs")
        for i, row in enumerate(fastq_rows):
            fastqs[row["sample_id"]].append(row["s3_folder"])
            related[row["sample_id"]].extend(fastq_touched.get(i, ()))
        claims, summary_corrections = defaultdict(list), defaultdict(list)
        summary_rows, summary_touched = self.curation.records("timepoint_summary")
        for i, row in enumerate(summary_rows):
            claims[row["timepoint"]].append(("timepoint_summary", row["date"], row["location"]))
            summary_corrections[row["timepoint"]].extend(summary_touched.get(i, ()))
        for event in self.timeline.select(lane="Time points"):
            claims[event.timepoint].append(("events", event.date, None))
        result = []
        for i, row in enumerate(rows):
            names = [n.strip() for n in re.split(r"[;|]", row.get("bam_display_names", "")) if n.strip()]
            disagreements = []
            if row.get("timepoint") and normalize_tissue(row.get("tissue")) == "tumor":
                for source, day, location in claims.get(row["timepoint"], ()):
                    if day != row["collection_date"]:
                        disagreements.append(dict(source=source, field="date", value=day))
                    site = (location or "").split(" (")[0]
                    if site and site != row.get("collection_site"):
                        disagreements.append(dict(source=source, field="site", value=location))
            result.append(dict(
                sample_id=row["sample_id"], timepoint=row.get("timepoint") or None,
                date=row.get("collection_date") or None, tissue=normalize_tissue(row.get("tissue")),
                site=row.get("collection_site") or None, specimen=row.get("tissue_source") or None,
                vendors=tuple(v.strip() for v in row.get("vendors_involved", "").split(";") if v.strip()),
                assays=tuple(a.strip() for a in row.get("assays_run", "").split(";") if a.strip()),
                assets=tuple(self._object_key(metadata[n]["s3_path"]) for n in names
                             if metadata.get(n, {}).get("s3_path")),
                unmatched_bams=tuple(n for n in names if not metadata.get(n, {}).get("s3_path")),
                fastq_folders=tuple(fastqs.get(row["sample_id"], ())),
                disagreements=disagreements,
                corrections=tuple(dict.fromkeys((*touched.get(i, ()), *related[row["sample_id"]],
                                                 *(summary_corrections[row["timepoint"]]
                                                   if row.get("timepoint") else ())))),
                raw=dict(row)))
        return Table(result, source=self.manifest["sources"]["specimens"])

    def _object_key(self, path):
        """Bucket key of a relative key or a URL under this snapshot's download base."""
        if not urlsplit(path).scheme:
            return path
        return object_key(path, self._download_header.get("download_base", BUCKET))

    def describe_samples(self, *, timepoint=None, tissue=None, width=None):
        """Readable specimen/sequence-type overview; no data acquisition."""
        from .explore import samples_view
        return samples_view(self, timepoint=timepoint, tissue=tissue, width=width)

    def assets_for_sample(self, sample_id, **filters):
        """Registry-linked alignments and files under the specimen's FASTQ folders."""
        row = next((r for r in self.specimens if r["sample_id"] == sample_id), None)
        if row is None:
            raise KeyError(f"Unknown sample: {sample_id}")
        keys = set(row["assets"])
        prefixes = tuple(self._object_key(p).rstrip("/") + "/" for p in row["fastq_folders"])
        return self.assets.where(lambda a: a.key in keys or a.key.startswith(prefixes)).select(**filters)

    @cached_property
    def measurements(self):
        """MRD, lab, and cytometry values with raw strings and an explicit kind."""
        self._require_timeline()
        rows = []
        for event in self.timeline.select(source="mrd"):
            rows.append(dict(date=event.date, source="mrd", category="MRD",
                             measurement=event.details["assay_name"], value=event.value,
                             unit=event.details["unit"], kind=event.details["value_kind"],
                             reference_low="", reference_high="", out_of_range="",
                             corrections=";".join(event.corrections)))
        for name in ("labs", "cytometry"):
            source_rows, touched = self.curation.records(name)
            for i, row in enumerate(source_rows):
                value = row.get("value", "")
                kind = "missing" if value == "" else "numeric" if re.fullmatch(r"-?[\d.]+", value) else "text"
                rows.append(dict(date=row["date"], source=name, category=row.get("category", ""),
                                 measurement=row.get("measurement", ""), value=value, unit=row.get("unit", ""),
                                 kind=kind, reference_low=row.get("reference_low", ""),
                                 reference_high=row.get("reference_high", ""),
                                 out_of_range=row.get("out_of_range", ""),
                                 corrections=";".join(touched.get(i, ()))))
        return Table(sorted(rows, key=lambda r: (r["date"], r["source"], r["measurement"])))

    @property
    def corrections(self):
        """Every correction's status (applied, fixed_upstream, stale, disabled) and evidence."""
        return Table(self.curation.report())

    @property
    def unrecognized(self):
        """Source labels outside osteosarc.curation's vocabulary (possible upstream drift)."""
        return Table(unrecognized_values(
            bams=self._json("bams"), metadata=self.curation.records("bam_metadata")[0],
            vafs=self.vafs, source_variants=self.curation.records("source_variants")[0],
            data_page=[values for _, values in data_page_rows(read_text(self.source_path("data_page")))]),
            columns=("source", "field", "value", "records"))

    @cached_property
    def vafs(self):
        """Count/annotation rows; no inference from missing counts.

        With corrections enabled, a final 'corrections' column lists the IDs of
        corrections that edited or flagged each row.
        """
        source = self.manifest["sources"]["vafs"]
        rows, touched = self.curation.records("vafs")
        columns = self._columns["vafs"]
        diagnostics = copy.deepcopy(self._table_diagnostics["vafs"])
        if not self.curation.enabled:
            return Table((dict(row) for row in rows), columns=columns, source=source, diagnostics=diagnostics)
        return Table((dict(row, corrections=";".join(touched.get(i, ()))) for i, row in enumerate(rows)),
                     columns=(*columns, "corrections"), source=source, diagnostics=diagnostics)

    @cached_property
    def assets(self):
        """All listed files plus metadata-only catalog objects and site tables."""
        bams = self._json("bams")
        files, bucket_touched = self.curation.records("bucket")
        listing = self._download_header
        rows, bams_touched = self.curation.records("bams")
        metadata, metadata_touched = self.curation.records("bam_metadata")
        rows = iter(rows)
        bams = dict(bams, categories=[dict(c, bams=[next(rows) for _ in c["bams"]]) for c in bams["categories"]])
        listing = dict(listing, files=[[f["key"], f["size"], f["modified"]] for f in files])
        # Site tables pinned by this snapshot keep the snapshot's URLs (e.g. a mirror).
        tables = {name: (self.manifest["sources"].get(name, {}).get("url", url), format)
                  for name, (url, format) in TABLE_SOURCES.items()}
        assets = build_assets(listing, bams, metadata, self.vafs,
                              parse_data_paths(read_text(self.source_path("data_page"))), tables=tables)
        touched = {}
        base = listing.get("download_base", BUCKET)
        for i, ids in bucket_touched.items():
            touched.setdefault(files[i]["key"], []).extend(ids)
        flat = [row for category in bams["categories"] for row in category["bams"]]
        for i, ids in bams_touched.items():
            touched.setdefault(object_key(flat[i]["url"], bams["baseUrl"]), []).extend(ids)
        for i, ids in metadata_touched.items():
            touched.setdefault(object_key(metadata[i]["s3_path"], base), []).extend(ids)
        for asset in assets:
            if asset.key in touched:
                asset.metadata["corrections"] = tuple(dict.fromkeys(touched[asset.key]))
        # The listing is large; evaluated corrections keep their matches if it is reloaded.
        self.curation.release("bucket")
        return assets

    @cached_property
    def _download_header(self):
        """The bucket listing's fields other than its file rows (download_base, dates)."""
        if self._bucket_header is None:
            self.curation.raw("bucket")
            self.curation.release("bucket")
        return self._bucket_header

    @cached_property
    def _variants(self):
        index, index_touched = self.curation.records("variant_index")
        records, record_touched = self.curation.records("source_variants")
        mutations, overlap_touched = self.curation.records("vaccine_overlap")
        counts, count_touched = self.curation.records("vafs")
        # Variants get their own copies, so edits to their annotations never reach the sources.
        index, records, mutations = copy.deepcopy(index), copy.deepcopy(records), copy.deepcopy(mutations)
        variants = parse_variants(index, Table(copy.deepcopy(counts), columns=self._columns["vafs"],
                                              diagnostics=copy.deepcopy(self._table_diagnostics["vafs"])),
                                  source_variants=records,
                                  vaccine_overlap=dict(self._json("vaccine_overlap"), mutations=mutations),
                                  source={"snapshot_id": self.id, "receipts": self.manifest["sources"],
                                          "corrections": self.curation.applied()})
        # A correction belongs to a variant when it touches the variant's own records
        # or all of its count rows. One that touches only some count rows (e.g. one
        # BAM's counts) is listed separately as a count correction.
        own, partial = defaultdict(list), defaultdict(list)
        for rows, marks in ((index, index_touched), (records, record_touched)):
            for i, ids in marks.items():
                own[rows[i]["id"]].extend(ids)
        overlap_ids = {id(mutations[i]): ids for i, ids in overlap_touched.items()}
        rows_per_variant, touched_rows = Counter(), defaultdict(Counter)
        for i, row in enumerate(counts):
            rows_per_variant[row["variant_id"]] += 1
            for correction in count_touched.get(i, ()):
                touched_rows[row["variant_id"]][correction] += 1
        for variant in variants:
            ids = own[variant.id]
            for record in variant.annotations.get("vaccine_overlap_records", ()):
                ids.extend(overlap_ids.get(id(record), ()))
            for correction, n in touched_rows[variant.id].items():
                (ids if n == rows_per_variant[variant.id] else partial[variant.id]).append(correction)
            variant.annotations["corrections"] = tuple(dict.fromkeys(ids))
            variant.annotations["count_corrections"] = tuple(dict.fromkeys(partial[variant.id]))
        return variants

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

    def _traced(self, name):
        """Copies of a JSON source's records, each naming the corrections that touched it."""
        rows, touched = self.curation.records(name)
        rows = copy.deepcopy(rows)
        if self.curation.enabled:
            for i, row in enumerate(rows):
                row["corrections"] = touched.get(i, ())
        return rows

    @property
    def vaccines(self):
        """Vaccine-overlap rows including unmodified ELISPOT states (copies; edits do not persist)."""
        return Table(self._traced("vaccine_overlap"), source=self.manifest["sources"]["vaccine_overlap"])

    @property
    def vaccine_names(self):
        return tuple(self._json("vaccine_overlap")["vaccine_names"])

    def vaccine_peptides(self, vaccine=None):
        """Published peptide sequences and experiments, keyed by exact variant ID."""
        rows = []
        for record in self.curation.records("source_variants")[0]:
            for peptide in record.get("vaccine_peptides", []):
                if vaccine is None or vaccine in peptide.get("in_vaccines", []):
                    rows.append(dict(copy.deepcopy(peptide), variant_id=record["id"], gene=record["gene"]))
        return Table(rows, source=self.manifest["sources"]["source_variants"])

    @property
    def annotations(self):
        """Variant annotation records, including peptides and validation (copies)."""
        return Table(self._traced("source_variants"), source=self.manifest["sources"]["source_variants"])

    @property
    def pipeline_names(self):
        return tuple(sorted({p for record in self.curation.records("source_variants")[0]
                             for p in record.get("detection", {})}))

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

    @cached_property
    def _asset_index(self):
        index = defaultdict(list)
        for asset in self.assets:
            for name in dict.fromkeys((asset.id, asset.key, asset.url, asset.metadata.get("resource"))):
                if name:
                    index[name].append(asset)
        return index

    def asset(self, key_or_id):
        """Resolve an exact bucket key, URL, resource name, or asset ID."""
        matches = self._asset_index.get(key_or_id, [])
        if len(matches) != 1:
            raise KeyError(f"Expected one asset for {key_or_id!r}, found {len(matches)}")
        return matches[0]

    def download(self, asset, *, refresh=False, verify_size=True):
        """Fetch one full object and bind its receipt to this metadata snapshot.

        After the first acquisition, URL refreshes elsewhere cannot change this
        snapshot's bytes. Use a new snapshot to acquire a newer object. The
        binding is populated lazily, so metadata-only use downloads no data.
        """
        return self._download(asset, cache=self.cache, refresh=refresh, verify_size=verify_size)

    def _download(self, asset, *, cache, refresh=False, verify_size=True):
        """Acquire a snapshot-bound asset with the supplied cache's network policy."""
        asset = asset if isinstance(asset, Asset) else self.asset(asset)
        if asset.id != stable_id(asset.url):
            raise ValueError("Asset ID must be the stable hash of its URL")
        # Metadata tables already pinned by this snapshot never drift to latest.
        for receipt in self.manifest["sources"].values():
            if receipt["url"] == asset.url:
                if refresh:
                    raise ValueError("Create a new snapshot to refresh pinned metadata")
                return cache.path(receipt)
        binding = cache.workspace / "bindings" / self.id / (asset.id + ".json")
        with file_lock(binding.with_suffix(".lock")):
            if binding.exists():
                if refresh:
                    raise ValueError("Create a new snapshot to refresh an acquired object")
                receipt = Receipt(**json.loads(binding.read_text()))
                if receipt.url != asset.url:
                    raise IntegrityError("Snapshot object binding has the wrong URL")
                try:
                    return cache.path(receipt)
                except OfflineError:
                    if cache.offline:
                        raise
                    # A pruned object is restored only if the server still has the same bytes.
                    return cache.path(cache.fetch(asset.url, sha256=receipt.sha256, size=receipt.size))
            md5s = {r.get("md5sum") for r in asset.metadata.get("metadata_rows", []) if r.get("md5sum")}
            if len(md5s) > 1:
                raise IntegrityError("Conflicting published MD5 claims")
            receipt = cache.fetch(asset.url, refresh=refresh,
                                  size=asset.size if verify_size else None,
                                  md5=next(iter(md5s), None))
            _check_inventory_time(asset, receipt)
            write_json(binding, receipt.to_dict())
            return cache.path(receipt)

    def parse(self, asset):
        """Download explicitly selected data and parse it with original columns."""
        asset = asset if isinstance(asset, Asset) else self.asset(asset)
        if asset.url == self.manifest["sources"]["vafs"]["url"]:
            return self.vafs  # the pinned count table, with corrections when enabled
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

    def extract_reads(self, asset, regions=None, *, variants=None, padding=0, **kwargs):
        """Extract an indexed region union; see osteosarc.extract_reads.

        Supply regions or selected osteosarc variants with optional padding.
        The listed index is downloaded once, bound to this snapshot like any
        other object, and reused for every query.
        """
        from .reads import extract_reads, require_samtools
        if variants is not None:
            if regions is not None:
                raise ValueError("Supply either regions or variants, not both")
            regions = tuple(v.region(padding=padding) for v in variants)
        elif padding:
            raise ValueError("padding requires variants; pad explicit regions when constructing them")
        regions = tuple(regions) if regions is not None else ()
        if not regions or not all(isinstance(r, Region) for r in regions):
            raise CoordinateError("Provide a nonempty sequence of regions or ready variants")
        asset = asset if isinstance(asset, Asset) else self.asset(asset)
        if kwargs.get("index") is None and asset.index_urls:
            # Try the pinned/local index first so cached reads work without
            # SAMtools. Check capabilities before an index needs downloading.
            offline = Cache(self.cache.root, offline=True, timeout=self.cache.timeout)
            try:
                index_path = self._download(asset.index_urls[0], cache=offline)
            except OfflineError:
                if self.cache.offline:
                    raise
                require_samtools(fetch_pairs=kwargs.get("fetch_pairs", False), filters=kwargs.get("filters"))
                index_path = self.download(asset.index_urls[0])
            kwargs["index"] = str(index_path)
        return extract_reads(asset, regions, cache=self.cache, snapshot_id=self.id, **kwargs)

    def generate_bundle(self, recipe, destination, **kwargs):
        """Generate a fixture bundle with this snapshot's verified source assets."""
        from .bundles import generate_bundle
        return generate_bundle(recipe, destination, dataset=self, **kwargs)

    def select_fixtures(self, recipe, sources):
        """Execute a pinned fixture recipe on explicitly supplied local inputs.

        Membership and reasons are identical to osteosarc.select_fixtures.
        The snapshot does not override the recipe's historical source identity.
        """
        from .fixtures import select_fixtures
        return select_fixtures(recipe, sources)

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
