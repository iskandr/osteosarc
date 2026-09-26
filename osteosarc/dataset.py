"""One consistent entry point for a reproducible snapshot of osteosarc data."""

from __future__ import annotations

import copy
import json
import re
import warnings
from collections import Counter, defaultdict
from datetime import datetime, timezone
from email.utils import parsedate_to_datetime
from functools import cached_property
from pathlib import Path, PurePosixPath
from urllib.parse import urlsplit

from .cache import Cache, Receipt, file_lock, place, stable_id, write_json
from .catalog import (
    BUCKET,
    SNAPSHOT_SOURCES,
    TABLE_SOURCES,
    TIMELINE_SOURCES,
    build_files,
    data_page_rows,
    object_key,
    parse_data_paths,
)
from .curation import (
    CORRECTIONS,
    Curation,
    normalize_provider,
    normalize_tissue,
    sequencing_pairs,
    unrecognized_values,
)
from .errors import CoordinateError, IntegrityError, NoSnapshotsError, OfflineError, SchemaError
from .models import File, Region
from .parsing import (
    PARSE_FORMATS,
    Table,
    parse_file,
    parse_table,
    parse_variant_index,
    parse_variants,
    read_text,
)


def _check_inventory_time(file, receipt):
    """Refuse bytes newer (or older) than the object this snapshot's inventory listed."""
    if not isinstance(file.modified, (int, float)) or not receipt.last_modified:
        return
    served = parsedate_to_datetime(receipt.last_modified).timestamp()
    if abs(served - file.modified) > 1:
        raise IntegrityError(
            f"{file.key} was modified at {receipt.last_modified}, not at the time this snapshot's "
            f"inventory lists ({datetime.fromtimestamp(file.modified, timezone.utc).isoformat()}); "
            "create a new snapshot to use the current object")


GUIDE = """What's here:
  data.samples          tumor, organoid and blood samples; data.samples["T1_tumor"] shows one
  data.files            every file in the bucket, and the site's tables
  data.variants()       the variant catalogue, with alleles
  data.vaccines         vaccine targets and ELISPOT results
  data.timeline         treatments, procedures, scans, MRD and lab draws
  data.corrections      known problems in the website's data, and their fixes

Get data:
  data.download(key)                       download a whole file; returns its local path
  data.extract_reads(key, variants=...)    reads around variants, as a small local BAM
  data.downloads()                         what's already on this computer, and where"""

DATE_SELECTOR = re.compile(r"\d{4}(-\d{2}(-\d{2})?)?")
DATED_NAME = re.compile(r"\d{4}-\d{2}-\d{2}(\.\d+)?")


def _asked_for(variants, regions):
    """A short file-name part for an extraction: its first variant or region, and how many more."""
    first = variants[0].id if variants is not None else f"{regions[0].contig}_{regions[0].start + 1}-{regions[0].end}"
    count = len(variants) if variants is not None else len(regions)
    return first + (f"+{count - 1}" if count > 1 else "")


def _read_receipt(path):
    """A receipt stored as JSON, or None if it's missing or unreadable."""
    try:
        return Receipt(**json.loads(Path(path).read_text()))
    except (OSError, ValueError, TypeError):
        return None


def downloaded_at(manifest):
    """UTC time of a snapshot's latest source download (its creation time if it has none)."""
    return max((r["retrieved_at"] for r in manifest["sources"].values()), default=manifest.get("created_at"))


def choose_snapshot(rows, name=None, *, date=None, root=None):
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
        raise NoSnapshotsError(root)
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
        dataset.files
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
            name = choose_snapshot(cls.snapshots(cache=cache), name, date=date, root=cache.root)
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
        samples, MRD, flow draws, and lab/cytometry draw dates."""
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
    def _sample_sources(self):
        """Registry rows, with each sample's BAM keys, missing BAM names and FASTQ table rows."""
        self._require_timeline()
        rows, touched = self.curation.records("specimens")
        if any(not row.get("sample_id") for row in rows):
            raise SchemaError("Every row of the sample registry needs a sample_id")
        metadata = {r["display_name"]: r for r in self.curation.records("bam_metadata")[0]}
        fastqs, related = defaultdict(list), defaultdict(list)
        fastq_rows, fastq_touched = self.curation.records("fastqs")
        for i, row in enumerate(fastq_rows):
            if row.get("sample_id") and row.get("s3_folder"):
                fastqs[row["sample_id"]].append(row)
            related[row.get("sample_id")].extend(fastq_touched.get(i, ()))
        result = []
        for i, row in enumerate(rows):
            names = [n.strip() for n in re.split(r"[;|]", row.get("bam_display_names", "")) if n.strip()]
            result.append(dict(
                row=row, corrections=(*touched.get(i, ()), *related[row["sample_id"]]),
                bams=tuple(dict.fromkeys(self._object_key(metadata[n]["s3_path"]) for n in names
                                         if metadata.get(n, {}).get("s3_path"))),
                missing=tuple(n for n in names if not metadata.get(n, {}).get("s3_path")),
                fastqs=tuple(fastqs.get(row["sample_id"], ()))))
        return result

    @cached_property
    def samples(self):
        """What was collected: tumor, organoid and blood samples, each with its
        sequencing, BAMs and FASTQ folders. Index by ID: data.samples["T1_tumor"].

        Sequencing and providers combine the site's sample registry with its
        FASTQ table, which names each folder's assay even where the registry
        doesn't. Dates or sites that other sources state differently are
        listed in disagreements.
        """
        from .models import Sample, Samples
        claims, summary_corrections = defaultdict(list), defaultdict(list)
        summary_rows, summary_touched = self.curation.records("timepoint_summary")
        for i, row in enumerate(summary_rows):
            claims[row["timepoint"]].append(("timepoint_summary", row["date"], row["location"]))
            summary_corrections[row["timepoint"]].extend(summary_touched.get(i, ()))
        for event in self.timeline.select(lane="Time points"):
            claims[event.timepoint].append(("events", event.date, None))
        result = []
        for source in self._sample_sources:
            row = source["row"]
            tissue = normalize_tissue(row.get("tissue"))
            disagreements = []
            if row.get("timepoint") and tissue == "tumor":
                for name, day, location in claims.get(row["timepoint"], ()):
                    if day != row["collection_date"]:
                        disagreements.append(dict(source=name, field="date", value=day))
                    site = (location or "").split(" (")[0]
                    if site and site != row.get("collection_site"):
                        disagreements.append(dict(source=name, field="site", value=location))
            labels = [a.strip() for a in row.get("assays_run", "").split(";") if a.strip()]
            labels += [f["assay"] for f in source["fastqs"] if f.get("assay")]
            providers = [v.strip() for v in row.get("vendors_involved", "").split(";") if v.strip()]
            providers += [f["provider"] for f in source["fastqs"] if f.get("provider")]
            sample = Sample(
                id=row["sample_id"], timepoint=row.get("timepoint") or None,
                date=row.get("collection_date") or None, tissue=tissue,
                site=row.get("collection_site") or None, description=row.get("tissue_source") or None,
                providers=tuple(sorted({normalize_provider(p) for p in providers})),
                sequencing=sequencing_pairs(labels),
                bams=source["bams"], missing_bams=source["missing"],
                fastq_folders=tuple(dict.fromkeys(self._object_key(f["s3_folder"]).rstrip("/")
                                                  for f in source["fastqs"])),
                notes=row.get("notes", ""), disagreements=tuple(disagreements),
                corrections=tuple(dict.fromkeys((*source["corrections"], *(
                    summary_corrections[row["timepoint"]] if row.get("timepoint") else ())))),
                details=dict(registry=dict(row), fastqs=[
                    dict(f, folder=self._object_key(f["s3_folder"]).rstrip("/")) for f in source["fastqs"]]))
            object.__setattr__(sample, "_dataset", self)
            result.append(sample)
        return Samples(result, source={"snapshot_id": self.id})

    def _tag_samples(self, files):
        """Record on each file the samples whose BAM it is or whose FASTQ folder holds it."""
        try:
            sources = self._sample_sources
        except SchemaError as error:
            # An older snapshot has no sample registry; any other problem is worth saying.
            if set(TIMELINE_SOURCES) <= set(self.manifest["sources"]):
                warnings.warn(f"Files aren't linked to samples: {error}", stacklevel=3)
            return
        owners, folders = defaultdict(list), defaultdict(list)
        for source in sources:
            sample = source["row"]["sample_id"]
            for key in source["bams"]:
                owners[key].append(sample)
            for row in source["fastqs"]:
                folders[self._object_key(row["s3_folder"]).rstrip("/")].append(sample)
        tops = {folder.split("/", 1)[0] for folder in folders}
        for file in files:
            found = list(owners.get(file.key, ()))
            if file.key.split("/", 1)[0] in tops:  # most of the bucket is under other folders
                parts = file.key.split("/")
                for depth in range(1, len(parts)):
                    found.extend(folders.get("/".join(parts[:depth]), ()))
            if found:
                file.metadata["samples"] = tuple(dict.fromkeys(found))

    def _object_key(self, path):
        """Bucket key of a relative key or a URL under this snapshot's download base."""
        if not urlsplit(path).scheme:
            return path
        return object_key(path, self._download_header.get("download_base", BUCKET))

    def summary(self):
        """Counts of everything in this snapshot, and what to try next."""
        from .display import Text
        from .views import summary_view
        return Text(summary_view(self) + "\n\n" + GUIDE)

    def __repr__(self):
        from .views import snapshot_line
        return f"Osteosarc {snapshot_line(self)}\n\n{GUIDE}"

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
    def files(self):
        """All listed files plus metadata-only catalog objects and site tables."""
        bams = self._json("bams")
        objects, bucket_touched = self.curation.records("bucket")
        listing = self._download_header
        rows, bams_touched = self.curation.records("bams")
        metadata, metadata_touched = self.curation.records("bam_metadata")
        rows = iter(rows)
        bams = dict(bams, categories=[dict(c, bams=[next(rows) for _ in c["bams"]]) for c in bams["categories"]])
        listing = dict(listing, files=[[f["key"], f["size"], f["modified"]] for f in objects])
        # Site tables pinned by this snapshot keep the snapshot's URLs (e.g. a mirror).
        tables = {name: (self.manifest["sources"].get(name, {}).get("url", url), format)
                  for name, (url, format) in TABLE_SOURCES.items()}
        files = build_files(listing, bams, metadata, self.vafs,
                            parse_data_paths(read_text(self.source_path("data_page"))), tables=tables)
        touched = {}
        base = listing.get("download_base", BUCKET)
        for i, ids in bucket_touched.items():
            touched.setdefault(objects[i]["key"], []).extend(ids)
        flat = [row for category in bams["categories"] for row in category["bams"]]
        for i, ids in bams_touched.items():
            touched.setdefault(object_key(flat[i]["url"], bams["baseUrl"]), []).extend(ids)
        for i, ids in metadata_touched.items():
            touched.setdefault(object_key(metadata[i]["s3_path"], base), []).extend(ids)
        for file in files:
            if file.key in touched:
                file.metadata["corrections"] = tuple(dict.fromkeys(touched[file.key]))
        self._tag_samples(files)
        # The listing is large; evaluated corrections keep their matches if it is reloaded.
        self.curation.release("bucket")
        return files

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

    def vaccine_peptides(self, vaccine=None):
        """Published peptide sequences and experiments, keyed by exact variant ID."""
        rows = []
        for record in self.curation.records("source_variants")[0]:
            for peptide in record.get("vaccine_peptides", []):
                if vaccine is None or vaccine in peptide.get("in_vaccines", []):
                    rows.append(dict(copy.deepcopy(peptide), variant_id=record["id"], gene=record["gene"]))
        return Table(rows, source=self.manifest["sources"]["source_variants"])

    @cached_property
    def _file_index(self):
        index = defaultdict(list)
        for file in self.files:
            for name in dict.fromkeys((file.id, file.key, file.url, file.metadata.get("resource"))):
                if name:
                    index[name].append(file)
        return index

    def file(self, key_or_id):
        """Resolve an exact bucket key, URL, resource name, or file ID."""
        matches = self._file_index.get(key_or_id, [])
        if len(matches) != 1:
            raise KeyError(f"Expected one file for {key_or_id!r}, found {len(matches)}")
        return matches[0]

    def download(self, file, *, to=None, refresh=False, verify_size=True):
        """Download one whole file into the cache and return its local path.

        to names a directory to put the file in under its own name, and the
        path returned is then that one; a BAM's or VCF's index goes there too.
        Each is a read-only hard link to the cached copy where possible (so
        it takes no space and can't corrupt the cache), else a copy.

        The download is bound to this snapshot: after the first, changes on
        the server can't alter this snapshot's bytes (use a new snapshot for
        newer ones). Metadata-only use downloads nothing.
        """
        file = file if isinstance(file, File) else self.file(file)
        path = self._download(file, cache=self.cache, refresh=refresh, verify_size=verify_size)
        if to is None:
            return path
        directory = Path(to).expanduser()
        if file.index_urls:
            index = self.file(file.index_urls[0])
            place(self._download(index, cache=self.cache), directory / PurePosixPath(index.key).name)
        return place(path, directory / PurePosixPath(file.key).name)

    def local_path(self, file):
        """The downloaded copy of a file, or None; never uses the network.

        A copy downloaded for another snapshot counts when its size matches
        this snapshot's listing.
        """
        file = file if isinstance(file, File) else self.file(file)
        for name, receipt in self.manifest["sources"].items():
            if receipt["url"] == file.url:
                return self.source_path(name)
        pointers = (self.cache.workspace / "bindings" / self.id / (file.id + ".json"),
                    self.cache.workspace / "urls" / (stable_id(file.url) + ".json"))
        for pointer in pointers:
            path = self._local_copy(file, _read_receipt(pointer))
            if path is not None:
                return path
        return None

    def _local_copy(self, file, receipt):
        """The cached object behind a receipt, if it's this file's bytes (same URL and size)."""
        if receipt is None or receipt.url != file.url or (file.size is not None and receipt.size != file.size):
            return None
        try:
            return self.cache.path(receipt, verify=False)
        except (OfflineError, IntegrityError):
            return None

    def downloads(self):
        """Files downloaded into the cache and reads extracted from BAMs, with local paths.

        Rows have kind ('file' or 'reads'), key (the file's key, or its URL if
        this snapshot doesn't list it), url, size, path and, for extracted
        reads, regions. Only this dataset's bucket is listed (the cache is
        shared with other tools), and not the snapshot's own metadata. A file
        counts as downloaded exactly when local_path finds it.
        """
        from .views import regions_text
        by_url = {f.url: f for f in self.files}
        sources = {r["url"] for r in self.manifest["sources"].values()}
        base = self._download_header.get("download_base", BUCKET)
        receipts = {}
        # This snapshot's own downloads first, then what the cache has for other snapshots.
        for folder in (self.cache.workspace / "bindings" / self.id, self.cache.workspace / "urls"):
            for pointer in sorted(folder.glob("*.json")):
                receipt = _read_receipt(pointer)
                if receipt is not None:
                    receipts.setdefault(receipt.url, []).append(receipt)
        rows = []
        for url, candidates in sorted(receipts.items()):
            file = by_url.get(url)
            if url in sources or not (file or url.startswith(base)):
                continue
            for receipt in candidates:
                path = self._local_copy(file or File(stable_id(url), url, url, "", ""), receipt)
                if path is not None:
                    rows.append(dict(kind="file", key=file.key if file else url, url=url, size=receipt.size,
                                     path=str(path), regions="", downloaded=receipt.retrieved_at))
                    break
        for receipt_path in sorted((self.cache.workspace / "derived").glob("*/receipt.json")):
            try:
                receipt = json.loads(receipt_path.read_text())
                source = receipt["request"]["source"]
                regions = receipt["request"].get("regions", ())
            except (OSError, ValueError, KeyError, TypeError):
                continue
            bam = receipt_path.parent / "reads.bam"
            file = by_url.get(source)
            if not bam.is_file() or not (file or str(source).startswith(base)):
                continue
            rows.append(dict(kind="reads", key=file.key if file else source, url=source, size=bam.stat().st_size,
                             path=str(bam), regions=regions_text(regions),
                             downloaded=datetime.fromtimestamp(bam.stat().st_mtime, timezone.utc).isoformat()))
        return Table(rows, columns=("kind", "key", "url", "size", "path", "regions", "downloaded"))

    def local_urls(self):
        """URLs of every file with a local copy: downloads and the snapshot's own tables."""
        return ({r["url"] for r in self.downloads() if r["kind"] == "file"}
                | {r["url"] for r in self.manifest["sources"].values()})

    def _download(self, file, *, cache, refresh=False, verify_size=True):
        """Acquire a snapshot-bound file with the supplied cache's network policy."""
        file = file if isinstance(file, File) else self.file(file)
        if file.id != stable_id(file.url):
            raise ValueError("File ID must be the stable hash of its URL")
        # Metadata tables already pinned by this snapshot never drift to latest.
        for receipt in self.manifest["sources"].values():
            if receipt["url"] == file.url:
                if refresh:
                    raise ValueError("Create a new snapshot to refresh pinned metadata")
                return cache.path(receipt)
        binding = cache.workspace / "bindings" / self.id / (file.id + ".json")
        with file_lock(binding.with_suffix(".lock")):
            if binding.exists():
                if refresh:
                    raise ValueError("Create a new snapshot to refresh an acquired object")
                receipt = Receipt(**json.loads(binding.read_text()))
                if receipt.url != file.url:
                    raise IntegrityError("Snapshot object binding has the wrong URL")
                try:
                    return cache.path(receipt)
                except OfflineError:
                    if cache.offline:
                        raise
                    # A pruned object is restored only if the server still has the same bytes.
                    return cache.path(cache.fetch(file.url, sha256=receipt.sha256, size=receipt.size))
            md5s = {r.get("md5sum") for r in file.metadata.get("metadata_rows", []) if r.get("md5sum")}
            if len(md5s) > 1:
                raise IntegrityError("Conflicting published MD5 claims")
            try:
                receipt = cache.fetch(file.url, refresh=refresh,
                                      size=file.size if verify_size else None,
                                      md5=next(iter(md5s), None))
            except OfflineError as error:
                raise OfflineError(f"{error}. This snapshot was opened offline; reopen it with "
                                   "Dataset.open(offline=False), or use osteosarc download, "
                                   "to download it") from None
            _check_inventory_time(file, receipt)
            write_json(binding, receipt.to_dict())
            return cache.path(receipt)

    def parse(self, file):
        """Download explicitly selected data and parse it with original columns."""
        file = file if isinstance(file, File) else self.file(file)
        if file.url == self.manifest["sources"]["vafs"]["url"]:
            return self.vafs  # the pinned count table, with corrections when enabled
        if file.format not in PARSE_FORMATS:
            raise ValueError(f"No built-in parser for {file.format!r}; download the original file")
        return parse_file(self.download(file), format=file.format,
                          source=dict(url=file.url, snapshot_id=self.id))

    def inspect_alignment(self, file, **kwargs):
        """Cache the alignment header to inspect assembly before choosing regions."""
        from .reads import inspect_alignment
        file = file if isinstance(file, File) else self.file(file)
        return inspect_alignment(file, cache=self.cache, snapshot_id=self.id, **kwargs)

    def extract_reads(self, file, regions=None, *, variants=None, padding=0, to=None, name=None, **kwargs):
        """Stream just the reads in some regions, or around variants, into a small indexed BAM.

        Give regions, or variants (from data.variants()) with optional padding.
        Only the needed parts of the remote BAM are read, and asking again
        reuses the result. to names a directory to also put the BAM and its
        index in, as NAME.bam: by default the source file's name and what was
        asked for, such as BG003082.MAP2-chr2-209694768.bam. Other options:
        see osteosarc.extract_reads.
        """
        from .reads import ReadSubset, extract_reads, require_samtools
        if variants is not None:
            if regions is not None:
                raise ValueError("Supply either regions or variants, not both")
            regions = tuple(v.region(padding=padding) for v in variants)
        elif padding:
            raise ValueError("padding requires variants; pad explicit regions when constructing them")
        regions = tuple(regions) if regions is not None else ()
        if not regions or not all(isinstance(r, Region) for r in regions):
            raise CoordinateError("Provide a nonempty sequence of regions or ready variants")
        file = file if isinstance(file, File) else self.file(file)
        if kwargs.get("index") is None and file.index_urls:
            # Try the pinned/local index first so cached reads work without
            # SAMtools. Check capabilities before an index needs downloading.
            offline = Cache(self.cache.root, offline=True, timeout=self.cache.timeout)
            try:
                index_path = self._download(file.index_urls[0], cache=offline)
            except OfflineError:
                if self.cache.offline:
                    raise
                require_samtools(fetch_pairs=kwargs.get("fetch_pairs", False), filters=kwargs.get("filters"))
                index_path = self.download(file.index_urls[0])
            kwargs["index"] = str(index_path)
        subset = extract_reads(file, regions, cache=self.cache, snapshot_id=self.id, **kwargs)
        if to is None:
            return subset
        name = name or f"{PurePosixPath(file.key).name.split('.')[0]}.{_asked_for(variants, regions)}"
        directory = Path(to).expanduser()
        return ReadSubset(place(subset.path, directory / f"{name}.bam"),
                          place(subset.index_path, directory / f"{name}.bam.bai"), subset.receipt)
