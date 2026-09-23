"""Indexed BAM/CRAM extraction and explicit, separate fixture downsampling."""

from __future__ import annotations

import json
import os
import re
import subprocess
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from urllib.parse import urlsplit

from .cache import (
    Cache,
    digest,
    file_identity,
    file_lock,
    http_identity,
    share,
    stable_id,
    write_json,
)
from .errors import CoordinateError, IntegrityError, OfflineError, OsteosarcError
from .models import Asset, Region

ASSEMBLY_LENGTHS = {
    "GRCh38": {"1": 248956422, "2": 242193529, "3": 198295559, "X": 156040895},
    "GRCh37": {"1": 249250621, "2": 243199373, "3": 198022430, "X": 155270560},
}


def normalize_assembly(value):
    """Normalize known assembly family names; mitochondrial sequence still matters."""
    return {"hg38": "GRCh38", "GRCh38": "GRCh38", "hg19": "GRCh37",
            "hs37d5": "GRCh37", "GRCh37": "GRCh37", "b37": "GRCh37"}.get(value, value)


def assembly_from_header(header):
    """Require at least two distinguishing contigs and no conflicting lengths."""
    sequences = header.get("SQ", [])
    matches = []
    for assembly, expected in ASSEMBLY_LENGTHS.items():
        relevant = [(row["SN"].removeprefix("chr"), row["LN"]) for row in sequences
                    if row["SN"].removeprefix("chr") in expected]
        if len({name for name, _ in relevant}) >= 2 and all(
                length == expected[name] for name, length in relevant):
            matches.append(assembly)
    return matches[0] if len(matches) == 1 else None


def resolve_regions(regions, header):
    """Check assembly, aliases, bounds and mitochondrial length, then merge union.

    Exact contig matches win. Aliases are limited to chr prefixes and M/MT;
    multiple possible aliases are an error. Returns zero-based Region objects
    in header order so indexed extraction preserves coordinate sorting.
    """
    regions = tuple(regions)
    if not regions:
        raise CoordinateError("At least one region is required; refusing a whole-file query")
    observed = assembly_from_header(header)
    if observed is None:
        raise CoordinateError("Cannot establish GRCh37/GRCh38 from the alignment header")
    sequences = header.get("SQ", [])
    lengths = {row["SN"]: row["LN"] for row in sequences}
    if len(lengths) != len(sequences):
        raise CoordinateError("Duplicate contig names in alignment header")
    order = {name: i for i, name in enumerate(lengths)}
    resolved = []
    for region in regions:
        if normalize_assembly(region.assembly) != observed:
            raise CoordinateError(f"Requested {region.assembly}; header establishes {observed}")
        canonical = region.contig.removeprefix("chr")
        aliases = {"M", "MT", "chrM", "chrMT"} if canonical in ("M", "MT") else {canonical, "chr" + canonical}
        candidates = {region.contig} if region.contig in lengths else aliases & lengths.keys()
        if len(candidates) != 1:
            raise CoordinateError(f"Missing or ambiguous contig {region.contig}: {sorted(candidates)}")
        contig = next(iter(candidates))
        if region.end > lengths[contig]:
            raise CoordinateError(f"Region extends beyond {contig}")
        expected_length = region.reference_length
        if canonical in ("M", "MT"):
            if observed == "GRCh37" and expected_length is None:
                raise CoordinateError("GRCh37 mitochondrial regions require an explicit reference_length")
            expected_length = expected_length or 16569
        if expected_length is not None and expected_length != lengths[contig]:
            raise CoordinateError(f"Reference length mismatch for {contig}")
        resolved.append(Region(contig, region.start, region.end, observed, expected_length))
    merged = []
    for region in sorted(resolved, key=lambda r: (order[r.contig], r.start, r.end)):
        if merged and merged[-1].contig == region.contig and region.start <= merged[-1].end:
            previous = merged[-1]
            merged[-1] = Region(region.contig, previous.start, max(previous.end, region.end),
                                observed, previous.reference_length or region.reference_length)
        else:
            merged.append(region)
    return tuple(merged)


@dataclass(frozen=True)
class ReadFilter:
    """Explicit acquisition filters; defaults preserve all overlapping records."""

    min_mapq: int = 0
    exclude_flags: int = 0
    require_flags: int = 0
    barcodes: tuple[str, ...] = ()
    barcode_tag: str = "CB"

    def __post_init__(self):
        # One barcode may be given as a string; never split it into characters.
        barcodes = (self.barcodes,) if isinstance(self.barcodes, str) else tuple(self.barcodes)
        object.__setattr__(self, "barcodes", barcodes)
        if not 0 <= self.min_mapq <= 255 or self.exclude_flags < 0 or self.require_flags < 0:
            raise ValueError("Invalid MAPQ or SAM flags")
        if len(self.barcode_tag) != 2 or not self.barcode_tag.isalnum():
            raise ValueError("barcode_tag must be a two-character SAM tag")
        if any(not value or "\n" in value or "\r" in value for value in self.barcodes):
            raise ValueError("Invalid barcode")


@dataclass(frozen=True)
class ReadSubset:
    path: Path
    index_path: Path
    receipt: dict

    def open(self):
        """Return a pysam AlignmentFile context manager for this indexed subset."""
        import pysam
        return pysam.AlignmentFile(str(self.path), index_filename=str(self.index_path))


@dataclass(frozen=True)
class AlignmentInfo:
    """Original header and acquisition evidence; assembly can be unresolved."""

    path: Path
    header: dict
    receipt: dict

    @property
    def assembly(self):
        return assembly_from_header(self.header)


def _run(command, timeout):
    return subprocess.run(command, check=True, capture_output=True, timeout=timeout)


def _samtools_version():
    # Distribution build flags can contain non-UTF-8 bytes after the version.
    return _run(["samtools", "--version"], 30).stdout.splitlines()[0].decode("utf-8", errors="replace")


def _remote_identity(url, timeout):
    return http_identity(url, timeout)


def require_samtools(*, header_only=False, fetch_pairs=False, filters=None):
    """Fail before acquisition when the installed binary lacks required options."""
    try:
        result = _run(["samtools", "view", "--help"], 30)
    except FileNotFoundError as error:
        raise OsteosarcError("Read extraction requires samtools on PATH; install SAMtools 1.21 or newer.") from error
    except subprocess.CalledProcessError as error:
        # Some versions print usable help but exit nonzero.
        result = error
    help_text = ((result.stdout or b"") + (result.stderr or b"")).decode(errors="replace")
    needed = ["--no-PG"]
    if not header_only:
        needed += ["-M", "-X"]
        if fetch_pairs:
            needed.append("--fetch-pairs")
        if filters is not None and filters.barcodes:
            needed.append("-D")
    missing = [flag for flag in needed if not re.search(r"(?<![\w-])" + re.escape(flag) + r"(?![\w-])", help_text)]
    if missing:
        raise OsteosarcError("samtools view lacks required options: " + ", ".join(missing)
                            + ". Install SAMtools 1.21 or newer; requested read filters and mate recovery cannot be omitted.")


def _verified_receipt(directory, request):
    path = directory / "receipt.json"
    if not path.exists():
        return None
    receipt = json.loads(path.read_text())
    if receipt["request"] != request:
        raise IntegrityError("Request differs from its cached receipt")
    for filename, expected in receipt["files"].items():
        if (Path(filename).name != filename or not (directory / filename).is_file()
                or digest(directory / filename) != expected):
            raise IntegrityError(f"Cached derivative was modified: {filename}")
    return receipt


def _cached_subset(directory, request):
    receipt = _verified_receipt(directory, request)
    if receipt is None:
        return None
    return ReadSubset(directory / "reads.bam", directory / "reads.bam.bai", receipt)


def _alignment_source(source):
    asset = source if isinstance(source, Asset) else None
    location = asset.url if asset else str(Path(source).resolve()) if not urlsplit(str(source)).scheme else str(source)
    remote = urlsplit(location).scheme in ("https", "http")
    if urlsplit(location).scheme and not remote:
        raise ValueError("Source must be a local path or HTTP(S) URL")
    if asset is not None and asset.format not in ("bam", "cram"):
        raise ValueError("Alignment access requires a BAM or CRAM asset")
    return asset, location, remote


def inspect_alignment(source, *, cache=None, snapshot_id=None, timeout=600):
    """Inspect and cache a BAM/CRAM header without requiring an index.

    Accepts the same sources as extract_reads. Remote headers are pinned on
    first inspection for a snapshot, checked against inventory size, and reused
    offline. This establishes assembly from contig lengths, not viewer labels.
    Extraction rejects a remote object that changed since header inspection.
    """
    import pysam
    cache = cache if isinstance(cache, Cache) else Cache(cache)
    asset, location, remote = _alignment_source(source)
    request = dict(schema_version=1, operation="inspect_alignment", source=location,
                   source_sha256=None if remote else cache.file_digest(location),
                   source_size=asset.size if asset else None,
                   source_modified=asset.modified if asset else None, snapshot_id=snapshot_id)
    directory = cache.workspace / "headers" / stable_id(request)
    with file_lock(cache.workspace / "locks" / (directory.name + ".lock")):
        receipt = _verified_receipt(directory, request)
        if receipt is None:
            if remote and cache.offline:
                raise OfflineError("Alignment header is not cached")
            require_samtools(header_only=True)
            before = _remote_identity(location, min(timeout, 60)) if remote else None
            if before and asset and asset.size is not None and before["content-length"] is not None:
                if int(before["content-length"]) != asset.size:
                    raise IntegrityError("Remote alignment size differs from the pinned inventory")
            identity = None if remote else file_identity(location)
            command = ["samtools", "view", "--no-PG", "-H", location]
            header_text = _run(command, timeout).stdout.decode()
            pysam.AlignmentHeader.from_text(header_text)
            after = _remote_identity(location, min(timeout, 60)) if remote else None
            if before != after or (not remote and file_identity(location) != identity):
                raise IntegrityError("Alignment changed during header inspection")
            directory.parent.mkdir(parents=True, exist_ok=True)
            with tempfile.TemporaryDirectory(dir=directory.parent, prefix=".header-") as temporary:
                work = Path(temporary)
                (work / "header.sam").write_text(header_text)
                receipt = dict(request=request, files={"header.sam": digest(work / "header.sam")},
                               remote_identity=before, command=command,
                               samtools_version=_samtools_version())
                write_json(work / "receipt.json", receipt)
                share(work)
                os.replace(work, directory)
        path = directory / "header.sam"
        return AlignmentInfo(path, pysam.AlignmentHeader.from_text(path.read_text()).to_dict(), receipt)


def extract_reads(source, regions, *, cache=None, index=None, filters=None, reference=None,
                  fetch_pairs=False, snapshot_id=None, timeout=600, recovery=None):
    """Fetch the indexed union of regions, retaining original record multiplicity.

    source can be an Asset, local BAM/CRAM, or HTTP(S) alignment URL. An index
    must be listed, supplied explicitly, or exist next to a local alignment.
    Remote requests never fall back to whole-file scans. Uses samtools -M -X;
    pysam validates the resulting BAM and creates its index. CRAM requires a
    local indexed reference FASTA, recorded by SHA256.

    The default retains only overlapping records. fetch_pairs=True additionally
    retrieves paired mates, not every supplementary alignment of each template.
    Cached results are immutable snapshot derivatives; they are verified offline
    without rechecking the remote object. New requests check remote identity
    before and after extraction and retain that identity in their receipt.
    """
    if recovery is not None:
        from .recovery import recover_reads
        return recover_reads(source, regions, policy=recovery, cache=cache, index=index,
                             filters=filters, reference=reference, fetch_pairs=fetch_pairs,
                             snapshot_id=snapshot_id, timeout=timeout)
    import pysam
    cache = cache if isinstance(cache, Cache) else Cache(cache)
    filters = filters or ReadFilter()
    regions = tuple(regions)
    if not regions or not all(isinstance(r, Region) for r in regions):
        raise CoordinateError("Provide a nonempty sequence of Region objects")
    asset, location, remote = _alignment_source(source)
    if index is None:
        if asset and asset.index_urls:
            index = asset.index_urls[0]
        elif not remote:
            suffix = Path(location).suffix.lower()
            candidates = [location + ".bai", str(Path(location).with_suffix(".bai")), location + ".csi"] if suffix == ".bam" else [location + ".crai", str(Path(location).with_suffix(".crai"))]
            index = next((p for p in candidates if Path(p).is_file()), None)
    if index is None:
        raise ValueError("No known index; full-alignment fallback is disabled")
    index = str(index)
    remote_index = urlsplit(index).scheme in ("http", "https")
    if not remote_index:
        index = str(Path(index).resolve())
    is_cram = (asset.format if asset else Path(urlsplit(location).path).suffix.lstrip(".")) == "cram"
    if is_cram and reference is None:
        raise CoordinateError("CRAM extraction requires an explicit local reference FASTA")
    reference = Path(reference).resolve() if reference is not None else None
    if reference is not None and not Path(str(reference) + ".fai").is_file():
        raise CoordinateError("Reference FASTA must already have a .fai index")
    local_identities = [file_identity(p) for p in ([] if remote else [location])
                        + ([] if remote_index else [index])]
    request = dict(schema_version=1, operation="extract_reads", source=location,
                   source_sha256=None if remote else cache.file_digest(location),
                   source_size=asset.size if asset else None,
                   source_modified=asset.modified if asset else None,
                   index=index, index_sha256=None if remote_index else cache.file_digest(index),
                   snapshot_id=snapshot_id,
                   regions=sorted((asdict(r) for r in regions), key=lambda r: json.dumps(r, sort_keys=True)),
                   filters=asdict(filters), fetch_pairs=fetch_pairs,
                   reference_sha256=cache.file_digest(reference) if reference else None,
                   reference_index_sha256=cache.file_digest(str(reference) + ".fai") if reference else None)
    # JSON normalization makes tuples and serialized lists compare identically.
    request = json.loads(json.dumps(request))
    directory = cache.workspace / "derived" / stable_id(request)
    with file_lock(cache.workspace / "locks" / (directory.name + ".lock")):
        cached = _cached_subset(directory, request)
        if cached is not None:
            return cached
        if remote and cache.offline:
            raise OfflineError("Regional reads are not cached")
        require_samtools(fetch_pairs=fetch_pairs, filters=filters)
        directory.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=directory.parent, prefix=".reads-") as temporary:
            work = Path(temporary)
            before = _remote_identity(location, min(timeout, 60)) if remote else None
            if before and asset and asset.size is not None and before["content-length"] is not None:
                if int(before["content-length"]) != asset.size:
                    raise IntegrityError("Remote alignment size differs from the pinned inventory")
            info = inspect_alignment(source, cache=cache, snapshot_id=snapshot_id, timeout=timeout)
            if remote and before != info.receipt["remote_identity"]:
                raise IntegrityError("Remote alignment changed since header inspection; use a new snapshot")
            header_text = info.path.read_text()
            resolved = resolve_regions(regions, info.header)
            index_receipt = None
            if remote_index:
                index_receipt = cache.fetch(index, refresh=not cache.offline)
                local_index = cache.path(index_receipt)
            else:
                local_index = Path(index)
            bed = work / "regions.bed"
            bed.write_text("".join(f"{r.contig}\t{r.start}\t{r.end}\n" for r in resolved))
            output = work / "reads.bam"
            command = ["samtools", "view", "--no-PG", "-b", "-M", "-X", "-L", str(bed),
                       "-q", str(filters.min_mapq), "-F", str(filters.exclude_flags),
                       "-f", str(filters.require_flags), "-o", str(output)]
            if reference:
                command += ["-T", str(reference)]
            if fetch_pairs:
                command += ["--fetch-pairs"]
            if filters.barcodes:
                barcodes = work / "barcodes.txt"
                barcodes.write_text("\n".join(sorted(set(filters.barcodes))) + "\n")
                command += ["-D", filters.barcode_tag + ":" + str(barcodes)]
            command += [location, str(local_index)]
            _run(command, timeout)
            _run(["samtools", "quickcheck", "-v", str(output)], min(timeout, 60))
            with pysam.AlignmentFile(output) as bam:
                count = sum(1 for _ in bam)
            pysam.index(str(output))
            after = _remote_identity(location, min(timeout, 60)) if remote else None
            if before != after:
                raise IntegrityError("Remote alignment changed during extraction")
            if [file_identity(p) for p in ([] if remote else [location])
                    + ([] if remote_index else [index])] != local_identities:
                raise IntegrityError("Local alignment or index changed during extraction")
            (work / "header.sam").write_text(header_text)
            files = {name: digest(work / name) for name in ("reads.bam", "reads.bam.bai", "header.sam", "regions.bed")}
            receipt = dict(request=request, files=files, records=count,
                           resolved_regions=[asdict(r) for r in resolved],
                           remote_identity=before,
                           header_receipt=info.receipt,
                           index_receipt=index_receipt.to_dict() if index_receipt else None,
                           samtools_version=_samtools_version(),
                           pysam_version=pysam.__version__, command=command,
                           scope="regional_records_and_paired_mates" if fetch_pairs else "regional_records")
            write_json(work / "receipt.json", receipt)
            # Only a complete directory becomes visible. No receipt means no cache hit.
            share(work)
            os.replace(work, directory)
        return _cached_subset(directory, request)


def subset_templates(source, count, *, cache=None, seed="0"):
    """Select up to count templates deterministically for a small local fixture.

    Template identity is (read group, query name). All available regional records
    for selected templates survive, including exact repeats. This is sampling
    from the supplied BAM, not remote mate recovery or evidence-based selection.
    """
    import pysam
    if not isinstance(count, int) or count < 0:
        raise ValueError("count must be a nonnegative integer")
    cache = cache if isinstance(cache, Cache) else Cache(cache)
    source = Path(source.path if isinstance(source, ReadSubset) else source).resolve()
    request = dict(schema_version=1, operation="subset_templates", source_sha256=cache.file_digest(source),
                   count=count, seed=str(seed))
    directory = cache.workspace / "derived" / stable_id(request)
    with file_lock(cache.workspace / "locks" / (directory.name + ".lock")):
        cached = _cached_subset(directory, request)
        if cached is not None:
            return cached
        def identity(read):
            return (read.get_tag("RG") if read.has_tag("RG") else "", read.query_name)
        with pysam.AlignmentFile(source) as bam:
            names = {identity(read) for read in bam}
        selected = set(sorted(names, key=lambda name: stable_id([str(seed), name]))[:count])
        directory.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=directory.parent, prefix=".subset-") as temporary:
            work = Path(temporary)
            output = work / "reads.bam"
            records = 0
            with pysam.AlignmentFile(source) as bam, pysam.AlignmentFile(output, "wb", template=bam) as out:
                for read in bam:
                    if identity(read) in selected:
                        out.write(read)
                        records += 1
            pysam.index(str(output))
            if cache.file_digest(source) != request["source_sha256"]:
                raise IntegrityError("Input BAM changed during template selection")
            receipt = dict(request=request, records=records, templates=len(selected),
                           files={name: digest(work / name) for name in ("reads.bam", "reads.bam.bai")},
                           selected_templates=sorted(selected), scope="sampled_fixture",
                           suitable_for_vaf=False, pysam_version=pysam.__version__)
            write_json(work / "receipt.json", receipt)
            share(work)
            os.replace(work, directory)
        return _cached_subset(directory, request)
