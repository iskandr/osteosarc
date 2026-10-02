"""Indexed BAM/CRAM extraction and explicit, separate fixture downsampling."""

from __future__ import annotations

import json
import os
import re
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict, dataclass
from pathlib import Path
from threading import Event, Timer
from urllib.parse import urlsplit

from .cache import (
    ATTEMPTS,
    Cache,
    digest,
    file_identity,
    file_lock,
    http_identity,
    share,
    stable_id,
    write_json,
)
from .errors import CoordinateError, IntegrityError, OfflineError, OsteosarcError, RecordLimitError
from .models import File, Region
from .read_receipts import query_names_asset, read_read_receipt, write_read_receipt

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
        if region.assembly is not None and normalize_assembly(region.assembly) != observed:
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
    query_names: tuple[str, ...] = ()

    def __post_init__(self):
        # One barcode may be given as a string; never split it into characters.
        barcodes = (self.barcodes,) if isinstance(self.barcodes, str) else tuple(self.barcodes)
        object.__setattr__(self, "barcodes", barcodes)
        names = (self.query_names,) if isinstance(self.query_names, str) else tuple(self.query_names)
        if any(not isinstance(name, str) or not re.fullmatch(r"[!-?A-~]{1,254}", name) for name in names):
            raise ValueError("Invalid SAM query name")
        object.__setattr__(self, "query_names", tuple(sorted(set(names))))
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

    @property
    def receipt_path(self):
        """On-disk provenance; use read_receipt_files to pin its shared assets."""
        return self.path.parent / "receipt.json"

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


def _recorded(path, cache):
    """A local file as receipts and cache keys record it: its path within the cache
    (objects/sha256/<checksum>..., osteosarc/...), or else its file name. Its checksum
    is recorded beside it, so neither depends on where the cache or the file lives."""
    path = Path(path).resolve()
    try:
        return path.relative_to(cache.root).as_posix()
    except ValueError:
        return path.name


def _recorded_command(command, work, cache):
    """A command as receipts record it: files in the work folder by their names there,
    other local files as _recorded gives them (also after a samtools TAG: prefix)."""
    def record(arg):
        tag, rest = (arg[:3], arg[3:]) if re.fullmatch(r"[A-Za-z][A-Za-z0-9]:/.+", arg) else ("", arg)
        if rest.startswith(str(work) + os.sep):
            rest = Path(rest).relative_to(work).as_posix()
        elif rest.startswith("/"):
            rest = _recorded(rest, cache)
        return tag + rest
    return [record(arg) for arg in command]


def merge_spans(spans):
    """Overlapping or touching (contig, start, end) spans merged, sorted: the same bases, fewer spans."""
    merged = []
    for contig, start, end in sorted(spans):
        if merged and merged[-1][0] == contig and start <= merged[-1][2]:
            merged[-1][2] = max(merged[-1][2], end)
        else:
            merged.append([contig, start, end])
    return merged


def overlap_test(spans):
    """A test of whether contig, start, end overlaps any of these (contig, start, end) spans."""
    from bisect import bisect_left
    by_contig = {}
    for contig, start, end in merge_spans(spans):
        starts, ends = by_contig.setdefault(contig, ([], []))
        starts.append(start)
        ends.append(end)

    def overlaps(contig, start, end):
        starts, ends = by_contig.get(contig, ((), ()))
        last = bisect_left(starts, end) - 1  # the last span starting before end
        return last >= 0 and ends[last] > start
    return overlaps


def _add_mates(regional, mates, inside, output):
    """Write the regional records and, in coordinate order among them, the mates'
    records that don't overlap the regions inside tests for (the others are regional
    already); returns how many of those. Records are copied byte for byte (samtools
    merge would move their RG tags); records starting at one position may come in
    another order than the source's. With none to add, regional is output as it is."""
    import heapq

    import pysam
    with pysam.AlignmentFile(str(mates)) as second:  # where samtools takes a record to end: bam_endpos
        outside = [r for r in second
                   if not inside(r.reference_name, r.reference_start, r.reference_end or r.reference_start + 1)]
    if not outside:
        os.replace(regional, output)
        return 0
    with pysam.AlignmentFile(str(regional)) as first, pysam.AlignmentFile(str(output), "wb", template=first) as out:
        for read in heapq.merge(first, outside, key=lambda r: (r.reference_id, r.reference_start)):
            out.write(read)
    regional.unlink()
    return len(outside)


def _run_bounded(command, output, max_records, timeout):
    """Stop indexed acquisition on overflow; never publish a partial BAM."""
    import pysam
    command = list(command)
    command[command.index("-o") + 1] = "-"
    expired = Event()
    with tempfile.TemporaryFile() as stderr:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=stderr)

        def expire():
            expired.set()
            process.kill()

        timer = Timer(timeout, expire)
        timer.daemon = True
        timer.start()
        try:
            with pysam.AlignmentFile(process.stdout, "rb") as bam:
                with pysam.AlignmentFile(output, "wb", template=bam) as out:
                    for count, read in enumerate(bam, 1):
                        if count > max_records:
                            raise RecordLimitError("Acquisition exceeds record limit")
                        out.write(read)
            returncode = process.wait()
            if returncode:
                stderr.seek(0)
                raise subprocess.CalledProcessError(returncode, command, stderr=stderr.read())
        except Exception as error:
            if not isinstance(error, (IntegrityError, subprocess.CalledProcessError)):
                # Reap early command failures before choosing the error to report.
                # The timer still bounds a process blocked after malformed output.
                process.wait()
            if expired.is_set():
                raise subprocess.TimeoutExpired(command, timeout) from error
            if process.poll() not in (None, 0) and not isinstance(error, (IntegrityError, subprocess.CalledProcessError)):
                stderr.seek(0)
                raise subprocess.CalledProcessError(process.returncode, command, stderr=stderr.read()) from error
            raise
        finally:
            timer.cancel()
            if process.poll() is None:
                process.kill()
            process.wait()
            process.stdout.close()
        if expired.is_set():
            raise subprocess.TimeoutExpired(command, timeout)
    return command


def _samtools_version():
    # Distribution build flags can contain non-UTF-8 bytes after the version.
    return _run(["samtools", "--version"], 30).stdout.splitlines()[0].decode("utf-8", errors="replace")


def _remote_identity(url, timeout, attempts=ATTEMPTS):
    return http_identity(url, timeout, attempts=attempts)


def require_samtools(*, header_only=False, fetch_pairs=False, unplaced_mates=True, filters=None):
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
        if fetch_pairs and unplaced_mates:
            needed.append("--fetch-pairs")
        if filters is not None and filters.barcodes:
            needed.append("-D")
        if (fetch_pairs and not unplaced_mates) or (filters is not None and filters.query_names):
            needed.append("-N")
    missing = [flag for flag in needed if not re.search(r"(?<![\w-])" + re.escape(flag) + r"(?![\w-])", help_text)]
    if missing:
        raise OsteosarcError("samtools view lacks required options: " + ", ".join(missing)
                            + ". Install SAMtools 1.21 or newer; requested read filters and mate recovery cannot be omitted.")


def _verified_receipt(directory, request):
    path = directory / "receipt.json"
    if not path.exists():
        return None
    receipt = read_read_receipt(path, directory.parent.parent)
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
    file = source if isinstance(source, File) else None
    location = file.url if file else str(Path(source).resolve()) if not urlsplit(str(source)).scheme else str(source)
    remote = urlsplit(location).scheme in ("https", "http")
    if urlsplit(location).scheme and not remote:
        raise ValueError("Source must be a local path or HTTP(S) URL")
    if file is not None and file.format not in ("bam", "cram"):
        raise ValueError("Alignment access requires a BAM or CRAM file")
    return file, location, remote


def inspect_alignment(source, *, cache=None, snapshot_id=None, timeout=600):
    """Inspect and cache a BAM/CRAM header without requiring an index.

    Accepts the same sources as extract_reads. Remote headers are pinned on
    first inspection for a snapshot, checked against inventory size, and reused
    offline. This establishes assembly from contig lengths, not viewer labels.
    Extraction rejects a remote object that changed since header inspection.
    """
    import pysam
    cache = cache if isinstance(cache, Cache) else Cache(cache)
    file, location, remote = _alignment_source(source)
    request = dict(schema_version=1, operation="inspect_alignment",
                   source=location if remote else _recorded(location, cache),
                   source_sha256=None if remote else cache.file_digest(location),
                   source_size=file.size if file else None,
                   source_modified=file.modified if file else None, snapshot_id=snapshot_id)
    directory = cache.workspace / "headers" / stable_id(request)
    with file_lock(cache.workspace / "locks" / (directory.name + ".lock")):
        receipt = _verified_receipt(directory, request)
        if receipt is None:
            if remote and cache.offline:
                raise OfflineError("Alignment header is not cached")
            require_samtools(header_only=True)
            before = _remote_identity(location, min(timeout, 60)) if remote else None
            if before and file and file.size is not None and before["content-length"] is not None:
                if int(before["content-length"]) != file.size:
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
                               remote_identity=before, command=_recorded_command(command, work, cache),
                               samtools_version=_samtools_version())
                write_json(work / "receipt.json", receipt)
                share(work)
                os.replace(work, directory)
        path = directory / "header.sam"
        return AlignmentInfo(path, pysam.AlignmentHeader.from_text(path.read_text()).to_dict(), receipt)


def extract_reads(source, regions, *, cache=None, index=None, filters=None, reference=None,
                  fetch_pairs=False, unplaced_mates=True, snapshot_id=None, timeout=600, recovery=None,
                  max_records=None):
    """Fetch the indexed union of regions, retaining original record multiplicity.

    source can be an File, local BAM/CRAM, or HTTP(S) alignment URL. An index
    must be listed, supplied explicitly, or exist next to a local alignment.
    Remote requests never fall back to whole-file scans. Uses samtools -M -X;
    pysam validates the resulting BAM and creates its index. CRAM requires a
    local indexed reference FASTA, recorded by SHA256.

    The default retains only overlapping records. fetch_pairs=True additionally
    retrieves paired mates, not every supplementary alignment of each template.
    With unplaced_mates=False it leaves out mates that are unmapped and have no
    position of their own (as STAR writes them). They carry no alignment, and
    finding them means reading every unplaced read in the BAM, so leaving them out
    is quicker when the regions hold reads whose mates have none, as RNA-seq
    regions usually do; otherwise it's about a second slower, for a second read.
    Cached results are immutable snapshot derivatives; they are verified offline
    without rechecking the remote object. New requests check remote identity
    before and after extraction and retain that identity in their receipt.
    max_records stops acquisition on overflow, discards the partial output and
    raises RecordLimitError; asking again, even offline, raises it without reading.
    """
    if max_records is not None and (type(max_records) is not int or max_records < 1):
        raise ValueError("max_records must be a positive integer")
    if recovery is not None:
        if max_records is not None:
            raise ValueError("Use RecoveryPolicy.max_records with recovery")
        if not unplaced_mates:
            raise ValueError("Recovery follows mates itself; unplaced_mates=False applies only to fetch_pairs")
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
    file, location, remote = _alignment_source(source)
    if index is None:
        if file and file.index_urls:
            index = file.index_urls[0]
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
    is_cram = (file.format if file else Path(urlsplit(location).path).suffix.lstrip(".")) == "cram"
    if is_cram and reference is None:
        raise CoordinateError("CRAM extraction requires an explicit local reference FASTA")
    reference = Path(reference).resolve() if reference is not None else None
    if reference is not None and not Path(str(reference) + ".fai").is_file():
        raise CoordinateError("Reference FASTA must already have a .fai index")
    local_identities = [file_identity(p) for p in ([] if remote else [location])
                        + ([] if remote_index else [index])]
    request = dict(schema_version=1, operation="extract_reads",
                   source=location if remote else _recorded(location, cache),
                   source_sha256=None if remote else cache.file_digest(location),
                   source_size=file.size if file else None,
                   source_modified=file.modified if file else None,
                   index=index if remote_index or index is None else _recorded(index, cache),
                   index_sha256=None if remote_index else cache.file_digest(index),
                   snapshot_id=snapshot_id,
                   regions=sorted((asdict(r) for r in regions), key=lambda r: json.dumps(r, sort_keys=True)),
                   filters=asdict(filters), fetch_pairs=fetch_pairs,
                   # Only when set: requests for every mate keep the key they always had.
                   **({"unplaced_mates": False} if fetch_pairs and not unplaced_mates else {}),
                   reference_sha256=cache.file_digest(reference) if reference else None,
                   reference_index_sha256=cache.file_digest(str(reference) + ".fai") if reference else None)
    # JSON normalization makes tuples and serialized lists compare identically.
    if max_records is not None:
        request["max_records"] = max_records
    request = json.loads(json.dumps(request))
    directory = cache.workspace / "derived" / stable_id(request)
    # Reads past max_records aren't kept, but that they were too many is, so asking
    # again (even offline) gets the same answer without reading them again.
    over_limit = directory.with_name(directory.name + ".over-limit.json")
    with file_lock(cache.workspace / "locks" / (directory.name + ".lock")):
        cached = _cached_subset(directory, request)
        if cached is not None:
            return cached
        if over_limit.exists() and read_read_receipt(over_limit, cache.workspace).get("request") == request:
            raise RecordLimitError("Acquisition exceeds record limit")
        if remote and cache.offline:
            raise OfflineError("Regional reads are not cached")
        require_samtools(fetch_pairs=fetch_pairs, unplaced_mates=unplaced_mates, filters=filters)
        directory.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=directory.parent, prefix=".reads-") as temporary, \
                ThreadPoolExecutor(max_workers=1) as background:
            work = Path(temporary)
            # The remote object is checked before reading it (the request runs while the
            # header is found) and after.
            checking = background.submit(_remote_identity, location, min(timeout, 60)) if remote else None
            info = inspect_alignment(source, cache=cache, snapshot_id=snapshot_id, timeout=timeout)
            header_text = info.path.read_text()
            resolved = resolve_regions(regions, info.header)
            before = checking.result() if remote else None
            if before and file and file.size is not None and before["content-length"] is not None:
                if int(before["content-length"]) != file.size:
                    raise IntegrityError("Remote alignment size differs from the pinned inventory")
            if remote and before != info.receipt["remote_identity"]:
                raise IntegrityError("Remote alignment changed since header inspection; use a new snapshot")
            index_receipt = None
            if remote_index:
                index_receipt = cache.fetch(index, refresh=not cache.offline)
                local_index = cache.path(index_receipt)
            else:
                local_index = Path(index)
            deadline = time.monotonic() + timeout  # for every samtools run, however many
            bed = work / "regions.bed"
            bed.write_text("".join(f"{r.contig}\t{r.start}\t{r.end}\n" for r in resolved))
            output = work / "reads.bam"
            barcodes = names = None
            if filters.barcodes:
                barcodes = work / "barcodes.txt"
                barcodes.write_text("\n".join(sorted(set(filters.barcodes))) + "\n")
            if filters.query_names:
                names = query_names_asset(filters.query_names, cache.workspace)

            def left():
                return max(1, deadline - time.monotonic())

            def view(regions_bed, destination, *, name_list=None, pairs=False):
                """samtools view over these regions, with every filter; checks its output."""
                command = ["samtools", "view", "--no-PG", "-b", "-M", "-X", "-L", str(regions_bed),
                           "-q", str(filters.min_mapq), "-F", str(filters.exclude_flags),
                           "-f", str(filters.require_flags), "-o", str(destination)]
                if reference:
                    command += ["-T", str(reference)]
                if pairs:
                    command += ["--fetch-pairs"]
                if barcodes:
                    command += ["-D", filters.barcode_tag + ":" + str(barcodes)]
                if name_list:
                    command += ["-N", str(name_list)]
                command += [location, str(local_index)]
                if max_records is None:
                    _run(command, left())
                else:
                    command = _run_bounded(command, destination, max_records, left())
                _run(["samtools", "quickcheck", "-v", str(destination)], min(left(), 60))
                return command

            def acquire():
                """Read the records. Returns the receipt's commands, and the number of
                records if it's known."""
                if not fetch_pairs or unplaced_mates:
                    return dict(command=view(bed, output, name_list=names, pairs=fetch_pairs)), None
                # samtools --fetch-pairs reads the regions a second time, with the positions of
                # the mates it takes to lie outside them and all unplaced reads, keeping the
                # records of those mates' reads' names. This reads just the placed positions.
                command = view(bed, output, name_list=names)
                inside = overlap_test((r.contig, r.start, r.end) for r in resolved)
                templates, mates, count = set(), set(), 0
                with pysam.AlignmentFile(str(output)) as bam:
                    for count, read in enumerate(bam, 1):
                        if not read.is_paired:
                            continue
                        contig, position = read.next_reference_name, read.next_reference_start
                        if read.next_reference_id < 0 or position < 0:
                            templates.add(read.query_name)  # its mate, if any, has no position
                        elif not inside(contig, position, position):  # samtools: a mate starting
                            templates.add(read.query_name)  # at a region's start is outside it
                            mates.add((contig, position, position + 1))
                if not mates:  # nothing to read elsewhere
                    return dict(command=command), count
                regional = work / "regions.bam"
                os.replace(output, regional)
                (work / "templates.txt").write_text("\n".join(sorted(templates)) + "\n")
                (work / "mates.bed").write_text("".join(f"{c}\t{start}\t{end}\n" for c, start, end in merge_spans(mates)))
                mates_command = view(work / "mates.bed", work / "mates.bam", name_list=work / "templates.txt")
                count += _add_mates(regional, work / "mates.bam", inside, output)
                if max_records is not None and count > max_records:
                    raise RecordLimitError("Acquisition exceeds record limit")
                _run(["samtools", "quickcheck", "-v", str(output)], min(left(), 60))
                (work / "mates.bam").unlink()  # made again by mates_command
                return dict(command=[str(regional) if c == str(output) else c for c in command],
                            mates_command=mates_command), count

            try:
                commands, count = acquire()
            except RecordLimitError:
                try:  # remembered only if the remote object is still the one read
                    unchanged = not remote or _remote_identity(location, min(timeout, 60)) == before
                except OsteosarcError:
                    unchanged = False
                if unchanged:
                    write_read_receipt(over_limit, dict(request=request, max_records=max_records), cache.workspace)
                raise
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as error:
                if remote:  # samtools may have failed because the object changed: say so if it did
                    try:  # once: this only says why the read failed
                        changed = _remote_identity(location, min(timeout, 60), attempts=1) != before
                    except OsteosarcError:
                        changed = False
                    if changed:
                        raise IntegrityError("Remote alignment changed during extraction") from error
                raise
            if count is None:
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
            if not fetch_pairs:
                scope = "regional_records"
            elif unplaced_mates:
                scope = "regional_records_and_paired_mates"
            else:
                scope = "regional_records_and_placed_mates"
            files = {name: digest(work / name) for name in ("reads.bam", "reads.bam.bai", "header.sam", "regions.bed",
                                                            "mates.bed", "templates.txt") if (work / name).exists()}
            receipt = dict(request=request, files=files, records=count,
                           resolved_regions=[asdict(r) for r in resolved],
                           remote_identity=before,
                           header_receipt=info.receipt,
                           index_receipt=index_receipt.to_dict() if index_receipt else None,
                           samtools_version=_samtools_version(),
                           pysam_version=pysam.__version__,
                           # With placed mates only: the regions' records, and the mates'
                           # records outside the regions (by read name), in coordinate order.
                           # This follows samtools 1.21's --fetch-pairs (samtools_version).
                           **{name: _recorded_command(c, work, cache) for name, c in commands.items()},
                           scope=scope)
            write_read_receipt(work / "receipt.json", receipt, cache.workspace)
            # Only a complete directory becomes visible. No receipt means no cache hit.
            share(work)
            os.replace(work, directory)
        return _cached_subset(directory, request)
