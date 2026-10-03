import gzip
import random
import subprocess
import sys

import pysam
import pytest

from osteosarc import IntegrityError, RecordLimitError
from osteosarc.reads import _run_bounded
from osteosarc.records import record_multiset


@pytest.fixture
def large_header_bam(tmp_path):
    rng = random.Random(106)
    header = dict(SQ=[dict(SN="chr1", LN=10000)],
                  CO=["".join(rng.choices("ACGT0123456789", k=400)) for _ in range(1300)])
    path = tmp_path / "input.bam"
    with pysam.AlignmentFile(path, "wb", header=header) as out:
        for i in range(3):
            read = pysam.AlignedSegment(out.header)
            read.query_name = "record-%d" % i
            read.reference_id, read.reference_start = 0, i * 100
            read.query_sequence, read.cigarstring = "ACGTACGTAA", "10M"
            read.query_qualities = [i + 10] * 10 if i else None
            read.set_tag("xf", 1.23456789, value_type="f")
            read.set_tag("xi", 12, value_type="I")
            read.set_tag("CB", "cell")
            out.write(read)
    return path


def producer(path, script):
    return [sys.executable, "-c", script, str(path), "-o", "-"]


def test_mid_header_timeout_has_no_seek_cleanup_error(large_header_bam, tmp_path, monkeypatch):
    unraisable = []
    monkeypatch.setattr(sys, "unraisablehook", lambda error: unraisable.append(str(error.exc_value)))
    script = ("import pathlib,sys,time; "
              "sys.stdout.buffer.write(pathlib.Path(sys.argv[1]).read_bytes()[:17400]); "
              "sys.stdout.buffer.flush(); time.sleep(10)")
    output = tmp_path / "out.bam"
    with pytest.raises(subprocess.TimeoutExpired):
        _run_bounded(producer(large_header_bam, script), output, 10, 0.5)
    assert not unraisable
    # A failed stream must not interfere with the next acquisition.
    _run_bounded(["samtools", "view", "--no-PG", "-b", "-o", "-", str(large_header_bam)],
                 output, 3, 10)
    assert record_multiset(output) == record_multiset(large_header_bam)


@pytest.mark.parametrize("limit", [1, 2, 3])
def test_record_caps_preserve_stored_values(large_header_bam, tmp_path, monkeypatch, limit):
    unraisable = []
    monkeypatch.setattr(sys, "unraisablehook", lambda error: unraisable.append(str(error.exc_value)))
    command = ["samtools", "view", "--no-PG", "-b", "-o", "-", str(large_header_bam)]
    output = tmp_path / "out.bam"
    if limit < 3:
        with pytest.raises(RecordLimitError):
            _run_bounded(command, output, limit, 10)
    else:
        _run_bounded(command, output, limit, 10)
        assert record_multiset(output) == record_multiset(large_header_bam)
    assert not unraisable


@pytest.mark.parametrize("compressed", [False, True])
def test_truncated_blocks_are_integrity_errors(large_header_bam, tmp_path, compressed):
    damaged = tmp_path / "damaged.bam"
    if compressed:
        damaged.write_bytes(large_header_bam.read_bytes()[:17400])
    else:
        raw = gzip.decompress(large_header_bam.read_bytes())
        damaged.write_bytes(gzip.compress(raw[:-1]))
    script = "import pathlib,sys; sys.stdout.buffer.write(pathlib.Path(sys.argv[1]).read_bytes())"
    with pytest.raises(IntegrityError, match="Truncated"):
        _run_bounded(producer(damaged, script), tmp_path / "out.bam", 10, 10)


def test_producer_failure_keeps_stderr_after_partial_header(large_header_bam, tmp_path):
    script = ("import pathlib,sys; "
              "sys.stdout.buffer.write(pathlib.Path(sys.argv[1]).read_bytes()[:17400]); "
              "sys.stdout.buffer.flush(); sys.stderr.write('remote source failed'); sys.exit(7)")
    with pytest.raises(subprocess.CalledProcessError) as caught:
        _run_bounded(producer(large_header_bam, script), tmp_path / "out.bam", 10, 10)
    assert caught.value.returncode == 7
    assert caught.value.stderr == b"remote source failed"
