"""Tests for derived/ cache pruning, staging leftover sweeps, OPENVAX_OFFLINE, and cache CLI (#119)."""

from __future__ import annotations

import json
import os
import time
from pathlib import Path

from osteosarc.cache import Cache, prune_derived, sweep_staging_leftovers
from osteosarc.cli import main


def _set_mtime(path: Path, mtime: float) -> None:
    os.utime(path, (mtime, mtime))


def test_sweep_staging_leftovers_removes_old_and_keeps_fresh(tmp_path: Path):
    derived = tmp_path / "derived"
    derived.mkdir()
    now = time.time()

    old_reads = derived / ".reads-killed123"
    old_reads.mkdir()
    (old_reads / "reads.bam").write_bytes(b"x" * 100)
    _set_mtime(old_reads, now - 48 * 3600)

    fresh_reads = derived / ".reads-active456"
    fresh_reads.mkdir()
    (fresh_reads / "reads.bam").write_bytes(b"y" * 50)
    _set_mtime(fresh_reads, now - 60)

    old_json = derived / ".json-dead"
    old_json.write_text("{}", encoding="utf-8")
    _set_mtime(old_json, now - 48 * 3600)

    normal = derived / "abcd1234"
    normal.mkdir()
    (normal / "reads.bam").write_bytes(b"z" * 20)
    _set_mtime(normal, now - 48 * 3600)

    dry = sweep_staging_leftovers(derived, older_than_seconds=24 * 3600, dry_run=True, now=now)
    assert {item["name"] for item in dry} == {".reads-killed123", ".json-dead"}
    assert old_reads.exists() and old_json.exists()

    removed = sweep_staging_leftovers(derived, older_than_seconds=24 * 3600, dry_run=False, now=now)
    assert {item["name"] for item in removed} == {".reads-killed123", ".json-dead"}
    assert not old_reads.exists()
    assert not old_json.exists()
    assert fresh_reads.exists()
    assert normal.exists()


def test_prune_derived_age_and_max_bytes_dry_run_then_apply(tmp_path: Path):
    ws = tmp_path / "osteosarc"
    derived = ws / "derived"
    derived.mkdir(parents=True)
    now = time.time()

    old_entry = derived / "old_entry"
    old_entry.mkdir()
    (old_entry / "reads.bam").write_bytes(b"a" * 500)
    _set_mtime(old_entry, now - 45 * 86400)

    mid_entry = derived / "mid_entry"
    mid_entry.mkdir()
    (mid_entry / "reads.bam").write_bytes(b"b" * 400)
    _set_mtime(mid_entry, now - 10 * 86400)

    new_entry = derived / "new_entry"
    new_entry.mkdir()
    (new_entry / "reads.bam").write_bytes(b"c" * 300)
    _set_mtime(new_entry, now - 1 * 86400)

    leftover = derived / ".reads-crash"
    leftover.mkdir()
    (leftover / "partial.bam").write_bytes(b"d" * 100)
    _set_mtime(leftover, now - 48 * 3600)

    # Dry-run with days=30 & max_bytes=350 should select leftover + old_entry (age) + mid_entry (size cap).
    plan = prune_derived(ws, days=30, max_bytes=350, dry_run=True, now=now)
    assert plan["dry_run"] is True
    assert {c["name"] for c in plan["candidates"]} == {".reads-crash", "old_entry", "mid_entry"}
    assert plan["deleted"] == []
    assert plan["kept_count"] == 1
    assert plan["kept_bytes"] == 300
    assert plan["reclaimable_bytes"] == 1000
    assert old_entry.exists() and mid_entry.exists() and leftover.exists()

    # Apply removes candidates and keeps newest entry.
    applied = prune_derived(ws, days=30, max_bytes=350, dry_run=False, now=now)
    assert applied["dry_run"] is False
    assert {d["name"] for d in applied["deleted"]} == {".reads-crash", "old_entry", "mid_entry"}
    assert applied["freed_bytes"] == 1000
    assert not old_entry.exists()
    assert not mid_entry.exists()
    assert not leftover.exists()
    assert new_entry.exists()


def test_cli_cache_dirs_is_read_only_and_honors_openvax_offline(monkeypatch, tmp_path: Path, capsys):
    root = tmp_path / "shared-cache"
    monkeypatch.setenv("OPENVAX_OFFLINE", "1")
    monkeypatch.delenv("DATACACHE_OFFLINE", raising=False)

    assert Cache(root).offline is True
    assert Cache(root, offline=False).offline is False

    rc = main(["--cache", str(root), "cache", "dirs", "--json"])
    assert rc == 0
    assert not root.exists(), "cache dirs must remain strictly read-only"
    payload = json.loads(capsys.readouterr().out)
    assert payload["package"] == "osteosarc"
    assert payload["root"] == str(root.resolve())
    assert payload["exists"] is False
    assert payload["offline"] is True
    assert payload["offline_source"] == "OPENVAX_OFFLINE=1"

    rc2 = main(["--cache", str(root), "cache", "dirs"])
    assert rc2 == 0
    out2 = capsys.readouterr().out
    assert "derived" in out2
    assert "enabled" in out2
    assert not root.exists()


def test_cli_cache_prune_dry_run_and_yes(tmp_path: Path, capsys):
    root = tmp_path / "cache"
    derived = root / "osteosarc" / "derived"
    entry = derived / "abc123"
    entry.mkdir(parents=True)
    (entry / "reads.bam").write_bytes(b"x" * 2048)

    rc = main(["--cache", str(root), "cache", "prune", "--all", "--json"])
    assert rc == 0
    plan = json.loads(capsys.readouterr().out)
    assert plan["dry_run"] is True
    assert len(plan["candidates"]) == 1
    assert entry.exists()

    rc2 = main(["--cache", str(root), "cache", "prune", "--all", "--yes", "--json"])
    assert rc2 == 0
    applied = json.loads(capsys.readouterr().out)
    assert applied["dry_run"] is False
    assert len(applied["deleted"]) == 1
    assert applied["freed_bytes"] == 2048
    assert not entry.exists()
