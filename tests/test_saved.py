import os
import time
import warnings

import pytest

from osteosarc import saved


def test_a_catalogue_is_built_once_and_its_warnings_raised_every_time(tmp_path):
    builds = []

    def build():
        builds.append(1)
        warnings.warn("stale correction", UserWarning)
        return {"files": [1, 2, 3]}
    for _ in range(2):
        with pytest.warns(UserWarning, match="stale correction"):
            assert saved.load_or_build(tmp_path, "snapshot-files-c", build) == {"files": [1, 2, 3]}
    assert len(builds) == 1
    (path,) = tmp_path.glob("*.pickle.gz")
    assert path.name.endswith(f"-u{os.getuid()}.pickle.gz")


def test_damaged_or_foreign_copies_are_built_again(tmp_path):
    saved.load_or_build(tmp_path, "s", lambda: 1)
    (path,) = tmp_path.glob("*.pickle.gz")
    path.write_bytes(b"not a pickle")
    assert saved.load_or_build(tmp_path, "s", lambda: 2) == 2
    # A copy others can write could have been replaced: it's never loaded.
    path.chmod(0o666)
    assert saved.load_or_build(tmp_path, "s", lambda: 3) == 3


def test_a_read_only_cache_still_works(tmp_path):
    folder = tmp_path / "catalogs"
    folder.mkdir()
    folder.chmod(0o500)
    try:
        assert saved.load_or_build(folder, "s", lambda: 4) == 4
        assert not list(folder.iterdir())
    finally:
        folder.chmod(0o700)


def test_each_catalogue_keeps_a_few_versions(tmp_path, monkeypatch):
    for n in range(saved.KEEP + 2):
        monkeypatch.setattr(saved, "code_key", lambda n=n: f"code{n}")
        saved.load_or_build(tmp_path, "s", lambda: n)
        recent = time.time() - 60 * (saved.KEEP + 2 - n)  # used a minute apart, most recent last
        os.utime(next(tmp_path.glob(f"s-code{n}-*")), (recent, recent))
    kept = sorted(p.name.split("-")[1] for p in tmp_path.glob("s-*.pickle.gz"))
    assert kept == [f"code{n}" for n in range(2, saved.KEEP + 2)]
    # One unused for more than 30 days goes, whatever its catalogue.
    old = tmp_path / f"other-code0-u{os.getuid()}.pickle.gz"
    old.write_bytes(b"")
    os.utime(old, (1_000_000_000, 1_000_000_000))
    saved.load_or_build(tmp_path, "t", lambda: 0)
    assert not old.exists()
    monkeypatch.setattr(saved, "code_key", lambda: None)  # no code to fingerprint: nothing saved
    assert saved.load_or_build(tmp_path / "none", "s", lambda: 9) == 9
    assert not (tmp_path / "none").exists()


def test_opening_a_snapshot_twice_reuses_its_catalogues(dataset, monkeypatch):
    import osteosarc.dataset as ds
    from osteosarc import Dataset
    key = dataset.files[0].key
    reopened = Dataset.open(dataset.name, cache=dataset.cache)
    monkeypatch.setattr(ds, "build_files", lambda *a, **k: pytest.fail("rebuilt the file catalogue"))
    assert reopened.file(key) == dataset.file(key)
    assert list(reopened.files) == list(dataset.files)
    sample = next(s for s in dataset.samples if s.files)
    assert list(reopened.samples[sample.id].files) == list(dataset.files.select(sample=sample.id))


def test_indexed_correction_matching_agrees_with_a_full_scan(dataset, monkeypatch):
    from osteosarc import Dataset
    from osteosarc.curation import Curation
    fast = dataset.curation.report()
    monkeypatch.setattr(Curation, "_candidates", lambda self, source, match: range(len(self.raw(source))))
    slow = Dataset.open(dataset.name, cache=dataset.cache).curation.report()
    assert fast == slow
