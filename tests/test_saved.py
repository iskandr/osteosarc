import os
import time
import warnings

import pytest

from osteosarc import saved


def test_a_catalogue_is_built_once(tmp_path, monkeypatch):
    builds = []

    def build():
        builds.append(1)
        return {"files": [1, 2, 3]}
    for _ in range(2):
        assert saved.load_or_build(tmp_path, "snapshot-files-c", build) == {"files": [1, 2, 3]}
    assert len(builds) == 1
    (path,) = tmp_path.glob("*.pickle.gz")
    assert path.name.endswith(f"-u{os.getuid()}.pickle.gz")
    # A copy that can't be marked as used is still a good copy.
    monkeypatch.setattr(os, "utime", lambda *a: (_ for _ in ()).throw(PermissionError()))
    assert saved.load_or_build(tmp_path, "snapshot-files-c", build) == {"files": [1, 2, 3]}
    assert len(builds) == 1


def test_damaged_or_foreign_copies_are_built_again(tmp_path):
    saved.load_or_build(tmp_path, "s", lambda: 1)
    (path,) = tmp_path.glob("*.pickle.gz")
    good = path.read_bytes()
    path.write_bytes(b"not a pickle")
    assert saved.load_or_build(tmp_path, "s", lambda: 2) == 2
    # A copy others can write could have been replaced: it's never loaded.
    path.chmod(0o666)
    assert saved.load_or_build(tmp_path, "s", lambda: 3) == 3
    # Nor is a link, even to a good copy.
    elsewhere = tmp_path / "elsewhere"
    elsewhere.write_bytes(good)
    path.unlink()
    path.symlink_to(elsewhere)
    assert saved.load_or_build(tmp_path, "s", lambda: 4) == 4


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


def test_a_saved_catalogue_warns_as_a_fresh_build_does(dataset, monkeypatch):
    import osteosarc.dataset as ds
    from osteosarc import Dataset
    from osteosarc.curation import CurationWarning

    def warned(data):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            data.files
        return sorted(str(w.message) for w in caught if w.category is CurationWarning)
    # Syncing built the variants before the files; corrections they share are still the files'.
    real_build = ds.build_files
    monkeypatch.setattr(ds, "build_files", lambda *a, **k: pytest.fail("rebuilt the file catalogue"))
    loaded = warned(Dataset.open(dataset.name, cache=dataset.cache))
    with warnings.catch_warnings():
        warnings.simplefilter("error", CurationWarning)  # an error, from the saved copy
        with pytest.raises(CurationWarning):
            Dataset.open(dataset.name, cache=dataset.cache).files
    monkeypatch.setattr(ds, "build_files", real_build)
    monkeypatch.setattr(saved, "code_key", lambda: None)  # nothing saved: built
    built = warned(Dataset.open(dataset.name, cache=dataset.cache))
    assert loaded == built and loaded


def test_problems_found_building_a_catalogue_are_warned_of_when_it_is_loaded(dataset):
    from osteosarc import CORRECTIONS, Change, Correction, Dataset
    broken = Correction("broken", "simulate a registry row without an ID", (
        Change("specimens", {"sample_id": "T3_tumor"}, set={"sample_id": ""}),))
    for _ in range(2):  # built, then loaded
        with pytest.warns(UserWarning, match="aren't linked to samples"):
            Dataset.open(dataset.name, cache=dataset.cache, corrections=[*CORRECTIONS, broken]).files


def test_the_corrections_key_is_the_same_in_every_process():
    import subprocess
    import sys
    script = ("from osteosarc import CORRECTIONS, Change, Correction, Dataset\n"
              "from osteosarc.dataset import Dataset\n"
              "c = Correction('tags', 'x', (Change('bams', {'a': 1}, set={'tags': frozenset('abcdefgh')}),))\n"
              "d = Dataset.__new__(Dataset)\n"
              "d.curation = type('C', (), dict(corrections=(*CORRECTIONS, c), enabled=True))()\n"
              "print(d._corrections_key)")
    keys = {subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, check=True,
                           env={**os.environ, "PYTHONHASHSEED": seed}).stdout for seed in ("1", "2", "3")}
    assert len(keys) == 1
