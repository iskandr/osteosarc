import os
import stat
import time
import warnings
from pathlib import Path

import pytest

from osteosarc import saved


def test_a_catalogue_is_built_once(tmp_path, monkeypatch):
    builds = []

    def build():
        builds.append(1)
        return {"files": [1, 2, 3]}
    for _ in range(2):
        assert saved.load_or_build(tmp_path, "files", "snapshot", build) == {"files": [1, 2, 3]}
    assert len(builds) == 1
    (path,) = (tmp_path / f"u{os.getuid()}").glob("*.pickle.gz")
    assert stat.S_IMODE(path.parent.stat().st_mode) == 0o700 and path.name.startswith("files-snapshot-")
    # A copy that can't be marked as used is still a good copy.
    monkeypatch.setattr(os, "utime", lambda *a, **k: (_ for _ in ()).throw(PermissionError()))
    assert saved.load_or_build(tmp_path, "files", "snapshot", build) == {"files": [1, 2, 3]}
    assert len(builds) == 1
    # With no key for what it's built from, nothing is saved.
    assert saved.load_or_build(tmp_path, "variants", None, lambda: 5) == 5
    assert not list(path.parent.glob("variants-*"))


def test_damaged_linked_or_misnamed_copies_are_built_again(tmp_path):
    saved.load_or_build(tmp_path, "files", "a", lambda: 1)
    (path,) = (tmp_path / f"u{os.getuid()}").glob("*.pickle.gz")
    good = path.read_bytes()
    path.write_bytes(b"not a pickle")
    assert saved.load_or_build(tmp_path, "files", "a", lambda: 2) == 2
    # A link isn't loaded, even to a good copy, nor is another catalogue under this one's name.
    elsewhere = tmp_path / "elsewhere"
    elsewhere.write_bytes(good)
    path.unlink()
    path.symlink_to(elsewhere)
    assert saved.load_or_build(tmp_path, "files", "a", lambda: 4) == 4
    other = path.with_name(path.name.replace("files-a-", "files-b-"))
    path.rename(other)
    assert saved.load_or_build(tmp_path, "files", "b", lambda: 5) == 5
    # A private folder made too open is made private again, and still used.
    path.parent.chmod(0o777)
    assert saved.load_or_build(tmp_path, "files", "b", lambda: 6) == 5
    assert stat.S_IMODE(path.parent.stat().st_mode) == 0o700


def test_a_read_only_cache_still_works(tmp_path, monkeypatch):
    monkeypatch.setattr(saved, "write_own", lambda *a: (_ for _ in ()).throw(PermissionError()))
    assert saved.load_or_build(tmp_path, "files", "a", lambda: 4) == 4
    assert not list((tmp_path / f"u{os.getuid()}").iterdir())


def test_each_kind_of_catalogue_keeps_a_few_copies(tmp_path, monkeypatch):
    folder = tmp_path / f"u{os.getuid()}"
    for n in range(saved.KEEP + 2):  # other snapshots
        saved.load_or_build(tmp_path, "files", f"snapshot{n}", lambda n=n: n)
        recent = time.time() - 60 * (saved.KEEP + 2 - n)  # used a minute apart, most recent last
        os.utime(next(folder.glob(f"files-snapshot{n}-*")), (recent, recent))
    saved.load_or_build(tmp_path, "header", "snapshot0", lambda: 0)  # another kind
    kept = sorted(p.name.split("-")[1] for p in folder.glob("files-*.pickle.gz"))
    assert kept == [f"snapshot{n}" for n in range(2, saved.KEEP + 2)]
    # One unused for more than 30 days goes, whatever its kind, as do leftovers of interrupted saves.
    old, leftover = folder / "variants-x-code.pickle.gz", folder / ".own-x"
    for path in (old, leftover):
        path.write_bytes(b"")
        os.utime(path, (1_000_000_000, 1_000_000_000))
    saved.load_or_build(tmp_path, "header", "snapshot1", lambda: 0)
    assert not old.exists() and not leftover.exists() and len(list(folder.glob("header-*"))) == 2
    monkeypatch.setattr(saved, "code_key", lambda: None)  # no code to fingerprint: nothing saved
    assert saved.load_or_build(tmp_path / "none", "files", "s", lambda: 9) == 9
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
    # Every kind loads with only the catalogue classes, and the header, built from no
    # corrected source, serves any corrections.
    header, variants = dataset._download_header, list(dataset.variants("all"))
    monkeypatch.setattr(ds, "parse_variants", lambda *a, **k: pytest.fail("rebuilt the variant catalogue"))
    monkeypatch.setattr(ds.Dataset, "_raw_download_header", lambda self: pytest.fail("rebuilt the header"))
    assert list(reopened.variants("all")) == variants
    assert Dataset.open(dataset.name, cache=dataset.cache, corrections=False)._download_header == header


def test_a_copy_naming_any_other_class_is_never_loaded(tmp_path):
    import fractions
    import gzip
    import pickle
    folder = tmp_path / f"u{os.getuid()}"
    saved.load_or_build(tmp_path, "files", "a", lambda: 1)
    (path,) = folder.glob("files-a-*")
    stem = path.name[:-len(".pickle.gz")]
    path.write_bytes(gzip.compress(pickle.dumps((stem, fractions.Fraction(1, 2)))))
    assert saved.load_or_build(tmp_path, "files", "a", lambda: 2) == 2


def test_a_catalogue_that_could_not_be_loaded_is_not_kept(tmp_path):
    import fractions
    with pytest.warns(RuntimeWarning, match="can't keep its files catalogue"):
        assert saved.load_or_build(tmp_path, "files", "a", lambda: fractions.Fraction(1, 2)) == 0.5
    assert not list((tmp_path / f"u{os.getuid()}").iterdir())


def test_the_code_fingerprint_passes_over_editors_files(tmp_path, monkeypatch):
    (tmp_path / "a.py").write_text("x = 1")
    (tmp_path / ".#a.py").symlink_to(tmp_path / "gone")  # an editor's lock
    monkeypatch.setattr(saved, "__file__", str(tmp_path / "saved.py"))
    assert saved._fingerprint() is not None
    (tmp_path / "b.py").symlink_to(tmp_path / "gone")
    assert saved._fingerprint() is None  # unreadable: nothing is saved


def test_indexed_correction_matching_agrees_with_a_full_scan(dataset, monkeypatch):
    from osteosarc import Dataset
    from osteosarc.curation import Curation
    fast = dataset.curation.report()
    monkeypatch.setattr(Curation, "_candidates", lambda self, source, match: range(len(self.raw(source))))
    slow = Dataset.open(dataset.name, cache=dataset.cache).curation.report()
    assert fast == slow


def test_a_saved_catalogue_warns_as_a_fresh_build_does(dataset, monkeypatch):
    import osteosarc.dataset as ds
    from osteosarc import CORRECTIONS, Dataset
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
    monkeypatch.setattr(ds, "build_files", real_build)
    # Warnings made errors are raised on every try, whether the catalogue is loaded or built.
    for data in (Dataset.open(dataset.name, cache=dataset.cache), Dataset.open(dataset.name, cache=dataset.cache,
                                                                                corrections=list(CORRECTIONS))):
        with warnings.catch_warnings():
            warnings.simplefilter("error", CurationWarning)
            for _ in range(2):
                with pytest.raises(CurationWarning):
                    data.files
    monkeypatch.setattr(saved, "code_key", lambda: None)  # nothing saved: built
    built = warned(Dataset.open(dataset.name, cache=dataset.cache))
    assert loaded == built and loaded


def test_building_catalogues_warns_only_through_what_is_saved(dataset, monkeypatch):
    # Anything else a build warned of would be lost once the catalogue is saved:
    # builders return their problems instead (see Dataset._saved).
    from osteosarc import Dataset
    from osteosarc.curation import CurationWarning
    monkeypatch.setattr(saved, "code_key", lambda: None)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        data = Dataset.open(dataset.name, cache=dataset.cache)
        data.files, data.variants("all"), data._download_header
    assert caught and all(w.category is CurationWarning for w in caught)


def test_problems_found_building_a_catalogue_are_warned_of_when_it_is_loaded(dataset):
    from osteosarc import CORRECTIONS, Change, Correction, Dataset
    broken = Correction("broken", "simulate a registry row without an ID", (
        Change("specimens", {"sample_id": "T3_tumor"}, set={"sample_id": ""}),))
    for _ in range(2):  # built, then loaded
        with pytest.warns(UserWarning, match="aren't linked to samples"):
            Dataset.open(dataset.name, cache=dataset.cache, corrections=[*CORRECTIONS, broken]).files


def test_the_corrections_key_is_the_same_in_every_process(dataset):
    import subprocess
    import sys

    from osteosarc import CORRECTIONS, Change, Correction, Dataset
    script = ("from osteosarc import CORRECTIONS, Change, Correction\n"
              "from osteosarc.dataset import Dataset\n"
              "c = Correction('tags', 'x', (Change('bams', {'a': 1}, set={'tags': frozenset('abcdefgh')}),))\n"
              "d = Dataset.__new__(Dataset)\n"
              "d.curation = type('C', (), dict(corrections=(*CORRECTIONS, c), enabled=True))()\n"
              "print(d._corrections_key)")
    keys = {subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, check=True,
                           env={**os.environ, "PYTHONHASHSEED": seed}).stdout for seed in ("1", "2", "3")}
    assert len(keys) == 1 and "None" not in keys
    # Values JSON can't write have no such key, so nothing is saved for them.
    custom = Correction("custom", "x", (Change("bams", {"a": 1}, set={"parse": lambda v: v}),))
    assert Dataset.open(dataset.name, cache=dataset.cache,
                        corrections=[*CORRECTIONS, custom])._corrections_key is None


def test_catalogues_list_the_sources_their_builds_correct(dataset, monkeypatch):
    # A catalogue warns of stale corrections to the sources it names (Dataset._saved);
    # its build must use no others.
    import osteosarc.dataset as ds
    from osteosarc import Dataset
    from osteosarc.curation import Curation
    declared, used = {}, set()
    real_saved, real_records = ds.Dataset._saved, Curation.records

    def saved_noting(self, name, build, sources):
        declared[name] = set(sources)
        return real_saved(self, name, build, sources)
    monkeypatch.setattr(ds.Dataset, "_saved", saved_noting)
    monkeypatch.setattr(Curation, "records", lambda self, source: used.add(source) or real_records(self, source))
    monkeypatch.setattr(saved, "code_key", lambda: None)
    for build in (lambda d: d.files, lambda d: d.variants("all"), lambda d: d._download_header):
        used.clear()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            build(Dataset.open(dataset.name, cache=dataset.cache))
        name = {"files", "variants", "header"} & set(declared)
        assert used <= set().union(*(declared[n] for n in name)), used
        declared.clear()


def test_saved_warnings_come_from_osteosarc(dataset):
    from osteosarc import Dataset
    from osteosarc.curation import CurationWarning
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        Dataset.open(dataset.name, cache=dataset.cache).files  # loaded
    files = {Path(w.filename).name for w in caught if w.category is CurationWarning}
    assert files == {"dataset.py"}


def test_corrections_differing_only_in_type_have_different_keys(dataset):
    from osteosarc import CORRECTIONS, Change, Correction, Dataset, glob

    def key(value, match="x*"):
        custom = Correction("custom", "x", (Change("bams", {"a": match}, set={"tags": value}),))
        return Dataset.open(dataset.name, cache=dataset.cache, corrections=[*CORRECTIONS, custom])._corrections_key
    assert len({key(("a",)), key(["a"]), key({1: "a"}), key({"1": "a"}), key(1), key("1"),
                key({"a": 1, "b": 2}), key({"b": 2, "a": 1})}) == 8  # changes are made in order
    assert key(1, glob("x*")) != key(1, {"pattern": "x*"})


def test_reporting_corrections_does_not_keep_the_bucket_listing(dataset):
    from osteosarc import Dataset
    data = Dataset.open(dataset.name, cache=dataset.cache)
    data.files  # loaded, so the bucket's corrections weren't evaluated
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        assert len(data.corrections)
    assert "bucket" not in data.curation._raw
