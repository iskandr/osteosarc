"""Isovar's required records for the shared test data: every fixture made of original Sid reads.

python scripts/shared_test_data/required_isovar.py ~/code/isovar required-isovar.json.gz [--rev origin/master]

Reads the committed Isovar tree at --rev (git show, not the working tree) and
writes one subset per fixture, with its source BAM's URL and its SAM lines
(repeats kept), for osteosarc.shared.build_shared_recipe:

- every fixture in isovar/data/sid-read-recipe.json.gz (records from
  sid-reads.json.gz);
- the three-fusions fixtures (tests/data/fusions/three-fusions/*.json.gz), which
  the recipe doesn't list;
- the DLG5 supplementary partners (tests/data/osteosarc/figure_comparisons/
  corpus/dlg5.json.gz: observed_windows[].partner_sam and the split_paths
  observed_cigar pairs), which the recipe's compiler skips.

Every SAM line must come back unchanged from pysam.AlignedSegment.fromstring
and to_string with the stored header, so it matches the source BAM's own
record text. Afterwards every SAM record in Isovar's Sid fixture directories
is checked against the subsets, and records no subset covers are listed.
"""

import argparse
import gzip
import json
import subprocess
import tempfile
from collections import Counter
from hashlib import sha256
from pathlib import Path

import pysam

BUCKET = "https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/"
# Directories of fixtures made of Sid reads, audited for records no subset covers.
SID_DIRECTORIES = ("tests/data/osteosarc", "tests/data/fusions", "tests/data/chimeric", "tests/data/read_ends")
NOT_SID = ("tests/data/fusions/coding-corpus/BCR--ABL1.input.json.gz",)  # K562 external control


class Tree:
    """Files of one committed Isovar revision."""

    def __init__(self, checkout, rev):
        self.checkout, self.rev = checkout, rev
        self.commit = self._git("rev-parse", rev + "^{commit}").decode().strip()

    def _git(self, *args):
        return subprocess.run(["git", "-C", str(self.checkout), *args], capture_output=True, check=True).stdout

    def read(self, path):
        return self._git("show", f"{self.commit}:{path}")

    def json(self, path):
        raw = self.read(path)
        return json.loads(gzip.decompress(raw) if path.endswith(".gz") else raw)

    def files(self, *directories):
        return self._git("ls-tree", "-r", "--name-only", self.commit, "--", *directories).decode().split()


def round_trip(lines, header, where):
    """Require each stored line to be exactly pysam's text of that record."""
    header = pysam.AlignmentHeader.from_text(header) if isinstance(header, str) else pysam.AlignmentHeader.from_dict(header)
    for line in lines:
        if pysam.AlignedSegment.fromstring(line, header).to_string() != line:
            raise ValueError(f"{where}: SAM text of {line.split(chr(9), 1)[0]} doesn't round-trip through pysam")
    return lines


def subset(source, lines, description):
    if not source.startswith(BUCKET) or not source.endswith(".bam"):
        raise ValueError(f"Not a Sid source BAM: {source}")
    return dict(source=source, sam=list(lines), description=description)


def recipe_subsets(tree):
    reads = tree.json("isovar/data/sid-reads.json.gz")
    recipe = tree.json("isovar/data/sid-read-recipe.json.gz")
    subsets = {}
    for name, fixture in sorted(recipe["fixtures"].items()):
        source = reads["sources"][fixture["source"]]
        lines = round_trip([source["records"][digest] for digest in fixture["records"]], source["header"], name)
        subsets["isovar/" + name] = subset(source["asset"]["url"], lines, f"Isovar fixture {name} (tests: "
                                           + ", ".join(fixture["tests"]) + ")")
    return subsets


def three_fusion_subsets(tree):
    directory = "tests/data/fusions/three-fusions"
    manifest = tree.json(directory + "/manifest.json")
    subsets = {}
    for fixture in manifest["fixtures"]:
        path = f"{directory}/{fixture['file']}"
        raw = tree.read(path)
        if sha256(raw).hexdigest() != fixture["sha256"]:
            raise ValueError(f"{path} doesn't match its manifest checksum")
        data = json.loads(gzip.decompress(raw))
        lines = list(data["records"].values())
        if len(lines) != fixture["records"] or any(sha256(line.encode()).hexdigest() != digest
                                                   for digest, line in data["records"].items()):
            raise ValueError(f"{path}: records don't match their digests or the manifest count")
        name = "fusions/three-fusions/" + fixture["file"] + "#/records"
        subsets["isovar/" + name] = subset(
            manifest["sources"][fixture["source"]]["url"], round_trip(lines, data["header"], name),
            f"Isovar three-fusions fixture {fixture['file']}: {fixture['event']} in {fixture['source']}, every "
            f"record of templates touching both +/-{manifest['window']} breakpoint windows, identical lines "
            f"stored once in digest order ({fixture['segments']} segments; tests: test_sv_rna_three_fusions.py, "
            "test_source_distribution.py)")
    return subsets


def dlg5_partner_subsets(tree):
    """Supplementary partners of DLG5 split reads, kept beside the recipe's original_sam records."""
    path = "osteosarc/figure_comparisons/corpus/dlg5.json.gz"
    data = tree.json("tests/data/" + path)
    subsets = {}
    for i, product in enumerate(data["products"]):
        url, header = product["source"]["url"], product["header"]
        for j, window in enumerate(product.get("observed_windows", [])):
            name = f"{path}#/products/{i}/observed_windows/{j}/partner_sam"
            subsets["isovar/" + name] = subset(url, round_trip([window["partner_sam"]], header, name),
                                               f"Isovar fixture {name}: supplementary partner of the DLG5 split "
                                               "read in original_sam (tests: test_extended_footprints.py)")
        for k, pair in enumerate(product.get("split_paths", {}).get("observed_cigar", [])):
            name = f"{path}#/products/{i}/split_paths/observed_cigar/{k}"
            subsets["isovar/" + name] = subset(url, round_trip(list(pair), header, name),
                                               f"Isovar fixture {name}: both records of one observed DLG5 split "
                                               "path (tests: test_extended_footprints.py)")
    return subsets


def sam_strings(value):
    if isinstance(value, str):
        for line in value.split("\n"):
            fields = line.split("\t")
            if len(fields) >= 11 and fields[1].isdigit() and fields[3].isdigit():
                yield line
    elif isinstance(value, dict):
        for item in value.values():
            yield from sam_strings(item)
    elif isinstance(value, list):
        for item in value:
            yield from sam_strings(item)


def uncovered(tree, subsets):
    """SAM records in Isovar's Sid fixture files that no subset lists, by file."""
    covered = {line for s in subsets.values() for line in s["sam"]}
    missing = {}
    with tempfile.TemporaryDirectory() as scratch:
        for path in tree.files(*SID_DIRECTORIES):
            if path in NOT_SID:
                continue
            if path.endswith((".bam", ".sam", ".sam.gz")):
                local = Path(scratch) / Path(path).name
                local.write_bytes(tree.read(path))
                with pysam.AlignmentFile(str(local), check_sq=False) as handle:
                    lines = [record.to_string() for record in handle]
            elif path.endswith((".json", ".json.gz")):
                lines = list(sam_strings(tree.json(path)))
            else:
                continue
            absent = [line for line in lines if line not in covered]
            if absent:
                missing[path] = absent
    return missing


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("isovar", type=Path, help="An Isovar checkout")
    parser.add_argument("output", type=Path, help="Where to write the required records (.json.gz)")
    parser.add_argument("--rev", default="origin/master", help="The committed Isovar revision to read")
    args = parser.parse_args()
    tree = Tree(args.isovar, args.rev)
    subsets = recipe_subsets(tree)
    added = {**three_fusion_subsets(tree), **dlg5_partner_subsets(tree)}
    if set(added) & set(subsets):
        raise ValueError(f"Subsets named twice: {sorted(set(added) & set(subsets))[:3]}")
    subsets.update(added)
    text = json.dumps(dict(consumer="isovar", commit=tree.commit, subsets=dict(sorted(subsets.items()))))
    args.output.write_bytes(gzip.compress(text.encode(), mtime=0))

    def per_source(chosen):
        counts = Counter()
        for s in chosen.values():
            counts[s["source"][len(BUCKET):]] += len(s["sam"])
        return counts

    everything = per_source(subsets)
    print(f"isovar {tree.commit[:12]}: {len(subsets)} subsets, {sum(everything.values())} records, "
          f"{len(everything)} sources")
    for key, n in everything.most_common():
        print(f"  {n:6}  {key}")
    extra = per_source(added)
    print(f"added outside the recipe: {len(added)} subsets, {sum(extra.values())} records")
    for key, n in extra.most_common():
        print(f"  {n:6}  {key}")
    for name, s in added.items():
        print(f"  {len(s['sam']):4}  {name}")
    missing = uncovered(tree, subsets)
    print(f"Sid fixture records no subset covers: {sum(map(len, missing.values()))}")
    for path, lines in sorted(missing.items()):
        print(f"  {len(lines):4}  {path}  e.g. {lines[0].split(chr(9), 1)[0]}")


if __name__ == "__main__":
    main()
