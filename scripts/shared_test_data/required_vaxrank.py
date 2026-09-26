"""Vaxrank's required records: each Sid read cohort's SAM lines, with its source file.

python scripts/shared_test_data/required_vaxrank.py VAXRANK_CHECKOUT OUT.json.gz [REF]

Reads REF (default origin/main) of the checkout: the package bundle
vaxrank/data/sid-test-data.zip and its cohort recipe
examples/osteosarc_test_data/recipe/selection.json.gz. Every line is what
pysam's AlignedSegment.to_string() prints for the record. BAM cohorts are read
as stored; SAM text (the legacy SAMs and fusion inputs) is parsed back through
pysam. Each cohort's records must reproduce its recipe record digests exactly,
in order and with repeats, or nothing is written.
"""
import gzip
import io
import json
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import pysam

from osteosarc.cohort_bundle import extract_bundle, record_digest

S3 = "https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/"
root, output = sys.argv[1], Path(sys.argv[2])
ref = sys.argv[3] if len(sys.argv) > 3 else "origin/main"


def git(*args):
    return subprocess.run(["git", "-C", root, *args], capture_output=True, check=True).stdout


commit = git("rev-parse", ref + "^{commit}").decode().strip()
plan = json.loads(gzip.decompress(git("show", f"{commit}:examples/osteosarc_test_data/recipe/selection.json.gz")))
archive = io.BytesIO(git("show", f"{commit}:vaxrank/data/sid-test-data.zip"))


def role(path):
    if "/shared-v1/" in path:
        return "one reviewed template for indexed retrieval"
    if "/selection_validation/" in path:
        return "reviewed reconstruction/ranking cohort (Isovar 1.8.1 corpus bytes)"
    if path.endswith(".sam.gz"):
        return "legacy six-locus context cohort"
    return "reviewed fusion junction records"


def read_cohort(cohort, path, scratch):
    """Records as pysam segments, plus the SAM text lines they were parsed from (None for BAM)."""
    if cohort["format"] == "bam":
        with pysam.AlignmentFile(str(path)) as reads:
            return list(reads), None
    if cohort["format"] == "sam.gz":
        text = gzip.decompress(path.read_bytes()).decode()
        scratch.write_text(text)
        with pysam.AlignmentFile(str(scratch)) as reads:  # the cohort's own header
            return list(reads), [line for line in text.splitlines() if line and not line.startswith("@")]
    lines = [r["sam"] for r in json.loads(gzip.decompress(path.read_bytes()))["original_records"]]
    contigs = sorted({f for line in lines for f in (line.split("\t")[2], line.split("\t")[6])} - {"*", "="})
    header = pysam.AlignmentHeader.from_dict({"SQ": [{"SN": c, "LN": 2 ** 31 - 1} for c in contigs]})
    return [pysam.AlignedSegment.fromstring(line, header) for line in lines], lines


subsets, problems, notes = {}, [], []
with tempfile.TemporaryDirectory(prefix="vaxrank-required-") as temporary:
    bundle = Path(temporary) / "bundle"
    bundle.mkdir()
    extract_bundle(archive, bundle)  # verifies every member against bundle.json
    provenance = json.loads((bundle / "provenance.json").read_text())
    # provenance lists cohorts grouped by source; compare them by path
    if ({c["path"]: (c["source"], c["records"]) for c in provenance["cohorts"]}
            != {c["path"]: (c["source"], c["records"]) for c in plan["cohorts"]}):
        problems.append("bundle provenance and selection recipe list different records")
    for cohort in plan["cohorts"]:
        records, text = read_cohort(cohort, bundle / cohort["path"], Path(temporary) / "cohort.sam")
        sam = [r.to_string() for r in records]
        digests = [record_digest(r, text_only=cohort["format"] != "bam") for r in records]
        if digests != cohort["records"]:
            wrong = sum(a != b for a, b in zip(digests, cohort["records"])) + abs(len(digests) - len(cohort["records"]))
            problems.append(f"{cohort['path']}: {wrong} of {len(cohort['records'])} records differ from the recipe digests")
        if text is not None and text != sam:
            notes.append(f"{cohort['path']}: {sum(a != b for a, b in zip(text, sam))} lines print differently via pysam")
        if not cohort["source"].startswith(S3):
            problems.append(f"{cohort['path']}: source is not a Sid S3 object: {cohort['source']}")
        label = ", ".join(cohort["variants"]) or Path(cohort["path"]).name.removesuffix(".input.json.gz")
        subsets["vaxrank/" + cohort["path"]] = dict(
            source=cohort["source"], sam=sam, description=f"Vaxrank {role(cohort['path'])}: {label}")

print("\n".join("note: " + n for n in notes) or "SAM text cohorts print identically after the pysam round trip")
if problems:
    sys.exit("\n".join(problems))
output.parent.mkdir(parents=True, exist_ok=True)
with gzip.GzipFile(output, "wb", mtime=0) as handle:
    handle.write(json.dumps(dict(consumer="vaxrank", commit=commit, subsets=subsets), sort_keys=True).encode())
per_source = Counter()
for s in subsets.values():
    per_source[s["source"]] += len(s["sam"])
print(len(subsets), "subsets,", sum(per_source.values()), "records, all digests match; vaxrank", commit[:7])
for url, n in per_source.most_common():
    print(f"  {n:6} {url.removeprefix(S3)}")
