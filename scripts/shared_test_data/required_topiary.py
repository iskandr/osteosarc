"""Topiary's required records: every fixture of real Sid reads, with its source BAM.

python scripts/shared_test_data/required_topiary.py ~/code/topiary OUT.json.gz [REV]

Reads the committed tree of REV (default origin/master) through git, never the
working tree, and changes nothing in the checkout. Each SAM line is pysam's
AlignedSegment.to_string() for the record, repeats kept. Covers the recipe's
13 read files (tests/data/sid-fixtures.json), isovar_repeats/reads.bam and each
rearrangement path embedded in osteosarc_rearrangements/*.json.gz.
"""
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import pysam

root, output = Path(sys.argv[1]).expanduser(), sys.argv[2]
rev = sys.argv[3] if len(sys.argv) > 3 else "origin/master"
commit = subprocess.run(["git", "-C", str(root), "rev-parse", rev + "^{commit}"],
                        check=True, capture_output=True, text=True).stdout.strip()


def blob(path):
    """Bytes of tests/data/<path> at the pinned commit."""
    return subprocess.run(["git", "-C", str(root), "cat-file", "blob", f"{commit}:tests/data/{path}"],
                          check=True, capture_output=True).stdout


manifest = json.loads(blob("manifest.json"))
pinned = {a["filename"]: a["sha256"] for a in manifest["assets"]}


def checked(path, sha256=None):
    data = blob(path)
    expected = sha256 or pinned[path]
    if hashlib.sha256(data).hexdigest() != expected:
        raise ValueError(f"{path} differs from its pinned SHA-256")
    return data


def records(path, data, work):
    """to_string() of every record in a BAM, SAM or gzip-compressed SAM, in file order."""
    local = Path(work) / Path(path).name.removesuffix(".gz")
    local.write_bytes(gzip.decompress(data) if path.endswith(".sam.gz") else data)
    with pysam.AlignmentFile(str(local), check_sq=False) as handle:
        lines = [read.to_string() for read in handle.fetch(until_eof=True)]
        header = handle.header
    return lines, header


def multiset_sha256(lines):
    """The digest osteosarc.regional_corpus records as read_selections[...]["sam_records_sha256"]."""
    return hashlib.sha256(json.dumps(sorted(Counter(lines).items()), separators=(",", ":")).encode()).hexdigest()


DESCRIPTIONS = {
    "osteosarc/": "six-locus RNA: pinned template selection (Isovar 1.8.0) at six vaccine loci",
    "osteosarc_indels/": "GLIS3/KTN1 indel reads: every record within 2 bases of the allele",
    "osteosarc_all_variants/": "184-entry T2 audit: every record within 2 bases of each ready allele",
    "osteosarc_rna_overlay/": "pVAC RNA overlay: every record within 2 bases of the 20 panel alleles",
    "osteosarc_shared/": "shared NTF3 case from Isovar's vaccine-rna-v1 corpus",
}
recipe = json.loads(checked("sid-fixtures.json"))
subsets, notes = {}, []
with tempfile.TemporaryDirectory() as work:
    ntf3_header = None
    for selection in recipe["reads"]:
        path = selection["filename"]
        lines, header = records(path, checked(path), work)
        expected = manifest["read_selections"][path]
        if len(lines) != expected["selected_records"] or multiset_sha256(lines) != expected["sam_records_sha256"]:
            raise ValueError(f"{path}: records differ from tests/data/manifest.json read_selections")
        if path.startswith("osteosarc_shared/"):
            ntf3_header = header
        subsets["topiary/" + path] = dict(
            source=selection["original_alignment_url"], sam=lines,
            description=next(d for prefix, d in DESCRIPTIONS.items() if path.startswith(prefix)))
    urls = {s["filename"]: s["original_alignment_url"] for s in recipe["reads"]}

    # A reviewed subset of the T2 audit BAM by read name, so its source is the audit's.
    repeats = json.loads(checked("isovar_repeats/recipe.json", json.loads(blob("isovar_repeats/manifest.json"))["recipe_sha256"]))
    files = json.loads(blob("isovar_repeats/manifest.json"))["files"]
    lines, _ = records("isovar_repeats/reads.bam", checked("isovar_repeats/reads.bam", files["reads.bam"]), work)
    if pinned[repeats["source"]] != repeats["source_sha256"]:
        raise ValueError("isovar_repeats no longer names the pinned T2 audit BAM")
    subsets["topiary/isovar_repeats/reads.bam"] = dict(
        source=urls[repeats["source"]], sam=lines,
        description="MT-ND5 chrM:12994 repeated-translation regression: 419 read names from the T2 audit BAM")

    # Rearrangement paths: SAM text of a primary alignment and its supplementary partner.
    group = "osteosarc_rearrangements/"
    rearrangements = json.loads(checked(group + "manifest.json"))
    observed = {(p["source"], p["read_id"]): p["sam_sha256"] for p in rearrangements["observed_paths"]}
    used = set()
    for entry in rearrangements["inputs"]:
        path = group + entry["file"]
        document = json.loads(gzip.decompress(checked(path, entry["sha256"])))
        product = entry["sample"] + "-ONT-tagged"
        source = rearrangements["source_receipts"][product]["url"]
        if document["fusion"]["provenance"]["source"] != source:
            raise ValueError(f"{path}: provenance source differs from the manifest")
        if len(document["reads"]) != len(document["original_records"]):
            raise ValueError(f"{path}: reads and original_records differ in length")
        for read, original in zip(document["reads"], document["original_records"]):
            sam = [original["sam"], original["partner_sam"]]
            rendered = [pysam.AlignedSegment.fromstring(line, ntf3_header).to_string() for line in sam]
            if rendered != sam:
                raise ValueError(f"{path}#{read['read_id']}: stored SAM text isn't pysam's rendering")
            if read["source"] != source or {line.split("\t", 1)[0] for line in sam} != {read["read_id"]}:
                raise ValueError(f"{path}#{read['read_id']}: path source or read name disagrees")
            # The manifest lists a path's records in junction order, not primary first.
            if sorted(hashlib.sha256(line.encode()).hexdigest() for line in sam) != sorted(observed[product, read["read_id"]]):
                raise ValueError(f"{path}#{read['read_id']}: records differ from the manifest's observed path")
            used.add((product, read["read_id"]))
            name = f"topiary/{path}#{read['read_id']}"
            if name in subsets:
                raise ValueError(f"{name} is named twice")
            subsets[name] = dict(source=source, sam=sam, description=(
                f"{entry['event']} {entry['sample']} junction path: primary record and supplementary partner"))
    for product, read_id in sorted(set(observed) - used):
        notes.append(f"{group}manifest.json observed path {product} {read_id}: only SHA-256 digests, no SAM text")

with gzip.open(output, "wt") as handle:
    json.dump(dict(consumer="topiary", commit=commit, subsets=subsets), handle)
per_source = Counter()
for s in subsets.values():
    per_source[s["source"]] += len(s["sam"])
print(f"topiary {commit[:12]}: {len(subsets)} subsets, {sum(per_source.values())} records")
for url, n in per_source.most_common():
    print(f"  {n:6} {url}")
for note in notes:
    print("  not included:", note)
