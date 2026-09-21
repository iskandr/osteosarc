"""Recheck issue #5's whole-source checksums and UCSC mapping evidence.

    python scripts/check_allele_sources.py --cache .cache/allele-resolution

Downloads through datacache; --offline reuses cached originals. The offline
pytest suite separately checks the VCF excerpts, RefSeq bases and haplotypes.
The UCSC chain is fetched for verification, not redistributed with this repo.
"""

import argparse
import gzip
import json
from pathlib import Path

from osteosarc import Cache

ROOT = Path(__file__).resolve().parents[1]


def chain_blocks(text):
    """Yield hg19 -> hg38 blocks in forward genomic coordinates (zero-based)."""
    header = None
    for line in text.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "chain":
            header = fields
            if header[4] != "+":
                raise ValueError("Expected forward hg19 chain coordinates")
            source, target = int(fields[5]), int(fields[10])
            continue
        size = int(fields[0])
        start, end = target, target + size
        if header[9] == "-":
            start, end = int(header[8]) - end, int(header[8]) - start
        yield dict(chain_id=header[12], strand=header[9],
                   source=[header[2], source, source + size],
                   target=[header[7], start, end])
        source += size
        target += size
        if len(fields) == 3:
            source += int(fields[1])
            target += int(fields[2])


def check_mapping(mapping, blocks):
    for wanted in mapping["blocks"]:
        matching = [b for b in blocks if b["chain_id"] == wanted["chain_id"]
                    and b["strand"] == wanted["strand"] == "+"
                    and b["source"][0] == wanted["source"][0]
                    and b["source"][1] <= wanted["source"][1]
                    and wanted["source"][2] <= b["source"][2]]
        actual, = matching
        offset = actual["target"][1] - actual["source"][1]
        assert wanted["target"] == [actual["target"][0], wanted["source"][1] + offset,
                                    wanted["source"][2] + offset]
    if "source_gap" in mapping:
        # The claimed gap must be exactly bounded by consecutive aligned blocks.
        left, right = mapping["blocks"]
        assert mapping["source_gap"] == [left["source"][0], left["source"][2], right["source"][1]]
        assert mapping["target_gap"] == [left["target"][0], left["target"][2], right["target"][1]]
        chrom, start, end = mapping["source_gap"]
        assert not any(b["source"][0] == chrom and b["source"][1] < end
                       and start < b["source"][2] for b in blocks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache", default=".cache/allele-resolution")
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    cache = Cache(args.cache, offline=args.offline)
    manifest = json.loads((ROOT / "tests/data/allele_resolution/provenance.json").read_text())
    for source in manifest["files"].values():
        cache.fetch(source["url"], sha256=source["sha256"])
        print("Verified", source["fixture"])
    chain = manifest["chain"]
    receipt = cache.fetch(chain["url"], sha256=chain["sha256"])
    with gzip.open(cache.path(receipt), "rt") as handle:
        blocks = list(chain_blocks(handle.read()))
    review = json.loads((ROOT / "osteosarc/data/allele_resolutions.json").read_text())
    for vid, entry in review.items():
        resolution = entry["resolution"]
        if "relationship" in resolution:
            resolution = resolution["relationship"]["evidence"]
        if "mapping" in resolution:
            check_mapping(resolution["mapping"], blocks)
            print("Verified mapping", vid)


if __name__ == "__main__":
    main()
