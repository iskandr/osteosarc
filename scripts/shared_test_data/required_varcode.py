"""Varcode's required records for the shared test data: the ONT reads behind its observed junctions.

python scripts/shared_test_data/required_varcode.py ~/code/varcode required-varcode.json.gz

tests/data/osteosarc_observed_junctions.json keeps each junction's read name
and sequence, not the read. This names each read, with the windows around both
of its junction ends, so the shared test data keeps the whole read from its ONT BAM.
"""

import argparse
import gzip
import json
import subprocess
from pathlib import Path

BUCKET = "https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/ONT/IPISRC044_ONT_upload/IPISRC044_ONT/processed/"
#: The junction file names its sources by Isovar's source IDs; its source summary
#: (tests/data/osteosarc/expansion/audit/SOURCE_SUMMARY.md) gives these labels.
SOURCES = {
    "7010dd389c50467d": BUCKET + "IPISRC044_T1_sclrs_ONT/IPISRC044_T1_sclrs_ONT/IPISRC044_T1_sclrs_ONT_dedup/"
                                 "IPISRC044_T1_sclrs_ONT_dedup.bam",  # T1 UCSF Tumor ONT 2024-06
    "2bc1fc291308debb": BUCKET + "IPISRC044_T2_sclrs_ONT/IPISRC044_T2_sclrs_ONT/IPISRC044_T2_sclrs_ONT_dedup/"
                                 "IPISRC044_T2_sclrs_ONT_dedup.bam",  # T2 UCSF Tumor ONT 2025-01
}
WINDOW = 1000


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("varcode", type=Path, help="A Varcode checkout")
    parser.add_argument("output", type=Path, help="Where to write the required records (.json.gz)")
    parser.add_argument("revision", nargs="?", default="origin/main")
    args = parser.parse_args()

    def git(*command):
        return subprocess.run(["git", "-C", str(args.varcode), *command], capture_output=True, text=True,
                              check=True).stdout
    commit = git("rev-parse", args.revision).strip()
    junctions = json.loads(git("show", f"{commit}:tests/data/osteosarc_observed_junctions.json"))
    subsets = {}
    for record in junctions["records"]:
        ends = [("chr" + record["contig"], int(record["start"])), ("chr" + record["mate_contig"], int(record["mate_start"]))]
        subsets[f"varcode/osteosarc_observed_junctions.json#{record['label']}"] = dict(
            source=SOURCES[record["source_id"]], names=[record["read_name"]],
            regions=[[contig, max(0, position - WINDOW), position + WINDOW] for contig, position in ends],
            description=f"The ONT read behind Varcode's {record['label']} junction fragment")
    text = json.dumps(dict(consumer="varcode", commit=commit, subsets=subsets))
    args.output.write_bytes(gzip.compress(text.encode(), mtime=0))
    print(f"{len(subsets)} subsets: " + ", ".join(subsets))


if __name__ == "__main__":
    main()
