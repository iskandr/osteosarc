"""Drive one pinned recipe through the libraries' own bundle builders, offline.

python -m scripts.check_fixture_consumers --topiary /checkout --vaxrank /checkout [--isovar /checkout]
An explicit adoption check, not a dependency of ordinary Osteosarc tests, for the
libraries that still build their own bundles (give at least two). A library that
has moved its test reads onto openvax-v1, as Isovar does from isovar#398, has no
builder to compare.
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import pysam

from osteosarc import verify_bundle
from osteosarc.records import record_multiset


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("isovar", "topiary", "vaxrank"):
        parser.add_argument("--" + name, type=Path)
    args = parser.parse_args()
    given = [name for name in ("isovar", "topiary", "vaxrank") if getattr(args, name)]
    if len(given) < 2:
        parser.error("give at least two library checkouts to compare")
    with tempfile.TemporaryDirectory() as temporary:
        work = Path(temporary)
        path = work / "input.bam"
        header = dict(HD=dict(SO="coordinate"), SQ=[dict(SN="chr1", LN=248956422), dict(SN="chr2", LN=242193529)])
        with pysam.AlignmentFile(path, "wb", header=header) as bam:
            r = pysam.AlignedSegment.fromstring("witness\t0\tchr1\t101\t60\t4M\t*\t0\t0\tACGT\t*", bam.header)
            bam.write(r)
            bam.write(r)
        pysam.index(str(path))
        recipe = dict(schema_version=1, id="cross-consumer-v1", targets={"v": dict(kind="small_variant", assembly="GRCh38",
                      reference=dict(id="GRCh38"), coordinates="one-based", contig="chr1", position=101, ref="A", alt="C")},
                      sources={"rna": dict(identity=dict(id="source"), assembly="GRCh38", sample="sample", library="library", product="original")},
                      members={"case": dict(target="v", source="rna", policy=dict(kind="exact", version=1, records=dict(record_multiset(path))))})
        (work / "recipe.json").write_text(json.dumps(recipe))
        # All subprocess Python network calls fail, including accidental imports.
        (work / "sitecustomize.py").write_text("import socket\ndef blocked(*a,**k): raise RuntimeError('network disabled')\nsocket.socket.connect=blocked\nsocket.create_connection=blocked\n")
        env = dict(os.environ, PYTHONPATH=os.pathsep.join([str(work), str(Path(__file__).resolve().parents[1]), os.environ.get("PYTHONPATH", "")]))
        commands = {
            "isovar": ["-m", "isovar.sid_data", "generate"],
            "topiary": ["-m", "scripts.generate_sid_fixtures"],
            "vaxrank": ["examples/osteosarc_test_data/build.py", "--cache", str(work / "cache")],
        }
        manifests = []
        for consumer, command in ((name, commands[name]) for name in given):
            destination = work / consumer
            subprocess.run([sys.executable, *command, "--panel-recipe", str(work / "recipe.json"),
                            "--panel-source", f"rna={path}", "--output", str(destination), "--offline"],
                           cwd=getattr(args, consumer), env=env, check=True)
            manifest = verify_bundle(destination)
            manifests.append((manifest["recipe_sha256"], manifest["members"], manifest["sources"]))
        assert all(manifest == manifests[0] for manifest in manifests)
        print(f"{', '.join(given)} builders agree on recipe, complete record multiset, reasons and "
              "source/header provenance")


if __name__ == "__main__":
    main()
