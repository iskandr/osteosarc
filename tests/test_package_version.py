"""The published wheel and import must use the same single version source."""

import os
import subprocess
import sys
from pathlib import Path

from osteosarc._version import __version__


def test_installed_wheel_reports_its_distribution_version(tmp_path):
    root = Path(__file__).resolve().parents[1]
    artifacts, target = tmp_path / "dist", tmp_path / "installed"
    subprocess.run([sys.executable, "-m", "build", "--wheel", "--no-isolation",
                    "--outdir", str(artifacts), str(root)], check=True, capture_output=True)
    wheel, = artifacts.glob("*.whl")
    subprocess.run([sys.executable, "-m", "pip", "install", "--no-deps", "--no-compile",
                    "--target", str(target), str(wheel)], check=True, capture_output=True)
    env = dict(os.environ, PYTHONPATH=str(target))
    check = (
        "import importlib.metadata, pathlib, sys, osteosarc;"
        "assert pathlib.Path(osteosarc.__file__).is_relative_to(pathlib.Path(sys.argv[1]));"
        "assert osteosarc.__version__ == importlib.metadata.version('osteosarc') == sys.argv[2]"
    )
    subprocess.run([sys.executable, "-c", check, str(target), __version__],
                   cwd=tmp_path, env=env, check=True, capture_output=True)

