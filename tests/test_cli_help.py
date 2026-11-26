"""Tests for shared datacache.cli_help integration, argcomplete, and lazy imports."""

from __future__ import annotations

import io
import os
import subprocess
import sys
from pathlib import Path

import pytest
from datacache.cli_help import (
    ANSI_BOLD,
    ANSI_BOLD_CYAN,
    ANSI_CYAN,
    ANSI_RESET,
    HelpArgumentParser,
    HelpFormatter,
)

import osteosarc
import osteosarc.cli as cli_mod
from osteosarc.cli import _DEFAULT_DERIVED_UNUSED_DAYS, main, parser


class _Stream(io.StringIO):
    def __init__(self, tty: bool):
        super().__init__()
        self._tty = tty

    def isatty(self) -> bool:
        return self._tty


def test_parser_uses_shared_help_argument_parser():
    root = parser()
    assert isinstance(root, HelpArgumentParser)
    assert issubclass(root.formatter_class, HelpFormatter)

    subparsers_action = next(
        a for a in root._actions if hasattr(a, "choices") and isinstance(a.choices, dict)
    )
    for sub in subparsers_action.choices.values():
        assert isinstance(sub, HelpArgumentParser)
        assert issubclass(sub.formatter_class, HelpFormatter)


def test_help_plain_vs_tty_and_no_color(monkeypatch):
    monkeypatch.delenv("NO_COLOR", raising=False)
    monkeypatch.delenv("TERM", raising=False)

    root = parser()
    plain = root.format_help()
    assert "\033[" not in plain
    assert "Browse:" in plain

    tty_stream = _Stream(True)
    root.print_help(file=tty_stream)
    styled = tty_stream.getvalue()
    assert ANSI_BOLD + "usage:" + ANSI_RESET in styled
    assert ANSI_CYAN + "--version" + ANSI_RESET in styled
    assert ANSI_BOLD_CYAN + "Browse:" + ANSI_RESET in styled

    non_tty_stream = _Stream(False)
    root.print_help(file=non_tty_stream)
    assert non_tty_stream.getvalue() == plain

    monkeypatch.setenv("NO_COLOR", "1")
    no_color_stream = _Stream(True)
    root.print_help(file=no_color_stream)
    assert no_color_stream.getvalue() == plain


def test_argcomplete_marker_and_registration(monkeypatch):
    pkg_dir = Path(osteosarc.__file__).parent
    for fname in ("cli.py", "__main__.py"):
        head = (pkg_dir / fname).read_text(encoding="utf-8")[:1024]
        assert "PYTHON_ARGCOMPLETE_OK" in head, f"Missing PYTHON_ARGCOMPLETE_OK in {fname}"

    called: list[object] = []
    monkeypatch.setattr(cli_mod, "enable_argcomplete", lambda p, **kw: called.append(p) or True)
    with pytest.raises(SystemExit) as exc:
        main(["--help"])
    assert exc.value.code == 0
    assert len(called) == 1
    assert isinstance(called[0], HelpArgumentParser)


def test_package_and_cli_help_stay_import_light():
    code = (
        "import sys; "
        "import osteosarc; "
        "import osteosarc.cli; "
        "assert osteosarc.__version__; "
        "p = osteosarc.cli.parser(); "
        "p.format_help(); "
        "heavy = [m for m in ('pysam', 'requests', 'bs4', 'beautifulsoup4', 'urllib3') "
        "if m in sys.modules]; "
        "assert not heavy, f'Unexpected heavy imports on CLI help: {heavy}'"
    )
    proc = subprocess.run(
        [sys.executable, "-c", code],
        env=os.environ.copy(),
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr


def test_top_level_lazy_exports_resolve():
    helper_exports = (
        "__version__",
        "verify_digest",
        "verify_gzip_digests",
        "verify_manifest_files",
        "digest",
        "parse_variants",
        "read_receipt_files",
        "ReadCount",
        "assembly_from_header",
        "count_reads",
        "inspect_alignment",
    )
    for name in (*osteosarc.__all__, *helper_exports):
        assert getattr(osteosarc, name) is not None
    with pytest.raises(AttributeError):
        getattr(osteosarc, "_definitely_not_an_export")


def test_cli_derived_unused_days_matches_cache():
    from osteosarc.cache import DERIVED_UNUSED_DAYS

    assert _DEFAULT_DERIVED_UNUSED_DAYS == DERIVED_UNUSED_DAYS
