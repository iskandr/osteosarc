"""Run every example in README.md and docs/*.md against the live public sources.

    python scripts/check_docs.py                      # all pages, fresh temporary cache
    python scripts/check_docs.py --cache DIR docs/tour.md

Python blocks on one page share a namespace, in order, as a reader would run
them. Shell blocks run their `osteosarc` lines (interactive commands receive
"quit"). Install, clone and development commands are listed but not run; CI
covers the development commands. A block preceded by
`<!-- docs-check: skip (reason) -->` is reported as skipped.
"""

import argparse
import json
import os
import re
import shlex
import subprocess
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAGES = ["README.md", "docs/index.md", "docs/tour.md", "docs/explore.md", "docs/variants.md",
         "docs/timeline.md", "docs/reads.md", "docs/consumers.md", "docs/curation.md", "docs/api.md",
         "docs/migration.md", "docs/design.md", "docs/validation.md"]
FENCE = re.compile(r"^```(\w*)\s*$")
SKIP = re.compile(r"<!--\s*docs-check:\s*skip\b(.*?)-->")

DRIVER = r"""
import contextlib, io, json, sys, time, traceback
blocks = json.load(open(sys.argv[1]))
namespace = {"__name__": "__docs__"}
for block in blocks:
    started, out = time.time(), io.StringIO()
    try:
        with contextlib.redirect_stdout(out):
            exec(compile(block["code"], f"{block['page']}:{block['line']}", "exec"), namespace)
        status, error = "ok", ""
    except BaseException:
        status, error = "failed", traceback.format_exc(limit=-3)
    print(json.dumps(dict(line=block["line"], status=status, seconds=round(time.time() - started, 1),
                          error=error, output=out.getvalue()[-2000:])), flush=True)
"""


def blocks(page):
    lines, result, current = (ROOT / page).read_text().splitlines(), [], None
    for number, line in enumerate(lines, 1):
        fence = FENCE.match(line.strip())
        if current is None and fence:
            previous = "\n".join(lines[max(0, number - 3):number - 1])
            match = SKIP.search(previous)
            current = dict(page=page, line=number + 1, language=fence[1] or "text", code=[],
                           skip=(match[1].strip(" ()") or "marked") if match else None)
        elif current is not None and line.strip() == "```":
            current["code"] = "\n".join(current["code"])
            result.append(current)
            current = None
        elif current is not None:
            current["code"].append(line)
    return result


def run_python(page, items, env, cwd, timeout):
    runnable = [b for b in items if not b["skip"]]
    if not runnable:
        return {}
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as handle:
        json.dump(runnable, handle)
    try:
        process = subprocess.run([sys.executable, "-c", DRIVER, handle.name], env=env, cwd=cwd,
                                 capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return {b["line"]: dict(status="failed", seconds=timeout, error="page timed out") for b in runnable}
    results = {r["line"]: r for r in map(json.loads, filter(None, process.stdout.splitlines()))}
    for block in runnable:
        results.setdefault(block["line"], dict(status="failed", seconds=0, error=process.stderr[-2000:]))
    return results


def run_shell(block, env, cwd, timeout):
    outcomes = []
    for raw in block["code"].splitlines():
        command = raw.split(" #")[0].strip()
        if not command or command.startswith("#"):
            continue
        if not command.startswith("osteosarc"):
            outcomes.append(dict(command=command, status="not run"))
            continue
        started = time.time()
        try:
            process = subprocess.run(shlex.split(command), env=env, cwd=cwd, input="quit\n",
                                     capture_output=True, text=True, timeout=timeout)
            status = "ok" if process.returncode == 0 else "failed"
            error = "" if status == "ok" else (process.stderr or process.stdout)[-1500:]
        except subprocess.TimeoutExpired:
            status, error = "failed", "timed out"
        outcomes.append(dict(command=command, status=status, seconds=round(time.time() - started, 1), error=error))
    return outcomes


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("pages", nargs="*", default=PAGES)
    parser.add_argument("--cache", help="OSTEOSARC_CACHE to use (default: a fresh temporary directory)")
    parser.add_argument("--timeout", type=int, default=1800, help="Seconds per page or command")
    parser.add_argument("--show-output", action="store_true", help="Print what each Python block printed")
    args = parser.parse_args(argv)
    work = Path(tempfile.mkdtemp(prefix="osteosarc-docs-"))
    cache = Path(args.cache).resolve() if args.cache else work / "cache"
    bin_dir = str(Path(sys.executable).parent)
    env = dict(os.environ, OSTEOSARC_CACHE=str(cache), PATH=bin_dir + os.pathsep + os.environ.get("PATH", ""))
    print(f"cache {cache}\nwork  {work}\n")
    failures = 0
    for page in args.pages:
        items = [b for b in blocks(page) if b["language"] in ("python", "sh", "bash")]
        python = run_python(page, [b for b in items if b["language"] == "python"], env, work, args.timeout)
        for block in items:
            where = f"{page}:{block['line']}"
            if block["skip"]:
                print(f"  skip    {where}  ({block['skip']})")
            elif block["language"] == "python":
                result = python[block["line"]]
                failures += result["status"] != "ok"
                print(f"  {result['status']:<7} {where}  {result['seconds']}s")
                if args.show_output and result.get("output"):
                    print("          | " + result["output"].rstrip().replace("\n", "\n          | "))
                if result["status"] != "ok":
                    print("          " + result["error"].strip().replace("\n", "\n          "))
            else:
                for outcome in run_shell(block, env, work, args.timeout):
                    failures += outcome["status"] == "failed"
                    print(f"  {outcome['status']:<7} {where}  $ {outcome['command']}"
                          + (f"  {outcome['seconds']}s" if "seconds" in outcome else ""))
                    if outcome["status"] == "failed":
                        print("          " + outcome["error"].strip().replace("\n", "\n          "))
    print(f"\n{failures} failing example(s)")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
