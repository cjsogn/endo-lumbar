"""Verify the analysis package's file allowlist and source portability.

Inputs: RELEASE_FILES.txt and either working files or the Git index.
Outputs: a pass message or an assertion identifying the failed check.
"""
from pathlib import Path
import argparse
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--staged",
    action="store_true",
    help="Inspect the Git index instead of working files",
)
args = parser.parse_args()


def read_working(name):
    return (ROOT / name).read_bytes()


def read_staged(name):
    return subprocess.check_output(["git", "show", ":" + name], cwd=ROOT)


read = read_staged if args.staged else read_working
allowed = set(read("RELEASE_FILES.txt").decode("utf-8").splitlines())
assert all(
    name and not name.startswith("/") and ".." not in Path(name).parts
    for name in allowed
)

if args.staged:
    tracked = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT)
    names = {name for name in tracked.decode().split("\0") if name}
else:
    names = {
        str(path.relative_to(ROOT))
        for path in ROOT.rglob("*")
        if path.is_file() and ".git" not in path.parts
    }

assert names == allowed, (
    f"Files differ from allowlist. Unexpected: {sorted(names - allowed)}. "
    f"Missing: {sorted(allowed - names)}"
)

for name in sorted(names):
    path = Path(name)
    assert path.suffix in (".R", ".py", ".md", ".txt") or name in (
        "LICENSE", ".gitignore"
    ), name
    assert not (ROOT / name).is_symlink(), name

    data = read(name)
    assert b"\x00" not in data, f"Binary content: {name}"
    text = data.decode("utf-8")

    if name.startswith("scripts/analysis/"):
        prohibited = r"/Users/|/home/|~/\.[^/]+/|TODO|FIXME|PLACEHOLDER"
        assert not re.search(prohibited, text, re.IGNORECASE), (
            f"Nonportable path or unfinished source: {name}"
        )

print(f"PASS: {len(names)} allowed source/documentation files.")
