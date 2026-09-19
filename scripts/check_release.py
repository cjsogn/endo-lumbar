"""Verify the explicit public-release allowlist before committing or archiving."""
from pathlib import Path
import argparse
import re
import subprocess

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--staged', action='store_true', help='Inspect the Git index instead of working files')
args = parser.parse_args()
allowed = set((ROOT / 'RELEASE_FILES.txt').read_text().splitlines())
assert all(x and not x.startswith('/') and '..' not in Path(x).parts for x in allowed)
if args.staged:
    names = subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).decode().split('\0')
    names = {x for x in names if x}
    read = lambda name: subprocess.check_output(['git', 'show', ':' + name], cwd=ROOT)
else:
    names = {str(p.relative_to(ROOT)) for p in ROOT.rglob('*') if p.is_file() and '.git' not in p.parts}
    read = lambda name: (ROOT / name).read_bytes()
assert names == allowed, f'Files differ from allowlist. Unexpected: {sorted(names-allowed)}. Missing: {sorted(allowed-names)}'
for name in sorted(names):
    p = Path(name)
    assert p.suffix in ('.R', '.py', '.md', '.txt') or name in ('LICENSE', '.gitignore'), name
    assert not (ROOT / name).is_symlink(), name
    data = read(name)
    assert b'\x00' not in data, f'Binary content: {name}'
    text = data.decode('utf-8')
    if name.startswith('scripts/analysis/'):
        prohibited = r'ggplot|ggsave|geom_|matplotlib|python-docx|library\(officer\)|\b(?:pdf|png|svg|plot)\s*\(|/Users/|/home/|\.claude|\.codex|TODO|FIXME|PLACEHOLDER'
        assert not re.search(prohibited, text, re.IGNORECASE), f'Excluded content: {name}'
print(f'PASS: {len(names)} allowed text/source files. No data, models, results, plotting or document generators.')
