"""Write MANIFEST-SHA256.txt: SHA-256 of every file of the release (sorted paths, `sha256sum -c` format).
Usage (from the release root): python3 tools/make_manifest.py"""
import os, hashlib
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SKIP_DIRS = {'__pycache__', '.git'}
SKIP_FILES = {'MANIFEST-SHA256.txt'}
SKIP_EXT = ('.aux', '.log', '.out')          # LaTeX by-products, skipped only under paper/ (evidence logs elsewhere are kept)
rows = []
for d, dirs, files in os.walk(ROOT):
    dirs[:] = sorted(x for x in dirs if x not in SKIP_DIRS)
    for f in sorted(files):
        rel = os.path.relpath(os.path.join(d, f), ROOT)
        if f in SKIP_FILES or f.endswith('.pyc') or (rel.startswith('paper' + os.sep) and f.endswith(SKIP_EXT)):
            continue
        rows.append((rel, hashlib.sha256(open(os.path.join(d, f), 'rb').read()).hexdigest()))
rows.sort()
with open(os.path.join(ROOT, 'MANIFEST-SHA256.txt'), 'w') as fh:
    for rel, h in rows:
        fh.write('%s  %s\n' % (h, rel))
print('MANIFEST-SHA256.txt:', len(rows), 'files')
