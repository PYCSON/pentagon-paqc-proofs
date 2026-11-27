"""Re-prove far-region leaves that tests/far_label_spot.py did not reproduce in one piece, by halving all six
coordinates (w in R^4, t, u) and applying the same exact test to every piece (a label that holds on every piece
holds on the union).  Only for labels T and TX of the E2/E3 charts (no table entries are used there).

usage: python3 -I label_subdivide.py TESTS_DIR OUT_DIR CHART BLOCK KEY(e.g. 16,5,8,0) J LABEL [max_depth=2]
"""
import sys, json, itertools
from fractions import Fraction as Fr
sys.path.insert(0, sys.argv[1])
import far_label_spot as F

out, chart, block = sys.argv[2], sys.argv[3], sys.argv[4]
key = [int(x) for x in sys.argv[5].split(',')]; J = int(sys.argv[6]); LABEL = sys.argv[7]
maxd = int(sys.argv[8]) if len(sys.argv) > 8 else 2
assert chart in ('E2', 'E3') and LABEL in ('T', 'TX')


def sub(lo, hi, d):
    """number of pieces used, or None if some piece fails at depth d"""
    if F.check(chart, key, J, lo, hi, LABEL, 0):
        return 1
    if d == 0:
        return None
    cnt = 0
    for bits in itertools.product((0, 1), repeat=6):
        l2, h2 = list(lo), list(hi)
        for i, b in enumerate(bits):
            m = (lo[i] + hi[i]) / 2
            if b:
                l2[i] = m
            else:
                h2[i] = m
        r = sub(l2, h2, d - 1)
        if r is None:
            return None
        cnt += r
    return cnt


n = bad = 0
for line in open('%s/%s/blocks/%s/leaves.jsonl' % (out, chart, block)):
    k, j, lo, hi, label = json.loads(line)
    if k != key or j != J or label != LABEL:
        continue
    lo = [Fr(x) for x in lo]; hi = [Fr(x) for x in hi]
    if F.check(chart, key, j, lo, hi, label, 0):
        continue
    n += 1
    print('leaf', [str(x) for x in lo], [str(x) for x in hi])
    r = sub(lo, hi, maxd)
    print('  ->', 'all pieces re-proved (%d pieces)' % r if r is not None else 'NOT re-proved up to depth %d' % maxd)
    bad += r is None
print('leaves needing subdivision: %d, still not re-proved: %d' % (n, bad))
