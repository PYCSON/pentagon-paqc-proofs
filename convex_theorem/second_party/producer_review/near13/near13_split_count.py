import sys, json, itertools, math
sys.path.insert(0, '<research tree>/probes/PENT/v2_run')
import near_lg as NL
from fractions import Fraction as Q
B = sys.argv[1]
IDS = '023112/23 023121/14 132012/151 132012/26 132012/47 132012/7 132022/157 132022/72 132102/10 133022/30 133022/40 232012/19 233021/2'.split()
def kids(lo, hi, k):
    # split every axis into k equal parts
    out = []
    for idx in itertools.product(range(k), repeat=6):
        a = [lo[i] + (hi[i] - lo[i]) * idx[i] / k for i in range(6)]; b = [lo[i] + (hi[i] - lo[i]) * (idx[i] + 1) / k for i in range(6)]
        out.append((a, b))
    return out
tot = {2: 0, 4: 0}
for tid in IDS:
    blk, ln = tid.split('/'); lo, hi, lab, N = json.loads(open(B + '/' + blk + '.jsonl').read().splitlines()[int(ln)])
    lo = [Q(x) for x in lo]; hi = [Q(x) for x in hi]
    a2, b2 = NL.d2_range(lo, hi)
    res = []
    for k in (2, 4):
        need = 0
        for a, b in kids(lo, hi, k):
            if Q(b[2]) <= Q(1, 2): continue
            if NL.d2_range(a, b)[0] >= NL.NEAR2: continue
            need += 1
        res.append('%d/%d' % (need, k ** 6)); tot[k] += need
    print(tid, 'dmin %.4f' % math.sqrt(a2), 'maxrad %.4f' % float(max((h - l) / 2 for l, h in zip(lo, hi))), 'need-LG after 2^6 split:', res[0], ' after 4^6 split:', res[1])
print('total', tot)
