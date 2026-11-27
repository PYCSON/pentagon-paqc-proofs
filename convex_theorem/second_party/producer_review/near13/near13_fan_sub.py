import sys, json, itertools, math
sys.path.insert(0, '<research tree>/probes/PENT/v2_run')
import near_lg as NL, lgbox as LB
B = sys.argv[1]
IDS = '023112/23 023121/14 132012/151 132012/26 132012/47 132012/7 132022/157 132022/72 132102/10 133022/30 133022/40 232012/19 233021/2'.split()
def fanmin(lo, hi):
    mid = [(a + b) / 2 for a, b in zip(lo, hi)]; V = LB.near_verts(mid); n = len(V)
    O = (sum(a for a, b in V) / n, sum(b for a, b in V) / n); m = 9
    cvx = 9
    for c in itertools.product((0, 1), repeat=6):
        z = [hi[i] if c[i] else lo[i] for i in range(6)]; V = LB.near_verts(z)
        for k in range(n):
            (ax, ay), (bx, by) = V[k], V[(k + 1) % n]; m = min(m, (ax - O[0]) * (by - O[1]) - (ay - O[1]) * (bx - O[0]))
            (px, py) = V[k - 1]; cvx = min(cvx, (ax - px) * (by - ay) - (ay - py) * (bx - ax))
    return m, cvx
for tid in IDS:
    blk, ln = tid.split('/'); lo, hi, lab, N = json.loads(open(B + '/' + blk + '.jsonl').read().splitlines()[int(ln)])
    m0, c0 = fanmin(lo, hi)
    ms = []; cs = []
    for k in (2, 4):
        worst = 9; worstc = 9; cnt = 0
        for idx in itertools.product(range(k), repeat=6):
            a = [lo[i] + (hi[i] - lo[i]) * idx[i] / k for i in range(6)]; b = [lo[i] + (hi[i] - lo[i]) * (idx[i] + 1) / k for i in range(6)]
            if NL.skip(a, b) is not None: continue
            cnt += 1; m, c = fanmin(a, b); worst = min(worst, m); worstc = min(worstc, c)
        ms.append(worst); cs.append(worstc)
    print('%-11s whole box: min fan det %.4f, min convexity %.3f | kept sub-boxes 2^6: fan %.4f cvx %.3f | 4^6: fan %.4f cvx %.3f' % (tid, m0, c0, ms[0], cs[0], ms[1], cs[1]), flush=True)
