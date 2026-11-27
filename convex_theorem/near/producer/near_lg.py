#!/usr/bin/env python3
"""Near region by LG box certificates in exact polynomial form, PROTOTYPE driver (Claude, 2026-11-16; not reviewed).

Claim proved for a box B of the near chart z = (x0..x5) (pentagon (0,0), (1,0), (x0,x1), (x2,x3), (x4,x5)):
    every z in B with  1/20 <= |z - C| < 3/5  and  x2 >= 1/2  has  F(P(z)) > F*_up   (x2 = 1/2 included),
C = plan_cells.C (the binary64 chart point of the regular pentagon used by all programs).  B is cut by bisection into
leaves, each closed by one of
    mirror   sup x2 <= 1/2 on the leaf             (the mirror image has x2 >= 1/2 and the same distance to C)
    far      inf |z - C|^2 >= 9/25                 (outside the near region)
    inner    sup |z - C|^2 <  1/400                (not covered here; covered by Theorem L, far/LOCAL_LEMMA.md and lglocal.py, radius 52429/2^20 > 1/20 around the exact regular pentagon)
    LG       lgpoly.check_box: for all z of the leaf the fan (O; vertices) is positively oriented and the Lemma-LG
             conditions hold with a trial family computed at the leaf's centre (exact rational arithmetic).
The three skip tests are exact rational comparisons.  Heuristics (not part of the proof): when to attempt, the axis to
bisect, the degree.  Usage:
    python3 near_lg.py OUT --box lo0,..,lo5:hi0,..,hi5   [--N 8] [--alpha 0.2] [--max-sec S]      (a sub-box, for probes)
    python3 near_lg.py OUT --probe R,SEED,HALF                                               (cube of half-width HALF at distance R)
    python3 near_lg.py OUT run [--jobs J] [--max-sec S]      the whole near region, block by block (cubes of side 5/16), resumable
    python3 near_lg.py OUT status"""
import sys, os, json, math, time, argparse
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import numpy as np
import lgfan as LF, lgbox as LB, lgpoly as LP

# C: the binary64 chart point of the regular pentagon used by all programs (v2_cover/plan_cells.C; tests/lgpoly_check.py
# compares): (1 + cos 72, sin 72, 1/2, sin 72 + sin 36, -cos 72, sin 72) rounded, with C[4] = 1 - C[0] exactly.
CF = [float.fromhex(v) for v in ('0x1.4f1bbcdcbfa53p+0', '0x1.e6f0e134454ffp-1', '0x1.0000000000000p-1', '0x1.89f188bdcd7aep+0', '-0x1.3c6ef372fe94cp-2', '0x1.e6f0e134454ffp-1')]
C = [Fraction(v) for v in CF]
FSTAR = 18.919134195733914
NEAR2 = Fraction(9, 25); INNER2 = Fraction(1, 400)
ASPECT = 0.3                 # a failed box is bisected along the axis with the largest loss among the axes of half-width >= ASPECT * the largest
MIN_RAD = 2.0 ** -11         # a box with all half-widths below this that still fails is left unresolved (label U)


def d2_range(lo, hi):
    """exact [inf, sup] of |z - C|^2 over the box (Fractions)"""
    a = Fraction(0); b = Fraction(0)
    for i in range(6):
        l = Fraction(lo[i]) - C[i]; h = Fraction(hi[i]) - C[i]
        b += max(l * l, h * h)
        if l > 0:
            a += l * l
        elif h < 0:
            a += h * h
    return a, b


def skip(lo, hi):
    if Fraction(hi[2]) <= Fraction(1, 2):                       # closed half space x2 <= 1/2: mirror images of points with x2 >= 1/2
        return 'mirror'
    a, b = d2_range(lo, hi)
    if a >= NEAR2:
        return 'far'
    if b < INNER2:
        return 'inner'
    return None


def cover(lo, hi, N=8, alpha=0.7, max_sec=0, log=None, order=1):
    """bisection cover of the box; returns (stats, leaves) with leaves = [lo, hi, label, info]"""
    t0 = time.time(); stack = [(list(lo), list(hi), 0, None)]; leaves = []
    st = dict(mirror=0, far=0, inner=0, LG=0, attempts=0, fail=0, open=0, sec_family=0.0, sec_check=0.0, why={})
    while stack:
        if max_sec and time.time() - t0 > max_sec:                 # out of time: the remaining boxes become leaves with the label 'open' (no claim)
            st['open'] = len(stack); leaves += [[a_, b_, 'open', None] for a_, b_, h_, w_ in stack]; break
        a, b, hold, w = stack.pop()
        s = skip(a, b)
        if s is not None:
            st[s] += 1; leaves.append([a, b, s, None]); continue
        mid = [(x + y) / 2 for x, y in zip(a, b)]; rad = [(y - x) / 2 for x, y in zip(a, b)]
        r = math.sqrt(sum((m - c) ** 2 for m, c in zip(mid, CF)))
        if hold == 0 and max(rad) <= alpha * max(r, 0.04):
            st['attempts'] += 1; t = time.time(); fam = LB.Family(LB.near_verts, mid, N if r >= 0.09 else N + 2, order=order); t1 = time.time(); st['sec_family'] += t1 - t
            ok = False; info = None
            if fam.ok:
                ok, info = LP.check_box(fam, LP.near_aff, mid, rad, FSTAR)
            st['sec_check'] += time.time() - t1
            if ok:
                st['LG'] += 1; leaves.append([a, b, 'LG', dict(N=fam.N, c0=info['c0'], qp=info['qp'], tail=info['tail'])]); continue
            st['fail'] += 1; why = (info or {}).get('stage', 'family') + '-' + (info or {}).get('why', ''); st['why'][why] = st['why'].get(why, 0) + 1
            # heuristics for the bisection: the loss of every axis in the failed test, and how many levels to bisect
            # before the next attempt (a tail several times larger than the quadratic bound needs more than one halving)
            w = None; hold = 1
            if info is not None and 'lin' in info:
                w = [info['lin'][j] + info['hi'][j] for j in range(6)]
                if info.get('why') == 'tail' and info['qp'] > 0:
                    ratio = info['tail'] / info['qp']; hold = 1 if ratio <= 2.5 else (2 if ratio <= 8 else 3)
        mx = max(rad)
        if mx < MIN_RAD:                                            # give up on this leaf: label U (unresolved), never a claim
            st['U'] = st.get('U', 0) + 1; leaves.append([a, b, 'U', None]); continue
        if w is not None:                                           # after a failed test: the axis with the largest loss within the aspect limit
            cand = [i for i in range(6) if rad[i] >= mx * ASPECT]; i = max(cand, key=lambda j: w[j]); w2 = list(w); w2[i] /= 4
        else:                                                       # not attempted yet: the longest axis, farthest from C first
            cand = [i for i in range(6) if rad[i] >= mx * 0.99]; i = max(cand, key=lambda j: abs(mid[j] - CF[j])); w2 = None
        h2 = max(hold - 1, 0)
        stack.append((a[:i] + [mid[i]] + a[i + 1:], b, h2, w2 if h2 else None)); stack.append((a, b[:i] + [mid[i]] + b[i + 1:], h2, w2 if h2 else None))
    st['sec'] = round(time.time() - t0, 1)
    return st, leaves


BLOCK_SEC = float(os.environ.get('NEAR_LG_BLOCK_SEC', '7200'))      # time limit of one block; what is left gets the label 'open'
ROOT_HALF = 0.625; BLOCK_LEVEL = 2            # root cube: centre C rounded to 2^-10, half-width 5/8 (contains the ball |z - C| < 3/5); blocks: side 5/16
                                              # (side 5/32 with BLOCK_LEVEL = 3 was used first: 31,296 blocks, most of them a single leaf)
ROOT_C = [round(v * 1024) / 1024 for v in CF]


def blocks():
    """the blocks (cubes of side 2 ROOT_HALF / 2^BLOCK_LEVEL of the root grid) that are not closed by a skip test: list of (id, lo, hi)"""
    n = 2 ** BLOCK_LEVEL; side = 2 * ROOT_HALF / n; out = []

    def rec(level, idx):
        m = 2 ** level; s = 2 * ROOT_HALF / m
        lo = [ROOT_C[i] - ROOT_HALF + idx[i] * s for i in range(6)]; hi = [l + s for l in lo]
        if skip(lo, hi) is not None:
            return
        if level == BLOCK_LEVEL:
            out.append((''.join(str(t) for t in idx), lo, hi)); return
        for c in range(64):
            rec(level + 1, [2 * idx[i] + ((c >> i) & 1) for i in range(6)])
    rec(0, [0] * 6)
    return out


def block_job(args):
    bid, lo, hi, out, N, alpha = args
    f = os.path.join(out, 'blocks', bid + '.jsonl')
    if os.path.exists(f):
        return bid, None
    st, leaves = cover(lo, hi, N, alpha, BLOCK_SEC)
    tmp = f + '.tmp'
    with open(tmp, 'w') as g:
        g.write(json.dumps(dict(block=bid, lo=lo, hi=hi, stats=st)) + '\n')
        for l in leaves:
            g.write(json.dumps([l[0], l[1], l[2], (l[3] or {}).get('N')]) + '\n')
        g.flush(); os.fsync(g.fileno())
    os.replace(tmp, f)
    return bid, st


def run(out, jobs, N, alpha, max_sec, only=''):
    from multiprocessing import Pool
    os.makedirs(os.path.join(out, 'blocks'), exist_ok=True); t0 = time.time()
    bl = blocks()
    if only:                                                        # testing: only these block ids
        keep = set(only.split(',')); bl = [b for b in bl if b[0] in keep]
    todo = [b for b in bl if not os.path.exists(os.path.join(out, 'blocks', b[0] + '.jsonl'))]
    # nearest blocks first (they are the expensive ones and the most informative)
    todo.sort(key=lambda b: sum(((b[1][i] + b[2][i]) / 2 - CF[i]) ** 2 for i in range(6)))
    print('blocks %d, to do %d, jobs %d' % (len(bl), len(todo), jobs), flush=True); n = 0; tot = dict(LG=0, attempts=0, fail=0, sec=0.0)
    with Pool(jobs) as pool:
        for bid, st in pool.imap_unordered(block_job, [(b[0], b[1], b[2], out, N, alpha) for b in todo], chunksize=1):
            n += 1
            if st is not None:
                for k in tot:
                    tot[k] += st[k]
            if n % 10 == 0 or n == len(todo):
                print('%s  blocks %d/%d  LG leaves %d  attempts %d  failed %d  CPU %.0f s  wall %.0f s' % (time.strftime('%H:%M:%S'), n, len(todo), tot['LG'], tot['attempts'], tot['fail'], tot['sec'], time.time() - t0), flush=True)
            if max_sec and time.time() - t0 > max_sec:
                pool.terminate(); print('stopped by --max-sec', flush=True); break
    print('done' if n == len(todo) else 'incomplete', flush=True)


def status(out):
    import glob
    bl = blocks(); done = set(os.path.basename(f)[:-6] for f in glob.glob(os.path.join(out, 'blocks', '*.jsonl'))); tot = {}
    for f in glob.glob(os.path.join(out, 'blocks', '*.jsonl')):
        st = json.loads(open(f).readline())['stats']
        for k in ('LG', 'attempts', 'fail', 'mirror', 'far', 'inner', 'open', 'sec', 'U'):
            tot[k] = tot.get(k, 0) + st.get(k, 0)
    print('blocks done %d / %d' % (len(done & set(b[0] for b in bl)), len(bl)), json.dumps(tot))


def main():
    if len(sys.argv) > 2 and sys.argv[2] in ('run', 'status', 'count'):
        ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('cmd'); ap.add_argument('--jobs', type=int, default=2); ap.add_argument('--N', type=int, default=8)
        ap.add_argument('--alpha', type=float, default=0.7); ap.add_argument('--max-sec', type=float, default=0); ap.add_argument('--only', default=''); a = ap.parse_args()
        if a.cmd == 'run':
            run(a.out, a.jobs, a.N, a.alpha, a.max_sec, a.only)
        elif a.cmd == 'status':
            status(a.out)
        else:
            t = time.time(); bl = blocks(); print('active blocks', len(bl), '%.0f s' % (time.time() - t))
        return
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--box', default=''); ap.add_argument('--probe', default='')
    ap.add_argument('--N', type=int, default=8); ap.add_argument('--alpha', type=float, default=0.7); ap.add_argument('--max-sec', type=float, default=0); ap.add_argument('--order', type=int, default=1)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    if a.probe:
        R, seed, half = a.probe.split(','); R = float(R); half = float(half); rng = np.random.default_rng(int(seed))
        while True:
            d = rng.normal(size=6); d /= np.linalg.norm(d); zc = [round((CF[i] + R * d[i]) * 2 ** 12) / 2 ** 12 for i in range(6)]
            if zc[2] - half >= 0.5:
                break
        lo = [v - half for v in zc]; hi = [v + half for v in zc]
    else:
        l, h = a.box.split(':'); lo = [float(v) for v in l.split(',')]; hi = [float(v) for v in h.split(',')]
    st, leaves = cover(lo, hi, a.N, a.alpha, a.max_sec, order=a.order)
    vol = 1.0
    for x, y in zip(lo, hi):
        vol *= y - x
    lgvol = sum(np.prod([y - x for x, y in zip(l[0], l[1])]) for l in leaves if l[2] == 'LG')
    st.update(box=[lo, hi], volume=vol, lg_volume=float(lgvol))
    json.dump(st, open(os.path.join(a.out, 'stats.json'), 'w'), indent=1)
    with open(os.path.join(a.out, 'leaves.jsonl'), 'w') as f:
        for l in leaves:
            f.write(json.dumps(l) + '\n')
    print(json.dumps(st))


if __name__ == '__main__':
    main()
