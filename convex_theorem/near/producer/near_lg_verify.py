#!/usr/bin/env python3
"""Producer-side check of a near_lg.py run (Claude, 2026-11-16; prototype).  For every block file:
  1. the header is a block of the root grid that near_lg.blocks() lists (recomputed here);
  2. the leaves tile the block exactly (recursive check of the bisection structure, exact binary64 comparisons);
  3. skip leaves: the exact rational test of their label (mirror / far / inner) is repeated;
  4. LG leaves (all, or a random sample with --sample K per block): re-proved from scratch with lgpoly.check_box
     (a new trial family; if that fails once more with a higher degree, then reported).
And once: the root cube contains the ball |z - C| < 3/5, and the mirror map fixes C.
Usage: python3 near_lg_verify.py OUT [--sample K] [--seed S] [--blocks id,id,...]
Exit code 0 and a last line 'VERIFY OK' (all blocks present) or 'VERIFY OK (partial: n of m blocks)'."""
import sys, os, json, glob, random, argparse
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import near_lg as NL, lgbox as LB, lgpoly as LP
LP.MAX_SUB = 256; LP.SUB_RATIO = float('inf')      # the re-proof may bisect a leaf more generously than the producer did


def tiles(lo, hi, leaves):
    """True if the boxes `leaves` (pairs lo, hi) are exactly the leaves of a bisection tree of the box"""
    if len(leaves) == 1:
        return list(leaves[0][0]) == list(lo) and list(leaves[0][1]) == list(hi)
    if not leaves:
        return False
    for i in range(6):
        m = (lo[i] + hi[i]) / 2
        left = [l for l in leaves if l[1][i] <= m]; right = [l for l in leaves if l[0][i] >= m]
        if left and right and len(left) + len(right) == len(leaves):
            return tiles(lo, hi[:i] + [m] + hi[i + 1:], left) and tiles(lo[:i] + [m] + lo[i + 1:], hi, right)
    return False


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--sample', type=int, default=-1); ap.add_argument('--seed', type=int, default=1); ap.add_argument('--blocks', default='')
    a = ap.parse_args(); rng = random.Random(a.seed); bad = []
    # the root cube and the mirror map
    if not all(abs(Fraction(NL.ROOT_C[i]) - NL.C[i]) + Fraction(3, 5) <= Fraction(NL.ROOT_HALF) for i in range(6)):
        bad.append('root cube does not contain the ball of radius 3/5')
    C = NL.C
    if not (1 - C[4] == C[0] and C[5] == C[1] and 1 - C[2] == C[2]):
        bad.append('the mirror map does not fix C')
    bl = {b[0]: (b[1], b[2]) for b in NL.blocks()}; files = sorted(glob.glob(os.path.join(a.out, 'blocks', '*.jsonl')))
    if a.blocks:
        keep = set(a.blocks.split(',')); files = [f for f in files if os.path.basename(f)[:-6] in keep]
    nb = 0; cnt = {}; nre = 0
    for f in files:
        bid = os.path.basename(f)[:-6]; L = [json.loads(l) for l in open(f)]; hd = L[0]; leaves = L[1:]; nb += 1
        if bid not in bl or hd['block'] != bid or hd['lo'] != bl[bid][0] or hd['hi'] != bl[bid][1]:
            bad.append('%s: not a block of the root grid' % bid); continue
        if hd['stats'].get('open'):
            bad.append('%s: unfinished' % bid); continue
        if not tiles(list(hd['lo']), list(hd['hi']), [(l[0], l[1]) for l in leaves]):
            bad.append('%s: leaves do not tile the block' % bid); continue
        lg = []
        for l in leaves:
            cnt[l[2]] = cnt.get(l[2], 0) + 1
            if l[2] == 'LG':
                lg.append(l)
            elif l[2] in ('U', 'open'):
                bad.append('%s: unresolved leaf %s (%s)' % (bid, l[0], l[2]))
            elif NL.skip(l[0], l[1]) != l[2]:
                bad.append('%s: skip label %s not reproduced' % (bid, l[2]))
        for l in (lg if a.sample < 0 else rng.sample(lg, min(a.sample, len(lg)))):
            mid = [(x + y) / 2 for x, y in zip(l[0], l[1])]; rad = [(y - x) / 2 for x, y in zip(l[0], l[1])]; ok = False
            for N in (l[3] or 8, (l[3] or 8) + 2):
                fam = LB.Family(LB.near_verts, mid, N, order=1)
                if fam.ok and LP.check_box(fam, LP.near_aff, mid, rad, NL.FSTAR)[0]:
                    ok = True; break
            nre += 1
            if not ok:
                bad.append('%s: LG leaf %s not re-proved' % (bid, l[0]))
    for b in bad[:20]:
        print(' ', b)
    print('blocks checked %d of %d, leaves %s, LG leaves re-proved %d' % (nb, len(bl), json.dumps(cnt, sort_keys=True), nre))
    if bad:
        print('VERIFY FAILED: %d problems' % len(bad)); sys.exit(1)
    print('VERIFY OK' if nb == len(bl) and not a.blocks else 'VERIFY OK (partial: %d of %d blocks)' % (nb, len(bl)))


if __name__ == '__main__':
    main()
