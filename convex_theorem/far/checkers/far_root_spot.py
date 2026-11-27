#!/usr/bin/env python3
"""far_root_spot.py RUN_OUT  (Claude, 2026-11-19)

Second implementation of the root ranges of the far cover (FAR_THEOREM.md section 4), from the note only (no production
module): enumerates every cell of side 1/16 of the three root ranges, decides in exact rational arithmetic whether it
can contain a strictly convex Q (all four invariants q2y, q3y, K, D positive somewhere: not "<= 0 at all 16 corners")
and, in the std chart, an admissible q (no edge 1, 2, 3 longer than 1 on the whole cell), and compares the result with
the root lists of the run (RUN_OUT/<chart>/roots.json; the lists of a complete run).  A cell that this program keeps
and the run does not is a gap of the cover; the converse is harmless (an unnecessary root).  Data only; python3 -I."""
import sys, os, json
from fractions import Fraction as Fr

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from far_label_spot import invariants, wcorners          # noqa: E402  (the chart maps of section 3, own implementation)

S = Fr(1, 16)


def box_min_norm2(lo, hi):
    """min of x^2 + y^2 over the rectangle [lo0, hi0] x [lo1, hi1]"""
    out = Fr(0)
    for a, b in zip(lo, hi):
        m = Fr(0) if a <= 0 <= b else min(abs(a), abs(b)); out += m * m
    return out


def keep(chart, lo, hi):
    vals = [invariants(chart, c) for c in wcorners(lo, hi)]
    if any(all(v[k] <= 0 for v in vals) for k in ('q2y', 'q3y', 'K', 'D')):
        return False
    if chart == 'std':                         # edges e1 = q2 - q1, e2 = q3 - q2, e3 = -q3 (each of length <= 1 if admissible)
        e1 = ((lo[0] - 1, lo[1]), (hi[0] - 1, hi[1]))
        e2 = ((lo[2] - hi[0], lo[3] - hi[1]), (hi[2] - lo[0], hi[3] - lo[1]))
        e3 = ((lo[2], lo[3]), (hi[2], hi[3]))
        if any(box_min_norm2(a, b) > 1 for a, b in (e1, e2, e3)):
            return False
    return True


def cells(chart):
    if chart == 'std':                         # [0,2] x [0,1] x [-1,1] x [0,1]
        rng = [range(0, 32), range(0, 16), range(-16, 16), range(0, 16)]
        org = [Fr(0)] * 4; side = [S] * 4
    elif chart == 'E2':                        # [-1,1] x [1/4,1] x [-4.5,4.5] x [0,1/8]  (cells of side 1/16, indices i/16)
        rng = [range(-16, 16), range(4, 16), range(-72, 72), range(0, 2)]
        org = [Fr(0)] * 4; side = [S] * 4
    else:                                      # E3: [1/4,2] x [0,1] x [-4,4] x [0,1/8]
        rng = [range(4, 32), range(0, 16), range(-64, 64), range(0, 2)]
        org = [Fr(0)] * 4; side = [S] * 4
    for a in rng[0]:
        for b in rng[1]:
            for c in rng[2]:
                for d in rng[3]:
                    i = (a, b, c, d)
                    yield i, [org[k] + i[k] * side[k] for k in range(4)], [org[k] + (i[k] + 1) * side[k] for k in range(4)]


def main():
    out = sys.argv[1]
    for chart in ('std', 'E2', 'E3'):
        f = os.path.join(out, chart, 'roots.json')
        if not os.path.exists(f):
            print(chart, ': no roots.json'); continue
        info = json.load(open(f)); run = {}
        for name, lst in info['blocks'].items():
            for r in lst:
                run[tuple(r[0][-4:])] = ([Fr(x) for x in r[1]], [Fr(x) for x in r[2]])
        mine = set(); n = 0; geom_bad = 0
        for i, lo, hi in cells(chart):
            n += 1
            if keep(chart, lo, hi):
                mine.add(i)
                if i in run and (run[i][0] != lo or run[i][1] != hi):
                    geom_bad += 1
        missing = sorted(mine - set(run)); extra = len(set(run) - mine)
        print('%s: cells %d, kept here %d, roots of the run %d (complete: %s); kept here but not a root: %d %s; '
              'roots not kept here: %d; roots with other bounds: %d' % (chart, n, len(mine), len(run), info.get('complete'),
                                                                          len(missing), missing[:10], extra, geom_bad))


if __name__ == '__main__':
    main()
