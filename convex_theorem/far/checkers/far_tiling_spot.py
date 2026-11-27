#!/usr/bin/env python3
"""far_tiling_spot.py RUN_OUT  (Claude, 2026-11-19)

Second implementation of the tiling check of the far cover (FAR_THEOREM.md section 7, appendix A), no production
module: for every block of the run, every root of the block and every j in {3, 0}, the leaves (leaves.jsonl) and the
certified boxes (cert.jsonl) lie in root x [0,1]^2, have total volume equal to its volume, and have pairwise disjoint
interiors (sweep along the first coordinate).  Exact rationals.  Data only; python3 -I."""
import sys, os, json, glob, collections
from fractions import Fraction as Fr
from math import prod


def main():
    out = sys.argv[1]; tot = collections.Counter(); bad = []
    for chart in ('std', 'E2', 'E3'):
        rf = os.path.join(out, chart, 'roots.json')
        if not os.path.exists(rf):
            continue
        blocks = json.load(open(rf))['blocks']
        for name in sorted(blocks):
            d = os.path.join(out, chart, 'blocks', name)
            if not os.path.exists(os.path.join(d, 'done.json')):
                tot['unfinished blocks'] += 1; continue
            roots = {tuple(r[0]): ([Fr(x) for x in r[1]] + [Fr(0), Fr(0)], [Fr(x) for x in r[2]] + [Fr(1), Fr(1)]) for r in blocks[name]}
            boxes = collections.defaultdict(list)
            for line in open(os.path.join(d, 'leaves.jsonl')):
                key, j, lo, hi, lab = json.loads(line); boxes[(tuple(key), j)].append(([Fr(x) for x in lo], [Fr(x) for x in hi]))
            cf = os.path.join(d, 'cert.jsonl')
            if os.path.exists(cf):
                for line in open(cf):
                    c = json.loads(line); boxes[(tuple(c['key']), c['j'])].append(([Fr(x) for x in c['lo']], [Fr(x) for x in c['hi']]))
            for key in roots:
                for j in (3, 0):
                    L = boxes.get((key, j), []); rlo, rhi = roots[key]; tot['columns'] += 1; tot['boxes'] += len(L)
                    ok = all(all(rlo[i] <= lo[i] < hi[i] <= rhi[i] for i in range(6)) for lo, hi in L)
                    ok = ok and sum(prod(hi[i] - lo[i] for i in range(6)) for lo, hi in L) == prod(rhi[i] - rlo[i] for i in range(6))
                    if ok:
                        act = []
                        for lo, hi in sorted(L, key=lambda b: b[0][0]):
                            act = [a for a in act if a[1][0] > lo[0]]
                            if any(all(max(a[0][i], lo[i]) < min(a[1][i], hi[i]) for i in range(6)) for a in act):
                                ok = False; break
                            act.append((lo, hi))
                    if not ok:
                        bad.append((chart, name, key, j))
            if any(k not in roots for k, _ in boxes):
                bad.append((chart, name, 'box of an unknown root'))
            tot['blocks'] += 1
    print(dict(tot), 'failures:', len(bad), bad[:10])


if __name__ == '__main__':
    main()
