import sys, json, glob, os, math
sys.path.insert(0, '<research tree>/probes/PENT/v2_run')
import near_lg as NL
B = sys.argv[1]; FAIL = set(sys.argv[2].split(','))
rows = []
for p in sorted(glob.glob(B + '/*.jsonl')):
    with open(p) as f:
        next(f)
        for i, line in enumerate(f, 1):
            lo, hi, lab, N = json.loads(line)
            if lab != 'LG': continue
            mid = [(a + b) / 2 for a, b in zip(lo, hi)]; rad = max((b - a) / 2 for a, b in zip(lo, hi))
            r = math.sqrt(sum((m - c) ** 2 for m, c in zip(mid, NL.CF)))
            rows.append((os.path.basename(p)[:-6] + '/' + str(i), r, rad, N))
print('LG near boxes', len(rows))
import collections
for rmin in (0.6, 0.65, 0.68, 0.7):
    sel = [x for x in rows if x[1] >= rmin]
    by = collections.Counter(round(x[2], 4) for x in sel)
    nf = sum(1 for x in sel if x[0] in FAIL)
    print('r >= %.2f: %d boxes, half-width counts %s, of which in the 13: %d' % (rmin, len(sel), dict(sorted(by.items())), nf))
sel = [x for x in rows if x[1] >= 0.68 and x[2] > 0.07]
sel.sort(key=lambda x: -x[1])
print('largest r among r>=0.68, halfwidth>0.07:', [(x[0], round(x[1], 3), x[2], x[0] in FAIL) for x in sel[:25]])
