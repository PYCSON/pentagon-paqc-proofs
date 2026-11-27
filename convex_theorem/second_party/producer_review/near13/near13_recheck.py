import sys, json, math, time
sys.path.insert(0, '<research tree>/probes/PENT/v2_run')
import lgbox as LB, lgpoly as LP, near_lg as NL
B = sys.argv[1]
ids = sys.argv[2:]
for tid in ids:
    blk, ln = tid.split('/')
    lines = open(B + '/' + blk + '.jsonl').read().splitlines()
    lo, hi, lab, N = json.loads(lines[int(ln)])
    mid = [(a + b) / 2 for a, b in zip(lo, hi)]; rad = [(b - a) / 2 for a, b in zip(lo, hi)]
    r = math.sqrt(sum((m - c) ** 2 for m, c in zip(mid, NL.CF)))
    a2, b2 = NL.d2_range(lo, hi)
    t = time.time(); fam = LB.Family(LB.near_verts, mid, N, order=1)
    ok, info = (LP.check_box(fam, LP.near_aff, mid, rad, NL.FSTAR) if fam.ok else (False, None))
    keep = {k: info[k] for k in ('c0', 'qp', 'tail', 'stage', 'why') if info and k in info}
    print(tid, lab, N, 'r=%.3f' % r, 'dmin=%.3f dmax=%.3f' % (math.sqrt(a2), math.sqrt(b2)), 'maxrad=%.4g' % max(rad), 'ok', ok, json.dumps(keep, default=float), '%.0fs' % (time.time() - t), flush=True)
