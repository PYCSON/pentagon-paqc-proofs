"""Quadrilateral table in a chart (Claude, 2026-11-15; FAR_SPEC 9.9).  DRAFT for review.
For every root box c of the chart (std: the cells of side 1/16 of t4table; E2 / E3: qchart.edge_roots) a rigorous
      LB_c <= F(Q)  for every convex quadrilateral Q of c
is the maximum of
   closed forms (qchart.quad_lb: floor 2 pi^2, triangle bound u3, strip bounds), and
   tv:  the tangent / Bernstein bound of tvn.minimize on the cell (n = 4, d = 4), from 16 corner certificates
        (t4table.corner_job: PolyMPS + FHM / HKS + Taylor sup), when all 16 corner quadrilaterals are strictly convex and
        certified.  The corner vertices must be exact in binary64 (Fraction check).
Corner store: key = t4table.corner_key(q) (the standard coordinates of the corner), shared by all charts.
Output OUT/cells.jsonl (key, LB, how, closed, tv), OUT/corners.jsonl (key, F, info), OUT/table_info.json.  Resumable.
Usage: python3 qtable.py OUT --chart std|E2|E3 [--std-keys T4_DIR] [--sample N] [--seed S] [--jobs 2] [--corners FILE ...]"""
import os, sys, json, time, argparse, random, itertools
from fractions import Fraction
from multiprocessing import Pool
HERE = os.path.dirname(os.path.abspath(__file__))
for p in ('.', '../far', '../scripts', '../v2_cover'):
    sys.path.insert(0, os.path.join(HERE, p))
import numpy as np
from flint import arb
import rig_paqc as R
import cutchart as CC
import qchart as QC
import t4table as T4
import jsonl_safe as JS

CORN4 = list(itertools.product([-1, 1], repeat=4))


def corner_ws(lo4, hi4):
    return [[hi4[i] if s[i] > 0 else lo4[i] for i in range(4)] for s in CORN4]


def std_q(w, chart):
    """standard coordinates of the corner w as binary64 numbers, exact (Fraction check)"""
    q = [float(v) for v in CC.to_std([float(v) for v in w] + [0.0, 0.0], chart)[:4]]
    qf = CC.to_std([Fraction(float(v)) for v in w] + [Fraction(0), Fraction(0)], chart)[:4]
    assert all(Fraction(a) == b for a, b in zip(q, qf)), 'corner not exact in binary64'
    return q


QLG = os.environ.get('PENT_QLG', '0') == '1'          # 2026-11-17: cell bounds by LG box certificates (lgpoly.box_lower) instead of the tangent bound
QLG_N = int(os.environ.get('PENT_QLG_N', '6')); QLG_DEPTH = int(os.environ.get('PENT_QLG_DEPTH', '1'))
QLG_TOP = 36.0                                         # the certificate cannot prove more than KS = 36.33
QLG_SAFE = 1 - 2.0 ** -10                              # the recorded bound is this fraction of the bound found (the verifier re-proves it with its own family)


def quad_aff(chart):
    """vertex map of the quadrilateral chart for polynomial arguments (lgpoly): (0,0), (1,0), q2, q3"""
    def f(x, const):
        q = CC.to_std(list(x) + [const(0), const(0)], chart)
        return [(const(0), const(0)), (const(1), const(0)), (q[0], q[1]), (q[2], q[3])]
    return f


def quad_verts(chart):
    def f(z):
        q = CC.to_std(list(z) + [z[0] * 0, z[0] * 0], chart)
        return [(z[0] * 0, z[0] * 0), (z[0] * 0 + 1, z[0] * 0), (q[0], q[1]), (q[2], q[3])]
    return f


def lg_lower(lo4, hi4, chart, L0, N=None):
    """a number L > L0 with F(Q) > L PROVED for every quadrilateral of the box with positive fan determinants (in
    particular every convex one), or None.  One family at the centre, exact polynomial test, no bisection."""
    import lgbox as LB, lgpoly as LP
    mid = [(a + b) / 2 for a, b in zip(lo4, hi4)]; rad = [(b - a) / 2 for a, b in zip(lo4, hi4)]
    try:
        fam = LB.Family(lambda z: [(float(a), float(b)) for a, b in quad_verts(chart)(list(z))], mid, N or QLG_N, order=1)
        if not fam.ok:
            return None
        L = LP.box_lower(fam, quad_aff(chart), mid, rad, float(L0), QLG_TOP)
    except Exception:
        return None
    return None if L is None else L * QLG_SAFE


def lg_check(lo4, hi4, chart, L, N=None):
    """True if F(Q) > L is PROVED on the box (verifier side of lg_lower: own family, degree N, then N + 2).
    The exact positivity test of P1 - L' P2 is a sufficient test and is NOT monotone in L' (2026-11-19: 10 of 8,288 table
    parts of the user's far5 run fail at the recorded L = QLG_SAFE * L_found but pass above it).  So the test is also tried
    at L' = L (1 + k 2^-14), k = 1..63 (the producer proved L / QLG_SAFE = L (1 + 2^-10 + ...)); F > L' >= L proves F > L."""
    import lgbox as LB, lgpoly as LP
    mid = [(a + b) / 2 for a, b in zip(lo4, hi4)]; rad = [float((b - a) / 2) for a, b in zip(lo4, hi4)]
    for n in ((N or QLG_N), (N or QLG_N) + 2):
        try:
            fam = LB.Family(lambda z: [(float(a), float(b)) for a, b in quad_verts(chart)(list(z))], mid, n, order=1)
            if not fam.ok:
                continue
            (P1, P2), PD, dets = LP.build(fam, quad_aff(chart), mid, None)
            if not all(LP.poly_positive(P, rad)[0] for P in list(dets) + [PD]):
                continue
            for k in range(64):
                if LP.poly_positive(P1 - LP.fq(L * (1 + k * 2.0 ** -14)) * P2, rad)[0]:
                    return True
        except Exception:
            pass
    return False


def closed_lb(lo4, hi4, chart):
    """(bound, how, no convex member): the closed-form bounds of the box"""
    inv = QC.invariants(lo4, hi4, chart)
    if QC.noncvx(inv):
        return float('inf'), 'noncvx', True
    closed, how = QC.quad_lb(inv); fb = QC.quad_fibre(lo4, hi4, chart, inv)
    return (fb, 'fibre', False) if fb > closed else (closed, how, False)


def lg_cell(lo4, hi4, chart, closed, depth=None):
    """(L, parts): the closed-form bound, or the LG bound of the cell, or (if the cell has none) of its 16 sub-cells.
    parts: list of [lo4, hi4, bound, 'lg' or the closed-form name]; the bound of the cell is the minimum over the parts."""
    depth = QLG_DEPTH if depth is None else depth
    L = lg_lower(lo4, hi4, chart, closed)
    best = (L, [[list(lo4), list(hi4), L, 'lg']]) if L is not None and L > closed else (closed, None)
    if depth > 0 and best[1] is None:                               # the 16 sub-cells only if the cell itself has no certificate
        import itertools
        mid = [(a + b) / 2 for a, b in zip(lo4, hi4)]; parts = []; mn = float('inf')
        for sgn in itertools.product((0, 1), repeat=4):
            l2 = [mid[i] if sgn[i] else lo4[i] for i in range(4)]; h2 = [hi4[i] if sgn[i] else mid[i] for i in range(4)]
            c2, how2, nc = closed_lb(l2, h2, chart)
            if nc:
                continue                                           # no convex member: nothing to bound
            L2 = lg_lower(l2, h2, chart, c2)
            if L2 is not None and L2 > c2:
                parts.append([l2, h2, L2, 'lg']); mn = min(mn, L2)
            else:
                parts.append([l2, h2, c2, how2]); mn = min(mn, c2)
            if mn <= best[0]:
                break                                               # cannot improve any more
        else:
            if parts and mn > best[0]:
                best = (mn, parts)
    return best


def cell_job(args):
    key, lo4, hi4, chart, Fl = args
    if Fl == 'lg':
        closed, how, nc = closed_lb(lo4, hi4, chart)
        rec = dict(key=list(key), chart=chart, closed=closed, LB=closed, how=how, tv='lg')
        if not nc:
            L, parts = lg_cell(lo4, hi4, chart, closed)
            if parts is not None and L > closed:
                rec['LB'] = L; rec['how'] = 'lg'; rec['lg'] = dict(N=QLG_N, parts=parts)
        return rec
    import tvn as TN
    inv = QC.invariants(lo4, hi4, chart); closed, how = QC.quad_lb(inv)
    fb = QC.quad_fibre(lo4, hi4, chart, inv)                        # fibre bound of the box (2026-11-15)
    if fb > closed:
        closed, how = fb, 'fibre'
    rec = dict(key=list(key), chart=chart, closed=closed, LB=closed, how=how)
    if Fl is None or any(f is None for f in Fl):
        rec['tv'] = 'corner'; return rec
    try:
        Vs = []
        for w in corner_ws(lo4, hi4):
            q = std_q(w, chart); Vs.append([[arb(0), arb(0)], [arb(1), arb(0)], [arb(q[0]), arb(q[1])], [arb(q[2]), arb(q[3])]])
        cn = TN.CellN(Vs, [arb(f) for f in Fl])
    except AssertionError as ex:
        rec['tv'] = 'gate'; return rec
    m = TN.minimize(cn)
    rec['tv'] = dict(lower=m['lower'], boxes=m['boxes'], point=m['point'], Fmin=min(Fl))
    if m['lower'] > rec['LB']:
        rec['LB'] = m['lower']; rec['how'] = 'tv'
    return rec


QFIB = 1.30          # a corner quadrilateral whose fibre bound is >= QFIB * F* gets no eigenvalue certificate


def trivial(lo4, hi4, chart):
    """float/driver heuristic: True if the closed-form bounds already close both columns of the root box (no convex
    member, handed over to an edge chart, or the column test TX for j = 3 and j = 0), so the tangent bound is not needed"""
    import cut_run as CR
    inv = QC.invariants(lo4, hi4, chart)
    if QC.noncvx(inv):
        return True
    if chart == 'std' and QC.defer(inv, lo4, hi4, CR.CFG['H2'], CR.CFG['H3']) is not None:
        return True
    lb = max(QC.quad_lb(inv)[0], QC.quad_fibre(lo4, hi4, chart, inv))
    if chart == 'std':
        lb = max([lb] + list(T4.fallback_bounds(lo4, hi4).values()))
    return all(QC.theta_T(inv, j, lb) > QC.x_threshold(inv, j) for j in (3, 0))


def fibre_corner(q):
    """fibre bound of the corner quadrilateral (0, (1,0), q2, q3) if it is strictly convex (exact test), else 0.0"""
    import fibre as FB
    V = [(0.0, 0.0), (1.0, 0.0), (q[0], q[1]), (q[2], q[3])]; Vf = [(Fraction(a), Fraction(b)) for a, b in V]
    for i in range(4):
        a, b, c = Vf[i], Vf[(i + 1) % 4], Vf[(i + 2) % 4]
        if (b[0] - a[0]) * (c[1] - b[1]) - (b[1] - a[1]) * (c[0] - b[0]) <= 0:
            return 0.0
    return FB.F_lower(V)[0]


FAN_N = int(os.environ.get('PENT_FAN_N', '5'))        # degree of the fan certificate (lgfan.py) of the corner quadrilaterals; 0: off


def qcorner_job(args):
    """value of one corner quadrilateral (2026-11-16): the fan certificate (lgfan.py, Lemma LG; 2 ms); without one the
    fibre bound if it reaches QFIB F*; otherwise the eigenvalue certificate t4table.corner_job as before.
    Returns (key, F lower bound or None, info)."""
    ck, q, settings = args
    if FAN_N:
        import lgfan
        try:
            r = lgfan.F_lower([(0.0, 0.0), (1.0, 0.0), (q[0], q[1]), (q[2], q[3])], FAN_N)
        except Exception:
            r = None
        if r is not None:
            return ck, r['F_lo'], dict(how='fan', lg=r)
        fb = fibre_corner(q)
        if fb >= QFIB * float(R.FSTAR_UP):
            return ck, fb, dict(how='fibre')
    return T4.corner_job(args)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('out'); ap.add_argument('--chart', default='std', choices=('std', 'E2', 'E3'))
    ap.add_argument('--std-keys', default='', help='std chart: take the root cells from this t4table directory'); ap.add_argument('--H', type=float, default=1 / 8)
    ap.add_argument('--sample', type=int, default=0); ap.add_argument('--seed', type=int, default=1); ap.add_argument('--jobs', type=int, default=2)
    ap.add_argument('--all-tv', action='store_true', help='tangent bound for every root (default: only where the closed forms do not close the root)')
    ap.add_argument('--blocks', default=''); ap.add_argument('--old-sample', action='store_true')
    ap.add_argument('--corners', nargs='*', default=[], help='extra corner stores to read (e.g. the t4table one)')
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True); t0 = time.time(); S0 = 1 / 16
    if a.chart == 'std':
        if a.old_sample:                                          # the cells of the old t4table that have a bound (surveys)
            T = T4.Table(a.std_keys, S0)
            roots = [(tuple(k), ) + tuple(T4.cell_box(tuple(k), S0)) for k, c in sorted(T.cells.items()) if c.get('LB') is not None and c['how'] != 'noncvx']
        else:
            roots = [r[:3] for r in QC.std_roots(a.std_keys, S0)]
    else:
        roots = QC.edge_roots(a.chart, a.H, S0)
    if a.blocks:                                                  # the root cells of these blocks (far_runner.block_of)
        keep = set(a.blocks.split(','))
        roots = [r for r in roots if '%d_%d' % ((r[0][1], r[0][2]) if a.chart == 'std' else (r[0][0], r[0][1])) in keep]
    if a.sample:
        cells = [(k, 0.0) for k, lo4, hi4 in roots]; random.Random(a.seed).shuffle(cells)          # the same draw as cut_run.py
        keep = set(k for k, _ in cells[:a.sample]); roots = [r for r in roots if r[0] in keep]
    fc = os.path.join(a.out, 'corners.jsonl'); fl = os.path.join(a.out, 'cells.jsonl'); store = {}
    for f in a.corners + [fc]:
        for d in JS.read_jsonl(f):
            store[tuple(d['key'])] = d['F']
    done = set(tuple(d['key']) for d in JS.read_jsonl(fl)); roots = [r for r in roots if r[0] not in done]
    print('chart %s: %d root cells to do, %d corners known' % (a.chart, len(roots), len(store)), flush=True)
    need = {}; cq = {}; nfib = 0; ntriv = 0; fstar = float(R.FSTAR_UP)
    with open(fc, 'a') as f:
        for k, lo4, hi4 in roots:
            if not a.all_tv and trivial(lo4, hi4, a.chart):
                cq[k] = None; ntriv += 1; continue                 # closed forms suffice: no corner values, no tangent bound
            if QLG:
                cq[k] = 'lg'; continue                             # LG box certificate of the cell: no corner values
            ks = []
            for w in corner_ws(lo4, hi4):
                try:
                    q = std_q(w, a.chart)
                except AssertionError:
                    ks = None; break
                ck = T4.corner_key(q); ks.append(ck)
                if ck not in store and ck not in need:
                    fb = 0.0 if FAN_N else fibre_corner(q)          # with the fan certificate everything is decided in qcorner_job
                    if fb >= QFIB * fstar:                          # elongated corner: the closed-form value is enough
                        store[ck] = fb; nfib += 1; f.write(json.dumps(dict(key=list(ck), F=fb, info=dict(how='fibre'))) + '\n')
                    else:
                        need[ck] = (ck, q, T4.SETTINGS)
            cq[k] = ks
    print('roots closed by closed forms: %d; corners with fibre value: %d; corner values to compute: %d' % (ntriv, nfib, len(need)), flush=True)
    with Pool(a.jobs) as pool:
        with open(fc, 'a') as f:
            for n, (ck, F, info) in enumerate(pool.imap_unordered(qcorner_job, list(need.values()), chunksize=16)):
                store[ck] = F; f.write(json.dumps(dict(key=list(ck), F=F, info=info)) + '\n')
                if n % 200 == 0:
                    f.flush(); print('  corners %d/%d %.0fs' % (n + 1, len(need), time.time() - t0), flush=True)
        jobs = [(k, lo4, hi4, a.chart, None if cq[k] is None else ('lg' if cq[k] == 'lg' else [store.get(ck) for ck in cq[k]])) for k, lo4, hi4 in roots]
        hows = {}
        with open(fl, 'a') as f:
            for rec in pool.imap_unordered(cell_job, jobs, chunksize=8):
                f.write(json.dumps(rec) + '\n'); hows[rec['how']] = hows.get(rec['how'], 0) + 1
    JS.atomic_json(os.path.join(a.out, 'table_info.json'), dict(chart=a.chart, S0=S0, cells=len(done) + len(roots), new_corners=len(need), how=hows, sec=round(time.time() - t0, 1),
                                                                 settings=T4.SETTINGS, method='qtable/2: max(closed forms, LG box certificates)' if QLG else 'qtable/1: max(closed forms, tvn.minimize)'))
    print('done', hows, '%.0fs' % (time.time() - t0), flush=True)


if __name__ == '__main__':
    main()
