"""Producer-side replay of all cell certificates in ../certificates from their saved corner certificates
(corner eigenvalue certificates are NOT recomputed here; see second_party/ for the independent corner checker).
For each cell: convexity gate, geometry (min-corner factor or the T1 branch-and-bound), final comparison with F*_up.
Usage (from producer/): python3 verify_certificates.py      (about 15 CPU minutes; needs python-flint 0.9, numpy)"""
import os, json, time
import numpy as np
from flint import arb
import rig_paqc as R, rig_t1 as T

CERT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'certificates')
JOBS = [('cell_liu_point_rho0.025_v6.json', 'mincorner', None), ('cell_r0.1_id71_rho0.012_v6.json', 'mincorner', None),
        ('t1_liu_point_rho0.03_rigorous_v6.json', 't1', 'es'), ('t1_r0.1_id71_rho0.015_rigorous_v6.json', 't1', 'es'),
        ('t1_id301_rho0.0239_rigorous_v6.json', 't1', 'both'), ('t1_r0.1_id71_rho0.012_rigorous_v6.json', 't1', 'es')]
print('F*_up =', float(R.FSTAR_UP))
allok = True
for fn, kind, em in JOBS:
    t0 = time.time()
    d = json.load(open(os.path.join(CERT, fn)))
    p = np.array(d['p']); rho = d['rho']
    lo, hi = R.cell_endpoints(p, rho)
    assert [repr(float(v)) for v in lo] == d['cell_lo'] and [repr(float(v)) for v in hi] == d['cell_hi'], 'cell endpoints'
    for c, x in zip(d['corner_certificates'], R.cell_corners(lo, hi)):
        assert [float(v) for v in c['vertices_x']] == [float(v) for v in x], 'corner order'
    if kind == 'mincorner':
        Fl = [arb(R.cert_F_lo(c)) for c in d['corner_certificates']]
        g = R.geometry_bound(p, rho)
        LB = g['factor'] * R.amin(Fl)
        ok = bool(LB > R.FSTAR_UP)
        print('%-42s min-corner  factor %.10f  LB %.12f  certified %s  (%.0f s)' % (fn, float(g['factor'].lower()), float(LB.lower()), ok, time.time() - t0), flush=True)
    else:
        tmp = '/tmp/_corners_' + fn
        json.dump(d, open(tmp, 'w'))
        res = T.certify_cell_T1(p, rho, corners_from=tmp, emode=em)
        ok = bool(res['certified'])
        print('%-42s T1 (%s)  boxes %d  leaves %d  min leaf lower %.12f  certified %s  (%.0f s)' % (
            fn, em, res['T1']['boxes'], res['T1']['accepted'], res['T1'].get('min_lower', float('nan')), ok, time.time() - t0), flush=True)
    allok = allok and ok
print('ALL CELLS CERTIFIED' if allok else 'SOME CELL FAILED')
