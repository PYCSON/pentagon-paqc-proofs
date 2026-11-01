"""Re-certify all G1 certificates with the v6 code (outward rounding fixes from GPT PENT-02) from their SAVED corner
certificates (no corner eigenvalue recomputation).  Claude, 2026-10-31.  Usage: python3 recertify_g1.py OUTDIR"""
import sys, os, json, time
import numpy as np
from flint import arb
import rig_paqc as R, rig_t1 as T
G1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'G1')
if not os.path.isdir(G1):
    G1 = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'certificates')   # release layout
out = sys.argv[1]
jobs = [('cell_liu_point_rho0.025_v4.json', 'mincorner', None), ('cell_r0.1_id71_rho0.012_v4.json', 'mincorner', None),
        ('t1_r0.1_id71_rho0.012_rigorous.json', 't1', 'es'), ('t1_r0.1_id71_rho0.015_rigorous.json', 't1', 'es'),
        ('t1_liu_point_rho0.03_rigorous.json', 't1', 'es'), ('t1_id301_rho0.0239_rigorous.json', 't1', 'both')]
summary = []
for fn, kind, em in jobs:
    t0 = time.time()
    d = json.load(open(os.path.join(G1, fn)))
    p = np.array(d['p']); rho = d['rho']
    if kind == 'mincorner':
        lo, hi = R.cell_endpoints(p, rho)
        for c, x in zip(d['corner_certificates'], R.cell_corners(lo, hi)):
            assert [float(v) for v in c['vertices_x']] == [float(v) for v in x]
        Fl = [R.cert_F_lo(c) for c in d['corner_certificates']]
        g = R.geometry_bound(p, rho)
        LB = g['factor'] * R.amin([arb(v) for v in Fl])
        res = dict(d, geometry=g['info'], Fmin=min(Fl), LB=R.f_down(LB), Fstar_up=float(R.FSTAR_UP),
                   certified=bool(LB > R.FSTAR_UP), version='v6', corner_F_lo_v6=[repr(v) for v in Fl])
    else:
        tmp = os.path.join(out, 'corners_' + fn)
        json.dump(d, open(tmp, 'w'))
        res = T.certify_cell_T1(p, rho, corners_from=tmp, emode=em)
        res['version'] = 'v6'
    res['recertified'] = 'v6 (outward rounding: hull via exact arb endpoints, F_lo rounded down, F*_up rounded up), from saved corner certificates of ' + fn
    res['recert_sec'] = round(time.time() - t0, 1)
    name = fn.replace('_v4.json', '_v6.json').replace('_rigorous.json', '_rigorous_v6.json')
    json.dump(res, open(os.path.join(out, name), 'w'))
    row = dict(file=name, certified=res['certified'], sec=res['recert_sec'])
    if kind == 'mincorner':
        row.update(LB=res['LB'], margin=res['LB'] - float(R.FSTAR_UP))
    else:
        row.update(boxes=res['T1']['boxes'], min_lower=res['T1'].get('min_lower'), margin=(res['T1'].get('min_lower') or 0) - float(R.FSTAR_UP))
    summary.append(row); print(json.dumps(row), flush=True)
json.dump(dict(Fstar_up_outward=float(R.FSTAR_UP), rows=summary), open(os.path.join(out, 'recertify_summary.json'), 'w'), indent=1)
