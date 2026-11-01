"""Export the MPS data (vertices, K, lam_tilde, 40 coefficients) of EVERY corner of a certificate, for second-party
replay (Claude, 2026-10-31).  Only the MPS fit is recomputed (deterministic, seed 0); lam_tilde must reproduce bit-for-bit.
The sup / L2 / mu values are NOT recomputed here: the replayer computes them for this explicit function.
Usage: python3 export_all_mps.py CERT.json OUT.json"""
import sys, json
import numpy as np
import rig_paqc as R
import rig_box as B
from mps import PolyMPS, poly_area
src, out = sys.argv[1], sys.argv[2]
d = json.load(open(src))
recs = []
for idx, c in enumerate(d['corner_certificates']):
    x = np.array([float(v) for v in c['vertices_x']])
    V = np.array([[0, 0], [1, 0], [x[0], x[1]], [x[2], x[3]], [x[4], x[5]]], float)
    m = PolyMPS(V, K=c['K'])
    lam = m.solve(18.94 / poly_area(V), width=0.05, ngrid=11)
    assert repr(float(lam)) == c['lam_tilde'], (idx, lam, c['lam_tilde'])
    U = B.ArbU(m)
    recs.append(dict(corner_index=idx, vertices=[[repr(float(a)) for a in v] for v in V], K=c['K'], lam_tilde=repr(float(lam)),
                     coefficients=[repr(float(v)) for v in U.c], producer_F_lo_v6=repr(R.cert_F_lo(c)), producer_mu=c['mu']))
json.dump(dict(source=src, n=len(recs), basis='as in G1/replay_corners/*.json (export_corner_mps.py docstring)', corners=recs), open(out, 'w'))
print(src, len(recs), 'corners exported')
