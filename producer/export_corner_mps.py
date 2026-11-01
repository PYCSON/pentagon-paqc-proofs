"""Export a corner eigenvalue certificate in REPLAYABLE form (Claude, 2026-10-31): vertices, K, lam_tilde and the MPS
coefficients c, so that a second party can recompute mu = sqrt|P| sup_bdry|u| / ||u||_2 for the SAME explicit function
   u(x) = sum_{j=0..4} sum_{m=1..K} c[j*K+m-1] J_{nu_jm}(k r_j) sin(nu_jm th_j),   k = sqrt(lam_tilde),
   nu_jm = m pi / alpha_j,  alpha_j = interior angle at vertex j,  (r_j, th_j) polar coords at vertex j with th_j = 0 on
   the edge (V_j -> V_{j+1}) and th_j = alpha_j on the edge (V_j -> V_{j-1}).
u solves Delta u + lam_tilde u = 0 in P exactly (each term does).  Fox-Henrici-Moler: some Dirichlet eigenvalue lam_n has
|lam_tilde - lam_n| <= mu lam_n.  Any rescaling of the basis (e.g. Gamma(nu+1) J_nu) gives another valid certificate if the
SAME function is used for the sup and the L2 norm.  Usage: python3 export_corner_mps.py CERT.json INDEX OUT.json"""
import sys, json
import numpy as np
from flint import arb
import rig_paqc as R
import rig_box as B
from mps import PolyMPS, poly_area

src, idx, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
d = json.load(open(src))
c = d['corner_certificates'][idx]
x = np.array([float(v) for v in c['vertices_x']])
V = np.array([[0, 0], [1, 0], [x[0], x[1]], [x[2], x[3]], [x[4], x[5]]], float)
K, nseg = c['K'], c['nseg']
m = PolyMPS(V, K=K)
lam = m.solve(18.94 / poly_area(V), width=0.05, ngrid=11)
assert repr(float(lam)) == c['lam_tilde'], ('MPS not reproduced', lam, c['lam_tilde'])
U = B.ArbU(m)
sup = B.sup_boundary(U, nseg)
l2 = B.l2_lower(U)
xs = [B.A(v[0]) for v in V]; ys = [B.A(v[1]) for v in V]
area = abs(sum(xs[i] * ys[(i + 1) % 5] - xs[(i + 1) % 5] * ys[i] for i in range(5))) / 2
mu = area.sqrt() * sup / l2.sqrt()
lamt = arb(float(lam))
rec = dict(source=src, corner_index=idx, vertices=[[repr(float(a)) for a in v] for v in V], K=K, nseg_used=nseg,
           lam_tilde=repr(float(lam)), coefficients=[repr(float(v)) for v in U.c], coefficient_order='c[j*K + (m-1)], j = vertex 0..4, m = 1..K',
           basis='J_{nu}(k r_j) sin(nu th_j), nu = m pi/alpha_j, th_j measured from edge V_j->V_{j+1} towards V_j->V_{j-1}; Bessel J of the first kind (unnormalised)',
           producer_sup_boundary_upper=repr(R.f_up(sup)), producer_l2sq_lower=repr(R.f_down(l2)),
           producer_mu_upper=repr(R.f_up(mu)), area=str(area),
           producer_lam_lo=repr(R.f_down(lamt / (1 + mu))), producer_F_lo=repr(R.f_down(lamt / (1 + mu) * area)),
           hks_check='lam_tilde/(1-mu)*area < 2 pi j01^2 = %s : %s' % (str(B.HKS), bool(lamt / (1 - mu) * area < B.HKS)),
           l2_method='mean-value identity on the disk B(x0,R), x0 = vertex centroid, R = 0.999 * min distance to edge lines: int_B u = u(x0) 2 pi R J_1(kR)/k, ||u||^2 >= (int_B u)^2/(pi R^2)',
           sup_method='Taylor-2 on nseg segments per edge (value and slope at segment midpoint, interval 2nd derivative); the two vertex expansions of the edge endpoints vanish identically on that edge and are skipped')
json.dump(rec, open(out, 'w'), indent=1)
print(json.dumps({k: rec[k] for k in ('corner_index', 'lam_tilde', 'producer_mu_upper', 'producer_F_lo', 'hks_check')}))
print('cert_F_lo (v6) of the saved certificate:', repr(R.cert_F_lo(c)))
