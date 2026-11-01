"""G1: rigorous PAQC cell certificate (Claude, 2026-10-30).  Own code; ball arithmetic (python-flint arb).

Claim certified for a cell  C = { P : v0=(0,0), v1=(1,0), (v2,v3,v4) = x in the sup-norm box p +- rho } :
        lambda_1(P) * Area(P)  >  F*_up        for every P in C,
where F*_up is a rigorous upper bound for lambda_1 * Area of the regular pentagon.

Chain (all quantities enclosed in arb balls):
 (0) Reference pentagon P0 (regular, float vertices taken exactly), centroid fan T_k, E_k = [P0_k - c0, P0_{k+1} - c0].
     For a pentagon P with vertex centroid c, D_k = [V_k - c, V_{k+1} - c] E_k^{-1}, J_k = det D_k, G_k = J_k D_k^{-1} D_k^{-T}.
     v4 (after GPT Q7 audit): J_k > 0 alone does NOT give a homeomorphism (regular pentagram is a counterexample).
     Mandatory gate: the 15 half-plane determinants det(V_{i+1}-V_i, V_j-V_i), j not in {i,i+1}, are multiaffine in the
     6 raw coordinates, so their minimum over the cell is attained at the 64 corners; all > 0 => every P in the cell is a
     strictly convex positively oriented pentagon, its vertex centroid is interior and the fan map is a homeomorphism.
     Then lambda_1(P) = inf_v sum_k int_{T_k} grad v' G_k grad v / sum_k J_k int_{T_k} v^2.
 (1) Corner pentagons P_c (64), rigorous F_c^L = lambda_1(P_c)^L * Area(P_c) (MPS + Fox-Henrici-Moler + mean-value L2 bound
     + Hong-Krahn-Szego index, imported from MEAS-01 rig_box.py).  F_min = min_c F_c^L.
 (2) x is affine in u, the vertex map is affine in x, hence D_k(x) = D_k(p) + sum_i (x_i - p_i) B_{k,i} exactly, and the
     multilinear weights a_c(u) reproduce D exactly: D_P = sum_c a_c D_c.  Exact identity (adj is linear for 2x2):
        G_P,k = S_k - tbar_k V_k,   S_k = sum_c a_c G_c,k,  tbar_k = sum_c a_c J_c,k  (= J_P,k, see (3)),
        V_k = sum_c w_c (X_c - M)(X_c - M)' >= 0,  w_c = a_c J_c,k / tbar_k,  X_c = D_c^{-1},  M = sum_c w_c X_c.
     tbar V <= Jmax * Var_a(X) (weighted variance is minimal at its own mean), and the Loewner-order matrix Efron-Stein
     inequality for the product measure a gives  Var_a(X) <= sum_i p_i q_i E[dX_i dX_i'] <= 1/4 sum_i E[dX_i dX_i'],
     dX_i = edge difference along coordinate i.  v3: compare against Gl <= G_c (all corners), so S_k >= Gl; v5: Gl = Gcen - eta I or mu Gcen (lower_metric).
     v4: dX = R_i + d (R_i = mean edge difference), dX dX' <= (1+tau) R R' + (1+1/tau) d d', hence
        eps_k = Jmax [ (1+tau)/4 lam_max(Gl^{-1/2} (sum_i R_i R_i') Gl^{-1/2}) + (1+1/tau)/4 sum_i max_e lam_max(d_e' Gl^{-1} d_e) ]
     (best of a few fixed tau), and  tbar V <= eps Gl <= eps S,  so  G_P,k >= (1 - eps_k) S_k.
 (3) J_k is multiaffine in the raw coordinates (each B_{k,i} has rank 1, det B_{k,i} = 0, and J is quadratic with no
     pure squares), so J_P,k = tbar_k exactly; delta_k = sum_i (h_i/2)^2 |det B_{k,i}| is computed and is 0.
 (4) For every test function v:  N(G_P,v) >= (1 - eps_max) sum_c a_c N(G_c,v) >= (1 - eps_max) sum_c a_c lambda_c D(J_c,v),
        hence  F(P) >= m * F_min * min_k [ I(f_k)(u) / f_k(P) ],  m = min_k (1 - eps_k),  f_k = J_k / A,
        A = sum_j a0_j J_j  (a0_j = |T_j|), and  I(f_k)/f_k >= 1 - err_k / f_k^min,
        err_k = sum_i (h_i^2/8) sup |d^2 f_k / dx_i^2|  (tensor-linear interpolation, h_i = exact side length),
     sup over 64 sub-boxes whose arb enclosures provably cover the cell (exact binary64 endpoints, asserted).
 (5) Certified iff  m * F_min * (1 - max_k err_k/f_k^min) > F*_up, where F*_up bounds lambda_1*Area of the EXACT regular
     pentagon (fstar_exact: float pentagon shrunk by 1e-12, containment checked with arb trig, domain monotonicity).
     Every corner certificate (vertices, K, nseg, lam_tilde, mu, lam_lo/hi, area, HKS index, F_lo) is saved in the output.
 See rig_t1.py for the stronger T1 form (no min over corners, eps depending on u) checked by branch-and-bound.
Usage: python3 rig_paqc.py MODE ...   (see __main__)
"""
import os, sys, json, time, itertools
import numpy as np
from flint import arb, ctx
HERE = os.path.dirname(os.path.abspath(__file__))
C456 = os.path.join(HERE, '..', '..', 'MEAS-01', 'claude', 'scripts', 'C456')
if not os.path.isdir(C456):        # self-contained package layout: rig_box.py / mps.py next to this file
    C456 = HERE
sys.path.insert(0, C456)
import rig_box                      # certify_polygon (rigorous lambda_1 lower bound), HKS index check
from mps import regular

ctx.prec = 96
CORN = np.array(list(itertools.product([-1, 1], repeat=6)), float)
A = lambda x: arb(float(x))

# ---------- reference pentagon (exact float vertices) ----------
P0 = regular(5)
P0a = [[A(v[0]), A(v[1])] for v in P0]
c0 = [sum(v[0] for v in P0a) / 5, sum(v[1] for v in P0a) / 5]


def f_down(x):
    """largest-or-near binary64 that is <= every point of the ball x (outward rounding; GPT PENT-02 finding)."""
    import math
    f = float(x.lower())
    if not (arb(f) <= x.lower()):
        f = math.nextafter(f, -math.inf)
    assert arb(f) <= x.lower()
    return f


def f_up(x):
    """binary64 that is >= every point of the ball x (outward rounding)."""
    import math
    f = float(x.upper())
    if not (arb(f) >= x.upper()):
        f = math.nextafter(f, math.inf)
    assert arb(f) >= x.upper()
    return f


def cert_F_lo(cert):
    """rigorous binary64 lower bound of lambda_1 * area recomputed from a saved corner certificate (arb decimal strings
    are outward enclosures, so the product ball contains the true value).  2026-11-01: the index argument is re-checked
    from the saved lam_tilde, mu and area instead of trusting the saved flag: mu < 1 and lam_tilde/(1-mu)*area <
    2 pi j01^2 (Hong-Krahn-Szego), and lambda_1 >= lam_tilde/(1+mu) is recomputed; the saved lam_lo must agree."""
    # lam_tilde is repr() of the binary64 Helmholtz parameter actually used by the MPS function: parse it as a float
    # (arb('14.03...') would enclose the decimal string, which differs from the binary64 value by < 1/2 ulp)
    lt, mu, ar = arb(float(cert['lam_tilde'])), arb(cert['mu']), arb(cert['area'])
    assert mu < 1, 'FHM: mu >= 1'
    assert lt / (1 - mu) * ar < rig_box.HKS, 'HKS index check'
    lam_lo = lt / (1 + mu)
    assert lam_lo.overlaps(arb(cert['lam_lo'])), 'saved lam_lo inconsistent'
    return f_down(lam_lo * ar)


def inv2(M):
    a, b, c, d = M[0][0], M[0][1], M[1][0], M[1][1]
    det = a * d - b * c
    return [[d / det, -b / det], [-c / det, a / det]], det


def mul2(M, N):
    return [[M[0][0] * N[0][0] + M[0][1] * N[1][0], M[0][0] * N[0][1] + M[0][1] * N[1][1]],
            [M[1][0] * N[0][0] + M[1][1] * N[1][0], M[1][0] * N[0][1] + M[1][1] * N[1][1]]]


E = [[[P0a[k][0] - c0[0], P0a[(k + 1) % 5][0] - c0[0]], [P0a[k][1] - c0[1], P0a[(k + 1) % 5][1] - c0[1]]] for k in range(5)]
Einv = []
a0 = []
for k in range(5):
    Ei, dE = inv2(E[k]); Einv.append(Ei); a0.append(abs(dE) / 2)      # |T_k| = |det E_k| / 2


def verts(x):
    """vertex list (arb) for chart point x (6 arb/float)."""
    xs = [A(v) if not isinstance(v, arb) else v for v in x]
    return [[arb(0), arb(0)], [arb(1), arb(0)], [xs[0], xs[1]], [xs[2], xs[3]], [xs[4], xs[5]]]


def Dmats(V):
    c = [sum(v[0] for v in V) / 5, sum(v[1] for v in V) / 5]
    out = []
    for k in range(5):
        M = [[V[k][0] - c[0], V[(k + 1) % 5][0] - c[0]], [V[k][1] - c[1], V[(k + 1) % 5][1] - c[1]]]
        out.append(mul2(M, Einv[k]))
    return out


def Bmats():
    """B_{k,i} = dD_k/dx_i (constant, exact)."""
    B = [[None] * 6 for _ in range(5)]
    for i in range(6):
        dV = [[arb(0), arb(0)] for _ in range(5)]
        dV[2 + i // 2][i % 2] = arb(1)
        dc = [sum(v[0] for v in dV) / 5, sum(v[1] for v in dV) / 5]
        for k in range(5):
            M = [[dV[k][0] - dc[0], dV[(k + 1) % 5][0] - dc[0]], [dV[k][1] - dc[1], dV[(k + 1) % 5][1] - dc[1]]]
            B[k][i] = mul2(M, Einv[k])
    return B


BK = Bmats()


def det2(M):
    return M[0][0] * M[1][1] - M[0][1] * M[1][0]


def sq(x):
    return x * x          # NB: arb ** 2 returns nan for balls containing 0 (python-flint 0.9)


def fro2(M):
    return sq(M[0][0]) + sq(M[0][1]) + sq(M[1][0]) + sq(M[1][1])


def spec2(M):
    """squared spectral norm of a 2x2 matrix (Loewner-order matrix Efron-Stein: V <= Jmax * 1/4 sum_i max ||dX||_2^2 I)."""
    f = fro2(M); d = det2(M)
    disc = f * f - 4 * d * d
    disc = arb(max(0.0, float(disc.lower())), 0) if disc.lower() < 0 else disc
    r = (f + disc.sqrt()) / 2
    return arb(min(r.upper(), f.upper()))   # also <= Frobenius^2


def gen_lam_max(M, Gm):
    """largest root of det(M - lam*Gm) = 0 for symmetric 2x2 M >= 0 and Gm > 0, i.e. lam_max(Gm^{-1/2} M Gm^{-1/2})."""
    m11, m12, m22 = M[0][0], (M[0][1] + M[1][0]) / 2, M[1][1]
    g11, g12, g22 = Gm[0][0], (Gm[0][1] + Gm[1][0]) / 2, Gm[1][1]
    dG = g11 * g22 - g12 * g12
    dM = m11 * m22 - m12 * m12
    b = m11 * g22 + m22 * g11 - 2 * m12 * g12
    disc = b * b - 4 * dG * dM
    if disc.lower() < 0:
        disc = arb(0).union(disc)
        disc = arb(max(disc.upper(), 0.0))     # roots are real; widen to [0, upper]
    return (b + disc.sqrt()) / (2 * dG)


def lam_max_sym2(M):
    a, b, d = M[0][0], (M[0][1] + M[1][0]) / 2, M[1][1]
    return (a + d) / 2 + (sq((a - d) / 2) + b * b).sqrt()


def spec_sym2(M):
    """spectral norm of a symmetric 2x2 matrix."""
    a, b, d = M[0][0], (M[0][1] + M[1][0]) / 2, M[1][1]
    r = (sq((a - d) / 2) + b * b).sqrt()
    m = abs((a + d) / 2) + r
    return arb(m.upper())


def amax(xs):
    m = xs[0]
    for x in xs[1:]:
        m = arb(max(m.upper(), x.upper()))
    return arb(m.upper())


def amin(xs):
    m = xs[0]
    for x in xs[1:]:
        m = arb(min(m.lower(), x.lower()))
    return arb(m.lower())


def lam_min_sym2(M):
    a, b, d = M[0][0], (M[0][1] + M[1][0]) / 2, M[1][1]
    return (a + d) / 2 - (sq((a - d) / 2) + b * b).sqrt()


def lower_metric(Gc):
    """A matrix Gl with Gl <= G_c (Loewner) for every corner metric G_c, as large as practical (v5, 2026-10-31).
    Candidates: (a) Gcen - eta I, eta = max_c ||G_c - Gcen||_2 (absolute; bad for anisotropic G);
                (b) mu Gcen, mu = min_c lam_min(Gcen^{-1} G_c) (relative; always PD when all G_c are PD).
    Both satisfy Gl <= G_c.  Pick the one with the larger relative size lam_min(Gcen^{-1} Gl)  (a choice, not a bound).
    Returns (Gl, which) or (None, why)."""
    Gcen = [[sum(G[a][b] for G in Gc) / len(Gc) for b in range(2)] for a in range(2)]
    g11, g12, g22 = Gcen[0][0], (Gcen[0][1] + Gcen[1][0]) / 2, Gcen[1][1]
    dC = g11 * g22 - g12 * g12
    lmin_cen = lam_min_sym2(Gcen)
    cands = []
    eta = amax([spec_sym2([[G[a][b] - Gcen[a][b] for b in range(2)] for a in range(2)]) for G in Gc])
    if lmin_cen - eta > 0:
        cands.append((1 - eta / lmin_cen, [[Gcen[0][0] - eta, Gcen[0][1]], [Gcen[1][0], Gcen[1][1] - eta]], 'abs'))
    mus = []
    for G in Gc:                       # smallest root of det(G - mu Gcen) = 0
        h11, h12, h22 = G[0][0], (G[0][1] + G[1][0]) / 2, G[1][1]
        b = h11 * g22 + h22 * g11 - 2 * h12 * g12
        dG = h11 * h22 - h12 * h12
        disc = b * b - 4 * dC * dG
        if disc.lower() < 0:
            disc = arb(0).union(disc); disc = arb(max(disc.upper(), 0.0))
        mus.append((b - disc.sqrt()) / (2 * dC))
    mu = amin(mus)
    if mu > 0:
        mu = arb(mu.lower())
        cands.append((mu, [[mu * Gcen[0][0], mu * Gcen[0][1]], [mu * Gcen[1][0], mu * Gcen[1][1]]], 'rel'))
    if not cands:
        return None, 'no PD lower metric'
    best = max(cands, key=lambda t: float(t[0].lower()))
    return best[1], best[2]


def cell_endpoints(p, rho):
    """The cell is the box with EXACT binary64 endpoints lo_i = fl(p_i - rho), hi_i = fl(p_i + rho); corners use them."""
    lo = np.array([float(np.float64(p[i]) - np.float64(rho)) for i in range(6)])
    hi = np.array([float(np.float64(p[i]) + np.float64(rho)) for i in range(6)])
    return lo, hi


def cell_corners(lo, hi):
    return [np.where(s > 0, hi, lo) for s in CORN]


def convexity_gate(lo, hi):
    """15 half-plane determinants at the 64 corners (multiaffine => corner minimum is the cell minimum)."""
    mn = None
    for x in cell_corners(lo, hi):
        V = verts(x)
        for i in range(5):
            a, b = V[i], V[(i + 1) % 5]
            for j in range(5):
                if j in (i, (i + 1) % 5):
                    continue
                c = V[j]
                d = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
                mn = d if mn is None else arb(min(mn.lower(), d.lower()))
    return arb(mn.lower())


def geometry_bound(p, rho):
    """Everything in steps (2)-(4) except the corner eigenvalues.  Returns factor such that F(P) >= factor * F_min."""
    lo, hi = cell_endpoints(p, rho)
    cg = convexity_gate(lo, hi)
    if not (cg > 0):
        return dict(ok=False, why='convexity gate failed')
    corners = cell_corners(lo, hi)
    Ds = [Dmats(verts(x)) for x in corners]
    out = {'convexity_gate_min': float(cg.lower())}
    m_list = []; err_ratio = []
    for k in range(5):
        Jc = [det2(D[k]) for D in Ds]
        if not all(j > 0 for j in Jc):
            return dict(ok=False, why='corner J<=0')
        Xc = [inv2(D[k])[0] for D in Ds]
        Gc = [[[Jc[c] * (sq(Xc[c][0][0]) + sq(Xc[c][0][1])), Jc[c] * (Xc[c][0][0] * Xc[c][1][0] + Xc[c][0][1] * Xc[c][1][1])],
               [Jc[c] * (Xc[c][0][0] * Xc[c][1][0] + Xc[c][0][1] * Xc[c][1][1]), Jc[c] * (sq(Xc[c][1][0]) + sq(Xc[c][1][1]))]]
              for c in range(64)]
        gmin = amin([lam_min_sym2(G) for G in Gc])
        Jmax = amax(Jc); tmin = amin(Jc)
        # v3: anisotropy-aware relative bound.  Gl = Gcen - eta*I <= G_c for all corners (eta = max_c ||G_c - Gcen||_2),
        # so S_k >= Gl, and  lambda_max(Gl^{-1/2} V Gl^{-1/2}) <= Jmax * 1/4 sum_i max_e lambda_max(dX_e' Gl^{-1} dX_e).
        Gl, which = lower_metric(Gc)                 # v5: Gcen - eta I or mu Gcen, whichever is relatively larger
        if Gl is None or not (lam_min_sym2(Gl) > 0):
            return dict(ok=False, why='Gl not PD')
        Gli = inv2(Gl)[0]
        # v4: joint (Loewner) sum.  For each direction i let R_i be the mean edge difference and d_e = dX_e - R_i.
        # dX dX' <= (1+tau) R R' + (1+1/tau) d d'  (any tau>0), so after congruence by Gl^{-1/2}
        # eps_k <= Jmax * [ (1+tau)/4 lam_max(sum_i R_i' Gl^{-1} R_i)  +  (1+1/tau)/4 sum_i max_e lam_max(d_e' Gl^{-1} d_e) ].
        Rsum = [[arb(0), arb(0)], [arb(0), arb(0)]]
        dev = arb(0)
        for i in range(6):
            bit = 1 << (5 - i)                 # CORN is itertools.product order: coordinate 0 is the most significant bit
            dXs = [[[Xc[c][a][b] - Xc[c ^ bit][a][b] for b in range(2)] for a in range(2)] for c in range(64) if not (c & bit)]
            Rm = [[sum(d[a][b] for d in dXs) / 32 for b in range(2)] for a in range(2)]
            Rt = [[Rm[0][0], Rm[1][0]], [Rm[0][1], Rm[1][1]]]
            T = mul2(Rm, Rt)                   # R R'  (V is built from (X-M)(X-M)'; orientation matters for the SUM)
            Rsum = [[Rsum[a][b] + T[a][b] for b in range(2)] for a in range(2)]
            vals = []
            for d in dXs:
                dd = [[d[a][b] - Rm[a][b] for b in range(2)] for a in range(2)]
                ddt = [[dd[0][0], dd[1][0]], [dd[0][1], dd[1][1]]]
                vals.append(lam_max_sym2(mul2(mul2(ddt, Gli), dd)))
            dev += amax(vals)
        main = gen_lam_max(Rsum, Gl)           # lam_max(Gl^{-1/2} (sum_i R_i R_i') Gl^{-1/2})
        best = None
        for tau in (0.003, 0.01, 0.03, 0.1, 0.3, 1.0):
            v = (1 + arb(tau)) / 4 * main + (1 + 1 / arb(tau)) / 4 * dev
            best = v if best is None else arb(min(best.upper(), v.upper()))
        eps = Jmax * best
        delta = sum(sq(hw_ / 2) * abs(det2(BK[k][i])) for i, hw_ in enumerate([A(hi[j]) - A(lo[j]) for j in range(6)]))
        if not (tmin - delta > 0):
            return dict(ok=False, why='J may vanish')
        m_k = (1 - eps) * tmin / (tmin + delta)
        m_list.append(m_k)
        out[f'k{k}'] = dict(eps=float(eps.upper()), delta_over_t=float((delta / tmin).upper()), gmin=float(gmin.lower()))
    m = amin(m_list)
    # step (4): interpolation error of f_k = J_k / Aw over the cell (interval evaluation on the whole box)
    errmax = arb(0)
    mid = [float((lo[i] + hi[i]) / 2) for i in range(6)]
    hw = [A(hi[i]) - A(lo[i]) for i in range(6)]            # exact side lengths h_i
    for sub in CORN:                           # sup over 64 sub-boxes built from EXACT endpoints (v4: arb unions)
        xbox = [(A(lo[i]).union(A(mid[i])) if sub[i] < 0 else A(mid[i]).union(A(hi[i]))) for i in range(6)]
        for i in range(6):                     # assert coverage of the child interval
            a_, b_ = (lo[i], mid[i]) if sub[i] < 0 else (mid[i], hi[i])
            assert xbox[i].lower() <= a_ and xbox[i].upper() >= b_
        r_ = _interp_err_box(xbox, hw)
        if r_ is None:
            return dict(ok=False, why='J box not positive')
        errmax = arb(max(errmax.upper(), r_.upper()))
    factor = m * (1 - errmax)
    out.update(m=float(m.lower()), interp_err=float(errmax.upper()), factor=float(factor.lower()))
    return dict(ok=True, factor=arb(factor.lower()), info=out)


def _interp_err_box(xbox, hw):
    Dbox = Dmats(verts(xbox))
    Jb = [det2(Dbox[k]) for k in range(5)]
    Aw = sum(a0[j] * Jb[j] for j in range(5))
    for k in range(5):
        if not (Jb[k] > 0):
            return None
    errmax = arb(0)
    for k in range(5):
        fmin = arb((Jb[k] / Aw).lower())
        err = arb(0)
        for i in range(6):
            # derivatives along x_i: J_j(t) = J_j + t*tr(adj(D_j)B_ji) + t^2 det(B_ji)
            Ji = [Dbox[j][1][1] * BK[j][i][0][0] - Dbox[j][1][0] * BK[j][i][0][1] - Dbox[j][0][1] * BK[j][i][1][0] + Dbox[j][0][0] * BK[j][i][1][1]
                  for j in range(5)]
            Jii = [2 * det2(BK[j][i]) for j in range(5)]
            Ai = sum(a0[j] * Ji[j] for j in range(5)); Aii = sum(a0[j] * Jii[j] for j in range(5))
            f2 = Jii[k] / Aw - 2 * Ji[k] * Ai / sq(Aw) - Jb[k] * Aii / sq(Aw) + 2 * Jb[k] * sq(Ai) / (sq(Aw) * Aw)
            err += sq(hw[i]) / 8 * abs(f2)             # tensor-linear interpolation: h_i^2/8 sup|f_ii|
        errmax = arb(max(errmax.upper(), (err / fmin).upper()))
    return errmax


def corner_F(x, K=8, nseg=500):
    V = np.array([[0, 0], [1, 0], [x[0], x[1]], [x[2], x[3]], [x[4], x[5]]], float)
    from mps import poly_area
    guess = 18.94 / poly_area(V)
    r = rig_box.certify_polygon(V, K, guess, nseg=nseg, width=0.05)
    if not r['ok_index']:
        return None, r
    return r['lam_lo'] * r['area'], r


def certify_cell(p, rho, K=8, nseg=1000, workers=1):
    t0 = time.time()
    g = geometry_bound(p, rho)
    res = dict(rho=rho, p=[float(v) for v in p])
    if not g['ok']:
        res.update(status='fail_geometry', why=g['why']); return res
    res['geometry'] = g['info']
    lo, hi = cell_endpoints(p, rho)
    corners = cell_corners(lo, hi)
    if workers > 1:
        from multiprocessing import Pool
        with Pool(workers) as pool:
            outs = pool.starmap(_cf, [(list(x), K, nseg) for x in corners])
    else:
        outs = [_cf(list(x), K, nseg) for x in corners]
    if any(o is None for o in outs):
        res.update(status='fail_corner_index'); return res
    res['cell_lo'] = [repr(float(v)) for v in lo]; res['cell_hi'] = [repr(float(v)) for v in hi]
    res['corner_certificates'] = [o[2] for o in outs]
    Fs = [arb(o[0], o[1]) for o in outs]
    Fmin = amin(Fs)
    LB = g['factor'] * Fmin
    res.update(Fmin=float(Fmin.lower()), LB=float(LB.lower()), Fstar_up=float(FSTAR_UP.upper()),
               certified=bool(LB > FSTAR_UP), sec=round(time.time() - t0, 1))
    return res


def _cf(x, K, nseg):
    F, r = corner_F(np.array(x), K, nseg)
    if F is None:
        return None
    lo_ = f_down(F)                          # v6: rigorous LOWER endpoint, rounded DOWN to binary64 (was round-to-nearest)
    cert = dict(vertices_x=[repr(float(v)) for v in x], K=K, nseg=nseg, lam_tilde=repr(float(r['lam'])),
                mu=str(r['mu']), lam_lo=str(r['lam_lo']), lam_hi=str(r['lam_hi']), area=str(r['area']),
                hks_index_ok=bool(r['ok_index']), F_lo=repr(lo_),
                coefficients=r.get('coefficients'), sup_bdry_upper=repr(f_up(r['sup_bdry'])) if 'sup_bdry' in r else None,
                l2sq_lower=repr(f_down(r['l2sq_lower'])) if 'l2sq_lower' in r else None)   # v6: replayable MPS data
    return (lo_, 0.0, cert)


def fstar_exact(K=8, nseg=500, shrink="1e-12"):
    """Rigorous upper bound of lambda_1*Area for the EXACT regular pentagon R (unit circumradius, vertex angle 2 pi j/5 + rot).
    Pt = float vertices from mps.regular(5).  s*Pt (s = 1 - shrink, about the origin) is checked to lie inside R with exact
    arb trigonometry; domain monotonicity gives lambda_1(R) <= lambda_1(s Pt) = lambda_1(Pt)/s^2 <= lam_hi(Pt)/s^2."""
    from mps import poly_area
    Pt = regular(5)
    rot = float(np.arctan2(Pt[0][1], Pt[0][0]))
    th0 = arb.atan2(A(Pt[0][1]), A(Pt[0][0]))
    Rv = [[(th0 + 2 * arb.pi() * j / 5).cos(), (th0 + 2 * arb.pi() * j / 5).sin()] for j in range(5)]
    s_ = 1 - arb(shrink)
    for i in range(5):
        a, b = Rv[i], Rv[(i + 1) % 5]
        for v in Pt:
            x, y = s_ * A(v[0]), s_ * A(v[1])
            cr = (b[0] - a[0]) * (y - a[1]) - (b[1] - a[1]) * (x - a[0])
            assert cr > 0, 'containment failed'
    r = rig_box.certify_polygon(Pt, K, 18.9191042862 / poly_area(Pt), nseg=nseg)
    assert r['ok_index']
    AR = arb(5) / 2 * (2 * arb.pi() / 5).sin()
    Fup = r['lam_hi'] / (s_ * s_) * AR
    return Fup, dict(K=K, nseg=nseg, shrink=shrink, theta0=str(th0), lam_tilde=repr(float(r['lam'])), mu=str(r['mu']),
                     lam_hi=str(r['lam_hi']), area_exact=str(AR), F_up=str(Fup), F_up_upper=repr(f_up(Fup)))


_FS = os.path.join(HERE, '..', 'G1', 'fstar_exact.json')
for _alt in (os.path.join(HERE, 'fstar_exact.json'), os.path.join(HERE, '..', 'certificates', 'fstar_exact.json')):
    if not os.path.exists(_FS):
        _FS = _alt                                          # package / release layouts
if os.path.exists(_FS):
    _fj = json.load(open(_FS))
    if 'F_up_upper_outward' in _fj:
        FSTAR_UP = arb(_fj['F_up_upper_outward'])
    else:                                   # v6: the stored float was round-to-nearest; derive an outward one from the ball
        FSTAR_UP = arb(f_up(arb(_fj['F_up'])))
        import math as _m
        FSTAR_UP = arb(max(float(FSTAR_UP), _m.nextafter(float(_fj['F_up_upper']), _m.inf)))
else:
    FSTAR_UP = None

if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'fstar':
        K = int(sys.argv[2]) if len(sys.argv) > 2 else 8; nseg = int(sys.argv[3]) if len(sys.argv) > 3 else 500
        Fup, info = fstar_exact(K, nseg)
        print(json.dumps(info))
        if '--save' in sys.argv:
            json.dump(info, open(_FS, 'w'), indent=1); print('saved', _FS)
    elif mode == 'geom':          # quick: geometry factor only, at chart point p (json list) for several rho
        p = np.array(json.loads(sys.argv[2]))
        for rho in [float(x) for x in sys.argv[3:]]:
            g = geometry_bound(p, rho)
            print(rho, g.get('why'), g.get('info'), flush=True)
    elif mode == 'cell':
        p = np.array(json.loads(sys.argv[2])); rho = float(sys.argv[3]); w = int(sys.argv[4]) if len(sys.argv) > 4 else 1
        print(json.dumps(certify_cell(p, rho, workers=w)), flush=True)
