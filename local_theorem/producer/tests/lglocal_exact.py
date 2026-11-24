#!/usr/bin/env python3
"""exact tests of lglocal.Local against the formulas of Proposition R (far/LOCAL_LEMMA.md).  Run from v2_run:
    python3 tests/lglocal_exact.py [N]
1. With F0 fixed to a point (debug option), the polynomial K D and the linear functional applied to an arbitrary test
   function agree with the pointwise formulas evaluated in exact rational arithmetic.
2. Cauchy-Schwarz: for the test function e (an arbitrary H^1_0 polynomial) every coefficient of degree <= 3 of the linear
   part satisfies  |coefficient| <= |grad e| NL_alpha + (1/2) |e| NT_alpha  with the dual norms computed by Local."""
import sys
from fractions import Fraction as Fr
import os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import lglocal as LL, lgbox as LB, lgpoly as LP, near_lg as NL
from flint import fmpq
N = int(sys.argv[1]) if len(sys.argv) > 1 else 5
fam, ce, rng = LL._test_family(N, 1); X = LL.Local(N, test_list=ce, verbose=False, debug_point_F0=True); q = X.q; qq = X.qq; NV = LL.NV; ND = LL.ND
fr = lambda x: Fr(int(x.p), int(x.q)); mid = lambda b: LL.ball(b)[0]
G = [[[fr(v) for v in row] for row in g] for g in X.GR]; Mm = [[[fr(X.tb[k]['m'][a][b]) for b in range(qq)] for a in range(qq)] for k in range(NV)]
KS = fr(LL.KSQ); F0 = mid(X.F0); lam0 = mid(X.lam0); Cq = [Fr(v) for v in NL.CF]; O = (Fr(fam.O[0]), Fr(fam.O[1]))
gc = [tuple(mid(x) for x in g) for g in X.gc]
def geo(z):
    V = [(Fr(0), Fr(0)), (Fr(1), Fr(0)), (z[0], z[1]), (z[2], z[3]), (z[4], z[5])]; out = []
    for k in range(NV):
        e1 = (V[k][0] - O[0], V[k][1] - O[1]); e2 = (V[(k + 1) % NV][0] - O[0], V[(k + 1) % NV][1] - O[1])
        out.append((e1[0] ** 2 + e1[1] ** 2, e1[0] * e2[0] + e1[1] * e2[1], e2[0] ** 2 + e2[1] ** 2, e1[0] * e2[1] - e1[1] * e2[0]))
    return out
Z = [Fr(0)] * (6 * qq)
def vec(I, w):
    v = list(Z)
    for b in range(qq):
        v[I * qq + b] = w[b]
    return v
lin2 = lambda c1, v1, c2, v2: [c1 * a + c2 * b for a, b in zip(v1, v2)]
def comb(*terms):
    v = list(Z)
    for c, u in terms:
        v = [a + c * b for a, b in zip(v, u)]
    return v
def ip(Gk, x, y):
    nzx = [(i, a) for i, a in enumerate(x) if a]; nzy = [(j, b) for j, b in enumerate(y) if b]
    return sum(a * Gk[i][j] * b for i, a in nzx for j, b in nzy)
def run(d):
    z = [c + x for c, x in zip(Cq, d)]; gz = geo(z); area = sum(g[3] for g in gz) / 2; lh = F0 / area; DD = Fr(1)
    for g in gz:
        DD *= g[3]
    w = [Fr(1)] + list(d) + [Fr(0)]; w0 = [Fr(1)] + [Fr(0)] * (qq - 1); we = [Fr(0)] * q + [Fr(1)]
    psi = k2 = k3 = lin = lth = Fr(0)
    for k in range(NV):
        Gk = G[k]; dot = lambda X_, Y_: ip(Gk, X_[0], Y_[0]) + ip(Gk, X_[1], Y_[1])
        def fields(g, x, l):
            n11, n12, n22, det = g; ah, bh, xa, xb, ua, ub = [vec(I, x) for I in range(6)]
            sa = comb((-det, ah), (Fr(1), xb)); sb = comb((-det, bh), (Fr(-1), xa))
            return (sa, sb), (ua, ub), (comb((n11 / det, sa), (n12 / det, sb), (-1 / l, ua)), comb((n12 / det, sa), (n22 / det, sb), (-1 / l, ub)))
        Mof = lambda g: ((g[2] / g[3], -g[1] / g[3]), (-g[1] / g[3], g[0] / g[3])); Nof = lambda g: ((g[0] / g[3], g[1] / g[3]), (g[1] / g[3], g[2] / g[3]))
        ap = lambda A_, x: (comb((A_[0][0], x[0]), (A_[0][1], x[1])), comb((A_[1][0], x[0]), (A_[1][1], x[1])))
        g = gz[k]; det = g[3]; Mz, Nz, Mc = Mof(g), Nof(g), Mof(gc[k]); ddet = det - gc[k][3]
        s, gv, rho = fields(g, w, lh); _, _, rhoc = fields(gc[k], w0, lam0); mrc = ap(Mc, rhoc)
        bz = det * sum(w[a] * Mm[k][a][b] * w[b] for a in range(qq) for b in range(qq))
        psi += area * (KS + F0) * bz - area ** 2 * dot(gv, ap(Mz, gv)) - F0 * KS * dot(s, ap(Nz, s))
        k2 += 2 * F0 * KS * dot(rho, mrc); k3 += F0 * KS * dot(mrc, ap(Nz, mrc))
        ge = (vec(4, we), vec(5, we)); ye = (vec(0, we), vec(1, we)); mr = ap(Mc, rho); mz = ap(Mz, rho); mnm = ap(Mc, ap(Nz, mrc)); nm = ap(Nz, mrc)
        Lam = tuple(comb((2 * F0 * KS / lam0, mr[i]), (-2 * F0 * area, mz[i]), (-2 * F0 * KS / lam0, mnm[i])) for i in range(2))
        Th = tuple(comb((-2 * F0 * KS * ddet, rho[i]), (2 * F0 * KS * ddet, nm[i])) for i in range(2))
        lin += dot(ge, Lam); lth += dot(ye, Th)
    return DD * (psi + k2 - k3), DD * (lin + lth)
def pev(P, d):
    return sum(fr(c) * Fr(1).__mul__(1) * __import__('math').prod(Fr(x) ** e for x, e in zip(d, m)) for m, c in LL.pdict(P).items())
def poly_side(d):
    kp = pev(X.KD.m, d); lp = Fr(0); w = [Fr(1)] + list(d); ix = lambda I, b: I * qq + b
    for k in range(NV):
        La, Lb, Ta, Tb = X.fields[k]; Gk = G[k]
        for Fd, I2 in ((La, 4), (Lb, 5), (Ta, 0), (Tb, 1)):
            for I in range(6):
                if Fd.W[I] is not None:
                    lp += pev(Fd.W[I].m, d) * sum(w[b] * Gk[ix(I, b)][ix(I2, q)] for b in range(q))
                if Fd.V[I] is not None:
                    lp += pev(Fd.V[I].m, d) * Gk[ix(I, 0)][ix(I2, q)]
    return kp, lp
import random; rnd = random.Random(2); WORST = [0.0, 0.0]
for size in (0, 1e-3, 2e-2):
    d = [Fr(round(rnd.uniform(-1, 1) * size * 2 ** 40), 2 ** 40) for _ in range(ND)]
    kf, lf = run(d); kp, lp = poly_side(d)
    rk = float(abs(kf - kp) / abs(X.scale)); rl = float(abs(lf - lp) / abs(lf)); WORST[0] = max(WORST[0], rk); WORST[1] = max(WORST[1], rl)
    print('|d| ~ %.0e: K D pointwise %.12e polynomial %.12e rel diff %.2e | linear part pointwise %.12e polynomial %.12e rel diff %.2e'
          % (size, float(kf), float(kp), float(abs(kf - kp) / abs(X.scale)), float(lf), float(lp), float(abs(lf - lp) / abs(lf))), flush=True)
ok1 = WORST[0] < 1e-14 and WORST[1] < 1e-9
print('1. polynomials equal the pointwise formulas (exact arithmetic, F0 a point): %s (K D: %.1e of the scale, linear part: %.1e relative)' % ('ok' if ok1 else 'FAIL', WORST[0], WORST[1]))
# 2. Cauchy-Schwarz with the dual norms of Local, for the test function
from flint import arb
A = LL.A; ix = lambda I, b: I * qq + b
ge2 = sum((g[2] * A(X.tb[k]['aa'][q][q]) - 2 * g[1] * A(X.tb[k]['ab'][q][q]) + g[0] * A(X.tb[k]['bb'][q][q])) / g[3] for k, g in enumerate(X.gc))      # |grad e|^2 on P(c)
e2 = sum(g[3] * A(X.tb[k]['m'][q][q]) for k, g in enumerate(X.gc))
worst = 0.0; n = 0
for a in X.NLs:
    val = Fr(0)
    for k in range(NV):
        Gk = G[k]
        for Fd, I2 in zip(X.fields[k], (4, 5, 0, 1)):
            for I in range(6):
                if Fd.W[I] is not None:
                    Wd = LL.pdict(Fd.W[I].m)
                    for b in range(q):
                        if b == 0:
                            c = Wd.get(a)
                        elif a[b - 1] >= 1:
                            c = Wd.get(tuple(a[i] - (1 if i == b - 1 else 0) for i in range(ND)))
                        else:
                            c = None
                        if c is not None:
                            val += fr(c) * Gk[ix(I, b)][ix(I2, q)]
                if Fd.V[I] is not None:
                    c = LL.pdict(Fd.V[I].m).get(a)
                    if c is not None:
                        val += fr(c) * Gk[ix(I, 0)][ix(I2, q)]
    bound = ge2.sqrt() * X.NLs[a] + e2.sqrt() * X.NTs[a] / 2; r = abs(float(val)) / float(LL.up(bound)); worst = max(worst, r); n += 1
ok2 = worst <= 1.0
print('2. Cauchy-Schwarz for the test function, %d coefficients of degree <= 3: largest |value| / bound = %.3f  %s' % (n, worst, 'ok' if ok2 else 'FAIL'))
print('PASS' if ok1 and ok2 else 'FAILED'); sys.exit(0 if ok1 and ok2 else 1)
