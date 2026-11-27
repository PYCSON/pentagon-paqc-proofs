"""Checks of lgbox.py / lgpoly.py (LG box certificate in exact polynomial form).  python3 tests/lgpoly_check.py
 1. the polynomial identity, by another route: at random rational points z = mid + d the trial pair u(z), chi(z) is
    rebuilt on the polygon P(z) as polynomials in (x, y) on every triangle, A0, A1, A2 are integrated by Green's theorem
    (lgcert.integrate), and  PC(d) = C D,  PD(d) = Dc D  must hold EXACTLY (rational equality), for orders 0 and 1.
 2. poly_positive: never accepts a polynomial that is <= 0 at a sample point of the box (random polynomials), and
    accepts clearly positive ones.
 3. certified boxes of the near chart: at random points of the box the pointwise bound |P| Dc / (KS A2 - |P| A1) is
    above F*, and below a Rayleigh-Ritz value of the same polygon from another trial space (upper bound of F).
 4. a box that contains the regular pentagon is NOT certified; nor is a box in which a fan triangle degenerates.
 5. near_lg: the centre point equals plan_cells.C; the three skip labels on hand-made boxes.
Exit code 0 and the last line PASS if all checks pass."""
import sys, os, math, random
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..')); sys.path.insert(0, HERE)
import numpy as np
import lgfan as LF, lgcert as LG, lgbox as LB, lgpoly as LP
from flint import fmpq
bad = 0; FSTAR = 18.919134195733914


def check(name, cond, info=''):
    global bad
    bad += not cond; print('%-74s %s %s' % (name, 'ok' if cond else 'FAIL', info), flush=True)


def fq(v):
    v = Fraction(v); return fmpq(v.numerator, v.denominator)


X, Y, CTX = LG.X, LG.Y, LG.CTX


def exact_A(V, O, N, n, cu, cx):
    """exact A0, A1, A2 (fmpq) of u = sum cu B (degree N, zero on the boundary), w = -(p - O)(2 + E)^-1 u + curl chi,
    chi = sum cx B (degree n), on the fan (O; V); cu, cx: rational coefficient lists in the global numbering"""
    nv = len(V); G, nd = LF.dofs(nv, N, True); Gc, nc = LF.dofs(nv, n, False); A = [fmpq(0)] * 3; ox, oy = fq(O[0]), fq(O[1])
    one = CTX.from_dict({(0, 0): fmpq(1)})

    def bern(coef, deg, lam, a, b):
        f = CTX.from_dict({})
        for m, (i, j, l) in enumerate(LF._loc(deg)):
            if coef[m]:
                f = f + fq(Fraction(coef[m]) * (math.factorial(deg) // (math.factorial(i) * math.factorial(j) * math.factorial(l)))) * (lam ** i if i else one) * (a ** j if j else one) * (b ** l if l else one)
        return f
    for k in range(nv):
        v1, v2 = V[k], V[(k + 1) % nv]; e1 = (v1[0] - O[0], v1[1] - O[1]); e2 = (v2[0] - O[0], v2[1] - O[1]); det = e1[0] * e2[1] - e1[1] * e2[0]
        a = fq(e2[1] / det) * (X - ox) + fq(-e2[0] / det) * (Y - oy); b = fq(-e1[1] / det) * (X - ox) + fq(e1[0] / det) * (Y - oy); lam = 1 - a - b
        u = bern([cu[g] if g >= 0 else 0 for g in G[k]], N, lam, a, b); chi = bern([cx[g] for g in Gc[k]], n, lam, a, b)
        us = u.compose(X + ox, Y + oy, ctx=CTX); h = CTX.from_dict({e: cf / (e[0] + e[1] + 2) for e, cf in us.to_dict().items()}).compose(X - ox, Y - oy, ctx=CTX)
        w1 = -(X - ox) * h + chi.derivative(1); w2 = -(Y - oy) * h - chi.derivative(0); tri = [O, v1, v2]; ux = u.derivative(0); uy = u.derivative(1)
        A = [A[0] + LG.integrate(ux * ux + uy * uy, tri), A[1] + LG.integrate(u * u, tri), A[2] + LG.integrate(w1 * w1 + w2 * w2, tri)]
    return A


def verts_q(z):
    return [(Fraction(0), Fraction(0)), (Fraction(1), Fraction(0)), (z[0], z[1]), (z[2], z[3]), (z[4], z[5])]


c72, s72 = math.cos(2 * math.pi / 5), math.sin(2 * math.pi / 5); XREG = [1 + c72, s72, 0.5, s72 + math.sin(math.pi / 5), -c72, s72]
rng = random.Random(11); nprng = np.random.default_rng(11)
# 1
for order, N in ((0, 3), (1, 3), (1, 4)):
    dv = nprng.normal(size=6); dv /= np.linalg.norm(dv); mid = [round((XREG[i] + 0.2 * dv[i]) * 64) / 64 for i in range(6)]
    fam = LB.Family(LB.near_verts, mid, N, order=order); PC, PD, dets = LP.build(fam, LP.near_aff, mid, FSTAR); ok = True
    for trial in range(2):
        d = [Fraction(rng.randint(-8, 8), 256) for _ in range(6)]; z = [Fraction(mid[i]) + d[i] for i in range(6)]
        Dv = [Fraction(1)] + (d if order else []); su = Fraction(1, 2 ** LF.CB); sx = Fraction(2) ** fam.xe
        cu = [su * sum(Dv[a] * fam.CI[a][g] for a in range(len(Dv))) for g in range(len(fam.CI[0]))]
        cx = [sx * sum(Dv[a] * fam.XI[a][g] for a in range(len(Dv))) for g in range(len(fam.XI[0]))]
        V = verts_q(z); O = (Fraction(fam.O[0]), Fraction(fam.O[1])); A0, A1, A2 = exact_A(V, O, fam.N, fam.n, cu, cx)
        dk = [(V[k][0] - O[0]) * (V[(k + 1) % 5][1] - O[1]) - (V[k][1] - O[1]) * (V[(k + 1) % 5][0] - O[0]) for k in range(5)]; area = sum(dk) / 2; D = Fraction(1)
        for t in dk:
            D *= t
        KS = Fraction(int(LF.KS_LO.p), int(LF.KS_LO.q)); fr = lambda t: Fraction(int(t.p), int(t.q)); A0, A1, A2 = fr(A0), fr(A1), fr(A2)
        Dc = KS * A1 - area * A0; Cc = area * Dc - Fraction(FSTAR) * (KS * A2 - area * A1)
        ev = lambda P: fr(P(*[fq(t) for t in d]))
        ok = ok and ev(PD) == Dc * D and ev(PC) == Cc * D and all(ev(dets[k]) == dk[k] for k in range(5))
    check('polynomial identity, order %d, degree %d (2 rational points, exact)' % (order, N), ok, 'terms %d, degree %d' % (len(PC.to_dict()), PC.total_degree()))
# 2
C2, DV2 = LP.ctx(2); wrong = 0; acc = 0; tot = 0
for trial in range(400):
    cf = {(i, j): fq(Fraction(rng.randint(-20, 20), 8)) for i in range(4) for j in range(4 - i)}; cf[(0, 0)] = fq(Fraction(rng.randint(1, 60), 8))
    P = C2.from_dict(cf); rho = [2.0 ** -rng.randint(0, 3), 2.0 ** -rng.randint(0, 3)]; ok, info = LP.poly_positive(P, rho); tot += 1; acc += ok
    if ok:
        for s in range(200):
            pt = [fq(Fraction(rng.randint(-64, 64), 64) * Fraction(rho[0])), fq(Fraction(rng.randint(-64, 64), 64) * Fraction(rho[1]))]
            if not (P(*pt) > 0):
                wrong += 1; break
check('poly_positive on 400 random cubics in 2 variables: no false acceptance', wrong == 0 and acc > 20, 'accepted %d of %d' % (acc, tot))
P = C2.from_dict({(0, 0): fmpq(1), (2, 0): fmpq(1), (0, 2): fmpq(1), (1, 1): fmpq(-1), (1, 0): fmpq(-1)})
check('poly_positive accepts 1 - x + x^2 - x y + y^2 on [-1, 1]^2 (minimum 2/3; convexity needed)', LP.poly_positive(P, [1.0, 1.0])[0])
check('poly_positive refuses x^2 + y^2 - x y - x + 1/3 on [-1, 1]^2 (minimum 0)', not LP.poly_positive(P - fmpq(2, 3), [1.0, 1.0])[0])
# 3
sys.path.insert(0, HERE)
import lgfan_check_ritz as RZ
low = 9.0; high = 0.0; n = 0; viol = 0
for r in (0.06, 0.15, 0.35):
    dv = nprng.normal(size=6); dv /= np.linalg.norm(dv); mid = [round((XREG[i] + r * dv[i]) * 4096) / 4096 for i in range(6)]
    fam = LB.Family(LB.near_verts, mid, 8, order=1); rho = 2.0 ** math.floor(math.log2(0.15 * r)); ok, info = LP.check_box(fam, LP.near_aff, mid, [rho] * 6, FSTAR)
    if not ok:
        rho /= 2; ok, info = LP.check_box(fam, LP.near_aff, mid, [rho] * 6, FSTAR)
    check('near box at distance %.2f, half-width %.5f certified' % (r, rho), ok, str({k: info[k] for k in ('stage', 'c0', 'qp', 'tail') if k in info}))
    for s in range(6):
        z = [mid[i] + rho * (rng.choice((-1, 1)) if s < 3 else rng.uniform(-1, 1)) for i in range(6)]
        b = LB.bound_at(fam, LB.near_verts, z); V = [(0.0, 0.0), (1.0, 0.0), (z[0], z[1]), (z[2], z[3]), (z[4], z[5])]; lam, A = RZ.ritz(V); F = lam * A
        if b is None:
            viol += 1; continue
        bl = float(b.lower()); n += 1; viol += not (FSTAR < bl <= F * (1 + 1e-9)); low = min(low, bl / F); high = max(high, bl / F)
check('pointwise bound at 18 points of the certified boxes: F* < bound <= Ritz value', viol == 0 and n == 18, 'bound / Ritz in [%.6f, %.6f]' % (low, high))
# 4
mid = [round(v * 1024) / 1024 for v in XREG]; fam = LB.Family(LB.near_verts, mid, 8, order=1)
check('box around the regular pentagon (half-width 1/64): not certified', not LP.check_box(fam, LP.near_aff, mid, [1 / 64] * 6, FSTAR)[0])
ok, info = LP.check_box(fam, LP.near_aff, mid, [2.0] * 6, FSTAR)
check('box of half-width 2 (fan triangles degenerate): refused', not ok, info['stage'])
dk = LP.build(fam, LP.near_aff, mid, FSTAR)[2]
check('  ... and a fan determinant is not proved positive on it', not all(LP.poly_positive(t, [2.0] * 6)[0] for t in dk))
# 5
import near_lg as NL
try:
    sys.path.insert(0, os.path.join(HERE, '..', '..', 'v2_cover'))
    import plan_cells as PC
    check('near_lg.C equals plan_cells.C (binary64, exact)', [float(v) for v in PC.C] == NL.CF)
except ImportError:
    print('plan_cells not available: C not compared')
check('skip tests: a box inside the ball of radius 1/20 is inner, one beyond 3/5 is far, one with x2 <= 1/2 is mirror',
      NL.skip([c - 0.01 for c in NL.CF], [c + (0.0 if i == 2 else 0.01) for i, c in enumerate(NL.CF)]) == 'mirror' and NL.skip([c + (0.001 if i == 2 else -0.01) for i, c in enumerate(NL.CF)], [c + 0.01 for c in NL.CF]) == 'inner'
      and NL.skip([c + (0.5 if i == 2 else 0.3) for i, c in enumerate(NL.CF)], [c + 0.6 for c in NL.CF]) == 'far' and NL.skip([c + (0.001 if i == 2 else -0.03) for i, c in enumerate(NL.CF)], [c + 0.03 for c in NL.CF]) is None)
print('PASS' if bad == 0 else 'FAILED %d' % bad); sys.exit(1 if bad else 0)
