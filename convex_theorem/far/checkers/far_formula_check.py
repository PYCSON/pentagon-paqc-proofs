#!/usr/bin/env python3
"""far_formula_check.py  (Claude, 2026-11-18): exact (symbolic) check of every algebraic identity used in
far/FAR_PROOFS.md, written from the geometric definitions and not from the label code; at the end the production
formulas (cutchart.verts, to_std, STUB, the invariants) are compared with the ones derived here.
Sections refer to FAR_PROOFS.md.  Needs sympy; the last part needs python-flint (import of cutchart).
Exit code 0 and a final line 'ALL IDENTITIES OK' if everything holds."""
import os, sys, random
from fractions import Fraction
import sympy as sp

OK = []


def chk(name, expr):
    e = sp.expand(expr) if not isinstance(expr, bool) else (0 if expr else 1)
    assert e == 0, (name, e)
    OK.append(name)


def cr(u, v):
    return u[0] * v[1] - u[1] * v[0]


def sub(u, v):
    return (u[0] - v[0], u[1] - v[1])


def area2(V):                                              # twice the signed area
    return sum(cr(V[i], V[(i + 1) % len(V)]) for i in range(len(V)))


def abc(V, k):                                             # A_k, B_k, C_k from the vertices (section 1)
    n = len(V); e = [sub(V[(i + 1) % n], V[i]) for i in range(n)]
    return cr(e[k - 1], e[k]), cr(e[k], e[(k + 1) % n]), cr(e[k - 1], e[(k + 1) % n])


# ---------------------------------------------------------------- section 1: the ear of an edge (generic polygon)
ex = sp.symbols('a0:6')                                    # three consecutive edge vectors e_{k-1}, e_k, e_{k+1}
em, ek, ep = (ex[0], ex[1]), (ex[2], ex[3]), (ex[4], ex[5])
A, B, C = cr(em, ek), cr(ek, ep), cr(em, ep)
s, r = B / C, A / C                                        # W = V_k + s e_{k-1} = V_{k+1} - r e_{k+1}
chk('1: W lies on both lines', sp.simplify(s * em[0] + r * ep[0] - ek[0]))
chk('1: W lies on both lines (y)', sp.simplify(s * em[1] + r * ep[1] - ek[1]))
chk('1: ear area', sp.simplify(sp.Rational(1, 2) * cr((s * em[0], s * em[1]), ek) - A * B / (2 * C)))

# ---------------------------------------------------------------- section 2: quadrilateral, invariants, vertex maps
q2x, q2y, q3x, q3y, t, u = sp.symbols('q2x q2y q3x q3y t u')
q0, q1, q2, q3 = (0, 0), (1, 0), (q2x, q2y), (q3x, q3y)
e0, e1, e2, e3 = sub(q1, q0), sub(q2, q1), sub(q3, q2), sub(q0, q3)
D = q2x * q3y - q2y * q3x; K = D - q3y + q2y; M = D - q3y; L = q3y - q2y
chk('2: q3y = e3 x e0', cr(e3, e0) - q3y); chk('2: q2y = e0 x e1', cr(e0, e1) - q2y)
chk('2: K = e1 x e2', cr(e1, e2) - K); chk('2: D = e2 x e3', cr(e2, e3) - D)
chk('2: M = e3 x e1', cr(e3, e1) - M); chk('2: L = e0 x e2', cr(e0, e2) - L)
chk('2: 2|Q| = q2y + D', area2([q0, q1, q2, q3]) - (q2y + D))


def verts(j):
    if j == 3:
        return [q0, q1, q2, (q3x + t * (q2x - q3x), q3y + t * (q2y - q3y)), ((1 - u) * q3x, (1 - u) * q3y)]
    return [(t, 0), q1, q2, q3, (u * q3x, u * q3y)]


CUT = {3: 3, 0: 4}
tau = {3: D / (q2y + D), 0: q3y / (q2y + D)}
earcut2 = {3: t * u * D, 0: t * u * q3y}                    # twice the ear of the cut edge
for j in (3, 0):
    V = verts(j)
    chk('2: |P| = |Q| (1 - t u tau), j=%d' % j, sp.simplify(area2(V) - (q2y + D) * (1 - t * u * tau[j])))
    a_, b_, c_ = abc(V, CUT[j])
    chk('2: ear of the cut edge, j=%d' % j, sp.simplify(a_ * b_ / c_ - earcut2[j]))
    chk('2: C of the cut edge, j=%d' % j, c_ - (1 - t) * (1 - u) * (D if j == 3 else q3y))

# ---------------------------------------------------------------- section 5 (X): the eight ear comparisons
# for every other edge k':  A B, C in closed form, and   C * (2 ear_cut) - A B  =  (positive factor) * (-g)
# so that, where C > 0 and the factor > 0:   ear_k' < ear_cut  <=>  g < 0.
X = {
    (3, 0): ((1 - u) * M, (1 - u), q2y * q3y - t * u * D * M),                 # untouched: t u > q2y q3y / (D M), M > 0
    (3, 1): ((1 - t) * L, (1 - t), q2y * K - t * u * L * D),                   # untouched: t u > q2y K / (L D), L > 0
    (3, 2): (t * K - u * M, u * D, K * (1 - 2 * t) + t * u * M),               # stub
    (3, 4): (u * q3y + t * (q2y - q3y), t * D, q3y * (1 - 2 * u) - t * u * (q2y - q3y)),      # stub
    (0, 1): ((1 - t) * L, (1 - t), q2y * K - t * u * q3y * L),                 # untouched: t u > q2y K / (q3y L), L > 0
    (0, 2): (-(1 - u) * M, (1 - u), D * K + t * u * q3y * M),                  # untouched: t u > D K / (q3y (-M)), M < 0
    (0, 0): (u * M + t * q2y, u * q3y, -(t * u * M + q2y * (2 * t - 1))),      # stub
    (0, 3): (u * D + t * (q2y - q3y), t * q3y, -(D * (2 * u - 1) + t * u * (q2y - q3y))),     # stub
}
for (j, k), (Cf, fac, g) in X.items():
    V = verts(j); a_, b_, c_ = abc(V, k)
    chk('5X: C_k closed form, j=%d k=%d' % (j, k), c_ - Cf)
    chk('5X: comparison, j=%d k=%d' % (j, k), c_ * earcut2[j] - a_ * b_ + fac * g)

# ---------------------------------------------------------------- section 3: charts and the corner principle
q3x_, q3y_, a, h2, q2x_, q2y_, b, h3 = sp.symbols('q3x q3y a h2 q2x q2y b h3')
CH = {
    'std': ((q2x, q2y, q3x, q3y), {}),
    'E2': ((q3x, q3y, a, h2), {q2x: 1 + a * (q3x - 1), q2y: a * q3y + h2}),
    'E3': ((q2x, q2y, b, h3), {q3x: b * q2x, q3y: b * q2y + h3}),
}
chk('3: E2, K = h (1 - q3x)', K.subs(CH['E2'][1]) - h2 * (1 - q3x))
chk('3: E3, D = h q2x', D.subs(CH['E3'][1]) - h3 * q2x)
chk('3: E2, q2 = q1 + a (q3 - q1) + h (0,1)', (q2x.subs(CH['E2'][1]) - (1 + a * (q3x - 1))) ** 2 + (q2y.subs(CH['E2'][1]) - (a * q3y + h2)) ** 2)
for name, (w, rep) in CH.items():
    for j in (3, 0):
        V = [(sp.expand(sp.sympify(x).subs(rep)), sp.expand(sp.sympify(y).subs(rep))) for x, y in verts(j)]
        z = list(w) + [t, u]
        for m in z:
            dv = [(sp.diff(x, m), sp.diff(y, m)) for x, y in V]
            for x, y in V:                                  # affine in each coordinate
                chk('3: %s j=%d affine in %s' % (name, j, m), sp.diff(x, m, 2)); chk('3: %s j=%d affine in %s (y)' % (name, j, m), sp.diff(y, m, 2))
            for i in range(5):                              # one common direction per coordinate
                for k in range(i + 1, 5):
                    chk('3: %s j=%d common direction of %s' % (name, j, m), cr(dv[i], dv[k]))
        for (jj, k), (Cf, fac, g) in X.items():             # the stub polynomials are multi-affine in the chart variables
            if jj == j and k in ((2, 4) if j == 3 else (0, 3)):
                for f in (Cf, g):
                    p = sp.Poly(sp.expand(f.subs(rep)), *z)
                    chk('5X: stub polynomial multi-affine, %s j=%d k=%d' % (name, j, k), max(max(mon) for mon in p.monoms()) <= 1)
        for f, nm in ((D, 'D'), (K, 'K'), (M, 'M'), (L, 'L'), (q2y, 'q2y'), (q3y, 'q3y'), (q2y + D, 'A2')):
            p = sp.Poly(sp.expand(f.subs(rep)), *w)
            chk('3: invariant %s multi-affine in %s' % (nm, name), max(max(mon) for mon in p.monoms()) <= 1)

# ---------------------------------------------------------------- section 4: hand-over
aa = (1 - q2x) / (1 - q3x); hh = K / (1 - q3x)
chk('4: E2 coordinates of q', sp.simplify(q2y - (aa * q3y + hh)))
bb = q3x / q2x; hb = D / q2x
chk('4: E3 coordinates of q', sp.simplify(q3y - (bb * q2y + hb)))

# ---------------------------------------------------------------- section 5 (X iii): homogeneous form at the E3 face
rr, x = sp.symbols('r x')
E3h = {q2x: x, q3x: b * x, q3y: b * q2y + rr * q2y}         # h = r q2y
chk('5X3: K', sp.expand(K.subs(E3h)) - sp.expand(q2y * (rr * (x - 1) + (1 - b)))); chk('5X3: L', sp.expand(L.subs(E3h) - q2y * (rr - (1 - b))))
chk('5X3: M', sp.expand(M.subs(E3h) - q2y * (rr * (x - 1) - b))); chk('5X3: D', sp.expand(D.subs(E3h) - q2y * x * rr))
Kr, Lr, Mr, q3r, Dr = rr * (x - 1) + (1 - b), rr - (1 - b), rr * (x - 1) - b, b + rr, x * rr
tu = sp.Symbol('tu')
F3 = {(0, 1): tu * q3r * Lr - Kr, (0, 2): tu * q3r * (-Mr) - Dr * Kr, (3, 0): tu * Dr * Mr - q3r, (3, 1): tu * Lr * Dr - Kr}
for (j, k), f in F3.items():
    g = X[(j, k)][2].subs(E3h).subs(t * u, tu)
    chk('5X3: homogeneous comparison j=%d k=%d' % (j, k), sp.expand(g.subs({t: tu, u: 1}) + q2y ** 2 * f) if (t in g.free_symbols or u in g.free_symbols) else sp.expand(g + q2y ** 2 * f))

# ---------------------------------------------------------------- section 5 (R): |T| = |P| + ear_m + ear_{m+2}, exact, random pentagons
random.seed(1); nR = 0
while nR < 500:
    ang = sorted(random.sample(range(1, 3600), 5)); rad = [Fraction(random.randint(60, 100), 100) for _ in range(5)]
    import math
    V = [(Fraction(round(1000 * rad[i] * math.cos(ang[i] * math.pi / 1800)), 1000), Fraction(round(1000 * rad[i] * math.sin(ang[i] * math.pi / 1800)), 1000)) for i in range(5)]
    if any(abc(V, k)[0] <= 0 for k in range(5)):
        continue
    for m in range(5):
        m2 = (m + 2) % 5
        Am, Bm, Cm = abc(V, m); A2_, B2_, C2_ = abc(V, m2)
        if Cm <= 0 or C2_ <= 0:
            continue
        lines = [(V[i], sub(V[(i + 1) % 5], V[i])) for i in ((m + 1) % 5, (m + 3) % 5, (m + 4) % 5)]

        def meet(l1, l2):
            (p, d), (p2, d2) = l1, l2; s_ = cr(sub(p2, p), d2) / cr(d, d2); return (p[0] + s_ * d[0], p[1] + s_ * d[1])
        T = [meet(lines[0], lines[1]), meet(lines[1], lines[2]), meet(lines[2], lines[0])]
        assert abs(area2(T)) == area2(V) + Am * Bm / Cm + A2_ * B2_ / C2_, (V, m)
        for (p, d) in lines:                                # P lies in T: its vertices are on the inner side of the three lines
            assert all(cr(d, sub(v, p)) >= 0 for v in V)
        nR += 1
OK.append('5R: |T| = |P| + ear_m + ear_{m+2} and P inside T on %d exact random cases' % nR)

# ---------------------------------------------------------------- production code against the formulas above
if '--no-code' not in sys.argv:
    HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..'))
    try:
        import cutchart as CC
    except ImportError as ex:
        print(len(OK), 'checks'); print('IDENTITIES OK; the production formulas were NOT compared (%s); run from the repository tree or with --no-code' % ex); sys.exit(0)
    z = [q2x, q2y, q3x, q3y, t, u]
    for j in (3, 0):
        for (x1, y1), (x2, y2) in zip(CC.verts(z, j), verts(j)):
            chk('code: verts j=%d' % j, (x1 - x2)); chk('code: verts j=%d (y)' % j, (y1 - y2))
        for k, Cc, gc in CC.STUB[j]:
            chk('code: STUB C j=%d k=%d' % (j, k), Cc(z) - X[(j, k)][0]); chk('code: STUB g j=%d k=%d' % (j, k), gc(z) - X[(j, k)][2])
        chk('code: tau j=%d' % j, sp.simplify(CC.tau(z, j) - tau[j]))
    chk('code: D', CC._D(z) - D); chk('code: M', CC._M(z) - M); chk('code: K', CC._K(z) - K); chk('code: quad_area', 2 * CC.quad_area(z) - (q2y + D))
    for name, (w, rep) in CH.items():
        got = CC.to_std(list(w) + [t, u], name); want = [sp.sympify(v).subs(rep) for v in (q2x, q2y, q3x, q3y)]
        for g_, w_ in zip(got[:4], want):
            chk('code: to_std %s' % name, g_ - w_)
print(len(OK), 'checks'); print('ALL IDENTITIES OK')
