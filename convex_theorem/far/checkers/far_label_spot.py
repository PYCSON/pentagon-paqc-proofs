#!/usr/bin/env python3
"""far_label_spot.py RUN_OUT [--sample K] [--seed S] [--chart std|E2|E3] [--part i/n]   (Claude, 2026-11-19)

A second implementation of the closed-form far labels, written from FAR_THEOREM.md / FAR_PROOFS.md only (it imports
no production module and re-derives the eight ear comparisons; see the table in FAR_PROOFS.md section 5).  Exact
rational arithmetic throughout.  It re-proves, for a random sample of leaves (or all with --sample 0):
  V          an invariant q2y, q3y, K, D <= 0 at the 16 corners of the w-box
  E2 / E3    the hand-over conditions of section 4 (std chart)
  T          sup t * sup u < theta_T = (1 - F*/LB)/tau_bar, with LB = max(2 pi^2, the two strip bounds, the triangle
             bound for the four edges of Q, Makai's bound, the root's table value from roots.json, or
             the smallest value of the table parts (qtable/cells.jsonl) that meet the w-box; if that is not enough, the
             fibre bound of FIBRE_LEMMA.md for the hull of the corner quadrilaterals, own implementation with python-flint)  [the table itself is checked by far_verify --table and by the second party]
  X          one of the eight rows: (i) untouched edges via inf t inf u > g_bar, (ii) stub edges via the signs of two
             multi-affine polynomials at the 64 corners, (iii) the E3 corner face: homogeneous forms in r = h/q2y,
             adaptive bisection of r in (0, 1] and s = 1/r in [0, 1] with exact interval arithmetic
  TX         theta_T > g_bar for an untouched edge whose sign condition holds on the w-box
Labels S, R, O, N and the certified boxes C are not re-checked here (they are counted).  A leaf that is not re-proved
is reported as 'not reproduced' (this checker may be weaker than the producer's criterion; it is not a failure of the
proof by itself).  Data only: the run directory is read, nothing in it is executed.  Run with python3 -I."""
import sys, os, json, random, itertools, collections, glob
from fractions import Fraction as Fr

PI2_LO = Fr(9869604401089358, 10 ** 15)         # < pi^2 = 9.86960440108935861883...
FSTAR = Fr(18.919134195733914)
CF = [float.fromhex(v) for v in ('0x1.4f1bbcdcbfa53p+0', '0x1.e6f0e134454ffp-1', '0x1.0000000000000p-1', '0x1.89f188bdcd7aep+0',
                                  '-0x1.3c6ef372fe94cp-2', '0x1.e6f0e134454ffp-1')]   # C, the binary64 chart point (FAR_THEOREM 0)
TRI_LO = 4 * PI2_LO / Fr(17320508075688773, 10 ** 16)   # < 4 pi^2/sqrt 3 (equilateral triangle); sqrt 3 = 1.73205080756887729...
H2 = H3 = Fr(1, 8)


def q_of(chart, w):
    """the two free vertices q2, q3 of Q = (0,0), (1,0), q2, q3 from the chart coordinates w (section 3)"""
    if chart == 'std':
        return (w[0], w[1]), (w[2], w[3])
    if chart == 'E2':                          # w = (q3x, q3y, a, h): q2 = q1 + a (q3 - q1) + h (0, 1)
        q3x, q3y, a, h = w
        return (1 + a * (q3x - 1), a * q3y + h), (q3x, q3y)
    q2x, q2y, b, h = w                         # E3: w = (q2x, q2y, b, h): q3 = b q2 + h (0, 1)
    return (q2x, q2y), (b * q2x, b * q2y + h)


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def invariants(chart, w):
    q2, q3 = q_of(chart, w)
    e1 = (q2[0] - 1, q2[1]); e2 = (q3[0] - q2[0], q3[1] - q2[1]); e3 = (-q3[0], -q3[1])
    K = cross(e1, e2); D = cross(e2, e3); q2y = q2[1]; q3y = q3[1]
    return dict(q2y=q2y, q3y=q3y, K=K, D=D, M=D - q3y, L=q3y - q2y, q3x=q3[0], q2x=q2[0],
                W=max(q2y, q3y), A2=q2y + D)                      # A2 = 2|Q|


def wcorners(lo, hi):
    return [tuple(hi[i] if (m >> i) & 1 else lo[i] for i in range(4)) for m in range(16)]


def per_lower(chart, lo, hi):
    """lower bound of the perimeter of Q over the w-box: each edge vector is multi-affine in w, so each component
    ranges between its corner values; |e| >= sqrt(inf|ex|^2 + inf|ey|^2) >= a rational lower bound"""
    cs = [q_of(chart, c) for c in wcorners(lo, hi)]
    edges = [lambda q2, q3: (Fr(1), Fr(0)), lambda q2, q3: (q2[0] - 1, q2[1]),
             lambda q2, q3: (q3[0] - q2[0], q3[1] - q2[1]), lambda q2, q3: (-q3[0], -q3[1])]
    tot = Fr(0)
    for e in edges:
        v = [e(*c) for c in cs]; sq = Fr(0)
        for k in (0, 1):
            a, b = min(x[k] for x in v), max(x[k] for x in v)
            m = Fr(0) if a <= 0 <= b else min(abs(a), abs(b)); sq += m * m
        tot += isqrt_lo(sq)
    return tot


def isqrt_lo(x):
    """a rational lower bound of sqrt(x), x >= 0"""
    if x <= 0:
        return Fr(0)
    import math
    r = Fr(math.sqrt(float(x)) * (1 - 1e-12))
    while r * r > x:
        r *= Fr(999999, 1000000)
    return r


def hull(pts):
    """convex hull, counter-clockwise, exact (Andrew's monotone chain)"""
    P = sorted(set(pts))
    if len(P) <= 2:
        return P
    def half(seq):
        h = []
        for p in seq:
            while len(h) >= 2 and cross((h[-1][0] - h[-2][0], h[-1][1] - h[-2][1]), (p[0] - h[-2][0], p[1] - h[-2][1])) <= 0:
                h.pop()
            h.append(p)
        return h
    lo, up = half(P), half(P[::-1])
    return lo[:-1] + up[:-1]


def fibre_lower(H, d, nseg=48):
    """lower bound of lambda_1 of the convex polygon H (ccw) by FIBRE_LEMMA.md: rotate-and-scale by d, step minorant of
    pi^2/w^2 on nseg segments, positive solution in ball arithmetic.  Returns a Fraction (0 if it fails)."""
    import flint, numpy as np
    D2 = d[0] * d[0] + d[1] * d[1]
    X = [p[0] * d[0] + p[1] * d[1] for p in H]; Y = [d[0] * p[1] - d[1] * p[0] for p in H]
    n = len(H); up = []; low = []
    for i in range(n):
        k = (i + 1) % n
        if X[k] < X[i]:
            up.append((X[i], Y[i], X[k], Y[k]))
        elif X[k] > X[i]:
            low.append((X[i], Y[i], X[k], Y[k]))
    if not up or not low:
        return Fr(0)
    a, b = min(X), max(X); h = (b - a) / nseg
    def line(e, x):
        x0, y0, x1, y1 = e
        return y0 + (y1 - y0) * (x - x0) / (x1 - x0)
    c = []
    for kk in range(nseg):
        s0, s1 = a + kk * h, a + (kk + 1) * h
        W = min(max(line(U, s0) - line(L, s0), line(U, s1) - line(L, s1)) for U in up for L in low)
        if W <= 0:
            return Fr(0)
        c.append(PI2_LO / (W * W))
    # floating estimate of mu: finite differences
    M = 40; N = nseg * M; hh = float(h) / M
    V = np.repeat([float(x) for x in c], M)[1:]
    main = 2 / hh ** 2 + V; off = -np.ones(N - 2) / hh ** 2
    import scipy.linalg
    mu = scipy.linalg.eigh_tridiagonal(main, off, select='i', select_range=(0, 0), eigvals_only=True)[0]
    flint.ctx.prec = 96
    for f in (0.9985, 0.994, 0.98, 0.95):
        m = Fr(mu * f).limit_denominator(10 ** 12)
        mm = flint.arb(m.numerator) / m.denominator; ha = flint.arb(h.numerator) / h.denominator
        phi, dphi, ok = flint.arb(0), flint.arb(1), True
        for ck in c:
            ca = flint.arb(ck.numerator) / ck.denominator; q = ca - mm
            if ck > m:
                kap = q.sqrt(); ch, sh = (kap * ha).cosh(), (kap * ha).sinh()
                phi, dphi = phi * ch + dphi * sh / kap, phi * kap * sh + dphi * ch
            elif ck < m:
                kap = (-q).sqrt()
                if not (kap * ha < flint.arb.pi()):
                    ok = False; break
                co, si = (kap * ha).cos(), (kap * ha).sin()
                phi, dphi = phi * co + dphi * si / kap, -phi * kap * si + dphi * co
            else:
                phi, dphi = phi + dphi * ha, dphi
            if not (phi > 0):
                ok = False; break
        if ok:
            return m * D2
    return Fr(0)


def quad_fibre(chart, lo, hi, need=None):
    """F(Q) >= lambda_1(hull of the corner quadrilaterals) * min corner area, best of a few directions"""
    cs = wcorners(lo, hi); pts = []
    for c in cs:
        q2, q3 = q_of(chart, c); pts += [(Fr(0), Fr(0)), (Fr(1), Fr(0)), q2, q3]
    H = hull(pts)
    if len(H) < 3:
        return Fr(0)
    amin = min(pos(invariants(chart, c)['A2']) for c in cs) / 2
    best = Fr(0)
    import math
    for step, nseg in ((10, 48), (2, 128)):   # any direction is valid; rational approximations of the angles
        for ang in range(0, 180, step):
            d = (Fr(math.cos(math.radians(ang))).limit_denominator(1000), Fr(math.sin(math.radians(ang))).limit_denominator(1000))
            best = max(best, fibre_lower(H, d, nseg) * amin)
            if need is not None and best > need:
                return best
        if need is None:
            break
    return best


def verts(chart, w, t, u, j):
    """the five vertices of P (section 2): j = 3: q0, q1, q2, q3 + t(q2 - q3), (1 - u) q3; j = 0: (t,0), q1, q2, q3, u q3"""
    q2, q3 = q_of(chart, w); O = (Fr(0), Fr(0)); q1 = (Fr(1), Fr(0))
    if j == 3:
        return [O, q1, q2, (q3[0] + t * (q2[0] - q3[0]), q3[1] + t * (q2[1] - q3[1])), ((1 - u) * q3[0], (1 - u) * q3[1])]
    return [(t, Fr(0)), q1, q2, q3, (u * q3[0], u * q3[1])]


def area2(V):
    return sum(cross(V[i], V[(i + 1) % len(V)]) for i in range(len(V)))


def zcorners(lo, hi):
    return [tuple(hi[i] if (m >> i) & 1 else lo[i] for i in range(6)) for m in range(64)]


def sqrt_up(x):
    import math
    if x <= 0:
        return Fr(0)
    r = Fr(math.sqrt(float(x)) * (1 + 1e-12))
    while r * r < x:
        r *= Fr(1000001, 1000000)
    return r


def pent_labels(chart, lo, hi, j, label):
    """S, R, O, N from FAR_PROOFS.md section 5, own implementation; True if the criterion is reproduced"""
    Z = zcorners(lo, hi); VV = [verts(chart, z[:4], z[4], z[5], j) for z in Z]
    amin = min(area2(V) for V in VV) / 2      # the area is multi-affine: its minimum is a corner value
    if amin <= 0 and label != 'S':
        return False
    E = [[(V[(i + 1) % 5][0] - V[i][0], V[(i + 1) % 5][1] - V[i][1]) for i in range(5)] for V in VV]
    if label == 'O':
        H = hull([p for V in VV for p in V])
        need = FSTAR / amin
        import math
        for step, nseg in ((10, 48), (2, 128), (1, 320)):
            for ang in range(0, 180, step):
                d = (Fr(math.cos(math.radians(ang))).limit_denominator(1000), Fr(math.sin(math.radians(ang))).limit_denominator(1000))
                if fibre_lower(H, d, nseg) > need:
                    return True
        return False
    if label == 'R':
        for m in range(5):
            tot = Fr(0); ok = True
            for k in (m, (m + 2) % 5):
                A = [cross(e[(k - 1) % 5], e[k]) for e in E]; B = [cross(e[k], e[(k + 1) % 5]) for e in E]; C = [cross(e[(k - 1) % 5], e[(k + 1) % 5]) for e in E]
                if min(C) <= 0:
                    ok = False; break
                tot += pos(max(A)) * pos(max(B)) / (2 * min(C))
            if ok and TRI_LO / (1 + tot / amin) > FSTAR:
                return True
        return False
    if label == 'S':
        for a in range(5):                     # S(b): edge a -> a+1
            b = (a + 1) % 5
            dx = [V[b][0] - V[a][0] for V in VV]; dy = [V[b][1] - V[a][1] for V in VV]
            l2 = sum((Fr(0) if min(v) <= 0 <= max(v) else min(abs(min(v)), abs(max(v)))) ** 2 for v in (dx, dy))
            cmax = max(cross((V[b][0] - V[a][0], V[b][1] - V[a][1]), (V[k][0] - V[a][0], V[k][1] - V[a][1])) for V in VV for k in range(5))
            if cmax > 0 and PI2_LO * l2 / (2 * cmax) > FSTAR:
                return True
        import math
        cen = verts(chart, [(lo[i] + hi[i]) / 2 for i in range(4)], (lo[4] + hi[4]) / 2, (lo[5] + hi[5]) / 2, j)
        dirs = [(Fr(math.cos(math.radians(a))).limit_denominator(1000), Fr(math.sin(math.radians(a))).limit_denominator(1000)) for a in range(0, 180, 5)]
        for a in range(5):                     # and the normals of the edges and diagonals at the centre
            for b in range(5):
                if a < b:
                    dirs.append((-(cen[b][1] - cen[a][1]), cen[b][0] - cen[a][0]))
        if amin <= 0:                          # S(a) needs a positive area at the corners
            return False
        for n in dirs:                         # S(a): direction n
            n2 = n[0] ** 2 + n[1] ** 2; tau = None
            for V in VV:
                p = [x * n[0] + y * n[1] for x, y in V]; W = max(p) - min(p)
                r = n2 * area2(V) / 2 / (W * W); tau = r if tau is None else min(tau, r)
            if PI2_LO * tau > FSTAR:
                return True
        return False
    if label == 'N':
        import cmath
        C6 = [Fr(x) for x in CF]
        Cp = [complex(float(C6[0]), float(C6[1])), complex(float(C6[2]), float(C6[3])), complex(float(C6[4]), float(C6[5]))]
        for rev in (False, True):
            for s0 in range(5):
                order = [(s0 + i) % 5 for i in range(5)] if not rev else [(s0 - i) % 5 for i in range(5)]
                def cpl(V, i):
                    x, y = V[order[i]]
                    return (x, -y) if rev else (x, y)                  # reversed order: reflect (conjugate)
                Vc = [Fr(sum(VV[k][i][0] for k in range(64)), 64) for i in range(5)]
                cen = verts(chart, [(lo[i] + hi[i]) / 2 for i in range(4)], (lo[4] + hi[4]) / 2, (lo[5] + hi[5]) / 2, j)
                def cc(V, i):
                    x, y = cpl(V, i); return complex(float(x), float(y))
                D0 = cc(cen, 1) - cc(cen, 0)
                if abs(D0) == 0:
                    continue
                r = [(cc(cen, p) - cc(cen, 0)) / D0 for p in (2, 3, 4)]
                if sum(abs(r[i] - Cp[i]) ** 2 for i in range(3)) > 0.3:
                    continue
                rq = [(Fr(z.real).limit_denominator(10 ** 9), Fr(z.imag).limit_denominator(10 ** 9)) for z in r]
                eps = (Fr(D0.real / abs(D0)).limit_denominator(10 ** 6), Fr(D0.imag / abs(D0)).limit_denominator(10 ** 6))
                e2 = eps[0] ** 2 + eps[1] ** 2
                delta = None; rho2 = [Fr(0)] * 3
                for V in VV:
                    W0 = cpl(V, 0); W1 = cpl(V, 1); Dl = (W1[0] - W0[0], W1[1] - W0[1])
                    re = (Dl[0] * eps[0] + Dl[1] * eps[1])                # Re(Delta conj(eps))
                    delta = re if delta is None else min(delta, re)
                    for i, pidx in enumerate((2, 3, 4)):
                        Wp = cpl(V, pidx); Np = (Wp[0] - W0[0], Wp[1] - W0[1])
                        rr = rq[i]; rD = (rr[0] * Dl[0] - rr[1] * Dl[1], rr[0] * Dl[1] + rr[1] * Dl[0])
                        rho2[i] = max(rho2[i], (Np[0] - rD[0]) ** 2 + (Np[1] - rD[1]) ** 2)
                # |Delta| >= Re(Delta conj eps)/|eps|; |eps| <= sqrt_up(e2)
                if delta is None or delta <= 0:
                    continue
                dl = delta / sqrt_up(e2)
                tot = Fr(0)
                for i in range(3):
                    off = sqrt_up((rq[i][0] - C6[2 * i]) ** 2 + (rq[i][1] - C6[2 * i + 1]) ** 2)
                    tot += (off + sqrt_up(rho2[i]) / dl) ** 2
                if tot < Fr(9, 25):
                    return True
        return False
    return None


def ranges(chart, lo, hi):
    vals = [invariants(chart, c) for c in wcorners(lo, hi)]
    return {k: (min(v[k] for v in vals), max(v[k] for v in vals)) for k in vals[0]}, vals


def pos(x):
    return x if x > 0 else Fr(0)


def tau_bar(R, j):
    if j == 3:                                 # tau3 = D/(q2y + D): increasing in D, decreasing in q2y
        d, q = pos(R['D'][1]), pos(R['q2y'][0])
        return Fr(1) if d + q == 0 else min(Fr(1), d / (q + d))
    den = R['A2'][0]                           # tau0 = q3y/(q2y + D)
    return Fr(1) if den <= 0 else min(Fr(1), pos(R['q3y'][1]) / den)


def lb_quad(R, table, per_lo=0):
    """F(Q) >= LB for every strictly convex Q of the w-box (section 6, the parts re-implemented here)"""
    out = [2 * PI2_LO]
    W = R['W'][1]
    if W > 0:
        A = pos(R['A2'][0]) / 2                # |Q| >= min corner area (multi-affine)
        out.append(PI2_LO * A / W ** 2)        # strip of width W = max(q2y, q3y)
        out.append(PI2_LO / (2 * W))           # |Q| >= W/2 (triangle over the base)
    area = pos(R['A2'][0]) / 2
    if area > 0:                               # triangle bound: Q is the triangle Q + ear_k minus the ear (Lemma 1)
        for A, B, C in (('q3y', 'q2y', 'M'), ('q2y', 'K', 'L'), ('K', 'D', '-M'), ('D', 'q3y', '-L')):
            cinf = -R[C[1:]][1] if C[0] == '-' else R[C][0]
            if cinf > 0:
                ear = pos(R[A][1]) * pos(R[B][1]) / (2 * cinf)
                out.append(TRI_LO / (1 + ear / area))
    if per_lo and pos(R['A2'][1]) > 0:        # Makai / Hersch: F >= pi^2 Per^2/(16 |Q|) for convex domains
        out.append(PI2_LO * per_lo ** 2 / (8 * R['A2'][1]))
    if table:
        out.append(Fr(table))
    return max(out)


def g_rows(R, j):
    """untouched edges: (sign condition holds on the w-box, g_bar) for the rows of FAR_PROOFS section 5"""
    s = lambda k: pos(R[k][1])                 # sup of a factor that is positive for members
    i = lambda k: R[k][0]                      # inf
    rows = []
    if j == 3:
        if i('M') > 0 and i('D') > 0:
            rows.append(s('q2y') * s('q3y') / (i('D') * i('M')))          # k' = 0
        if i('L') > 0 and i('D') > 0:
            rows.append(s('q2y') * s('K') / (i('L') * i('D')))            # k' = 1
    else:
        if i('L') > 0 and i('q3y') > 0:
            rows.append(s('q2y') * s('K') / (i('q3y') * i('L')))          # k' = 1
        if R['M'][1] < 0 and i('q3y') > 0:
            rows.append(s('D') * s('K') / (i('q3y') * (-R['M'][1])))      # k' = 2
    return rows


def stub_ok(chart, lo, hi, j):
    """rows of the stub edges: C_k' > 0 and the comparison polynomial < 0 at all 64 corners"""
    polys = []
    if j == 3:
        polys = [(lambda v, t, u: t * v['K'] - u * v['M'], lambda v, t, u: v['K'] * (1 - 2 * t) + t * u * v['M']),
                 (lambda v, t, u: u * v['q3y'] - t * v['L'], lambda v, t, u: v['q3y'] * (1 - 2 * u) + t * u * v['L'])]
    else:
        polys = [(lambda v, t, u: t * v['q2y'] + u * v['M'], lambda v, t, u: v['q2y'] * (1 - 2 * t) - t * u * v['M']),
                 (lambda v, t, u: u * v['D'] - t * v['L'], lambda v, t, u: v['D'] * (1 - 2 * u) + t * u * v['L'])]
    vals = [invariants(chart, c) for c in wcorners(lo[:4], hi[:4])]
    for C, P in polys:
        if all(C(v, t, u) > 0 and P(v, t, u) < 0 for v in vals for t in (lo[4], hi[4]) for u in (lo[5], hi[5])):
            return True
    return False


class I:
    """closed rational intervals (exact)"""
    def __init__(s, a, b=None):
        s.a = Fr(a); s.b = Fr(a if b is None else b)
    def __add__(s, o):
        o = o if isinstance(o, I) else I(o); return I(s.a + o.a, s.b + o.b)
    __radd__ = __add__
    def __neg__(s):
        return I(-s.b, -s.a)
    def __sub__(s, o):
        return s + (-(o if isinstance(o, I) else I(o)))
    def __rsub__(s, o):
        return I(o) - s
    def __mul__(s, o):
        o = o if isinstance(o, I) else I(o); v = [s.a * o.a, s.a * o.b, s.b * o.a, s.b * o.b]; return I(min(v), max(v))
    __rmul__ = __mul__
    def __truediv__(s, o):
        o = o if isinstance(o, I) else I(o)
        if o.a <= 0 <= o.b:
            return None
        return s * I(1 / o.b, 1 / o.a)


def face_rows(x, b, r=None, s=None):
    """E3 corner face (FAR_PROOFS X(iii)): invariants / q2y as functions of x = q2x, b and r = h/q2y (or s = 1/r, after
    multiplication by s^2); returns K/q2y (sign only) and, per (j, row), (sign quantity that must be > 0, g)"""
    if s is None:
        K = r * (x - 1) + (1 - b)
        rows = {3: [(r * (x - 1) - b, lambda: (b + r) / (x * r * (r * (x - 1) - b))),
                    (r - (1 - b), lambda: K / ((r - (1 - b)) * x * r))],
                0: [(r - (1 - b), lambda: K / ((b + r) * (r - (1 - b)))),
                    (b - r * (x - 1), lambda: K * x * r / ((b + r) * (b - r * (x - 1))))]}
    else:
        K = (x - 1) + (1 - b) * s
        rows = {3: [((x - 1) - b * s, lambda: s * (b * s + 1) / (x * ((x - 1) - b * s))),
                    (1 - (1 - b) * s, lambda: s * K / (x * (1 - (1 - b) * s)))],
                0: [(1 - (1 - b) * s, lambda: s * K / ((b * s + 1) * (1 - (1 - b) * s))),
                    (b * s - (x - 1), lambda: x * K / ((b * s + 1) * (b * s - (x - 1))))]}
    return K, rows


def face_X(lo, hi, j, depth=26):
    """X(iii): for all r in (0, inf): no member (K <= 0) or one untouched-edge row holds with t u >= inf t inf u"""
    if not (lo[1] == 0 and lo[3] == 0 and lo[2] >= 0):
        return False
    x, b, tu = I(lo[0], hi[0]), I(lo[2], hi[2]), lo[4] * lo[5]
    def piece(a, c, use_s, d):
        v = I(a, c)
        K, rows = face_rows(x, b, s=v) if use_s else face_rows(x, b, r=v)
        if K.b <= 0:
            return True
        for sg, g in rows[j]:
            if sg.a > 0:
                G = g()
                if G is not None and G.b < tu:
                    return True
        if d == 0:
            return False
        m = (a + c) / 2
        return piece(a, m, use_s, d - 1) and piece(m, c, use_s, d - 1)
    return piece(Fr(0), Fr(1), False, depth) and piece(Fr(0), Fr(1), True, depth)


def check(chart, key, j, lo, hi, label, table):
    lo4, hi4 = lo[:4], hi[:4]
    R, vals = ranges(chart, lo4, hi4)
    if label == 'V':
        return any(R[k][1] <= 0 for k in ('q2y', 'q3y', 'K', 'D'))
    if label == 'E2':
        return chart == 'std' and R['q3y'][0] >= Fr(1, 4) and all(v['K'] < H2 * (1 - v['q3x']) for v in vals)
    if label == 'E3':
        return chart == 'std' and R['q2x'][0] >= Fr(1, 4) and all(v['D'] < H3 * v['q2x'] for v in vals)
    if label in ('T', 'TX'):
        LB = lb_quad(R, table, per_lower(chart, lo4, hi4))
        def ok(LB):
            if LB <= FSTAR:
                return False
            th = (1 - FSTAR / LB) / tau_bar(R, j)
            return hi[4] * hi[5] < th if label == 'T' else any(th > g for g in g_rows(R, j))
        if ok(LB):
            return True
        tb = tau_bar(R, j)
        x = hi[4] * hi[5] if label == 'T' else min(g_rows(R, j), default=Fr(1))   # T: sup t u; TX: the smallest g_bar
        need = FSTAR / (1 - tb * x) if tb * x < 1 else None
        return ok(max(LB, quad_fibre(chart, lo4, hi4, need)))
    if label == 'X':
        if any(lo[4] * lo[5] > g for g in g_rows(R, j)):
            return True
        if stub_ok(chart, lo, hi, j):
            return True
        return chart == 'E3' and face_X(lo, hi, j)
    if label in ('S', 'R', 'O', 'N'):
        return pent_labels(chart, lo, hi, j, label)
    return None


def main():
    out = sys.argv[1]; K = 200; seed = 1
    if '--sample' in sys.argv:
        K = int(sys.argv[sys.argv.index('--sample') + 1])
    if '--seed' in sys.argv:
        seed = int(sys.argv[sys.argv.index('--seed') + 1])
    charts = ('std', 'E2', 'E3'); part = (0, 1)
    if '--chart' in sys.argv:
        charts = (sys.argv[sys.argv.index('--chart') + 1],)
    if '--part' in sys.argv:                   # --part i/n: only the blocks with index = i mod n (for parallel runs)
        a, b = sys.argv[sys.argv.index('--part') + 1].split('/'); part = (int(a), int(b))
    rnd = random.Random(seed); stats = collections.defaultdict(collections.Counter); bad = []
    for chart in charts:
        rf = os.path.join(out, chart, 'roots.json')
        if not os.path.exists(rf):
            continue
        roots = {}
        for name, lst in json.load(open(rf))['blocks'].items():
            for r in lst:
                roots[tuple(r[0])] = r[3] if len(r) > 3 else 0
        parts = {}
        qf = os.path.join(out, 'qtable', 'cells.jsonl')
        if chart == 'std' and os.path.exists(qf):
            for line in open(qf):
                c = json.loads(line)
                if c.get('lg') and c['lg'].get('parts'):
                    parts[tuple(c['key'])] = [([Fr(x) for x in a], [Fr(x) for x in b], Fr(v)) for a, b, v, _ in c['lg']['parts']]
        for bi, lf in enumerate(sorted(glob.glob(os.path.join(out, chart, 'blocks', '*', 'leaves.jsonl')))):
            if bi % part[1] != part[0]:
                continue
            lines = open(lf).read().splitlines()
            pick = lines if K == 0 or len(lines) <= K else rnd.sample(lines, K)
            for line in pick:
                key, j, lo, hi, label = json.loads(line)
                lo = [Fr(x) for x in lo]; hi = [Fr(x) for x in hi]
                tab = roots.get(tuple(key), 0) if chart == 'std' else 0
                pl = parts.get(tuple(key))
                if pl:                         # the table parts that meet the w-box in positive volume (each closed part:
                    # F(Q) > value on it; every point of the box lies in one of them, since they cover the root)
                    hit = [v for a, b, v in pl if all(a[i] < hi[i] and lo[i] < b[i] for i in range(4))]
                    if hit:
                        tab = max(Fr(tab), min(hit))
                r = check(chart, key, j, lo, hi, label, tab or 0)
                st = stats[(chart, label)]
                if r is None:
                    st['not re-checked here'] += 1
                elif r:
                    st['re-proved'] += 1
                else:
                    st['not reproduced'] += 1
                    if len(bad) < 30:
                        bad.append((chart, os.path.basename(os.path.dirname(lf)), key, j, label))
    for (chart, label), st in sorted(stats.items()):
        print('%-4s %-3s %s' % (chart, label, dict(st)))
    print('not reproduced (first 30):', bad)


if __name__ == '__main__':
    main()
