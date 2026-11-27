#!/usr/bin/env python3
"""LG box certificate, PROTOTYPE (Claude, 2026-11-16; not used by any pipeline, not reviewed).

A trial FAMILY for a whole box of polygons.  lgfan.py (far/LG_LEMMA.md, computation B) writes u and chi in Bernstein
form on the fan (O; v_0..v_nv-1).  The construction of the flux w = -(p - O)(2 + E)^-1 u + curl chi and the proof that
u is in H^1_0, w in H(div), -div w = u use only the barycentric coordinates of the triangles; they hold for EVERY fan
with the same coefficient lists as long as all triangles (O, v_k, v_k+1) are positively oriented (the fan is then a
polygon P, star-shaped with respect to O; convexity is not needed).  Here the vertices are functions of parameters z in
a box X with centre zc, and the coefficient lists are affine in z:

      u(z) = u_0 + sum_i (z_i - zc_i) u_i,       chi(z) = chi_0 + sum_i (z_i - zc_i) chi_i        (order 1; order 0: u_0, chi_0 only).

For every z this is an admissible pair on P(z).  A0 = |grad u|^2, A1 = |u|^2, A2 = |w|^2 are sums over the triangles
of (quadratic polynomial in z - zc with exact rational coefficients, lgfan.tri_bilinear) x (rational function of the
vertices, see lgfan.tri_forms).  With nu(P) = KS / |P| <= lambda_2(P) (Krahn-Szego), Lemma LG gives for every z

      Dc := KS A1 - |P| A0 > 0     =>     F(P) = lambda_1 |P| >= |P| Dc / (KS A2 - |P| A1),

so F(P(z)) > F* follows from  Dc(z) > 0  and  C(z) := |P| Dc - F* (KS A2 - |P| A1) > 0.  Both are verified on X in ball
arithmetic: mean value form f(z) in f(mid) + grad f(X) . (X - mid), gradient by forward differentiation, bisection.
The lists come from floating point (Rayleigh-Ritz vector and least-squares flux at zc; central differences for
u_i, chi_i) and are rounded to integers times a power of two; the proof does not depend on their quality.
Order 0: the bound differs from F(z) by O(|z - zc|^2); order 1: by O(|z - zc|^4)."""
import sys, os, math, time
from fractions import Fraction
import numpy as np
import scipy.linalg as sla
from flint import arb, fmpq
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import lgfan as LF

NAMES = ('m', 'aa', 'ab', 'bb', 'qaa', 'qab', 'qbb', 'xa1', 'xb1', 'xa2', 'xb2', 'raa', 'rab', 'rbb')


def _arb(q):
    """ball enclosing a rational (fmpq, Fraction or int)"""
    q = Fraction(int(q.p), int(q.q)) if isinstance(q, fmpq) else Fraction(q)
    return arb(q.numerator) / arb(q.denominator)


class Dual:
    """value and gradient (balls): forward differentiation"""
    __slots__ = ('v', 'g')

    def __init__(self, v, g):
        self.v = v; self.g = g

    def __add__(self, o):
        if isinstance(o, Dual):
            return Dual(self.v + o.v, [a + b for a, b in zip(self.g, o.g)])
        return Dual(self.v + o, self.g)
    __radd__ = __add__

    def __sub__(self, o):
        if isinstance(o, Dual):
            return Dual(self.v - o.v, [a - b for a, b in zip(self.g, o.g)])
        return Dual(self.v - o, self.g)

    def __rsub__(self, o):
        return Dual(o - self.v, [-a for a in self.g])

    def __neg__(self):
        return Dual(-self.v, [-a for a in self.g])

    def __mul__(self, o):
        if isinstance(o, Dual):
            return Dual(self.v * o.v, [a * o.v + self.v * b for a, b in zip(self.g, o.g)])
        return Dual(self.v * o, [a * o for a in self.g])
    __rmul__ = __mul__

    def __truediv__(self, o):
        if isinstance(o, Dual):
            q = self.v / o.v
            return Dual(q, [(a - q * b) / o.v for a, b in zip(self.g, o.g)])
        return Dual(self.v / o, [a / o for a in self.g])


class J2:
    """value, gradient and Hessian (balls): second-order forward differentiation.  h[IDX[i][j]] = d^2/dz_i dz_j, i <= j"""
    __slots__ = ('v', 'g', 'h')

    def __init__(self, v, g, h):
        self.v = v; self.g = g; self.h = h

    @staticmethod
    def const(v, nd):
        z = arb(0); return J2(v, [z] * nd, [z] * (nd * (nd + 1) // 2))

    @staticmethod
    def var(v, i, nd):
        z = arb(0); return J2(v, [arb(1) if j == i else z for j in range(nd)], [z] * (nd * (nd + 1) // 2))

    def __add__(self, o):
        if isinstance(o, J2):
            return J2(self.v + o.v, [a + b for a, b in zip(self.g, o.g)], [a + b for a, b in zip(self.h, o.h)])
        return J2(self.v + o, self.g, self.h)
    __radd__ = __add__

    def __sub__(self, o):
        if isinstance(o, J2):
            return J2(self.v - o.v, [a - b for a, b in zip(self.g, o.g)], [a - b for a, b in zip(self.h, o.h)])
        return J2(self.v - o, self.g, self.h)

    def __neg__(self):
        return J2(-self.v, [-a for a in self.g], [-a for a in self.h])

    def __mul__(self, o):
        if isinstance(o, J2):
            nd = len(self.g); g1 = self.g; g2 = o.g; v1 = self.v; v2 = o.v; h = []; k = 0
            for i in range(nd):
                for j in range(i, nd):
                    h.append(self.h[k] * v2 + g1[i] * g2[j] + g1[j] * g2[i] + v1 * o.h[k]); k += 1
            return J2(v1 * v2, [a * v2 + v1 * b for a, b in zip(g1, g2)], h)
        return J2(self.v * o, [a * o for a in self.g], [a * o for a in self.h])
    __rmul__ = __mul__

    def inv(self):
        r = 1 / self.v; r2 = r * r; nd = len(self.g); g = self.g; h = []; k = 0
        for i in range(nd):
            for j in range(i, nd):
                h.append((2 * g[i] * g[j] * r - self.h[k]) * r2); k += 1
        return J2(r, [-a * r2 for a in g], h)

    def __truediv__(self, o):
        if isinstance(o, J2):
            return self * o.inv()
        return J2(self.v / o, [a / o for a in self.g], [a / o for a in self.h])


KS = _arb(LF.KS_LO)
HALF = arb(1) / 2
QP = True                    # second-order test: use the Lagrangian bound of the quadratic model (_qp_lower)


def conditions(V, O, T):
    """(area, A0, A1, A2, [det_k]) for the fan (O; V) with the triangle numbers T[k][name]; arb or Dual"""
    nv = len(V); A0 = A1 = A2 = area = None; dets = []
    for k in range(nv):
        e1x = V[k][0] - O[0]; e1y = V[k][1] - O[1]; e2x = V[(k + 1) % nv][0] - O[0]; e2y = V[(k + 1) % nv][1] - O[1]; t = T[k]
        n11 = e1x * e1x + e1y * e1y; n12 = e1x * e2x + e1y * e2y; n22 = e2x * e2x + e2y * e2y; det = e1x * e2y - e1y * e2x; dets.append(det)
        a1 = det * t['m']
        a0 = (n22 * t['aa'] - (n12 * t['ab']) * 2 + n11 * t['bb']) / det
        a2 = det * (n11 * t['qaa'] + (n12 * t['qab']) * 2 + n22 * t['qbb']) + (n12 * (t['xa1'] - t['xb2']) + n22 * t['xb1'] - n11 * t['xa2']) * 2 \
            + (n22 * t['raa'] - (n12 * t['rab']) * 2 + n11 * t['rbb']) / det
        A0 = a0 if A0 is None else A0 + a0; A1 = a1 if A1 is None else A1 + a1; A2 = a2 if A2 is None else A2 + a2
        area = det * HALF if area is None else area + det * HALF
    return area, A0, A1, A2, dets


def targets(area, A0, A1, A2, fstar):
    Dc = A1 * KS - area * A0
    return area * Dc - (A2 * KS - area * A1) * fstar, Dc


class Family:
    """trial family of a box.  verts_f: float vertex map z -> list of (x, y); zc: centre (binary64 numbers).
    After construction: ok, O (two binary64 numbers), q = 1 + (number of parameters if order 1), and the symmetric
    tables S[k][name][a][b] (balls enclosing exact rationals) of the 14 forms of every triangle."""

    def __init__(self, verts_f, zc, N=8, extra=1, order=1, hstep=None, O=None):
        self.zc = [float(v) for v in zc]; nd = len(self.zc); self.N = N; self.n = N + extra; self.order = order; self.ok = False
        Vc = verts_f(self.zc); nv = len(Vc); self.nv = nv; g = 2 ** 20
        self.O = (round(sum(float(a) for a, b in Vc) / nv * g) / g, round(sum(float(b) for a, b in Vc) / nv * g) / g)
        if O is not None:                                           # fan centre given by the caller (2026-11-17: non-convex polygons need a kernel point)
            self.O = (round(float(O[0]) * g) / g, round(float(O[1]) * g) / g)
        try:
            c0, x0 = self._pair(verts_f, self.zc); cs = [c0]; xs = [x0]
            if order >= 1:
                h = hstep or 1e-3 * max(1e-3, min(math.hypot(a - self.O[0], b - self.O[1]) for a, b in Vc))
                for i in range(nd):
                    zp = list(self.zc); zm = list(self.zc); zp[i] += h; zm[i] -= h
                    cp, xp = self._pair(verts_f, zp); cm, xm = self._pair(verts_f, zm)
                    cs.append((cp - cm) / (2 * h)); xs.append((xp - xm) / (2 * h))
        except Exception:
            return
        sx = max(np.abs(x).max() for x in xs)
        if not (sx > 0) or not all(np.isfinite(c).all() for c in cs + xs):
            return
        e = math.frexp(sx)[1]; self.xe = e - LF.CB
        self.CI = [[int(v) for v in np.rint(c * 2.0 ** LF.CB)] for c in cs]; self.XI = [[int(v) for v in np.rint(x * 2.0 ** (LF.CB - e))] for x in xs]
        self._tables(); self.ok = True

    def _pair(self, verts_f, z):
        V = verts_f(z); nv = len(V)
        geo = [((a - self.O[0]) ** 2 + (b - self.O[1]) ** 2, (a - self.O[0]) * (c - self.O[0]) + (b - self.O[1]) * (d - self.O[1]),
                (c - self.O[0]) ** 2 + (d - self.O[1]) ** 2, (a - self.O[0]) * (d - self.O[1]) - (b - self.O[1]) * (c - self.O[0]))
               for (a, b), (c, d) in ((V[k], V[(k + 1) % nv]) for k in range(nv))]
        K, M, G = LF._assemble(geo, self.N, True); c = sla.eigh(K, M, subset_by_index=[0, 0], check_finite=False)[1][:, 0]; c = c / c[0]
        H = LF._refh(self.N, self.n); Kc, _, Gc = LF._assemble(geo, self.n, False); Gd = LF.dofs(nv, self.N, True)[0]
        E = H['E'][2]; hs = np.array(H['Hn'], float) / H['L']; X = {k: v[2] for k, v in H['X'].items()}; r = np.zeros(Kc.shape[0]); loc = len(Gd[0])
        for k, (n11, n12, n22, det) in enumerate(geo):
            cT = np.zeros(loc); act = Gd[k] >= 0; cT[act] = c[Gd[k][act]]; h = hs * (E @ cT)
            r[Gc[k]] += -(-n12 * (h @ X[('a', 1)]) - n22 * (h @ X[('b', 1)]) + n11 * (h @ X[('a', 2)]) + n12 * (h @ X[('b', 2)]))
        x = np.zeros(Kc.shape[0]); x[1:] = sla.solve(Kc[1:, 1:], -r[1:], assume_a='pos', check_finite=False)
        return c, x

    def _tables(self):
        tb = LF.tri_bilinear(self.nv, self.N, self.n, self.CI, self.XI, self.xe); q = len(self.CI); self.q = q; self.S = []
        for t in tb:
            self.S.append({nm: [[_arb((t[nm][a][b] + t[nm][b][a]) / 2) for b in range(q)] for a in range(q)] for nm in NAMES})

    def forms(self, d, mode):
        """the 14 numbers of every triangle at z - zc = d (list of balls).  mode 0: values (balls); 1: Dual (value and
        gradient); 2: J2 (value, gradient, Hessian; the Hessian of a quadratic polynomial is constant)"""
        q = self.q; nd = len(d); out = []; zero = arb(0)
        for S in self.S:
            t = {}
            for nm in NAMES:
                F = S[nm]
                if q == 1:
                    t[nm] = F[0][0] if mode == 0 else (Dual(F[0][0], [zero] * nd) if mode == 1 else J2.const(F[0][0], nd)); continue
                row = [F[a][0] + sum((F[a][b + 1] * d[b] for b in range(nd)), zero) for a in range(q)]      # (F (1, d))_a
                val = row[0] + sum((row[i + 1] * d[i] for i in range(nd)), zero)
                if mode == 0:
                    t[nm] = val
                elif mode == 1:
                    t[nm] = Dual(val, [2 * row[i + 1] for i in range(nd)])
                else:
                    t[nm] = J2(val, [2 * row[i + 1] for i in range(nd)], [2 * F[i + 1][j + 1] for i in range(nd) for j in range(i, nd)])
            out.append(t)
        return out


def float_margins(fam, Vs, fstar):
    """HEURISTIC (floating point, not part of any proof; 2026-11-17): for the polygons Vs[c] = fam's vertex map at points z_c
    (array (m, nv, 2)) given together with D[c] = z_c - fam.zc (array (m, nd)) as the pair (Vs, D): the array of
    min(det_k) (must be > 0), Dc, and the quotient bound |P| Dc / (KS A2 - |P| A1) / fstar - 1 of the transplanted family
    at every point.  A box certificate can only succeed if all three are positive at every point of the box."""
    Vs, D = Vs; m = len(D)
    if not hasattr(fam, '_Sf'):
        fam._Sf = np.array([[[[float(S[nm][a][b].mid()) for b in range(fam.q)] for a in range(fam.q)] for nm in NAMES] for S in fam.S])
    dd = np.concatenate([np.ones((m, 1)), D], 1) if fam.q > 1 else np.ones((m, 1))
    t = np.einsum('ca,knab,cb->ckn', dd, fam._Sf, dd); T = {nm: t[:, :, i] for i, nm in enumerate(NAMES)}
    e1 = Vs - np.array(fam.O); e2 = np.roll(e1, -1, axis=1)
    n11 = (e1 * e1).sum(2); n12 = (e1 * e2).sum(2); n22 = (e2 * e2).sum(2); det = e1[:, :, 0] * e2[:, :, 1] - e1[:, :, 1] * e2[:, :, 0]
    with np.errstate(all='ignore'):
        A1 = (det * T['m']).sum(1); A0 = ((n22 * T['aa'] - 2 * n12 * T['ab'] + n11 * T['bb']) / det).sum(1)
        A2 = (det * (n11 * T['qaa'] + 2 * n12 * T['qab'] + n22 * T['qbb']) + 2 * (n12 * (T['xa1'] - T['xb2']) + n22 * T['xb1'] - n11 * T['xa2'])
              + (n22 * T['raa'] - 2 * n12 * T['rab'] + n11 * T['rbb']) / det).sum(1)
        area = det.sum(1) / 2; ks = float(LF.KS_LO); Dc = ks * A1 - area * A0
        return det.min(1), Dc, area * Dc / (ks * A2 - area * A1) / fstar - 1


def check_box(fam, verts, lo, hi, fstar, max_boxes=200):
    """True if F(P(z)) > fstar is PROVED for every z in the box [lo, hi] (binary64 numbers).  verts(z): vertex list for a
    list z of numbers (arb or Dual).  Returns (ok, number of boxes used, reason)."""
    nd = len(lo); O = (arb(fam.O[0]), arb(fam.O[1])); fs = arb(fstar); zero = arb(0); one = arb(1); zc = fam.zc
    stack = [(list(lo), list(hi))]; used = 0
    while stack:
        a, b = stack.pop(); used += 1
        if used > max_boxes:
            return False, used, 'budget'
        mid = [(x + y) / 2 for x, y in zip(a, b)]; rad = [(y - x) / 2 for x, y in zip(a, b)]
        X = [Dual(arb(m, r), [one if i == j else zero for j in range(nd)]) for i, (m, r) in enumerate(zip(mid, rad))]
        area, A0, A1, A2, dets = conditions(verts(X), O, fam.forms([arb(m, r) - arb(c) for m, r, c in zip(mid, rad, zc)], 1))
        if not all(d.v > 0 for d in dets):
            return False, used, 'fan'                               # a triangle of the fan may degenerate in this box
        C, Dc = targets(area, A0, A1, A2, fs)
        am, a0, a1, a2, _ = conditions(verts([arb(m) for m in mid]), O, fam.forms([arb(m) - arb(c) for m, c in zip(mid, zc)], 0))
        Cm, Dm = targets(am, a0, a1, a2, fs)
        wC = [g.abs_upper() * arb(r) for g, r in zip(C.g, rad)]; wD = [g.abs_upper() * arb(r) for g, r in zip(Dc.g, rad)]
        okC = (Cm - sum(wC, zero) > 0) or (C.v > 0); okD = (Dm - sum(wD, zero) > 0) or (Dc.v > 0)
        if okC and okD:
            continue
        if not (Cm > 0 and Dm > 0):
            return False, used, 'midpoint'                          # the condition is not provable at the midpoint itself
        w = wC if not okC else wD; i = max(range(nd), key=lambda j: float(w[j].mid()) if rad[j] > 0 else -1.0)
        if not (rad[i] > 0):
            return False, used, 'point'
        stack.append((a[:i] + [mid[i]] + a[i + 1:], b)); stack.append((a, b[:i] + [mid[i]] + b[i + 1:]))
    return True, used, 'ok'


def _taylor2_lower(fm, gm, H, rad):
    """lower bound (ball arithmetic, a float rounded down or None) of  fm + g.d + d^T H d / 2  over |d_i| <= rad_i,  H any
    symmetric matrix inside the ball matrix H (list over i <= j): for every i the one-variable minimum of
    g_i t + lo(H_ii) t^2 / 2 on [-rad_i, rad_i] (convexity is used when lo(H_ii) > 0), minus sum_{i<j} |H_ij| rad_i rad_j.
    Also returns the per-axis losses (floats, for the choice of the bisection axis)."""
    nd = len(rad); tot = fm; loss = [0.0] * nd; k = 0; idx = {}
    for i in range(nd):
        for j in range(i, nd):
            idx[(i, j)] = k; k += 1
    for i in range(nd):
        r = arb(rad[i]); g = gm[i].abs_upper(); hl = H[idx[(i, i)]].lower()          # hl <= H_ii
        if hl > 0 and g < hl * r:
            t = -(g * g) / (2 * hl)                                 # interior minimum of |g| t + hl t^2 / 2 ... >= -g^2 / (2 hl)
        else:
            t = -g * r + (hl * r * r) / 2                           # at the end point (also valid when hl <= 0)
        tot = tot + t; loss[i] += float((-t).upper()) if t.is_finite() else float('inf')
    for i in range(nd):
        for j in range(i + 1, nd):
            p = H[idx[(i, j)]].abs_upper() * arb(rad[i]) * arb(rad[j]); tot = tot - p
            pf = float(p.upper()) if p.is_finite() else float('inf'); loss[i] += pf / 2; loss[j] += pf / 2
    return tot, loss


def _qp_lower(fm, gm, H, rad):
    """lower bound (ball) of  fm + g.d + d^T H d / 2  over the box |d_i| <= rad_i  for EVERY symmetric H in the ball matrix
    H (list over i <= j), by Lagrangian duality:  H >= H0 := mid(H) - s I  (Loewner order; s = max row sum of the radii,
    Gershgorin), and for any mu >= 0 with A := H0 + 2 diag(mu) positive definite
         g.d + d^T H0 d / 2  >=  g.d + d^T A d / 2 - sum mu_i rad_i^2  >=  - g^T A^-1 g / 2 - sum mu_i rad_i^2      (|d_i| <= rad_i).
    mu comes from a floating-point solution of the box-constrained problem (KKT multipliers); positive definiteness and
    the solve are done in ball arithmetic (LDL^T with certainly positive pivots).  Returns the ball or None."""
    nd = len(rad); idx = {}; k = 0
    for i in range(nd):
        for j in range(i, nd):
            idx[(i, j)] = k; k += 1
    Hm = np.empty((nd, nd)); Rr = np.empty((nd, nd))
    for i in range(nd):
        for j in range(i, nd):
            e = H[idx[(i, j)]]
            if not e.is_finite():
                return None
            Hm[i, j] = Hm[j, i] = float(e.mid()); Rr[i, j] = Rr[j, i] = float(e.rad().upper())
    sarb = max((sum((H[idx[(min(i, j), max(i, j))]].rad() for j in range(nd)), arb(0)) for i in range(nd)), key=lambda t: float(t.upper()))
    sfl = float(sarb.upper()); g = np.array([float(t.mid()) for t in gm]); r = np.array(rad, float); H0 = Hm - sfl * np.eye(nd)
    # float: projected gradient + KKT multipliers of min g.d + d^T H0 d / 2 on the box
    d = np.clip(-g / np.maximum(np.diag(H0), 1e-300), -r, r) if np.all(np.diag(H0) > 0) else -np.sign(g) * r
    L = max(np.abs(H0).sum(1).max(), 1e-300)
    for _ in range(200):
        dn = np.clip(d - (g + H0 @ d) / L, -r, r)
        if np.abs(dn - d).max() <= 1e-12 * r.max():
            d = dn; break
        d = dn
    gr = g + H0 @ d; mu = np.zeros(nd)
    for i in range(nd):
        if abs(d[i]) >= r[i] * (1 - 1e-9) and r[i] > 0:
            mu[i] = max(0.0, -gr[i] / (2 * d[i]))
    lam = np.linalg.eigvalsh(H0 + 2 * np.diag(mu))[0]
    if lam <= 1e-9 * L:
        mu = mu + (1e-6 * L - lam) / 2                              # make A safely positive definite
    A = [[arb(float(Hm[i, j])) - (sarb if i == j else 0) + (2 * arb(float(mu[i])) if i == j else 0) for j in range(nd)] for i in range(nd)]
    # LDL^T in balls
    Lm = [[arb(0)] * nd for _ in range(nd)]; Dg = [None] * nd
    for j in range(nd):
        dj = A[j][j] - sum((Lm[j][k] * Lm[j][k] * Dg[k] for k in range(j)), arb(0))
        if not (dj > 0):
            return None
        Dg[j] = dj
        for i in range(j + 1, nd):
            Lm[i][j] = (A[i][j] - sum((Lm[i][k] * Lm[j][k] * Dg[k] for k in range(j)), arb(0))) / dj
    y = [None] * nd                                                 # L y = g
    for i in range(nd):
        y[i] = gm[i] - sum((Lm[i][k] * y[k] for k in range(i)), arb(0))
    quad = sum((y[i] * y[i] / Dg[i] for i in range(nd)), arb(0))    # g^T A^-1 g
    return fm - quad / 2 - sum((arb(float(mu[i])) * arb(rad[i]) * arb(rad[i]) for i in range(nd)), arb(0))


def check_box2(fam, verts2, lo, hi, fstar, max_boxes=200):
    """as check_box, with the second-order form  f(z) >= f(mid) + grad f(mid).d + d^T H(X) d / 2  (H(X): ball Hessian on the
    box by second-order forward differentiation).  verts2(z, kind): vertex list for z given as J2 (kind 2), Dual (1), arb (0)."""
    nd = len(lo); O = (arb(fam.O[0]), arb(fam.O[1])); fs = arb(fstar); zc = fam.zc; one = arb(1); zero = arb(0)
    stack = [(list(lo), list(hi))]; used = 0
    while stack:
        a, b = stack.pop(); used += 1
        if used > max_boxes:
            return False, used, 'budget'
        mid = [(x + y) / 2 for x, y in zip(a, b)]; rad = [(y - x) / 2 for x, y in zip(a, b)]
        X = [J2.var(arb(m, r), i, nd) for i, (m, r) in enumerate(zip(mid, rad))]
        area, A0, A1, A2, dets = conditions(verts2(X), O, fam.forms([arb(m, r) - arb(c) for m, r, c in zip(mid, rad, zc)], 2))
        if not all(d.v > 0 for d in dets):
            return False, used, 'fan'
        C, Dc = targets(area, A0, A1, A2, fs)
        Xm = [Dual(arb(m), [one if i == j else zero for j in range(nd)]) for i, m in enumerate(mid)]
        am, a0, a1, a2, _ = conditions(verts2(Xm), O, fam.forms([arb(m) - arb(c) for m, c in zip(mid, zc)], 1))
        Cm, Dm = targets(am, a0, a1, a2, fs)
        if not (Cm.v > 0 and Dm.v > 0):
            return False, used, 'midpoint'
        lC, lossC = _taylor2_lower(Cm.v, Cm.g, C.h, rad); lD, lossD = _taylor2_lower(Dm.v, Dm.g, Dc.h, rad)
        okC = (lC > 0) or (C.v > 0); okD = (lD > 0) or (Dc.v > 0)
        if not okC and QP:
            v = _qp_lower(Cm.v, Cm.g, C.h, rad); okC = v is not None and v > 0
        if not okD and QP:
            v = _qp_lower(Dm.v, Dm.g, Dc.h, rad); okD = v is not None and v > 0
        if okC and okD:
            continue
        w = lossC if not okC else lossD; i = max(range(nd), key=lambda j: w[j] if rad[j] > 0 else -1.0)
        if not (rad[i] > 0):
            return False, used, 'point'
        stack.append((a[:i] + [mid[i]] + a[i + 1:], b)); stack.append((a, b[:i] + [mid[i]] + b[i + 1:]))
    return True, used, 'ok'


def bound_at(fam, verts, z):
    """the certified lower bound of F(P(z)) at ONE point z (list of binary64 numbers), as a ball; None if Dc is not > 0"""
    O = (arb(fam.O[0]), arb(fam.O[1]))
    area, A0, A1, A2, dets = conditions(verts([arb(v) for v in z]), O, fam.forms([arb(v) - arb(c) for v, c in zip(z, fam.zc)], 0))
    Dc = A1 * KS - area * A0
    if not all(d > 0 for d in dets) or not (Dc > 0):
        return None
    return area * Dc / (A2 * KS - area * A1)


def near_verts(x):
    """the near chart: (0, 0), (1, 0), (x0, x1), (x2, x3), (x4, x5); numbers: float, arb or Dual"""
    if isinstance(x[0], J2):
        z = J2.const(arb(0), len(x[0].g)); o = J2.const(arb(1), len(x[0].g))
    elif isinstance(x[0], Dual):
        z = Dual(arb(0), [arb(0)] * len(x[0].g)); o = Dual(arb(1), [arb(0)] * len(x[0].g))
    elif isinstance(x[0], arb):
        z = arb(0); o = arb(1)
    else:
        z = 0.0; o = 1.0
    return [(z, z), (o, z), (x[0], x[1]), (x[2], x[3]), (x[4], x[5])]


if __name__ == '__main__':
    FSTAR = 18.919134195733914
    c72, s72 = math.cos(2 * math.pi / 5), math.sin(2 * math.pi / 5); XREG = np.array([1 + c72, s72, 0.5, s72 + math.sin(math.pi / 5), -c72, s72])
    rng = np.random.default_rng(1); N = int(sys.argv[1]) if len(sys.argv) > 1 else 8; mb = int(sys.argv[2]) if len(sys.argv) > 2 else 300; T2 = len(sys.argv) > 3 and sys.argv[3] == 't2'
    for r in (0.05, 0.1, 0.2, 0.4, 0.6):
        for trial in range(2):
            d = rng.normal(size=6); d /= np.linalg.norm(d); xc = [round(float(v) * 2 ** 16) / 2 ** 16 for v in XREG + r * d]
            for order in (0, 1):
                t = time.time(); P = Family(near_verts, xc, N, order=order); tf = time.time() - t; out = []
                for k in range(2, 11):                              # half-width 2^-k
                    rho = 2.0 ** -k; t = time.time()
                    ok, used, why = (check_box2 if T2 else check_box)(P, near_verts, [v - rho for v in xc], [v + rho for v in xc], FSTAR, mb)
                    out.append('2^-%d:%s/%d/%.1fs' % (k, why, used, time.time() - t))
                    if ok:
                        break
                print('r %.2f order %d (family %.2f s)  %s   rho / r = %.3f' % (r, order, tf, '  '.join(out[-3:]), rho / r if ok else 0), flush=True)
