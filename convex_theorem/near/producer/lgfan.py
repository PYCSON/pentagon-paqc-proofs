#!/usr/bin/env python3
"""Polynomial certificate on the fan triangulation (Claude, 2026-11-16): a rigorous lower bound of F(P) = lambda_1(P) |P|
for a convex polygon P, from Lemma LG (far/LG_LEMMA.md) with PIECEWISE polynomial trial functions.  Integer arithmetic only
in the proof part; no eigenvalue certificate, no Bessel functions.

  1. frame p -> (p - c)/s and a polygon P+ with vertices on the grid 2^-GRID that contains the image of P (exact test);
     lambda_1(P) >= s^-2 lambda_1(P+).  O: a grid point strictly inside P+.  Fan: the triangles T_k = (O, v_k, v_k+1).
  2. u: continuous, a polynomial of degree N on every T_k, zero on the boundary of P+.  Bernstein form: on T_k with the
     barycentric coordinates (lam, a, b) of (O, v_k, v_k+1), u = sum c_m N!/(i! j! l!) lam^i a^j b^l; coefficients on a
     common spoke are shared (so u is continuous), coefficients with i = 0 vanish (so u = 0 on the boundary).
  3. flux w = -(p - O) h + curl chi with h = (2 + E)^-1 u on every T_k (E the Euler operator about O; in the coordinates
     (a, b), which are linear in p - O: a^q b^r -> a^q b^r / (q + r + 2)) and chi continuous, piecewise of degree n.
     Then -div w = u on every T_k, the first part is continuous across the spokes (h on a spoke depends only on u on that
     spoke) and curl chi has continuous normal components: w is in H(div; P+) and -div w = u.
  4. A0 = |grad u|^2, A1 = |u|^2, A2 = |w|^2 exactly (quadratic forms with fixed rational matrices of the reference
     triangle, integral of lam^p a^q b^r over T = 2|T| p! q! r! / (p + q + r + 2)!).
  5. nu <= lambda_2(P+) (bounding rectangle or Krahn-Szego), D = A0 - nu A1 < 0, G = A0 - 2 nu A1 + nu^2 A2:
     lambda_1(P+) >= nu (-D)/(G - D)      (Lemma LG).
The coefficients of u and chi come from floating point (Rayleigh-Ritz; least squares) and are rounded to integers times a
power of two; the proof does not depend on their quality.  API: F_lower(V, N) -> record or None; recheck(V, record)."""
import math
from fractions import Fraction
from math import factorial as _f
import numpy as np
import scipy.linalg as sla
from flint import fmpq, fmpz_mat

GRID = 24                    # P+ and O are on the grid 2^-GRID of the frame
EPS0 = 2.0 ** -21            # first inflation 1 + EPS0 of P about its centroid (costs about 2 EPS0 of the bound)
CB = 30                      # coefficients are integers times 2^-CB (u) or 2^(e - CB) (chi)
SLACK = 1e-7                 # stated F_lo = proved bound times (1 - SLACK): a re-check with another flux still passes
PI2_LO = fmpq(98696044, 10 ** 7)                                   # 9.8696044 < pi^2 = 9.8696044010...
KS_LO = fmpq(3633682907, 10 ** 8)                                # 36.33682907 < 2 pi j_{0,1}^2 = 36.336829071...
_cache = {}


def _I(p, q, r):
    """(integral of lam^p a^q b^r over a triangle) / (2 |T|)"""
    return Fraction(_f(p) * _f(q) * _f(r), _f(p + q + r + 2))


def _loc(N):
    return [(i, j, N - i - j) for i in range(N + 1) for j in range(N + 1 - i)]


def _D(e, beta):
    """D_beta = d/d beta - d/d lam (beta = 1: a, 2: b) of N!/(i!j!l!) lam^i a^j b^l: list of (coefficient, exponents)"""
    i, j, l = e; C = _f(i + j + l) // (_f(i) * _f(j) * _f(l)); out = []
    if beta == 1 and j:
        out.append((C * j, (i, j - 1, l)))
    if beta == 2 and l:
        out.append((C * l, (i, j, l - 1)))
    if i:
        out.append((-C * i, (i - 1, j, l)))
    return out


def _pack(rows):
    """a table of Fractions -> (fmpz_mat, common denominator, float array)"""
    den = 1
    for r in rows:
        for v in r:
            den = den * v.denominator // math.gcd(den, v.denominator)
    return fmpz_mat([[int(v * den) for v in r] for r in rows]), den, np.array([[float(v) for v in r] for r in rows])


def _ref(N):
    """reference matrices of degree N: mass M0 and the three gradient forms Saa, Sab, Sbb (all divided by 2|T|):
    integral of grad f . grad g over T = (|e2|^2 Saa - e1.e2 (Sab + Sab^T) + |e1|^2 Sbb) / det for the Bernstein
    coefficient vectors, e1 = v_k - O, e2 = v_k+1 - O, det = e1 x e2 = 2|T|;  integral of f g = det M0."""
    if ('ref', N) not in _cache:
        loc = _loc(N); C = [_f(N) // (_f(i) * _f(j) * _f(l)) for i, j, l in loc]
        M0 = [[C[m] * C[k] * _I(loc[m][0] + loc[k][0], loc[m][1] + loc[k][1], loc[m][2] + loc[k][2]) for k in range(len(loc))] for m in range(len(loc))]
        Da = [_D(e, 1) for e in loc]; Db = [_D(e, 2) for e in loc]

        def form(A, B):
            return [[sum(ca * cb * _I(ea[0] + eb[0], ea[1] + eb[1], ea[2] + eb[2]) for ca, ea in A[m] for cb, eb in B[k]) for k in range(len(loc))] for m in range(len(loc))]
        _cache[('ref', N)] = dict(loc=loc, M0=_pack(M0), Saa=_pack(form(Da, Da)), Sab=_pack(form(Da, Db)), Sbb=_pack(form(Db, Db)))
    return _cache[('ref', N)]


def _refh(N, n):
    """matrices for the flux: E (Bernstein coefficients of degree N -> coefficients of the monomials a^q b^r, q + r <= N),
    Hn (the integers L/(q + r + 2), L their common multiple: h = Hn E c / L), Qaa, Qab, Qbb (integrals of a^2, ab, b^2
    times products of two monomials, / 2|T|) and X[s][beta] (integral of s a^q b^r D_beta B_m, B_m of degree n, / 2|T|)."""
    if ('refh', N, n) not in _cache:
        loc = _loc(N); mon = [(q, r) for q in range(N + 1) for r in range(N + 1 - q)]; pos = {m: k for k, m in enumerate(mon)}
        E = [[Fraction(0)] * len(loc) for _ in mon]
        for m, (i, j, l) in enumerate(loc):
            C = _f(N) // (_f(i) * _f(j) * _f(l))
            for s in range(i + 1):
                for t in range(i + 1 - s):                          # (1 - a - b)^i = sum (-1)^(s+t) i!/(s! t! (i-s-t)!) a^s b^t
                    E[pos[(j + s, l + t)]][m] += C * (-1) ** (s + t) * (_f(i) // (_f(s) * _f(t) * _f(i - s - t)))
        L = 1
        for k in range(2, N + 3):
            L = L * k // math.gcd(L, k)
        Hn = [L // (q + r + 2) for q, r in mon]

        def Q(da, db):
            return [[_I(0, q + q2 + da, r + r2 + db) for q2, r2 in mon] for q, r in mon]
        locn = _loc(n); X = {}
        for s, (da, db) in (('a', (1, 0)), ('b', (0, 1))):
            for beta in (1, 2):
                Dn = [_D(e, beta) for e in locn]
                X[(s, beta)] = _pack([[sum(c * _I(e[0], e[1] + q + da, e[2] + r + db) for c, e in Dn[m]) for m in range(len(locn))] for q, r in mon])
        _cache[('refh', N, n)] = dict(mon=mon, E=_pack(E), Hn=Hn, L=L, Qaa=_pack(Q(2, 0)), Qab=_pack(Q(1, 1)), Qbb=_pack(Q(0, 2)), X=X)
    return _cache[('refh', N, n)]


def dofs(nv, N, dirichlet):
    """global numbering of the continuous piecewise polynomials of degree N on a fan with nv triangles.
    Order: the centre; the spokes k = 0..nv-1 (spoke k = O v_k), each with i = 1..N-1 (i the exponent of lam);
    [only if not dirichlet: the vertices v_0..v_nv-1, then the edges k = 0..nv-1, each with j = 1..N-1];
    the interiors of T_0..T_nv-1, each in the order of _loc(N) (i ascending, then j ascending).
    Returns (table G[k][m] = global index of the local coefficient m of T_k, or -1 for a boundary coefficient; count)."""
    key = ('dofs', nv, N, dirichlet)
    if key not in _cache:
        loc = _loc(N); cnt = 1; spoke = []
        for k in range(nv):
            spoke.append(list(range(cnt, cnt + N - 1))); cnt += N - 1
        vert = [-1] * nv; edge = [[-1] * (N - 1) for _ in range(nv)]
        if not dirichlet:
            vert = list(range(cnt, cnt + nv)); cnt += nv
            for k in range(nv):
                edge[k] = list(range(cnt, cnt + N - 1)); cnt += N - 1
        G = np.empty((nv, len(loc)), int)
        for k in range(nv):
            for m, (i, j, l) in enumerate(loc):
                if i == N:
                    g = 0
                elif i == 0 and l == 0:
                    g = vert[k]
                elif i == 0 and j == 0:
                    g = vert[(k + 1) % nv]
                elif i == 0:
                    g = edge[k][j - 1]
                elif l == 0:
                    g = spoke[k][i - 1]
                elif j == 0:
                    g = spoke[(k + 1) % nv][i - 1]
                else:
                    g = cnt; cnt += 1
                G[k, m] = g
        _cache[key] = (G, cnt)
    return _cache[key]


def inflate(V, grid=None, eps0=None):
    """frame, P+ and O.  V: exact binary64 vertices (counter-clockwise, strictly convex).  Returns (c, s, Uplus, O, area of P)
    with exact Fractions, or None.  Containment of the image of P in P+ and of O in the interior of P+ are checked exactly."""
    grid = GRID if grid is None else grid; eps = EPS0 if eps0 is None else eps0
    Vf = [(Fraction(float(a)), Fraction(float(b))) for a, b in V]; nv = len(Vf)
    for i in range(nv):
        p, q, r = Vf[i], Vf[(i + 1) % nv], Vf[(i + 2) % nv]
        if (q[0] - p[0]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[0] - q[0]) <= 0:
            return None
    xs = [float(a) for a, b in Vf]; ys = [float(b) for a, b in Vf]
    half = max(max(xs) - min(xs), max(ys) - min(ys)) / 2
    if not (half > 0):
        return None
    s = Fraction(2) ** math.ceil(math.log2(half * 1.02)); k = 1024 / s
    c = (Fraction(round(((max(xs) + min(xs)) / 2) * float(k))) / k, Fraction(round(((max(ys) + min(ys)) / 2) * float(k))) / k)
    U = [((a - c[0]) / s, (b - c[1]) / s) for a, b in Vf]
    area = sum(Vf[i][0] * Vf[(i + 1) % nv][1] - Vf[(i + 1) % nv][0] * Vf[i][1] for i in range(nv)) / 2
    cg = (sum(float(a) for a, b in U) / nv, sum(float(b) for a, b in U) / nv); g = 2 ** grid
    for _ in range(4):
        Up = [(Fraction(round((cg[0] + (1 + eps) * (float(a) - cg[0])) * g), g), Fraction(round((cg[1] + (1 + eps) * (float(b) - cg[1])) * g), g)) for a, b in U]
        O = (Fraction(round(cg[0] * g), g), Fraction(round(cg[1] * g), g))
        if _admissible(U, Up, O):
            return c, s, Up, O, area
        eps *= 2
    return None


def _admissible(U, Up, O):
    """P+ strictly convex and counter-clockwise, every vertex of U in P+, O strictly inside P+ (exact)"""
    nv = len(Up)
    for i in range(nv):
        p, q, r = Up[i], Up[(i + 1) % nv], Up[(i + 2) % nv]
        if (q[0] - p[0]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[0] - q[0]) <= 0:
            return False
        if any((q[0] - p[0]) * (v[1] - p[1]) - (q[1] - p[1]) * (v[0] - p[0]) < 0 for v in U):
            return False
        if (p[0] - O[0]) * (q[1] - O[1]) - (p[1] - O[1]) * (q[0] - O[0]) <= 0:
            return False
    return True


def _geom(Up, O):
    """per triangle: (|e1|^2, e1.e2, |e2|^2, det), e1 = v_k - O, e2 = v_k+1 - O"""
    nv = len(Up); out = []
    for k in range(nv):
        e1 = (Up[k][0] - O[0], Up[k][1] - O[1]); e2 = (Up[(k + 1) % nv][0] - O[0], Up[(k + 1) % nv][1] - O[1])
        out.append((e1[0] * e1[0] + e1[1] * e1[1], e1[0] * e2[0] + e1[1] * e2[1], e2[0] * e2[0] + e2[1] * e2[1], e1[0] * e2[1] - e1[1] * e2[0]))
    return out


def _assemble(geo, N, dirichlet):
    R = _ref(N); nv = len(geo); G, nd = dofs(nv, N, dirichlet); K = np.zeros((nd, nd)); M = np.zeros((nd, nd))
    Saa, Sab, Sbb, M0 = R['Saa'][2], R['Sab'][2], R['Sbb'][2], R['M0'][2]
    for k, (n11, n12, n22, det) in enumerate(geo):
        act = G[k] >= 0; idx = G[k][act]; ix = np.ix_(idx, idx); ia = np.ix_(act, act)
        K[ix] += ((n22 * Saa - n12 * (Sab + Sab.T) + n11 * Sbb) / det)[ia]; M[ix] += (det * M0)[ia]
    return K, M, G


def ritz(geo, N):
    """float Rayleigh-Ritz: the integer coefficients (times 2^-CB) of u, as a list over the global numbering, or None"""
    K, M, G = _assemble(geo, N, True)
    try:
        c = sla.eigh(K, M, subset_by_index=[0, 0], check_finite=False)[1][:, 0]
    except Exception:
        return None
    sc = np.abs(c).max()
    if not (sc > 0) or not np.isfinite(sc):
        return None
    if c[0] < 0:
        c = -c
    return [int(v) for v in np.rint(c / sc * 2.0 ** CB)]


def flux(geo, N, n, ci):
    """float least squares: chi (continuous, degree n) minimising |w|^2, w = -(p - O) h + curl chi.
    Returns (integer coefficients, exponent xe): chi = sum x_m 2^xe B_m; or None"""
    nv = len(geo); H = _refh(N, n); Kc, _, Gc = _assemble(geo, n, False); G, nd = dofs(nv, N, True)
    E = H['E'][2]; hs = np.array(H['Hn'], float) / H['L']; X = {k: v[2] for k, v in H['X'].items()}
    cf = np.array(ci, float) * 2.0 ** -CB; r = np.zeros(Kc.shape[0]); loc = len(G[0])
    for k, (n11, n12, n22, det) in enumerate(geo):
        cT = np.zeros(loc); act = G[k] >= 0; cT[act] = cf[G[k][act]]; h = hs * (E @ cT)
        r[Gc[k]] += -(-n12 * (h @ X[('a', 1)]) - n22 * (h @ X[('b', 1)]) + n11 * (h @ X[('a', 2)]) + n12 * (h @ X[('b', 2)]))
    try:
        x = np.zeros(Kc.shape[0]); x[1:] = sla.solve(Kc[1:, 1:], -r[1:], assume_a='pos', check_finite=False)
    except Exception:
        return None
    sc = np.abs(x).max()
    if not np.isfinite(sc):
        return None
    if sc == 0:
        return [0] * len(x), 0
    e = math.frexp(sc)[1]
    return [int(v) for v in np.rint(x * 2.0 ** (CB - e))], e - CB


def _col(v):
    return fmpz_mat(len(v), 1, [int(t) for t in v])


def _qf(x, A, y):
    """x^T A y for integer columns and a packed matrix (fmpz_mat, den, float): an exact fmpq"""
    return fmpq(int((x.transpose() * A[0] * y)[0, 0]), A[1])


def _fq(v):
    return fmpq(v.numerator, v.denominator)


def tri_forms(nv, N, n, ci, xi, xe):
    """the geometry-free part of the proof: for every triangle k of the fan the 14 exact numbers (fmpq)
        m;  aa, ab, bb;  qaa, qab, qbb;  xa1, xb1, xa2, xb2;  raa, rab, rbb
    of u = 2^-CB sum ci B (degree N) and chi = 2^xe sum xi B (degree n), such that for ANY fan geometry
    (n11, n12, n22, det) = (|e1|^2, e1.e2, |e2|^2, e1 x e2), e1 = v_k - O, e2 = v_k+1 - O, det > 0:
        int_T u^2        = det m
        int_T |grad u|^2 = (n22 aa - 2 n12 ab + n11 bb) / det
        int_T |w|^2      = det (n11 qaa + 2 n12 qab + n22 qbb) + 2 (n12 xa1 + n22 xb1 - n11 xa2 - n12 xb2)
                           + (n22 raa - 2 n12 rab + n11 rbb) / det,        w = -(p - O)(2 + E)^-1 u + curl chi.
    None if the coefficient lists do not have the lengths of dofs(nv, N, True) and dofs(nv, n, False)."""
    R = _ref(N); Rn = _ref(n); H = _refh(N, n); G, nd = dofs(nv, N, True); Gc, nc = dofs(nv, n, False)
    if len(ci) != nd or len(xi) != nc:
        return None
    su = fmpq(1, 2 ** CB); sx = fmpq(2 ** xe) if xe >= 0 else fmpq(1, 2 ** -xe); sh = su / H["L"]; Hn = H['Hn']; X = H['X']; out = []
    for k in range(nv):
        cT = _col([ci[g] if g >= 0 else 0 for g in G[k]]); xT = _col([xi[g] for g in Gc[k]])
        Ec = H['E'][0] * cT; h = _col([int(Ec[m, 0]) * Hn[m] for m in range(len(Hn))])                    # h = sh * (this column)
        out.append(dict(m=su * su * _qf(cT, R['M0'], cT), aa=su * su * _qf(cT, R['Saa'], cT), ab=su * su * _qf(cT, R['Sab'], cT), bb=su * su * _qf(cT, R['Sbb'], cT),
                        qaa=sh * sh * _qf(h, H['Qaa'], h), qab=sh * sh * _qf(h, H['Qab'], h), qbb=sh * sh * _qf(h, H['Qbb'], h),
                        xa1=sh * sx * _qf(h, X[('a', 1)], xT), xb1=sh * sx * _qf(h, X[('b', 1)], xT), xa2=sh * sx * _qf(h, X[('a', 2)], xT), xb2=sh * sx * _qf(h, X[('b', 2)], xT),
                        raa=sx * sx * _qf(xT, Rn['Saa'], xT), rab=sx * sx * _qf(xT, Rn['Sab'], xT), rbb=sx * sx * _qf(xT, Rn['Sbb'], xT)))
    return out


def tri_bilinear(nv, N, n, CI, XI, xe):
    """the forms of tri_forms as bilinear forms of several coefficient lists (for trial functions that depend on
    parameters, lgbox.py): CI = [ci_0, ..., ci_q] for u_a = 2^-CB sum ci_a B, XI = [xi_0, ..., xi_q] for
    chi_a = 2^xe sum xi_a B.  Returns for every triangle a dict name -> (q+1) x (q+1) table F with F[a][b] = the form of
    (u_a, u_b), (h_a, h_b), (h_a, chi_b) or (chi_a, chi_b) (fmpq; NOT symmetrised).  For u = sum d_a u_a,
    chi = sum d_a chi_a every quantity of tri_forms is sum_ab d_a d_b F[a][b]."""
    R = _ref(N); Rn = _ref(n); H = _refh(N, n); G, nd = dofs(nv, N, True); Gc, nc = dofs(nv, n, False); q = len(CI)
    if len(XI) != q or any(len(c) != nd for c in CI) or any(len(x) != nc for x in XI):
        return None
    su = fmpq(1, 2 ** CB); sx = fmpq(2 ** xe) if xe >= 0 else fmpq(1, 2 ** -xe); sh = su / H["L"]; Hn = H['Hn']; X = H['X']; out = []

    def table(L, A, Rt, sc):
        P = L.transpose() * A[0] * Rt
        return [[sc * fmpq(int(P[a, b]), A[1]) for b in range(q)] for a in range(q)]
    for k in range(nv):
        Cm = fmpz_mat([[int(CI[a][g]) if g >= 0 else 0 for a in range(q)] for g in G[k]]); Xm = fmpz_mat([[int(XI[a][g]) for a in range(q)] for g in Gc[k]])
        Ec = H['E'][0] * Cm; Hm = fmpz_mat([[int(Ec[m, a]) * Hn[m] for a in range(q)] for m in range(len(Hn))])
        out.append(dict(m=table(Cm, R['M0'], Cm, su * su), aa=table(Cm, R['Saa'], Cm, su * su), ab=table(Cm, R['Sab'], Cm, su * su), bb=table(Cm, R['Sbb'], Cm, su * su),
                        qaa=table(Hm, H['Qaa'], Hm, sh * sh), qab=table(Hm, H['Qab'], Hm, sh * sh), qbb=table(Hm, H['Qbb'], Hm, sh * sh),
                        xa1=table(Hm, X[('a', 1)], Xm, sh * sx), xb1=table(Hm, X[('b', 1)], Xm, sh * sx), xa2=table(Hm, X[('a', 2)], Xm, sh * sx), xb2=table(Hm, X[('b', 2)], Xm, sh * sx),
                        raa=table(Xm, Rn['Saa'], Xm, sx * sx), rab=table(Xm, Rn['Sab'], Xm, sx * sx), rbb=table(Xm, Rn['Sbb'], Xm, sx * sx)))
    return out


def exact_stage(Up, O, N, n, ci, xi, xe, detail=False):
    """the proof: exact A0, A1, A2 for u = 2^-CB sum ci B, chi = 2^xe sum xi B; returns (lower bound of lambda_1(P+), nu,
    Rayleigh quotient of u) as fmpq [and A0, A1, A2 if detail], or None if a hypothesis fails"""
    nv = len(Up); tf = tri_forms(nv, N, n, ci, xi, xe)
    if tf is None:
        return None
    geo = _geom(Up, O); A0 = fmpq(0); A1 = fmpq(0); A2 = fmpq(0)
    for (n11, n12, n22, det), t in zip(geo, tf):
        if not (det > 0):
            return None                                             # O must be strictly inside P+
        n11, n12, n22, det = _fq(n11), _fq(n12), _fq(n22), _fq(det)
        A1 += det * t['m']
        A0 += (n22 * t['aa'] - 2 * n12 * t['ab'] + n11 * t['bb']) / det
        A2 += det * (n11 * t['qaa'] + 2 * n12 * t['qab'] + n22 * t['qbb']) + 2 * (n12 * t['xa1'] + n22 * t['xb1'] - n11 * t['xa2'] - n12 * t['xb2']) + (n22 * t['raa'] - 2 * n12 * t['rab'] + n11 * t['rbb']) / det
    if not (A1 > 0):
        return None
    best = fmpq(0)                                                  # nu <= lambda_2(P+): bounding rectangles ...
    dirs = [(Fraction(1), Fraction(0))] + [(Up[(k + 1) % nv][0] - Up[k][0], Up[(k + 1) % nv][1] - Up[k][1]) for k in range(nv)]
    for dx, dy in dirs:
        a1 = [dx * p[0] + dy * p[1] for p in Up]; a2 = [-dy * p[0] + dx * p[1] for p in Up]; e1 = max(a1) - min(a1); e2 = max(a2) - min(a2); dd = dx * dx + dy * dy
        v = _fq(dd * min(4 / (e1 * e1) + 1 / (e2 * e2), 1 / (e1 * e1) + 4 / (e2 * e2)))
        if v > best:
            best = v
    nu = PI2_LO * best
    ar = sum(Up[i][0] * Up[(i + 1) % nv][1] - Up[(i + 1) % nv][0] * Up[i][1] for i in range(nv)) / 2
    ks = KS_LO / _fq(ar)                                            # ... or Krahn-Szego: lambda_2 >= 2 pi j_{0,1}^2 / |P+|
    if ks > nu:
        nu = ks
    D = A0 - nu * A1
    if not (D < 0):
        return None
    Gq = A0 - 2 * nu * A1 + nu * nu * A2
    return (nu * (-D) / (Gq - D), nu, A0 / A1) + ((A0, A1, A2) if detail else ())


def _fgeo(Up, O):
    return [tuple(float(v) for v in t) for t in _geom(Up, O)]


def F_lower(V, N=4, extra=1, keep_flux=False, grid=None):
    """Rigorous lower bound of F(P) = lambda_1(P) |P| for a strictly convex counter-clockwise polygon with exact binary64
    vertices.  Returns the record (dict; F_lo is a float rounded down) or None.
    grid: P+ on the grid 2^-grid with the inflation 2^-(grid - 3) (default GRID / EPS0; the near region uses 40)."""
    grid = GRID if grid is None else grid
    fr = inflate(V, grid, EPS0 if grid == GRID else 2.0 ** -(grid - 3))
    if fr is None:
        return None
    c, s, Up, O, area = fr; n = N + extra; geo = _fgeo(Up, O)
    ci = ritz(geo, N)
    if ci is None:
        return None
    fl = flux(geo, N, n, ci)
    if fl is None:
        return None
    ex = exact_stage(Up, O, N, n, ci, fl[0], fl[1])
    if ex is None:
        return None
    lam, nu, rho = ex; Flo = lam * _fq(area) / _fq(s * s)
    v = float(Fraction(int(Flo.p), int(Flo.q))) * (1 - (0.0 if keep_flux else SLACK)); v = math.nextafter(v, -math.inf)
    rec = dict(kind='fan', F_lo=v, N=N, n=n, grid=grid, cb=CB, c=[str(c[0]), str(c[1])], s=str(s), Uplus=[[str(a), str(b)] for a, b in Up],
               O=[str(O[0]), str(O[1])], u=ci, nu=float(Fraction(int(nu.p), int(nu.q))), rho_plus=float(Fraction(int(rho.p), int(rho.q))))
    if keep_flux:
        rec['x'] = fl[0]; rec['xe'] = fl[1]
    return rec


def recheck(V, rec):
    """re-proves a stored record (same P+, O and u; the stored flux if there is one, otherwise a recomputed one).
    True if the proved bound is >= rec['F_lo']."""
    try:
        if rec.get('kind') != 'fan' or rec['cb'] != CB:
            return False
        Vf = [(Fraction(float(a)), Fraction(float(b))) for a, b in V]; nv = len(Vf)
        for i in range(nv):
            p, q, r = Vf[i], Vf[(i + 1) % nv], Vf[(i + 2) % nv]
            if (q[0] - p[0]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[0] - q[0]) <= 0:
                return False
        c = (Fraction(rec['c'][0]), Fraction(rec['c'][1])); s = Fraction(rec['s'])
        if not (s > 0):
            return False
        U = [((a - c[0]) / s, (b - c[1]) / s) for a, b in Vf]; Up = [(Fraction(a), Fraction(b)) for a, b in rec['Uplus']]
        O = (Fraction(rec['O'][0]), Fraction(rec['O'][1]))
        if len(Up) != nv or not _admissible(U, Up, O):
            return False
        N = int(rec['N']); n = int(rec['n']); ci = [int(v) for v in rec['u']]
        if 'x' in rec:
            xi, xe = [int(v) for v in rec['x']], int(rec['xe'])
        else:
            fl = flux(_fgeo(Up, O), N, n, ci)
            if fl is None:
                return False
            xi, xe = fl
        ex = exact_stage(Up, O, N, n, ci, xi, xe)
        if ex is None:
            return False
        area = sum(Vf[i][0] * Vf[(i + 1) % nv][1] - Vf[(i + 1) % nv][0] * Vf[i][1] for i in range(nv)) / 2
        Flo = ex[0] * _fq(area) / _fq(s * s)
        return Fraction(int(Flo.p), int(Flo.q)) >= Fraction(float(rec['F_lo']))
    except Exception:
        return False
