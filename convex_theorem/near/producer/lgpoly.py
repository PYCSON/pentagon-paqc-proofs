#!/usr/bin/env python3
"""LG box certificate in exact polynomial form, PROTOTYPE (Claude, 2026-11-16; not used by any pipeline, not reviewed).

Setting (lgbox.py, far/LG_LEMMA.md): a fan (O; v_0..v_4) whose vertices are AFFINE functions of parameters z, a box
|z_i - mid_i| <= rho_i, and a trial family u(z) = u_0 + sum (z_i - mid_i) u_i, chi(z) likewise (integer Bernstein
coefficient lists).  For every z with all fan determinants det_k(z) > 0 the pair (u(z), w(z)) is admissible on the
polygon P(z), and with nu = KS / |P| <= lambda_2 (Krahn-Szego) Lemma LG gives

    Dc := KS A1 - |P| A0 > 0   and   C := |P| Dc - F* (KS A2 - |P| A1) > 0     =>     F(P(z)) = lambda_1 |P| > F*.

A0, A1, A2 are rational functions of d = z - mid with denominators det_k (lgfan.tri_forms / tri_bilinear), so
    PD := Dc * D   and   PC := C * D,     D = det_0 ... det_4,
are POLYNOMIALS in d with rational coefficients; they are expanded exactly (python-flint fmpq_mpoly).  Then
    det_k > 0 on the box (k = 0..4),   PD > 0 on the box,   PC > 0 on the box
prove F > F* for every z of the box.  Positivity of a polynomial p on the box from its exact coefficients:
    p(d) >= min over the box of (p_0 + g.d + d^T H d / 2)  -  sum_{|alpha| >= 3} |p_alpha| rho^alpha,
and the minimum of the quadratic part is bounded below by Lagrangian duality: for any mu >= 0 with A = H + 2 diag(mu)
positive definite,  p_0 + g.d + d^T H d / 2 >= p_0 - g^T A^-1 g / 2 - sum mu_i rho_i^2  on the box (exact rational
LDL^T; mu from floating point).  No interval arithmetic: every number in the proof is an exact rational, except that the
sum over the higher terms is evaluated in binary64 with an explicit allowance for the rounding errors (see poly_positive)."""
import sys, os, math, time
from fractions import Fraction
import numpy as np
from flint import fmpq, fmpz, fmpq_mpoly_ctx
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import lgfan as LF
import lgbox as LB

_ctx = {}


def ctx(nd):
    if nd not in _ctx:
        c = fmpq_mpoly_ctx.get(tuple('d%d' % i for i in range(nd))); _ctx[nd] = (c, c.gens())
    return _ctx[nd]


def fq(v):
    if isinstance(v, fmpq):
        return v
    v = Fraction(v); return fmpq(v.numerator, v.denominator)


def build(fam, vert_aff, mid, fstar):
    """(PC, PD, [det_k]) as polynomials in d = z - mid.  fam: lgbox.Family with centre fam.zc; vert_aff(x): vertex list for
    a list x of polynomials (the vertex map must be affine in x: the caller's responsibility, e.g. lgbox.near_verts)."""
    nd = len(mid); C, DV = ctx(nd)
    const = lambda v: C.from_dict({(0,) * nd: fq(v)})
    V = vert_aff([const(mid[i]) + DV[i] for i in range(nd)], const); O = (const(fam.O[0]), const(fam.O[1])); nv = len(V)
    tb = LF.tri_bilinear(nv, fam.N, fam.n, fam.CI, fam.XI, fam.xe); q = len(fam.CI)
    dd = [const(1)] + [const(Fraction(mid[i]) - Fraction(fam.zc[i])) + DV[i] for i in range(nd)]      # (1, z - zc)
    zero = C.from_dict({})

    def form(F):
        if q == 1:
            return const(F[0][0])
        s = zero
        for a in range(q):
            row = zero
            for b in range(q):
                if F[a][b] != 0:
                    row = row + F[a][b] * dd[b]
            s = s + dd[a] * row
        return s
    dets = []; N = []; Q = []; R = []; M = []
    for k in range(nv):
        e1 = (V[k][0] - O[0], V[k][1] - O[1]); e2 = (V[(k + 1) % nv][0] - O[0], V[(k + 1) % nv][1] - O[1]); t = {nm: form(tb[k][nm]) for nm in LB.NAMES}
        n11 = e1[0] * e1[0] + e1[1] * e1[1]; n12 = e1[0] * e2[0] + e1[1] * e2[1]; n22 = e2[0] * e2[0] + e2[1] * e2[1]; det = e1[0] * e2[1] - e1[1] * e2[0]
        dets.append(det); M.append(det * t['m']); N.append(n22 * t['aa'] - 2 * n12 * t['ab'] + n11 * t['bb'])
        Q.append(det * (n11 * t['qaa'] + 2 * n12 * t['qab'] + n22 * t['qbb']) + 2 * (n12 * (t['xa1'] - t['xb2']) + n22 * t['xb1'] - n11 * t['xa2']))
        R.append(n22 * t['raa'] - 2 * n12 * t['rab'] + n11 * t['rbb'])
    D = dets[0]
    for k in range(1, nv):
        D = D * dets[k]

    def over(L):                                                    # (sum_k L_k / det_k) * D, a polynomial
        s = zero
        for k in range(nv):
            p = L[k]
            for j in range(nv):
                if j != k:
                    p = p * dets[j]
            s = s + p
        return s
    area = sum(dets[1:], dets[0]) * fmpq(1, 2); A1 = sum(M[1:], M[0]); A0D = over(N); A2D = sum(Q[1:], Q[0]) * D + over(R)
    PD = LF.KS_LO * A1 * D - area * A0D
    P1 = area * PD; P2 = LF.KS_LO * A2D - area * A1 * D                # C = P1 - fstar * P2 (linear in the target value)
    if fstar is None:
        return (P1, P2), PD, dets
    PC = P1 - fq(fstar) * P2
    return PC, PD, dets


def box_lower(fam, vert_aff, mid, rho, L0, L1, steps=12):
    """A number L in [L0, L1] with F(P(z)) > L PROVED for all |z_i - mid_i| <= rho_i, or None (not even L0 is proved).
    Same certificate as check_box (no bisection of the box): the fan determinants and PD positive on the box, and
    P1 - L P2 > 0 on the box; the condition is linear in L, the largest provable L is found by bisection (float search;
    every accepted L is proved by the exact test).  L0, L1: binary64 numbers; the L returned is a binary64 number."""
    (P1, P2), PD, dets = build(fam, vert_aff, mid, None); r = [float(v) for v in rho]
    for P in list(dets) + [PD]:
        if not poly_positive(P, r)[0]:
            return None
    ok = lambda L: poly_positive(P1 - fq(L) * P2, r)[0]
    if not ok(L0):
        return None
    if ok(L1):
        return float(L1)
    a, b = float(L0), float(L1)
    for _ in range(steps):
        m = (a + b) / 2
        if ok(m):
            a = m
        else:
            b = m
    return a


def _up(x):
    return math.nextafter(x, math.inf)


def _ipow(rf, E):
    """prod_i rf[i]^E[k, i] for every row k, by repeated IEEE multiplication only (no library pow): each of the at most
    (row degree + number of variables) multiplications is correctly rounded"""
    w = np.ones(len(E))
    for i, x in enumerate(rf):
        emax = int(E[:, i].max()) if len(E) else 0; tab = np.ones(emax + 1)
        for e in range(1, emax + 1):
            tab[e] = tab[e - 1] * x
        w = w * tab[E[:, i]]
    return w


def poly_positive(P, rho, exact_tail=False):
    """True if the polynomial P (fmpq_mpoly in d) is PROVED positive on the box |d_i| <= rho_i (rho: positive binary64
    numbers, taken exactly).  Returns (ok, info): info = dict(c0, lin, tail, qp) as floats for diagnostics / splitting,
    and the witness of the quadratic bound: info['used'] = 'dual' (with the multipliers info['mu'], binary64) or 'crude'.
    exact_tail: the sum of the terms of degree >= 3 as an exact rational (no floating point anywhere in the decision)."""
    nd = len(rho); r = [Fraction(v) for v in rho]; rf = [float(v) for v in rho]
    c0 = Fraction(0); g = [Fraction(0)] * nd; H = [[Fraction(0)] * nd for _ in range(nd)]; lin = [0.0] * nd; hi_loss = [0.0] * nd; T = 0.0
    co = P.coeffs(); mo = P.monoms()
    if co:
        E = np.array(mo, dtype=np.int64).reshape(len(co), nd); deg = E.sum(1)
        for k in np.nonzero(deg <= 2)[0]:                           # constant, linear and quadratic terms: exact rationals
            e = [int(t) for t in E[k]]; v = Fraction(int(co[k].p), int(co[k].q)); dk = int(deg[k])
            if dk == 0:
                c0 = v
            elif dk == 1:
                g[e.index(1)] = v
            else:
                idx = [i for i in range(nd) if e[i]]
                if len(idx) == 1:
                    H[idx[0]][idx[0]] = 2 * v
                else:
                    H[idx[0]][idx[1]] = H[idx[1]][idx[0]] = v
        hk = np.nonzero(deg >= 3)[0]
        if len(hk):
            # tail: sum |p_alpha| rho^alpha in binary64, using only IEEE-754 multiplications and additions (2026-11-17, after
            # GPT PENT-21: no library pow).  Error analysis, n = number of terms < 10^6, m <= 64 multiplications per term:
            #   a_k = fl(fl(|p / q|) (1 + 2^-52)) >= |p_alpha| (1 - 2^-104)        (int / int is correctly rounded);
            #   w_k = _ipow >= rho^alpha (1 - 2^-53)^m                              (no underflow: w_k > 1e-200 is checked);
            #   fl(sum of n non-negative numbers) >= (exact sum) (1 - (n - 1) 2^-53) for ANY order of summation.
            # So the exact tail is at most fl(sum) (1 + (n + m + 2) 2^-53 (1 + 1e-6)) < fl(sum) (1 + 1.2e-10); the factor
            # 1 + 1e-9 and one step up cover it.  The 1e-290 is not needed any more and stays as extra slack.
            a = np.array([abs(int(co[k].p) / int(co[k].q)) for k in hk]) * (1 + 2.0 ** -52); Eh = E[hk]
            w = _ipow(rf, Eh); terms = a * w
            assert len(terms) < 10 ** 6 and int(Eh.sum(1).max()) + nd <= 64
            if not (w.min() > 1e-200 and np.isfinite(terms).all()):
                return False, dict(c0=float(c0), tail=float('inf'), why='range')
            T = _up(float(terms.sum()) * (1 + 1e-9) + 1e-290)
            if exact_tail:
                T = Fraction(0)
                for k in hk:
                    t = Fraction(abs(int(co[k].p)), int(co[k].q))
                    for i in range(nd):
                        if E[k][i]:
                            t *= r[i] ** int(E[k][i])
                    T += t
            hi_loss = [float(v) for v in (terms[:, None] * Eh / deg[hk][:, None]).sum(0)]
    info = dict(c0=float(c0), tail=float(T))
    if not (c0 > 0):
        info['why'] = 'centre'; return False, info
    # quadratic part: Lagrangian bound with float multipliers, exact rational evaluation
    gf = np.array([float(v) for v in g]); Hf = np.array([[float(v) for v in row] for row in H]); rr = np.array(rf)
    L = max(np.abs(Hf).sum(1).max(), 1e-300); d = -np.sign(gf) * rr
    for _ in range(300):
        dn = np.clip(d - (gf + Hf @ d) / L, -rr, rr)
        if np.abs(dn - d).max() <= 1e-13 * rr.max():
            d = dn; break
        d = dn
    gr = gf + Hf @ d; mu = np.zeros(nd)
    for i in range(nd):
        if abs(d[i]) >= rr[i] * (1 - 1e-9):
            mu[i] = max(0.0, -gr[i] / (2 * d[i]))
    lam = np.linalg.eigvalsh(Hf + 2 * np.diag(mu))[0]
    if lam <= 1e-9 * L:
        # indefinite quadratic part: add a uniform shift tau to mu; tau by a float search on the dual value
        # -g^T (H + 2 diag(mu) + 2 tau)^-1 g / 2 - sum (mu_i + tau) rho_i^2   (any tau that makes A positive definite is valid)
        base = (1e-6 * L - lam) / 2; r2 = float((rr * rr).sum()); ev, Q = np.linalg.eigh(Hf + 2 * np.diag(mu)); gq = Q.T @ gf

        def dual(tau):
            return -0.5 * float((gq * gq / (ev + 2 * tau)).sum()) - tau * r2
        best = base; bv = dual(base); tau = base
        for _ in range(60):
            tau *= 1.5; v = dual(tau)
            if v > bv:
                best, bv = tau, v
            elif tau > 100 * best + 1e-300:
                break
        mu = mu + best
    muq = [Fraction(float(m)) for m in mu]
    A = [[H[i][j] + (2 * muq[i] if i == j else 0) for j in range(nd)] for i in range(nd)]
    Lm = [[Fraction(0)] * nd for _ in range(nd)]; Dg = [None] * nd; pd = True
    for j in range(nd):
        dj = A[j][j] - sum(Lm[j][k] * Lm[j][k] * Dg[k] for k in range(j))
        if dj <= 0:
            pd = False; break
        Dg[j] = dj
        for i in range(j + 1, nd):
            Lm[i][j] = (A[i][j] - sum(Lm[i][k] * Lm[j][k] * Dg[k] for k in range(j))) / dj
    lin = [abs(float(g[i])) * rf[i] for i in range(nd)]
    info['lin'] = lin; info['hi'] = hi_loss
    if pd:
        y = [None] * nd
        for i in range(nd):
            y[i] = g[i] - sum(Lm[i][k] * y[k] for k in range(i))
        qp = c0 - sum(y[i] * y[i] / Dg[i] for i in range(nd)) / 2 - sum(muq[i] * r[i] * r[i] for i in range(nd))
    crude = c0 - sum(abs(g[i]) * r[i] for i in range(nd)) - sum(abs(H[i][j]) * r[i] * r[j] for i in range(nd) for j in range(nd)) / 2   # every term by its absolute value
    info['used'] = 'dual'; info['mu'] = [float(m) for m in muq]
    if not pd or crude > qp:
        qp = crude; info['used'] = 'crude'
    info['qp'] = float(qp)
    ok = qp > 0 and qp - Fraction(T) > 0
    if not ok:
        info['why'] = 'quadratic' if not (qp > 0) else 'tail'
    return ok, info


MAX_SUB = int(os.environ.get('LGPOLY_MAX_SUB', '64'))          # at most this many boxes per leaf (the leaf itself and its parts)
SUB_RATIO = float(os.environ.get('LGPOLY_SUB_RATIO', '15'))    # bisect inside a leaf only if the failed test was a tail <= SUB_RATIO * quadratic bound


def shift(P, i, h):
    """the polynomial d -> P(d + h e_i), exactly (h: a binary64 number)"""
    nd = P.context().nvars(); C, DV = ctx(nd)
    args = list(DV); args[i] = DV[i] + C.from_dict({(0,) * nd: fq(h)})
    return P.compose(*args, ctx=C)


def check_box(fam, vert_aff, mid, rho, fstar, max_sub=None):
    """(ok, info): F(P(z)) > fstar PROVED for all |z_i - mid_i| <= rho_i (exact rational arithmetic).
    The polynomials det_k, PD, PC are built once (same trial family); if a positivity test fails on the box, the box is
    bisected (at most max_sub sub-boxes in total) and the failing polynomials are re-expanded at the centres of the
    halves (exact Taylor shift); a polynomial proved positive on a box is not tested again on its parts.
    info: the diagnostics of the first failing test on the whole box (for the caller's heuristics), 'sub' = boxes used."""
    max_sub = MAX_SUB if max_sub is None else max_sub
    PC, PD, dets = build(fam, vert_aff, mid, fstar); nd = len(rho)
    stack = [([('det%d' % k, dk) for k, dk in enumerate(dets)] + [('D', PD), ('C', PC)], [float(v) for v in rho])]; used = 0; first = None; last = None
    while stack:
        pl, r = stack.pop(); used += 1; rest = []; bad = None
        for name, P in pl:
            ok, info = poly_positive(P, r)
            if ok:
                last = dict(info, stage=name)
            else:
                rest.append((name, P)); bad = dict(info, stage=name)
        if first is None:
            first = bad if bad is not None else last
        if not rest:
            continue
        if bad.get('why') == 'centre' or used + len(stack) + 2 > max_sub or 'lin' not in bad:
            return False, dict(first, sub=used)
        if used == 1 and (bad.get('why') != 'tail' or bad['tail'] > SUB_RATIO * bad['qp']):
            return False, dict(first, sub=used)                     # heuristic: bisecting inside the leaf pays only for a moderate tail excess
        i = max((jj for jj in range(nd) if r[jj] > 0), key=lambda jj: bad['lin'][jj] + bad['hi'][jj]); h = r[i] / 2; r2 = list(r); r2[i] = h
        for sg in (h, -h):
            stack.append(([(name, shift(P, i, sg)) for name, P in rest], r2))
    return True, dict(first, sub=used)


def near_aff(x, const):
    """the near chart for polynomial arguments: (0, 0), (1, 0), (x0, x1), (x2, x3), (x4, x5)"""
    return [(const(0), const(0)), (const(1), const(0)), (x[0], x[1]), (x[2], x[3]), (x[4], x[5])]


if __name__ == '__main__':
    FSTAR = 18.919134195733914
    c72, s72 = math.cos(2 * math.pi / 5), math.sin(2 * math.pi / 5); XREG = np.array([1 + c72, s72, 0.5, s72 + math.sin(math.pi / 5), -c72, s72])
    rng = np.random.default_rng(1); N = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    for r in (0.05, 0.1, 0.2, 0.4, 0.6):
        for trial in range(2):
            d = rng.normal(size=6); d /= np.linalg.norm(d); xc = [round(float(v) * 2 ** 16) / 2 ** 16 for v in XREG + r * d]
            for order in (0, 1):
                t = time.time(); fam = LB.Family(LB.near_verts, xc, N, order=order); tf = time.time() - t; out = []
                for k in range(1, 12):
                    rho = 2.0 ** -k; t = time.time(); ok, info = check_box(fam, near_aff, xc, [rho] * 6, FSTAR)
                    out.append('2^-%d:%s/%.2fs' % (k, 'ok' if ok else info.get('stage', '') + '-' + info.get('why', ''), time.time() - t))
                    if ok:
                        break
                print('r %.2f order %d (family %.2f s)  %s   rho / r = %.3f' % (r, order, tf, '  '.join(out[-3:]), rho / r if ok else 0), flush=True)
