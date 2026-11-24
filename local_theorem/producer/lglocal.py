#!/usr/bin/env python3
"""lglocal.py - the local theorem at the regular pentagon (far/LOCAL_PLAN.md; EXPERIMENTAL, 2026-11-17, not reviewed).

Claim to be proved by `certify`: F(P(z)) > F(P(c)) for every chart point z != c with |z_i - c_i| <= rho (c: the chart
point of the exact regular pentagon).  Method: Lemma LG with the trial family  u0 + sum d_i U_i  where u0 is the TRUE
first eigenfunction of the regular pentagon (u0 = computed - e, with a proved bound on the H^1 norm of e) and the flux
family whose member at c is the true flux grad u0 / lambda0.  Then Q(d) = (prod det_k) Psi(c + d) is a polynomial
without constant and linear terms (exactness at c; gradient zero by symmetry), and

    Psi_true(z) = K(z) + sum_k int grad_r e . Lambda_z + sum_k int e_h (a, b) . Theta_z + q_z(e, e)

(residual form, LOCAL_PLAN.md section 3).  K is computed; the coefficients of the linear part are bounded by dual norms
(Gram matrices of the computed functions), q by a majorant.

Parts of this file:
  tables      exact Gram matrices of the 6q basis functions of every fan triangle
              (a h_b, b h_b, d_a chi_b, d_b chi_b, d_a U_b, d_b U_b;  b = 0..q-1 the coefficient lists)
  selftest    float check of the algebra of the residual form with an arbitrary polynomial 'error' (exact identity up to
              the admissibility defect, which is computed), and of the flux identity  int s . grad_r f = sum det int u f
  Local       the rigorous computation (polynomials with rational centre and rational radius; constants as balls)

Reference coordinates on T_k: p = O + a e1 + b e2;  grad_r = (d_a, d_b);  for g = (n11, n12, n22, det):
  M_g = [[n22, -n12], [-n12, n11]] / det,  N_g = M_g^-1 = [[n11, n12], [n12, n22]] / det,
  int |grad v|^2 = int grad_r v . M_g grad_r v,   flux s = (-det a h + d_b chi, -det b h - d_a chi),  int |w|^2 = int s . N_g s.
"""
import os, sys, math, time, itertools
from fractions import Fraction
import numpy as np
from flint import fmpq, fmpz_mat, arb, ctx as fctx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'scripts'))
import lgfan as LF, lgbox as LB, lgpoly as LP, near_lg as NL          # noqa: E402

ND = 6; NV = 5
KSQ = LF.KS_LO; KS = int(KSQ.p) / int(KSQ.q)


def _f(x):
    return int(x.p) / int(x.q)


# ---------------------------------------------------------------------------------------------------------------------
# exact tables
# ---------------------------------------------------------------------------------------------------------------------
def tables(N, n, CI, XI, xe, nv=NV):
    """(tb, GR): tb = lgfan.tri_bilinear (the 14 forms as q x q tables); GR[k] = the exact Gram matrix (6q x 6q list of
    fmpq, reference integrals over a, b >= 0, a + b <= 1, i.e. physical integral / det) of the basis functions, in the
    order  I q + b,  I = 0: a h_b, 1: b h_b, 2: d_a chi_b, 3: d_b chi_b, 4: d_a U_b, 5: d_b U_b."""
    q = len(CI); tb = LF.tri_bilinear(nv, N, n, CI, XI, xe)
    H = LF._refh(N, n); HN = LF._refh(N, N); G_, _ = LF.dofs(nv, N, True); Gc, _ = LF.dofs(nv, n, False)
    su = fmpq(1, 2 ** LF.CB); sx = fmpq(2 ** xe) if xe >= 0 else fmpq(1, 2 ** -xe); sh = su / H['L']; Hn = H['Hn']
    locn = LF._loc(n); locN = LF._loc(N)
    key = ('mixed', n, N)
    if key not in LF._cache:
        def mixed(al, be):
            Dn = [LF._D(e, al) for e in locn]; DN = [LF._D(e, be) for e in locN]
            return LF._pack([[sum(ca * cb * LF._I(ea[0] + eb[0], ea[1] + eb[1], ea[2] + eb[2]) for ca, ea in Dn[m] for cb, eb in DN[k])
                              for k in range(len(locN))] for m in range(len(locn))])
        LF._cache[key] = {(al, be): mixed(al, be) for al in (1, 2) for be in (1, 2)}
    SM = LF._cache[key]

    def table(L, A, Rt, sc):
        P = L.transpose() * A[0] * Rt
        return [[sc * fmpq(int(P[a, b]), A[1]) for b in range(q)] for a in range(q)]
    tr = lambda M: [[M[b][a] for b in range(q)] for a in range(q)]
    GR = []
    for k in range(nv):
        Cm = fmpz_mat([[int(CI[a][g]) if g >= 0 else 0 for a in range(q)] for g in G_[k]]); Xm = fmpz_mat([[int(XI[a][g]) for a in range(q)] for g in Gc[k]])
        Ec = H['E'][0] * Cm; Hm = fmpz_mat([[int(Ec[m, a]) * Hn[m] for a in range(q)] for m in range(len(Hn))]); T = tb[k]
        xu = {(s, be): table(Hm, HN['X'][(s, be)], Cm, sh * su) for s in 'ab' for be in (1, 2)}            # int s h_a d_be U_b
        cu = {(al, be): table(Xm, SM[(al, be)], Cm, sx * su) for al in (1, 2) for be in (1, 2)}            # int d_al chi_a d_be U_b
        B = [[T['qaa'], T['qab'], T['xa1'], T['xa2'], xu[('a', 1)], xu[('a', 2)]],
             [tr(T['qab']), T['qbb'], T['xb1'], T['xb2'], xu[('b', 1)], xu[('b', 2)]],
             [tr(T['xa1']), tr(T['xb1']), T['raa'], T['rab'], cu[(1, 1)], cu[(1, 2)]],
             [tr(T['xa2']), tr(T['xb2']), tr(T['rab']), T['rbb'], cu[(2, 1)], cu[(2, 2)]],
             [tr(xu[('a', 1)]), tr(xu[('b', 1)]), tr(cu[(1, 1)]), tr(cu[(2, 1)]), T['aa'], T['ab']],
             [tr(xu[('a', 2)]), tr(xu[('b', 2)]), tr(cu[(1, 2)]), tr(cu[(2, 2)]), tr(T['ab']), T['bb']]]
        GR.append([[B[I][J][a][b] for J in range(6) for b in range(q)] for I in range(6) for a in range(q)])
    return tb, GR


def geo_float(V, O):
    out = []
    for k in range(len(V)):
        e1 = (V[k][0] - O[0], V[k][1] - O[1]); e2 = (V[(k + 1) % len(V)][0] - O[0], V[(k + 1) % len(V)][1] - O[1])
        out.append((e1[0] ** 2 + e1[1] ** 2, e1[0] * e2[0] + e1[1] * e2[1], e2[0] ** 2 + e2[1] ** 2, e1[0] * e2[1] - e1[1] * e2[0]))
    return out


# ---------------------------------------------------------------------------------------------------------------------
# float check of the algebra (not part of the proof)
# ---------------------------------------------------------------------------------------------------------------------
def pointwise(GR, Mm, gc, Of, F0, lam0, z, w, we):
    """float evaluation at the chart point z of the pieces of the residual form (sums over the triangles, NOT multiplied
    by prod det_k).  GR, Mm: float Gram and mass tables; gc: geometry at c; w: weights of the lists of the family
    ((1, d, 0...)); we: weights of the 'error' list (or None).  Returns a dict."""
    q = len(w); vec = lambda I, x: np.concatenate([x if J == I else np.zeros(q) for J in range(6)])
    w0 = np.zeros(q); w0[0] = 1.0; gz = geo_float(LB.near_verts(z), Of); area = sum(g[3] for g in gz) / 2; lh = F0 / area

    def fields(g, x, l):
        n11, n12, n22, det = g; ah, bh, xa, xb, ua, ub = [vec(I, x) for I in range(6)]
        sa = -det * ah + xb; sb = -det * bh - xa
        return (sa, sb), (ua, ub), ((n11 * sa + n12 * sb) / det - ua / l, (n12 * sa + n22 * sb) / det - ub / l)
    Mof = lambda g: np.array([[g[2], -g[1]], [-g[1], g[0]]]) / g[3]
    Nof = lambda g: np.array([[g[0], g[1]], [g[1], g[2]]]) / g[3]
    ap = lambda A, x: (A[0, 0] * x[0] + A[0, 1] * x[1], A[1, 0] * x[0] + A[1, 1] * x[1])
    out = dict(lhs=0.0, psi=0.0, k2=0.0, k3=0.0, lin=0.0, lth=0.0, quad=0.0, defect=0.0, flux=0.0, area=area, DD=float(np.prod([g[3] for g in gz])))
    for k in range(len(gz)):
        G = GR[k]; ip = lambda x, y: float(x @ G @ y); dot = lambda X, Y: ip(X[0], Y[0]) + ip(X[1], Y[1])
        g = gz[k]; det = g[3]; Mz, Nz, Mc = Mof(g), Nof(g), Mof(gc[k]); ddet = det - gc[k][3]
        s, gv, rho = fields(g, w, lh); _, _, rhoc = fields(gc[k], w0, lam0); mrc = ap(Mc, rhoc)
        bz = lambda u1, u2: det * float(u1 @ Mm[k] @ u2)
        Bz = lambda v1, g1, s1: area * (KS + F0) * bz(v1, v1) - area ** 2 * dot(g1, ap(Mz, g1)) - F0 * KS * dot(s1, ap(Nz, s1))
        out['psi'] += Bz(w, gv, s); out['k2'] += 2 * F0 * KS * dot(rho, mrc); out['k3'] += F0 * KS * dot(mrc, ap(Nz, mrc))
        if we is None:
            continue
        ge = (vec(4, we), vec(5, we)); ye = (vec(0, we), vec(1, we))                               # grad_r e, (a, b) e_h
        tz = tuple(mrc[i] + ap(Mc, ge)[i] / lam0 - ddet * ye[i] for i in range(2))
        st = (s[0] - tz[0], s[1] - tz[1]); gt = (gv[0] - ge[0], gv[1] - ge[1])
        out['lhs'] += Bz(w - we, gt, st)
        Lam = tuple(2 * F0 * ((KS / lam0) * ap(Mc, rho)[i] - area * ap(Mz, rho)[i] - (KS / lam0) * ap(Mc, ap(Nz, mrc))[i]) for i in range(2))
        Th = tuple(-2 * F0 * KS * ddet * (rho[i] - ap(Nz, mrc)[i]) for i in range(2))
        out['lin'] += dot(ge, Lam); out['lth'] += dot(ye, Th)
        tp = tuple(ap(Mc, ge)[i] / lam0 - ddet * ye[i] for i in range(2))
        out['quad'] += area * (KS + F0) * bz(we, we) - area ** 2 * dot(ge, ap(Mz, ge)) - F0 * KS * dot(tp, ap(Nz, tp))
        out['defect'] += dot(gv, tz) - bz(w, we); out['flux'] += dot(s, ge) - bz(w, we)
    return out


def _test_family(N, seed):
    rng = np.random.default_rng(seed); C = list(NL.CF)
    fam = LB.Family(LB.near_verts, C, N, order=1)
    ce = [int(v) for v in rng.integers(-2 ** 20, 2 ** 20, size=len(fam.CI[0]))]                     # an arbitrary H^1_0 function
    return fam, ce, rng


def selftest(N=6, seed=1, verbose=True):
    """Residual form of Psi for an ARBITRARY polynomial 'error' e (an extra coefficient list): with u0 := u~0 - e and the
    formal flux error t_z = M_c (rho_c + grad_r e / lam0) - (det_k(z) - det_k(c)) (a, b) e_h the identity
        B_z(V - E, V - E) = K + int grad_r e . Lambda + int e_h (a,b) . Theta + q(e, e) + 2 kappa |P| defect_z
    holds exactly, defect_z = sum_k int grad_r v~ . t_z - sum_k det_k int v~ e (zero for the true eigenfunction).
    Also checks  sum_k int s~ . grad_r e = sum_k det_k int v~ e.  Returns the largest relative discrepancy."""
    fam, ce, rng = _test_family(N, seed); C = list(NL.CF); q0 = len(fam.CI); q = q0 + 1; n = fam.n
    tb, GRq = tables(N, n, list(fam.CI) + [ce], list(fam.XI) + [[0] * len(fam.XI[0])], fam.xe); GR = [np.array([[_f(x) for x in r] for r in g]) for g in GRq]
    Mm = [np.array([[_f(tb[k]['m'][a][b]) for b in range(q)] for a in range(q)]) for k in range(NV)]
    Of = (float(fam.O[0]), float(fam.O[1])); gc = geo_float(LB.near_verts(C), Of); areac = sum(g[3] for g in gc) / 2
    we = np.zeros(q); we[q0] = 1.0
    T0 = [{nm: _f(tb[k][nm][0][0]) for nm in LB.NAMES} for k in range(NV)]
    A1 = sum(g[3] * t['m'] for g, t in zip(gc, T0)); A0 = sum((g[2] * t['aa'] - 2 * g[1] * t['ab'] + g[0] * t['bb']) / g[3] for g, t in zip(gc, T0))
    lam0 = A0 / A1 * (1 - 1e-3); F0 = lam0 * areac                                                 # any numbers: the identity does not need the true values
    worst = 0.0
    for trial in range(4):
        d = rng.normal(size=ND) * 0.02; z = [c + x for c, x in zip(C, d)]; w = np.concatenate([[1.0], d, [0.0]])
        o = pointwise(GR, Mm, gc, Of, F0, lam0, z, w, we); dt = 2 * KS * o['area'] * o['defect']
        rhs = o['psi'] + o['k2'] - o['k3'] + o['lin'] + o['lth'] + o['quad'] + dt
        sc = abs(o['k2']) + abs(o['k3']) + abs(o['lin']) + abs(o['lth']) + abs(o['quad']) + abs(dt) + 1e-300
        worst = max(worst, abs(o['lhs'] - rhs) / sc, abs(o['flux']) / (abs(float(w @ Mm[0] @ we)) + 1e-300))
        if verbose:
            print('  z = c + %s...: B(V-E) %.6e  Psi~ %.6e | K-corr %.3e  lin %.3e  theta %.3e  quad %.3e  defect term %.3e | lhs - rhs %.2e (rel %.1e), flux identity %.1e'
                  % (np.round(d[:2], 4), o['lhs'], o['psi'], o['k2'] - o['k3'], o['lin'], o['lth'], o['quad'], dt, o['lhs'] - rhs, abs(o['lhs'] - rhs) / sc, o['flux']))
    return worst


# ---------------------------------------------------------------------------------------------------------------------
# rigorous part
# ---------------------------------------------------------------------------------------------------------------------
CX, DV = LP.ctx(ND)
ZERO = CX.from_dict({})
PT = type(DV[0])
BITS = 140                      # constants are rounded to multiples of 2^-BITS (the rounding is added to their radius)
HB = Fraction(1, 2 ** 50)       # |c_i - C_q,i| <= HB (verified in Local.__init__)
HS = Fraction(1, 2 ** 10)       # auxiliary shift for the bound of the Taylor shift by c - C_q
GBITS = 160                     # the certificate coefficients Tm, Em are multiples of 2^-GBITS


def pconst(v):
    v = LP.fq(v)
    return CX.from_dict({(0,) * ND: v}) if v != 0 else ZERO


def pdict(P):
    return {tuple(int(t) for t in m): c for m, c in zip(P.monoms(), P.coeffs())}


def pabs(P):
    return CX.from_dict({m: abs(c) for m, c in pdict(P).items()})


def peval(P, d):
    co = P.coeffs()
    if not co:
        return 0.0
    E = np.array(P.monoms(), dtype=np.int64).reshape(len(co), ND); a = np.array([_f(c) for c in co])
    return float((a * np.prod(np.power(np.array(d, float)[None, :], E), axis=1)).sum())


def _frac(a):
    m, e = a.man_exp(); return Fraction(int(m)) * Fraction(2) ** int(e)


def ball(x):
    """arb -> (mid, rad): Fractions with |x - mid| <= rad for every number x of the ball; mid a multiple of 2^-BITS"""
    m = _frac(x.mid()); r = _frac(x.rad()); mm = Fraction(round(m * 2 ** BITS), 2 ** BITS)
    return mm, r + abs(m - mm)


def up(x):
    """a rational upper bound of the ball x (arb, fmpq or Fraction)"""
    if not isinstance(x, arb):
        return Fraction(int(LP.fq(x).p), int(LP.fq(x).q))
    m, r = ball(x); return m + r


def lo(x):
    m, r = ball(x); return m - r


def A(x):
    """arb from a Fraction / fmpq / int / float (exact)"""
    return x if isinstance(x, arb) else arb(LP.fq(x))


def _ra(a, b):
    return b if a is None else (a if b is None else a + b)


class PB:
    """polynomial ball: every coefficient of the polynomial it stands for lies within r (>= 0, coefficientwise) of m"""
    __slots__ = ('m', 'r')

    def __init__(s, m, r=None):
        s.m = m; s.r = r if (r is not None and not r.is_zero()) else None

    @staticmethod
    def of(x):
        if isinstance(x, PB):
            return x
        if isinstance(x, PT):
            return PB(x)
        if isinstance(x, arb):
            m, r = ball(x); return PB(pconst(m), pconst(r))
        return PB(pconst(x))

    def __add__(s, o):
        o = PB.of(o); return PB(s.m + o.m, _ra(s.r, o.r))
    __radd__ = __add__

    def __sub__(s, o):
        o = PB.of(o); return PB(s.m - o.m, _ra(s.r, o.r))

    def __rsub__(s, o):
        return PB.of(o) - s

    def __neg__(s):
        return PB(ZERO - s.m, s.r)

    def __mul__(s, o):
        o = PB.of(o); r = None
        if o.r is not None:
            r = _ra(r, pabs(s.m) * o.r)
        if s.r is not None:
            r = _ra(r, s.r * pabs(o.m))
            if o.r is not None:
                r = _ra(r, s.r * o.r)
        return PB(s.m * o.m, r)
    __rmul__ = __mul__

    def maj(s):
        return pabs(s.m) + s.r if s.r is not None else pabs(s.m)

    def at_c(s):
        """the value at d' = c - C_q (|d'_i| <= HB < 1) as a ball"""
        m = pdict(s.m); m0 = m.get((0,) * ND, fmpq(0)); rest = sum((abs(v) for k, v in m.items() if any(k)), fmpq(0)); r0 = fmpq(0)
        if s.r is not None:
            for k, v in pdict(s.r).items():
                if any(k):
                    rest += v
                else:
                    r0 += v
        return A(m0) + arb(0, 1) * A(LP.fq(HB) * rest + r0)


class Fld:
    """a function-valued polynomial on one triangle:  sum_I W[I](z) (sum_b d_b phi_{I,b}) + sum_I V[I](z) phi_{I,0}
    (phi_{I,b}: basis function I of the coefficient list b; W, V: PB or None)."""
    __slots__ = ('W', 'V')

    def __init__(s, W=None, V=None):
        s.W = W if W is not None else [None] * 6; s.V = V if V is not None else [None] * 6

    @staticmethod
    def unit(I):
        return Fld(W=[PB.of(1) if J == I else None for J in range(6)])

    def __add__(s, o):
        f = lambda a, b: b if a is None else (a if b is None else a + b)
        return Fld([f(a, b) for a, b in zip(s.W, o.W)], [f(a, b) for a, b in zip(s.V, o.V)])

    def __mul__(s, p):
        p = PB.of(p); return Fld([None if a is None else a * p for a in s.W], [None if a is None else a * p for a in s.V])

    def __sub__(s, o):
        return s + o * (-1)


def _eig2(a, b, d):
    """(smallest, largest) eigenvalue of [[a, b], [b, d]] (balls)"""
    t = (a + d) / 2; h = (a - d) / 2; r = (h * h + b * b).sqrt(); return t - r, t + r          # (no ** : arb power of a ball around 0 is nan)


class Local:
    """the data of the local theorem for the order-1 family of degree N centred at C_q = near_lg.CF.
    test_list: an extra coefficient list (an arbitrary H^1_0 function) for `debug` (not used by the proof)."""

    def __init__(self, N=16, extra=1, verbose=True, test_list=None, tdeg=3, debug_point_F0=False):
        fctx.prec = 400; t0 = time.time(); say = (lambda *a: print(*a, flush=True)) if verbose else (lambda *a: None)
        Cq = list(NL.CF); fam = LB.Family(LB.near_verts, Cq, N, extra=extra, order=1); assert fam.ok
        q = len(fam.CI); n = fam.n; self.fam = fam; self.q = q; self.N = N; self.tdeg = tdeg
        CI = list(fam.CI) + ([test_list] if test_list is not None else []); XI = list(fam.XI) + ([[0] * len(fam.XI[0])] if test_list is not None else [])
        tb, GR = tables(N, n, CI, XI, fam.xe); qq = len(CI); self.qq = qq; self.tb = tb; self.GR = GR; ix = lambda I, b: I * qq + b
        say('  family and exact tables %.1fs' % (time.time() - t0))
        # --- the exact regular pentagon and the constants that belong to it (balls)
        r5 = arb(5).sqrt(); c72 = (r5 - 1) / 4; s72 = (10 + 2 * r5).sqrt() / 4; s36 = (10 - 2 * r5).sqrt() / 4
        carb = [1 + c72, s72, arb(1) / 2, s72 + s36, -c72, s72]; self.carb = carb
        assert all(abs(carb[i] - A(Fraction(Cq[i]))) < A(HB) for i in range(ND)), 'C_q is not within 2^-50 of the regular pentagon'
        Oa = (A(Fraction(fam.O[0])), A(Fraction(fam.O[1]))); Vc = [(arb(0), arb(0)), (arb(1), arb(0)), (carb[0], carb[1]), (carb[2], carb[3]), (carb[4], carb[5])]
        gc = []
        for k in range(NV):
            e1 = (Vc[k][0] - Oa[0], Vc[k][1] - Oa[1]); e2 = (Vc[(k + 1) % NV][0] - Oa[0], Vc[(k + 1) % NV][1] - Oa[1])
            gc.append((e1[0] * e1[0] + e1[1] * e1[1], e1[0] * e2[0] + e1[1] * e2[1], e2[0] * e2[0] + e2[1] * e2[1], e1[0] * e2[1] - e1[1] * e2[0]))
        assert all(g[3] > 0 for g in gc)
        areac = sum((g[3] for g in gc[1:]), gc[0][3]) / 2; KSa = A(KSQ); self.gc = gc; self.areac = areac
        T0 = [{nm: A(tb[k][nm][0][0]) for nm in LB.NAMES} for k in range(NV)]
        A1 = sum(g[3] * t['m'] for g, t in zip(gc, T0)); A0 = sum((g[2] * t['aa'] - 2 * g[1] * t['ab'] + g[0] * t['bb']) / g[3] for g, t in zip(gc, T0))
        A2 = sum(g[3] * (g[0] * t['qaa'] + 2 * g[1] * t['qab'] + g[2] * t['qbb']) + 2 * (g[1] * (t['xa1'] - t['xb2']) + g[2] * t['xb1'] - g[0] * t['xa2'])
                 + (g[2] * t['raa'] - 2 * g[1] * t['rab'] + g[0] * t['rbb']) / g[3] for g, t in zip(gc, T0))
        nu = KSa / areac; assert A0 < nu * A1, 'Lemma LG not applicable at the centre'
        lam_up = A0 / A1; lam_lo = (nu * A1 - A0) / (nu * A2 - A1); Flo = lo(lam_lo * areac); Fup = up(lam_up * areac); assert Flo < Fup
        F0 = A((Flo + Fup) / 2) + arb(0, 1) * A((Fup - Flo) / 2)
        if debug_point_F0:                                         # DEBUG ONLY (comparison with the pointwise formulas): F0 := the midpoint, no radius; not a proof
            F0 = A((Flo + Fup) / 2)
        lam0 = F0 / areac                 # balls that contain the true F0 = F(regular pentagon), lambda0
        e1sq = A(up((A0 - lam_lo * A1) / (1 - lam_up / nu))); e0sq = A(up(e1sq / nu)); eps1 = A(up(e1sq.sqrt())); eps0 = A(up(e0sq.sqrt()))
        self.F0 = F0; self.Flo = Flo; self.Fup = Fup; self.eps1 = eps1; self.eps0 = eps0; self.lam0 = lam0
        self.wit = dict(A0=A0, A1=A1, A2=A2, lam_lo=lam_lo, lam_up=lam_up, nu=nu, area_c=areac, eps1_sq=e1sq, eps0_sq=e0sq); self.rhoc = []; self.CI = CI; self.XI = XI
        say('  F(regular pentagon) in [%.10f, %.10f]; |grad e| <= %.3e (relative %.2e), |e| <= %.3e' % (float(Flo), float(Fup), float(up(eps1)), float(up(eps1 / A0.sqrt())), float(up(eps0))))
        # --- polynomials in d' = z - C_q
        Vp = LP.near_aff([pconst(Fraction(Cq[i])) + DV[i] for i in range(ND)], pconst); Op = (pconst(Fraction(fam.O[0])), pconst(Fraction(fam.O[1])))
        dd = [pconst(1)] + [DV[i] for i in range(ND)]; geoP = []
        for k in range(NV):
            e1 = (Vp[k][0] - Op[0], Vp[k][1] - Op[1]); e2 = (Vp[(k + 1) % NV][0] - Op[0], Vp[(k + 1) % NV][1] - Op[1])
            geoP.append((e1[0] * e1[0] + e1[1] * e1[1], e1[0] * e2[0] + e1[1] * e2[1], e2[0] * e2[0] + e2[1] * e2[1], e1[0] * e2[1] - e1[1] * e2[0]))
        DDp = geoP[0][3]
        for k in range(1, NV):
            DDp = DDp * geoP[k][3]
        Dh = []
        for k in range(NV):
            pk = pconst(1)
            for j in range(NV):
                if j != k:
                    pk = pk * geoP[j][3]
            Dh.append(pk)
        areaP = sum((g[3] for g in geoP[1:]), geoP[0][3]) * fmpq(1, 2); self.geoP = geoP; self.DDp = DDp; self.areaP = areaP
        # the family is  U~_0 + sum_i (z_i - c_i) U_i:  in the variable d' = z - C_q the weight of the list i is d'_i - eta_i,
        # eta = c - C_q, |eta_i| <= HB: a polynomial ball with centre d'_i and radius HB  (2026-11-17)
        wts = [PB.of(1)] + [PB(DV[i], pconst(HB)) for i in range(ND)]

        def form(F):
            out = PB.of(0)
            for a_ in range(q):
                row = PB.of(0)
                for b_ in range(q):
                    if F[a_][b_] != 0:
                        row = row + wts[b_] * F[a_][b_]
                out = out + wts[a_] * row
            return out
        Ms = []; Ns = []; Qs = []; Rs = []
        for k in range(NV):
            n11, n12, n22, det = geoP[k]; t = {nm: form(tb[k][nm]) for nm in LB.NAMES}
            Ms.append(t['m'] * det); Ns.append(t['aa'] * n22 - t['ab'] * (2 * n12) + t['bb'] * n11)
            Qs.append((t['qaa'] * n11 + t['qab'] * (2 * n12) + t['qbb'] * n22) * det + ((t['xa1'] - t['xb2']) * n12 + t['xb1'] * n22 - t['xa2'] * n11) * 2)
            Rs.append(t['raa'] * n22 - t['rab'] * (2 * n12) + t['rbb'] * n11)
        over = lambda L_: sum((L_[k] * Dh[k] for k in range(1, NV)), L_[0] * Dh[0])
        A1p = sum(Ms[1:], Ms[0]); A0D = over(Ns); A2D = sum(Qs[1:], Qs[0]) * DDp + over(Rs)
        PD = A1p * DDp * KSQ - A0D * areaP; P1 = PD * areaP; P2 = A2D * KSQ - A1p * DDp * areaP          # as lgpoly.build, with the weights above
        (P1x, P2x), PDx, _ = LP.build(fam, LP.near_aff, Cq, None); assert P1.m == P1x and P2.m == P2x and PD.m == PDx, 'build mismatch'
        self.PD = PD; self._P2m = P2.m; self.scale = _f(pdict(P2.m)[(0,) * ND])
        KD = P1 - PB.of(F0) * P2
        pf = PB(areaP) * PB.of(1 / F0); kl = PB.of(KSa / lam0); lin_ = lambda v: sum((PB.of(v[b]) * wts[b] for b in range(1, q)), PB.of(v[0]))
        ah, bh, xa, xb, ua, ub = [Fld.unit(I) for I in range(6)]
        Mlin = ZERO; Mq = ZERO; self.fields = []; lowdeg = [a for m in range(0, tdeg + 1) for a in itertools.product(range(m + 1), repeat=ND) if sum(a) == m]
        tight = lowdeg; NL2 = {a: arb(0) for a in tight}; NT2 = {a: arb(0) for a in tight}
        from flint import fmpq_mat
        for k in range(NV):
            n11, n12, n22, det = [PB(x) for x in geoP[k]]; g11, g12, g22, gdet = gc[k]; G = GR[k]
            dra = ah * (-(det * n11)) + bh * (-(det * n12)) + xb * n11 + xa * (-n12) + ua * (-(det * pf))                 # det rho_z, component a
            drb = ah * (-(det * n12)) + bh * (-(det * n22)) + xb * n12 + xa * (-n22) + ub * (-(det * pf))                 # component b
            atc = lambda x: arb(0) if x is None else x.at_c(); rca = [atc(dra.W[I]) / gdet for I in range(6)]; rcb = [atc(drb.W[I]) / gdet for I in range(6)]             # rho_c: coefficients of phi_{I,0}
            Mc = [[g22 / gdet, -g12 / gdet], [-g12 / gdet, g11 / gdet]]
            ma = [Mc[0][0] * rca[I] + Mc[0][1] * rcb[I] for I in range(6)]; mb = [Mc[1][0] * rca[I] + Mc[1][1] * rcb[I] for I in range(6)]      # M_c rho_c
            mca = Fld(V=[PB.of(v) for v in ma]); mcb = Fld(V=[PB.of(v) for v in mb])
            self.rhoc.append(dict(rho_a=rca, rho_b=rcb, M_rho_a=ma, M_rho_b=mb))
            dma = ah * (-(det * det)) + xb * det + ua * (-(pf * n22)) + ub * (pf * n12)                                 # det M_z rho_z
            dmb = bh * (-(det * det)) + xa * (-det) + ua * (pf * n12) + ub * (-(pf * n11))
            dna = mca * n11 + mcb * n12; dnb = mca * n12 + mcb * n22                                                    # det N_z M_c rho_c
            t1a = dra * Mc[0][0] + drb * Mc[0][1]; t1b = dra * Mc[1][0] + drb * Mc[1][1]
            t3a = dna * Mc[0][0] + dnb * Mc[0][1]; t3b = dna * Mc[1][0] + dnb * Mc[1][1]
            pre = PB.of(2 * F0) * PB(Dh[k]); arP = PB(areaP)
            La = (t1a * kl - dma * arP - t3a * kl) * pre; Lb = (t1b * kl - dmb * arP - t3b * kl) * pre                  # D Lambda_z
            preT = PB.of(-2 * F0 * KSa) * PB(Dh[k]) * (det - PB.of(gdet))
            Ta = (dra - dna) * preT; Tb = (drb - dnb) * preT                                                            # D Theta_z
            self.fields.append((La, Lb, Ta, Tb))
            # computed part: + 2 F0 K int rho_z . M_c rho_c - F0 K int (M_c rho_c) . N_z (M_c rho_c)     (times D)
            ga = [[sum((A(G[ix(I, b)][ix(J, 0)]) * ma[J] for J in range(6)), arb(0)) for b in range(q)] for I in range(6)]
            gb = [[sum((A(G[ix(I, b)][ix(J, 0)]) * mb[J] for J in range(6)), arb(0)) for b in range(q)] for I in range(6)]
            t2 = PB.of(0)
            for I in range(6):
                if dra.W[I] is not None:
                    t2 = t2 + dra.W[I] * lin_(ga[I])
                if drb.W[I] is not None:
                    t2 = t2 + drb.W[I] * lin_(gb[I])
            ipc = lambda x, y: sum((x[I] * A(G[ix(I, 0)][ix(J, 0)]) * y[J] for I in range(6) for J in range(6)), arb(0))
            KD = KD + PB.of(2 * F0 * KSa) * PB(Dh[k]) * t2 - PB.of(F0 * KSa) * PB(Dh[k]) * (n11 * PB.of(ipc(ma, ma)) + n12 * PB.of(2 * ipc(ma, mb)) + n22 * PB.of(ipc(mb, mb)))
            # --- linear part: tight dual norms for the coefficients of degree 2..tdeg
            Gm = fmpq_mat(6 * q, 6 * q, [G[ix(I, a)][ix(J, b)] for I in range(6) for a in range(q) for J in range(6) for b in range(q)])
            sg = [[A(G[ix(I, b)][ix(I, b)]).sqrt() for b in range(q)] for I in range(6)]
            lN, _ = _eig2(g11 / gdet, g12 / gdet, g22 / gdet); lmaxN = _eig2(g11 / gdet, g12 / gdet, g22 / gdet)[1]
            lminM, lmaxM = _eig2(g22 / gdet, -g12 / gdet, g11 / gdet)

            def cvecs(X):
                """for every tight alpha: (mid vector as 1 x 6q fmpq_mat, sum_j rad_j sqrt(G_jj))"""
                W = [({} if X.W[I] is None else pdict(X.W[I].m), {} if X.W[I] is None or X.W[I].r is None else pdict(X.W[I].r)) for I in range(6)]
                Vv = [({} if X.V[I] is None else pdict(X.V[I].m), {} if X.V[I] is None or X.V[I].r is None else pdict(X.V[I].r)) for I in range(6)]
                out = {}
                for a in tight:
                    mid = []; rad = arb(0)
                    for I in range(6):
                        for b in range(q):
                            if b == 0:
                                m_ = W[I][0].get(a, fmpq(0)) + Vv[I][0].get(a, fmpq(0)); r_ = W[I][1].get(a, fmpq(0)) + Vv[I][1].get(a, fmpq(0))
                            else:                                   # coefficient of d'^a in W (d'_b - eta_b):  W_{a - e_b} - eta_b W_a
                                r_ = LP.fq(HB) * (abs(W[I][0].get(a, fmpq(0))) + W[I][1].get(a, fmpq(0))); m_ = fmpq(0)
                                if a[b - 1] >= 1:
                                    a2 = tuple(a[i] - (1 if i == b - 1 else 0) for i in range(ND)); m_ = W[I][0].get(a2, fmpq(0)); r_ += W[I][1].get(a2, fmpq(0))
                            mid.append(m_)
                            if r_ != 0:
                                rad += A(r_) * sg[I][b]
                    out[a] = (fmpq_mat(1, 6 * q, mid), rad)
                return out
            vLa, vLb, vTa, vTb = cvecs(La), cvecs(Lb), cvecs(Ta), cvecs(Tb)
            for a in tight:
                xa_, ra_ = vLa[a]; xb_, rb_ = vLb[a]; qf = lambda x, y: (x * Gm * y.transpose())[0, 0]
                v = (g11 * A(qf(xa_, xa_)) + 2 * g12 * A(qf(xa_, xb_)) + g22 * A(qf(xb_, xb_))) / gdet
                nl = A(max(up(v), 0)).sqrt() + lmaxN.sqrt() * (ra_ + rb_); NL2[a] += nl * nl
                ya_, sa_ = vTa[a]; yb_, sb_ = vTb[a]
                nt = (A(up(A(qf(ya_, ya_)) + A(qf(yb_, yb_)))).sqrt() + sa_ + sb_) / gdet.sqrt(); NT2[a] += nt * nt
            # --- linear part: crude majorant for all degrees
            S = [sum((pconst(up(sg[I][b])) * (dd[b] + pconst(HB)) for b in range(1, q)), pconst(up(sg[I][0]))) for I in range(6)]      # majorant of sum_b w_b |phi_{I,b}|

            def majf(X):
                out = ZERO
                for I in range(6):
                    if X.W[I] is not None:
                        out = out + X.W[I].maj() * S[I]
                    if X.V[I] is not None:
                        out = out + X.V[I].maj() * LP.fq(up(sg[I][0]))
                return out
            Mlin = Mlin + (majf(La) + majf(Lb)) * LP.fq(up(eps1 * lmaxN.sqrt())) + (majf(Ta) + majf(Tb)) * LP.fq(up(eps0 / (2 * gdet.sqrt())))
            # --- quadratic part q_z(e, e): majorant
            mk = up(e0sq / gdet); Gk = up(e1sq / lminM); Gam = up(lmaxM * e1sq / (lam0 * lam0)); Yk = mk / 4; cr = up((A(Gam) * A(Yk)).sqrt())
            a11, a12, a22, adet = [pabs(x) for x in geoP[k]]; aDh = pabs(Dh[k]); aA = pabs(areaP); ddm = (det - PB.of(gdet)).maj(); s3 = a11 + a22 + a12
            Mq = Mq + pabs(DDp) * aA * adet * LP.fq(up(KSa + F0) * mk) + aA * aA * aDh * s3 * LP.fq(Gk) \
                + aDh * LP.fq(up(F0 * KSa)) * (s3 * LP.fq(Gam) + ddm * (s3 + a12) * LP.fq(2 * cr) + ddm * ddm * s3 * LP.fq(Yk))
            say('  triangle %d done %.1fs' % (k, time.time() - t0))
        self.KD = KD; self.Mlin = Mlin; self.Mq = Mq
        self.Elin = {a: up(eps1 * NL2[a].sqrt() + eps0 * NT2[a].sqrt() / 2) for a in tight}
        self.NLs = {a: NL2[a].sqrt() for a in tight}; self.NTs = {a: NT2[a].sqrt() for a in tight}      # the two dual norms (balls), for tests
        # --- Taylor shift from C_q to c:  |s_alpha - t_alpha| <= (HB / HS) (B(e + HS) - B(e))_alpha,  B >= |t| coefficientwise
        B = KD.maj() + Mlin + Mq; Bs = B
        for i in range(ND):
            Bs = LP.shift(Bs, i, HS)
        slack = (Bs - B) * LP.fq(HB / HS)
        km = pdict(KD.m); kr = pdict(KD.r) if KD.r is not None else {}; ml = pdict(Mlin); mq = pdict(Mq); sl = pdict(slack)
        fr = lambda x: Fraction(int(x.p), int(x.q))
        self.Tm = {a: fr(v) for a, v in km.items() if sum(a) >= 2}
        keys = set(a for a in list(km) + list(kr) + list(ml) + list(mq) + list(sl) if sum(a) >= 2)
        self.Em = {a: fr(kr.get(a, fmpq(0))) + fr(mq.get(a, fmpq(0))) + fr(sl.get(a, fmpq(0))) + (self.Elin[a] if a in self.Elin else fr(ml.get(a, fmpq(0)))) for a in keys}
        # consistency check (not part of the proof): the true coefficients of degree 0 and 1 vanish (at c), so the computed ones must
        # lie within the error bounds:  |KD.m_alpha| <= radius + linear bound + quadratic bound + shift
        self.low = {a: (fr(km.get(a, fmpq(0))), fr(kr.get(a, fmpq(0))) + self.Elin[a] + fr(mq.get(a, fmpq(0))) + fr(sl.get(a, fmpq(0)))) for a in tight if sum(a) <= 1}
        self.parts = dict(rad={a: fr(v) for a, v in kr.items()}, quad={a: fr(v) for a, v in mq.items()}, shift={a: fr(v) for a, v in sl.items()}, crude={a: fr(v) for a, v in ml.items()})
        # the certificate numbers: centres rounded to multiples of 2^-GBITS, the rounding added to the (rounded up) error bounds
        U_ = 2 ** GBITS; Tr = {}; Er = {}
        for a in keys:
            v = self.Tm.get(a, Fraction(0)); vr = Fraction(round(v * U_), U_); er = self.Em[a] + abs(v - vr)
            Tr[a] = vr; Er[a] = Fraction(-((-er.numerator * U_) // er.denominator), U_)               # ceiling
        self.Tm_exact = self.Tm; self.Em_exact = self.Em; self.Tm = {a: v for a, v in Tr.items() if v != 0}; self.Em = Er
        say('  assembled %.1fs: %d coefficients, %d error bounds' % (time.time() - t0, len(self.Tm), len(self.Em)))

    # -----------------------------------------------------------------------------------------------------------------
    def low_check(self, verbose=True):
        """True if 0 lies in the enclosure of every coefficient of degree 0 and 1 (a necessary condition for the error
        bounds to be right; the largest ratio |computed| / bound is reported)"""
        worst = max((abs(v) / b if b > 0 else Fraction(10 ** 9)) for v, b in self.low.values()); sc = Fraction(self.scale)
        if verbose:
            print('  degree 0: computed %.3e, bound %.3e;  degree 1: largest |computed| %.3e, bounds about %.3e  (units of F);  largest ratio %.3f'
                  % (float(self.low[(0,) * ND][0] / sc), float(self.low[(0,) * ND][1] / sc), max(abs(float(v / sc)) for a, (v, b) in self.low.items() if sum(a) == 1),
                     max(float(b / sc) for a, (v, b) in self.low.items() if sum(a) == 1), float(worst)))
        return worst <= 1

    def quadratic(self):
        """(H, T3): the quadratic form (matrix H, Q2 = d^T H d / 2) and the symmetric cubic tensor of the computed part"""
        H = [[Fraction(0)] * ND for _ in range(ND)]; T3 = [[[Fraction(0)] * ND for _ in range(ND)] for _ in range(ND)]
        for a, v in self.Tm.items():
            idx = [i for i in range(ND) for _ in range(a[i])]
            if len(idx) == 2:
                if idx[0] == idx[1]:
                    H[idx[0]][idx[0]] = 2 * v
                else:
                    H[idx[0]][idx[1]] = H[idx[1]][idx[0]] = v
            elif len(idx) == 3:
                perms = set(itertools.permutations(idx))
                for p_ in perms:
                    T3[p_[0]][p_[1]][p_[2]] += v / len(perms)
        return H, T3

    def allowance(self, rho):
        """tau with |Q_true(e) - e^T (H/2 + sum e_i T_i) e| <= tau |e|^2 for |e_i| <= rho: (total, parts)"""
        rho = Fraction(rho); w = lambda a: Fraction(1) if max(a) == sum(a) else Fraction(1, 2)
        a2 = max(self.Em.get(tuple(2 if j == i else 0 for j in range(ND)), 0)
                 + sum(self.Em.get(tuple(1 if j in (i, l) else 0 for j in range(ND)), 0) for l in range(ND) if l != i) / 2 for i in range(ND))
        a3 = rho * sum(v * w(a) for a, v in self.Em.items() if sum(a) == 3)
        a4 = sum(rho ** (sum(a) - 2) * (abs(self.Tm.get(a, 0)) + self.Em.get(a, 0)) * w(a) for a in set(self.Tm) | set(self.Em) if sum(a) >= 4)
        return a2 + a3 + a4, (a2, a3, a4)

    def certify_inner(self, rho, verbose=True):
        """True if  e^T (H/2 + sum_i e_i T_i - tau) e > 0  is PROVED for all 0 != e, |e_i| <= rho (exact LDL^T at the 64
        vertices), i.e. Q_true > 0 on the punctured box around c."""
        rho = Fraction(rho); H, T3 = self.quadratic(); tau, parts = self.allowance(rho); sc = Fraction(self.scale); worst = None
        for sgn in itertools.product((-1, 1), repeat=ND):
            M = [[H[i][j] / 2 + rho * sum(sgn[l] * T3[l][i][j] for l in range(ND)) - (tau if i == j else 0) for j in range(ND)] for i in range(ND)]
            Lm = [[Fraction(0)] * ND for _ in range(ND)]; Dg = [None] * ND
            for j in range(ND):
                dj = M[j][j] - sum(Lm[j][l] * Lm[j][l] * Dg[l] for l in range(j))
                if dj <= 0:
                    if verbose:
                        print('  inner box %s: NOT positive definite at vertex %s (pivot %d = %.4f); tau = %.4f (errors of degree 2: %.4f, degree 3: %.4f, degree >= 4: %.4f), in units of F'
                              % (float(rho), sgn, j, float(dj / sc), float(tau / sc), float(parts[0] / sc), float(parts[1] / sc), float(parts[2] / sc)))
                    return False
                Dg[j] = dj
                for i in range(j + 1, ND):
                    Lm[i][j] = (M[i][j] - sum(Lm[i][l] * Lm[j][l] * Dg[l] for l in range(j))) / dj
            mn = min(Dg); worst = mn if worst is None else min(worst, mn)
        if verbose:
            print('  inner box %s: positive definite at all 64 vertices (smallest pivot %.4f); tau = %.4f (errors of degree 2: %.4f, degree 3: %.4f, degree >= 4: %.4f), in units of F'
                  % (float(rho), float(worst / sc), float(tau / sc), float(parts[0] / sc), float(parts[1] / sc), float(parts[2] / sc)))
        return True

    # -----------------------------------------------------------------------------------------------------------------
    def side_conditions(self, W, verbose=True):
        """the hypotheses of Lemma LG for the TRUE family at every z with |z_i - c_i| <= W:
        (a) the fan determinants are positive (exact polynomial test on the cube of half-width W + HB around C_q);
        (b) kappa b_z(v) > |P(z)| a_z(v) for v = v~(d) - e, from
              kappa b - |P| a >= (kappa b~ - |P| a~) - 2 kappa sqrt(b~ b_e) - |P| (2 sqrt(a~ a_e) + a_e),
            b_e <= beta eps0^2 (det_k(z) <= beta det_k(c)),  a_e <= theta eps1^2 (M_z <= theta M_c), beta = theta = 2,
            upper bounds of a~, b~, |P| by ball arithmetic on the cube, and the polynomial test
              PD~ - rhs D > 0 on the cube  (PD~ = D (kappa b~ - |P| a~), D = prod det_k).
        Returns True / False."""
        Wq = Fraction(W) + HB; r = [float(Wq)] * ND; assert Fraction(r[0]) >= Wq or Fraction(math.nextafter(r[0], 1)) >= Wq
        r = [math.nextafter(float(Wq), 1.0)] * ND; Cq = list(NL.CF); fam = self.fam
        for k in range(NV):
            ok, info = LP.poly_positive(self.geoP[k][3], r)
            if not ok:
                if verbose:
                    print('  side conditions: det %d not proved positive' % k)
                return False
        zb = [arb(Cq[i], r[i]) for i in range(ND)]; Oa = (A(Fraction(fam.O[0])), A(Fraction(fam.O[1])))
        area, A0, A1, A2, dets = LB.conditions(LB.near_verts(zb), Oa, fam.forms([zb[i] - A(Fraction(fam.zc[i])) for i in range(ND)], 0))
        Vz = LB.near_verts(zb); two = arb(2)
        for k in range(NV):
            e1 = (Vz[k][0] - Oa[0], Vz[k][1] - Oa[1]); e2 = (Vz[(k + 1) % NV][0] - Oa[0], Vz[(k + 1) % NV][1] - Oa[1])
            n11 = e1[0] * e1[0] + e1[1] * e1[1]; n12 = e1[0] * e2[0] + e1[1] * e2[1]; n22 = e2[0] * e2[0] + e2[1] * e2[1]; det = e1[0] * e2[1] - e1[1] * e2[0]
            g11, g12, g22, gdet = self.gc[k]
            if not (det > 0 and two * gdet - det > 0):
                return False
            # theta M_c - M_z positive definite:  [[x, y], [y, w]] with x, w > 0 and x w - y^2 > 0
            x = two * g22 / gdet - n22 / det; w = two * g11 / gdet - n11 / det; y = -two * g12 / gdet + n12 / det
            if not (x > 0 and w > 0 and x * w - y * y > 0):
                if verbose:
                    print('  side conditions: M_z <= 2 M_c not proved on triangle %d' % k)
                return False
        KSa = A(KSQ); Amax = A(up(A0)); Bmax = A(up(A1)); Pmax = A(up(area))
        rhs = up(2 * KSa * (Bmax * two).sqrt() * self.eps0 + Pmax * (2 * (Amax * two).sqrt() * self.eps1 + two * self.eps1 * self.eps1))
        rsup = sum((v * LP.fq(Wq) ** sum(k_) for k_, v in pdict(self.PD.r).items()), fmpq(0)) if self.PD.r is not None else fmpq(0)     # sup of the radius of PD on the cube
        ok, info = LP.poly_positive(self.PD.m - self.DDp * LP.fq(rhs) - pconst(rsup), r)
        if verbose:
            print('  side conditions on the cube of half-width %s: determinants positive, M_z <= 2 M_c, det_z <= 2 det_c; Lemma LG applicable to the true family: %s (allowance %.3e against %.3e at the centre)'
                  % (float(W), ok, float(rhs), _f(pdict(self.PD.m)[(0,) * ND]) / _f(pdict(self.DDp)[(0,) * ND])))
        return ok

    def symmetry_check(self, verbose=True):
        """the ten relabelling maps of the chart fix c and the mean of their differentials at c has norm < 1 (it is 0):
        no non-zero vector is fixed by all of them, hence grad F(c) = 0.  Ball arithmetic, forward differentiation."""
        z0 = arb(0); D = LB.Dual; cst = lambda v: D(A(v), [z0] * ND)
        z = [D(self.carb[i], [arb(1) if j == i else z0 for j in range(ND)]) for i in range(ND)]
        V = [(cst(0), cst(0)), (cst(1), cst(0)), (z[0], z[1]), (z[2], z[3]), (z[4], z[5])]
        M = [[z0] * ND for _ in range(ND)]; fixed = True
        for refl in (0, 1):
            for j in range(NV):
                Wv = [V[(j - i) % NV] if refl else V[(j + i) % NV] for i in range(NV)]
                bx = Wv[1][0] - Wv[0][0]; by = Wv[1][1] - Wv[0][1]; nn = bx * bx + by * by; out = []
                for i in (2, 3, 4):
                    px = Wv[i][0] - Wv[0][0]; py = Wv[i][1] - Wv[0][1]
                    re = (px * bx + py * by) / nn; im = (py * bx - px * by) / nn                   # (p - W0) / (W1 - W0)
                    out += [re, -im if refl else im]
                for a in range(ND):
                    fixed = fixed and bool(abs(out[a].v - self.carb[a]) < A(Fraction(1, 2 ** 300)))
                    for b in range(ND):
                        M[a][b] = M[a][b] + out[a].g[b] / 10
        nrm = max(up(sum((abs(M[a][b]) for b in range(ND)), z0)) for a in range(ND))
        if verbose:
            print('  symmetry: the ten relabelling maps fix c: %s; norm of the mean differential <= %.2e' % (fixed, float(nrm)))
        return fixed and nrm < 1

    def certify_shell(self, rho_in, R, W, keep=6, alpha=0.7, aspect=0.3, max_sec=0, log=None, verbose=True):
        """Q_true(e) > 0 PROVED for every e with |e_i| <= W, e_2 >= 0, |e|_2 < R, that is not in the inner cube
        |e_i| <= rho_in  (W, rho_in dyadic; R rational; the half cube is bisected).  A leaf box m +- r is
          I  inside the inner cube (|m_i| + r_i <= rho_in for all i),
          O  outside the ball (sum_i max(|m_i| - r_i, 0)^2 >= R^2),
          P  proved: with p(x) = sum_{2 <= |alpha| <= keep} Tm_alpha (m + x)^alpha  and
             err = sum_alpha E'_alpha prod_i (|m_i| + r_i)^alpha_i   (E' = Em plus |Tm_alpha| for |alpha| > keep),
             p(x) - err > 0 for all |x_i| <= r_i, by lgpoly.poly_positive with exact_tail.
        EXACT rational arithmetic in every decision (2026-11-17, after GPT PENT-21); floating point only chooses the
        multipliers and the bisection.  log: a file; one JSON line per leaf, in the order of the depth-first bisection:
        {"t": "I" | "O" | "P", "m": [...], "r": [...]} (fractions as strings) and for P the witness "w": "dual" with "mu"
        (six binary64 numbers as fractions) or "crude".  Returns (ok, stats)."""
        import json
        t0 = time.time(); rho_in = Fraction(rho_in); R2 = Fraction(R) ** 2; W = Fraction(W); F_ = Fraction
        P0 = CX.from_dict({a: LP.fq(v) for a, v in self.Tm.items() if sum(a) <= keep})
        Ek = dict(self.Em)
        for a, v in self.Tm.items():
            if sum(a) > keep:
                Ek[a] = Ek.get(a, 0) + abs(v)
        EP = CX.from_dict({a: LP.fq(v) for a, v in Ek.items() if v != 0})                      # non-negative coefficients

        def esup(m, r):
            return EP(*[LP.fq(abs(F_(m[i])) + F_(r[i])) for i in range(ND)])                      # exact value at the far corner
        fs = lambda v: [str(F_(x)) for x in v]

        def out(t, m, r, **kw):
            if log is not None:
                log.write(json.dumps(dict(t=t, m=fs(m), r=fs(r), **kw)) + '\n')
        m0 = [0.0] * ND; r0 = [float(W)] * ND; m0[2] = float(W / 2); r0[2] = float(W / 2)
        stack = [(LP.shift(P0, 2, m0[2]), m0, r0)]; st = dict(inner=0, out=0, proved=0, tests=0, failed=0, minrad=1.0, crude=0); ok_all = True
        while stack:
            if max_sec and time.time() - t0 > max_sec:
                ok_all = False; st['open'] = len(stack); break
            P, m, r = stack.pop()
            if all(abs(F_(m[i])) + F_(r[i]) <= rho_in for i in range(ND)):
                st['inner'] += 1; out('I', m, r); continue
            d2 = sum(max(abs(F_(m[i])) - F_(r[i]), 0) ** 2 for i in range(ND))
            if d2 >= R2:
                st['out'] += 1; out('O', m, r); continue
            i = None
            if d2 > 0 and max(r) <= alpha * math.sqrt(float(d2)):
                st['tests'] += 1; ok, info = LP.poly_positive(P - CX.from_dict({(0,) * ND: esup(m, r)}), r, exact_tail=True)
                if ok:
                    st['proved'] += 1; st['minrad'] = min(st['minrad'], min(r))
                    if info['used'] == 'dual':
                        out('P', m, r, w='dual', mu=fs(info['mu']))
                    else:
                        st['crude'] += 1; out('P', m, r, w='crude')
                    continue
                st['failed'] += 1
                if 'lin' in info:
                    cand = [j for j in range(ND) if r[j] >= aspect * max(r)]; i = max(cand, key=lambda j: info['lin'][j] + info['hi'][j])
            if i is None:
                i = max(range(ND), key=lambda j: r[j])
            if r[i] < 2.0 ** -30:
                ok_all = False; st['stuck'] = st.get('stuck', 0) + 1; out('U', m, r); continue
            h = r[i] / 2; r2 = list(r); r2[i] = h
            for sg in (h, -h):
                m2 = list(m); m2[i] = m[i] + sg; stack.append((LP.shift(P, i, sg), m2, r2))
            if verbose and (st['tests'] % 2000 == 0) and st['tests']:
                print('    %d tests, %d proved, stack %d, %.0fs' % (st['tests'], st['proved'], len(stack), time.time() - t0), flush=True)
        st['sec'] = round(time.time() - t0, 1)
        return ok_all, st

    # -----------------------------------------------------------------------------------------------------------------
    def export(self, out, rho_in):
        """write the certificate data for independent checkers into the directory out (format: far/LOCAL_CERT.md):
        witness.json (Tm, Em: everything the final positivity proof uses), inner.json (the 64 vertex matrices' pivots),
        inputs.json (the coefficient lists of the family), gram.json.xz (exact Gram matrices), constants.json, lowdeg.json
        (intermediate numbers of the coefficients of degree <= 3), rho_c.json."""
        import json, lzma
        os.makedirs(out, exist_ok=True); fs = lambda v: str(Fraction(v)); key = lambda a: ','.join(str(x) for x in a)
        bl = lambda x: dict(mid=str(ball(x)[0]), rad=str(ball(x)[1]), approx=float(x.mid())) if isinstance(x, arb) else dict(mid=fs(x), rad='0', approx=float(x))
        json.dump(dict(format='pent-local-cert/1', rho_in=fs(rho_in), gbits=GBITS, Tm={key(a): fs(v) for a, v in sorted(self.Tm.items())},
                       Em={key(a): fs(v) for a, v in sorted(self.Em.items())}), open(os.path.join(out, 'witness.json'), 'w'))
        rho = Fraction(rho_in); H, T3 = self.quadratic(); tau, parts = self.allowance(rho); piv = []
        for sgn in itertools.product((-1, 1), repeat=ND):
            M = [[H[i][j] / 2 + rho * sum(sgn[l] * T3[l][i][j] for l in range(ND)) - (tau if i == j else 0) for j in range(ND)] for i in range(ND)]
            Lm = [[Fraction(0)] * ND for _ in range(ND)]; Dg = []
            for j in range(ND):
                dj = M[j][j] - sum(Lm[j][l] * Lm[j][l] * Dg[l] for l in range(j)); Dg.append(dj)
                if dj <= 0:
                    break
                for i in range(j + 1, ND):
                    Lm[i][j] = (M[i][j] - sum(Lm[i][l] * Lm[j][l] * Dg[l] for l in range(j))) / dj
            piv.append(dict(vertex=list(sgn), pivots=[fs(d) for d in Dg], min_float=float(min(Dg))))
        json.dump(dict(rho_in=fs(rho), H=[[fs(v) for v in r] for r in H], T3=[[[fs(v) for v in r] for r in t] for t in T3], tau=fs(tau), tau2=fs(parts[0]), tau3=fs(parts[1]),
                       tau_ge4=fs(parts[2]), scale=self.scale, tau_F_units=float(tau / Fraction(self.scale)), min_pivot=fs(min(Fraction(d) for v in piv for d in v['pivots'])),
                       min_pivot_F_units=min(v['min_float'] for v in piv) / self.scale, vertices=piv), open(os.path.join(out, 'inner.json'), 'w'))
        fam = self.fam
        json.dump(dict(N=self.N, n=fam.n, xe=fam.xe, CB=LF.CB, O=[fs(fam.O[0]), fs(fam.O[1])], Cq=[fs(v) for v in NL.CF], Cq_hex=[float(v).hex() for v in NL.CF],
                       note='u_b = 2^-CB sum CI[b] B (degree N, Dirichlet numbering of lgfan.dofs), chi_b = 2^xe sum XI[b] B (degree n); b = 0: centre, b = 1..6: first-order lists',
                       CI=[[int(v) for v in c] for c in self.CI[:self.q]], XI=[[int(v) for v in c] for c in self.XI[:self.q]]), open(os.path.join(out, 'inputs.json'), 'w'))
        q = self.q; qq = self.qq
        with lzma.open(os.path.join(out, 'gram.json.xz'), 'wt') as f:
            json.dump(dict(order='index I * 7 + b; I = 0: a h_b, 1: b h_b, 2: d_a chi_b, 3: d_b chi_b, 4: d_a U_b, 5: d_b U_b; reference integrals over a, b >= 0, a + b <= 1',
                           gram=[[[str(Fraction(int(g[I * qq + a][J * qq + b].p), int(g[I * qq + a][J * qq + b].q))) for J in range(6) for b in range(q)] for I in range(6) for a in range(q)] for g in self.GR],
                           mass=[[[str(Fraction(int(self.tb[k]['m'][a][b].p), int(self.tb[k]['m'][a][b].q))) for b in range(q)] for a in range(q)] for k in range(NV)]), f)
        w = self.wit
        json.dump(dict(KS=fs(Fraction(int(KSQ.p), int(KSQ.q))), HB=fs(HB), HS=fs(HS), F_lo=fs(self.Flo), F_up=fs(self.Fup), eps1=fs(up(self.eps1)), eps0=fs(up(self.eps0)),
                       forms_at_c={k: bl(v) for k, v in w.items()}, geometry_at_c=[dict(n11=bl(g[0]), n12=bl(g[1]), n22=bl(g[2]), det=bl(g[3])) for g in self.gc],
                       c=[bl(v) for v in self.carb]), open(os.path.join(out, 'constants.json'), 'w'))
        P = self.parts
        json.dump({key(a): dict(K_mid=fs(self.Tm_exact.get(a, 0)) if sum(a) >= 2 else fs(self.low[a][0]), K_rad=fs(P['rad'].get(a, 0)), NL=bl(self.NLs[a]), NT=bl(self.NTs[a]), E_lin=fs(self.Elin[a]),
                            E_quad=fs(P['quad'].get(a, 0)), E_shift=fs(P['shift'].get(a, 0)), E_lin_crude=fs(P['crude'].get(a, 0))) for a in sorted(self.Elin)}, open(os.path.join(out, 'lowdeg.json'), 'w'))
        json.dump([{k: [bl(v) for v in r[k]] for k in r} for r in self.rhoc], open(os.path.join(out, 'rho_c.json'), 'w'))
        # every part of every error bound (2026-11-18, for an independent level-B check of all coefficients; diagnostic:
        # the proof uses Em only).  Integers in units of 2^-gbits, each part rounded up.
        U_ = 2 ** GBITS; cu = lambda v: int(-((-Fraction(v).numerator * U_) // Fraction(v).denominator))
        with lzma.open(os.path.join(out, 'parts.json.xz'), 'wt') as f:
            json.dump(dict(format='pent-local-parts/1', gbits=GBITS,
                           note='for every exponent a of degree >= 2: [K_rad, E_lin_used, E_quad, E_shift, tight], the four parts as integers in units of 2^-gbits, each rounded up; '
                                'E_lin_used = eps1 NL + eps0 NT / 2 if tight = 1 (degree <= 3), else the coefficient of the crude majorant. '
                                'With S = the sum of the four integers and E = Em 2^gbits of witness.json:  S - 3 <= E <= S + 1  '
                                '(Em = ceiling of the exact sum plus the rounding |K_mid - Tm| <= 2^-gbits / 2 of the centre).',
                           parts={key(a): [cu(P['rad'].get(a, 0)), cu(self.Elin[a]) if a in self.Elin else cu(P['crude'].get(a, 0)), cu(P['quad'].get(a, 0)), cu(P['shift'].get(a, 0)),
                                           1 if a in self.Elin else 0] for a in sorted(self.Em_exact)}), f)

    # -----------------------------------------------------------------------------------------------------------------
    def debug(self, seed=3, size=2e-3):
        """float comparison of the polynomials with the pointwise formulas (needs test_list): the computed part KD and
        the linear functional applied to the test function.  Returns the largest relative difference."""
        assert self.qq == self.q + 1, 'debug needs a test list'
        rng = np.random.default_rng(seed); q = self.q; qq = self.qq; fam = self.fam; Cq = list(NL.CF); ix = lambda I, b: I * qq + b
        GR = [np.array([[_f(x) for x in r] for r in g]) for g in self.GR]; Mm = [np.array([[_f(self.tb[k]['m'][a][b]) for b in range(qq)] for a in range(qq)]) for k in range(NV)]
        Of = (float(fam.O[0]), float(fam.O[1])); gcf = [tuple(float(x.mid()) for x in g) for g in self.gc]; F0 = float(self.F0.mid()); lam0 = float(self.lam0.mid())
        we = np.zeros(qq); we[q] = 1.0; worst = 0.0
        for trial in range(3):
            d = rng.normal(size=ND) * size; z = [c + x for c, x in zip(Cq, d)]; w = np.concatenate([[1.0], d, [0.0]])
            o = pointwise(GR, Mm, gcf, Of, F0, lam0, z, w, we); kf = o['DD'] * (o['psi'] + o['k2'] - o['k3']); kp = peval(self.KD.m, d)
            lf = o['DD'] * (o['lin'] + o['lth']); lp = 0.0
            for k in range(NV):
                La, Lb, Ta, Tb = self.fields[k]; G = GR[k]
                for X, I2 in ((La, 4), (Lb, 5), (Ta, 0), (Tb, 1)):                                  # pair with d_a e, d_b e, a e_h, b e_h
                    for I in range(6):
                        if X.W[I] is not None:
                            lp += peval(X.W[I].m, d) * sum(w[b] * G[ix(I, b), ix(I2, q)] for b in range(q))
                        if X.V[I] is not None:
                            lp += peval(X.V[I].m, d) * G[ix(I, 0), ix(I2, q)]
            # scale of the linear part for the comparison: the same pairing with absolute values is not available; use |lf| + tiny
            worst = max(worst, abs(kf - kp) / (abs(self.scale) * 1e-6), abs(lf - lp) / (abs(lf) + 1e-300))
            print('  d = %s...: K D pointwise %.9e polynomial %.9e (difference / scale %.1e) | linear part pointwise %.9e polynomial %.9e (rel %.1e)'
                  % (np.round(d[:2], 5), kf, kp, (kf - kp) / self.scale, lf, lp, abs(lf - lp) / (abs(lf) + 1e-300)))
        return worst


if __name__ == '__main__':
    if len(sys.argv) > 1 and sys.argv[1] == 'selftest':
        N = int(sys.argv[2]) if len(sys.argv) > 2 else 6
        w = selftest(N); print('selftest N = %d: largest relative discrepancy %.2e  %s' % (N, w, 'PASS' if w < 1e-7 else 'FAIL'))
    elif len(sys.argv) > 1 and sys.argv[1] == 'debug':
        N = int(sys.argv[2]) if len(sys.argv) > 2 else 6
        fam, ce, rng = _test_family(N, 1); L = Local(N, test_list=ce); w = L.debug(); print('debug N = %d: %.2e' % (N, w))
    elif len(sys.argv) > 1 and sys.argv[1] == 'run':
        # python3 lglocal.py run OUT [N]: the whole local theorem  |z - c|_2 < 1/20 + 2^-20 (strictly more than the ball left out by near_lg.py)
        import json, hashlib
        out = sys.argv[2]; N = int(sys.argv[3]) if len(sys.argv) > 3 else 16; os.makedirs(out, exist_ok=True); t0 = time.time()
        RIN = Fraction(1, 128); RR = Fraction(52429, 2 ** 20); WW = Fraction(13, 256); assert RR > Fraction(1, 20) + Fraction(1, 2 ** 30) and WW > RR
        L = Local(N); res = dict(N=N, n=L.fam.n, rho_in=str(RIN), R=str(RR), W=str(WW), F_lo_float=float(L.Flo), F_up_float=float(L.Fup),
                                 eps1=float(up(L.eps1)), eps0=float(up(L.eps0)), KS=str(Fraction(int(KSQ.p), int(KSQ.q))), HB=str(HB))
        res['low_check'] = bool(L.low_check()); res['symmetry'] = bool(L.symmetry_check()); res['inner'] = bool(L.certify_inner(RIN)); res['side'] = bool(L.side_conditions(WW))
        cert = os.path.join(out, 'cert'); L.export(cert, RIN)
        with open(os.path.join(cert, 'leaves.jsonl'), 'w') as fl:
            fl.write(json.dumps(dict(t='root', R=str(RR), W=str(WW), rho_in=str(RIN), keep=6)) + '\n'); ok, st = L.certify_shell(RIN, RR, WW, log=fl)
        res['shell'] = bool(ok); res['shell_stats'] = st
        res['cert_sha256'] = {f: hashlib.sha256(open(os.path.join(cert, f), 'rb').read()).hexdigest() for f in sorted(os.listdir(cert))}
        res['ok'] = bool(res['symmetry'] and res['inner'] and res['side'] and res['shell']); res['sec'] = round(time.time() - t0, 1)
        res['claim'] = 'F(P(z)) > F(P(c)) for every chart point z != c with |z - c|_2 < R (c: the regular pentagon), if ok'
        json.dump(res, open(os.path.join(out, 'local.json'), 'w'), indent=1); print(json.dumps(res, indent=1)); print('LOCAL THEOREM: %s' % ('PROVED (computation complete)' if res['ok'] else 'NOT PROVED'))
    elif len(sys.argv) > 1 and sys.argv[1] == 'shell':
        N = int(sys.argv[2]); rin = Fraction(sys.argv[3]); R = Fraction(sys.argv[4]); W = Fraction(sys.argv[5]); L = Local(N)
        print('symmetry', L.symmetry_check()); print('inner', L.certify_inner(rin)); print('side', L.side_conditions(W))
        ok, st = L.certify_shell(rin, R, W, max_sec=float(sys.argv[6]) if len(sys.argv) > 6 else 0); print('shell', ok, st)
    elif len(sys.argv) > 1 and sys.argv[1] == 'inner':
        N = int(sys.argv[2]); L = Local(N); H, T3 = L.quadratic()
        ev = np.linalg.eigvalsh(np.array([[float(x) for x in r] for r in H]) / L.scale); print('Hessian eigenvalues (F units)', np.round(ev, 4))
        for rho in sys.argv[3].split(','):
            L.certify_inner(Fraction(rho))
