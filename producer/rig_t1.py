"""Rigorous T1 with u-dependent matrix Efron-Stein eps, by second-order Taylor branch-and-bound.
Claude, 2026-10-30/31.  Own code; arb ball arithmetic.  Builds on rig_paqc.py (same cell, corners, geometry gate).

For u in the parameter cube [-1,1]^6 (x = mid + half*u, cell = exact binary64 box), with product weights a_c(u):
   J_k(u) = sum_c a_c J_{c,k}           (exact: J multiaffine in raw coordinates, see rig_paqc (3))
   A(u)   = sum_j a0_j J_j(u)          (exact area)
   N_k(u) = sum_c a_c lambda_c J_{c,k}, lambda_c = F_c^L / A_c   (F_c^L rigorous corner lower bounds)
   eps_k(u) = Jmax_k [ (1+tau) lam_max( Gl_k^{-1/2} (sum_i pq_i(u) R_ik R_ik') Gl_k^{-1/2} ) + (1+1/tau) sum_i pq_i(u) dev_ik ]
       pq_i(u) = (1-u_i^2)/4 ; R_ik = mean edge difference of X along i ; dev_ik = max_e lam_max(d_e' Gl^{-1} d_e)
       (Loewner matrix Efron-Stein with the TRUE Bernoulli parameters p_i = (1+u_i)/2, tau > 0 fixed per box)
   T1:  F(P(u)) >= B(u) := (1 - max_k eps_k(u)) * min_k Q_k(u),   Q_k = A(u) N_k(u) / J_k(u)   (needs eps_k < 1, J_k > 0).
   NB (2026-10-31): the corner Rayleigh bound N_c(v) >= lambda_c D_c(v) holds for the SUM over sectors only, so the
   per-sector pairing min_k (1-eps_k) Q_k stated as "Theorem A" in GPT v2.3 T1_SCALE_NEUTRAL.md is not justified; the
   global factor (1 - max_k eps_k), as in GPT [Q1-T1-RAW] and in run2, is what the proof gives.
   Proof: N_P(v) = sum_k N_{P,k} >= sum_k (1-eps_k) sum_c a_c N_{c,k} >= (1-eps_max) sum_c a_c N_c >= (1-eps_max) sum_c a_c lambda_c D_c(v)
          = (1-eps_max) sum_k N_k(u) int_{T_k} v^2,  and D_P(v) = sum_k J_k(u) int_{T_k} v^2;  ratio of positive sums >= min ratio.
Box test on U = prod [c_i - r_i, c_i + r_i]  (v3, 2026-10-31):
   B(u) = min_{k,k'} Q_k(u) (1 - eps_k'(u))   (Q_k > 0; accepted boxes need this > F*_up > 0, which forces 1 - eps_k' > 0),
   and each of the 25 products is bounded by its own second-order Taylor form
      f(u) >= f(c) - sum_i |df/du_i (c)| r_i - 1/2 sum_ij sup_U |d2f/du_i du_j| r_i r_j   (Lagrange remainder),
   centre value/gradient from point jets, Hessian from interval jets over the whole box.  Keeping the PRODUCT matters:
   Q is smallest at the cell corners where eps = 0; bounding min Q and max eps separately stalls the B&B.
   eps_k' majorants per box: the Efron-Stein form with lam_max replaced by a POLYNOMIAL majorant
      lam_max(Mt) = tr/2 + sqrt(q) <= tr/2 + (q + s^2)/(2s)  (AM-GM, constant s > 0),  or by the trace;
   several s and tau are tried and the smallest Taylor upper bound is used (each is a valid majorant on the box).
   (A square-root closed form was used first; at id71 the matrix is nearly isotropic, the radicand ~ 0, and branch
   handling either broke the Taylor remainder (one early version, retracted) or fell back to the 2x-loose trace.)
Multiaffine pieces (J_k, A, N_k): exact ranges of value / gradient / mixed second derivatives on a sub-box from the
sub-box corner values (finite differences), pure second derivatives = 0.
Accept U if the bound > F*_up; otherwise bisect the coordinate with the largest Taylor-loss contribution (heuristic only).
Accepted leaves of the halving tree are an exact partition of [-1,1]^6 and are all saved for replay.
"""
import sys, json, time, heapq, itertools
import numpy as np
from flint import arb
import rig_paqc as R

N = 6
A = R.A


class Jet:
    __slots__ = ('v', 'g', 'h')

    def __init__(self, v, g=None, h=None):
        self.v = v
        self.g = g if g is not None else [arb(0)] * N
        self.h = h if h is not None else [[arb(0)] * N for _ in range(N)]

    def __add__(self, o):
        if not isinstance(o, Jet):
            return Jet(self.v + o, self.g, self.h)
        return Jet(self.v + o.v, [a + b for a, b in zip(self.g, o.g)], [[self.h[i][j] + o.h[i][j] for j in range(N)] for i in range(N)])

    __radd__ = __add__

    def __neg__(self):
        return Jet(-self.v, [-a for a in self.g], [[-x for x in r] for r in self.h])

    def __sub__(self, o):
        return self + (-o if isinstance(o, Jet) else -o)

    def __rsub__(self, o):
        return (-self) + o

    def scale(self, s):
        return Jet(self.v * s, [a * s for a in self.g], [[x * s for x in r] for r in self.h])

    def __mul__(self, o):
        if not isinstance(o, Jet):
            return self.scale(o)
        a, b = self, o
        g = [a.v * b.g[i] + b.v * a.g[i] for i in range(N)]
        h = [[a.v * b.h[i][j] + b.v * a.h[i][j] + a.g[i] * b.g[j] + a.g[j] * b.g[i] for j in range(N)] for i in range(N)]
        return Jet(a.v * b.v, g, h)

    __rmul__ = __mul__

    def recip(self):
        x = self.v
        if not (x > 0 or x < 0):
            raise ZeroDivisionError
        iv = 1 / x; iv2 = iv * iv; iv3 = iv2 * iv
        g = [-self.g[i] * iv2 for i in range(N)]
        h = [[-self.h[i][j] * iv2 + 2 * self.g[i] * self.g[j] * iv3 for j in range(N)] for i in range(N)]
        return Jet(iv, g, h)

    def sqrt(self):
        x = self.v
        if not (x > 0):
            raise ValueError('sqrt of non-positive')
        s = x.sqrt(); is_ = 1 / s
        g = [self.g[i] * is_ / 2 for i in range(N)]
        h = [[self.h[i][j] * is_ / 2 - self.g[i] * self.g[j] * is_ / (4 * x) for j in range(N)] for i in range(N)]
        return Jet(s, g, h)


def alpha_jets(u):
    """product weights a_c(u) = prod_i (1 + s_ci u_i)/2 as jets (u: list of 6 arb, point or interval)."""
    out = []
    f = {(i, s): (1 + s * u[i]) / 2 for i in range(N) for s in (-1, 1)}
    for s in R.CORN:
        fs = [f[(i, int(s[i]))] for i in range(N)]
        v = arb(1)
        for x in fs:
            v = v * x
        g = []
        for i in range(N):
            p = arb(s[i]) / 2
            for j in range(N):
                if j != i:
                    p = p * fs[j]
            g.append(p)
        h = [[arb(0)] * N for _ in range(N)]
        for i in range(N):
            for j in range(i + 1, N):
                p = arb(s[i] * s[j]) / 4
                for l in range(N):
                    if l != i and l != j:
                        p = p * fs[l]
                h[i][j] = p; h[j][i] = p
        out.append(Jet(v, g, h))
    return out


def lin(al, vals):
    """multiaffine function sum_c vals_c * a_c as a jet."""
    v = arb(0); g = [arb(0)] * N; h = [[arb(0)] * N for _ in range(N)]
    for a, c in zip(al, vals):
        v += a.v * c
        g = [g[i] + a.g[i] * c for i in range(N)]
        for i in range(N):
            for j in range(i + 1, N):
                h[i][j] += a.h[i][j] * c
    for i in range(N):
        for j in range(i + 1, N):
            h[j][i] = h[i][j]
    return Jet(v, g, h)


def pq_jet(u, i):
    g = [arb(0)] * N; g[i] = -u[i] / 2
    h = [[arb(0)] * N for _ in range(N)]; h[i][i] = arb(-0.5)
    return Jet((1 - u[i] * u[i]) / 4, g, h)


class CellData:
    def __init__(self, p, rho, Fc_lower, emode='es', delta_rel=0.05):
        """Fc_lower: list of 64 arb (rigorous corner lower bounds of lambda_1*Area), corner order = R.CORN.
        emode='exact' (v2, 2026-10-31): e_k(u) = lam_max(S_k^{-1} H_k), H_k = S_k - G_P,k = tbar V >= 0 computed from the
           exact rational expressions S_k = sum a_c G_c,k (multiaffine) and G_P,k = adj(D_P,k) adj(D_P,k)' / J_P,k
           (adj(D_P) and J_P multiaffine), with the polynomial majorants lam_max <= tr/2 + (q + s^2)/(2s),
           q = tr^2/4 - det, s > 0 constant per box (AM-GM), or tr itself (all eigenvalues >= 0); best per box.
           No Efron-Stein, no Gl, no Jmax: G_P >= (1 - e_k) S_k is then the sharpest uniform-in-sector statement.
        emode='es' (default): matrix Efron-Stein majorant (cheap, good on large boxes).
        emode='both': per box and sector, the smaller of the two majorants (both are valid upper bounds of e_k)."""
        self.emode = emode; self.delta_rel = delta_rel
        self.lo, self.hi = R.cell_endpoints(p, rho)
        g = R.convexity_gate(self.lo, self.hi)
        assert g > 0, 'convexity gate'
        corners = R.cell_corners(self.lo, self.hi)
        Ds = [R.Dmats(R.verts(x)) for x in corners]
        self.J = [[R.det2(D[k]) for k in range(5)] for D in Ds]              # 64 x 5
        self.Ac = [sum(R.a0[k] * self.J[c][k] for k in range(5)) for c in range(64)]
        self.lam = [Fc_lower[c] / self.Ac[c] for c in range(64)]
        sh = (2,) * N
        self.Jt = [np.array([self.J[c][k] for c in range(64)], dtype=object).reshape(sh) for k in range(5)]
        self.Nt = [np.array([self.lam[c] * self.J[c][k] for c in range(64)], dtype=object).reshape(sh) for k in range(5)]
        self.At = np.array(self.Ac, dtype=object).reshape(sh)
        self.pre = []
        self.St = []; self.Adt = []
        for k in range(5):
            Jc = [self.J[c][k] for c in range(64)]
            assert all(j > 0 for j in Jc)
            Xc = [R.inv2(D[k])[0] for D in Ds]
            Gc = [[[Jc[c] * (R.sq(Xc[c][0][0]) + R.sq(Xc[c][0][1])), Jc[c] * (Xc[c][0][0] * Xc[c][1][0] + Xc[c][0][1] * Xc[c][1][1])],
                   [Jc[c] * (Xc[c][0][0] * Xc[c][1][0] + Xc[c][0][1] * Xc[c][1][1]), Jc[c] * (R.sq(Xc[c][1][0]) + R.sq(Xc[c][1][1]))]]
                  for c in range(64)]
            self.St.append([np.array([G[0][0] for G in Gc], dtype=object).reshape(sh),
                            np.array([(G[0][1] + G[1][0]) / 2 for G in Gc], dtype=object).reshape(sh),
                            np.array([G[1][1] for G in Gc], dtype=object).reshape(sh)])
            self.Adt.append([np.array([D[k][1][1] for D in Ds], dtype=object).reshape(sh),       # adj = [[d, -b], [-c, a]]
                             np.array([-D[k][0][1] for D in Ds], dtype=object).reshape(sh),
                             np.array([-D[k][1][0] for D in Ds], dtype=object).reshape(sh),
                             np.array([D[k][0][0] for D in Ds], dtype=object).reshape(sh)])
            if emode == 'exact':                       # 'both' and 'es' need the Efron-Stein data
                self.pre.append(None); continue
            Gl, _which = R.lower_metric(Gc)
            assert Gl is not None and R.lam_min_sym2(Gl) > 0, 'no PD lower metric' 
            Gli = R.inv2(Gl)[0]
            RR = []; dev = []
            for i in range(N):
                bit = 1 << (5 - i)
                dXs = [[[Xc[c][a][b] - Xc[c ^ bit][a][b] for b in range(2)] for a in range(2)] for c in range(64) if not (c & bit)]
                Rm = [[sum(d[a][b] for d in dXs) / 32 for b in range(2)] for a in range(2)]
                RR.append(R.mul2(Rm, [[Rm[0][0], Rm[1][0]], [Rm[0][1], Rm[1][1]]]))
                vals = []
                for d in dXs:
                    dd = [[d[a][b] - Rm[a][b] for b in range(2)] for a in range(2)]
                    vals.append(R.lam_max_sym2(R.mul2(R.mul2([[dd[0][0], dd[1][0]], [dd[0][1], dd[1][1]]], Gli), dd)))
                dev.append(R.amax(vals))
            g11 = Gl[0][0]; g12 = (Gl[0][1] + Gl[1][0]) / 2; g22 = Gl[1][1]       # Cholesky Gl = L L'
            l11 = g11.sqrt(); l21 = g12 / l11; l22 = (g22 - l21 * l21).sqrt()
            Li = [[1 / l11, arb(0)], [-l21 / (l11 * l22), 1 / l22]]
            Rt = []
            for RRi in RR:
                Mt = R.mul2(R.mul2(Li, RRi), [[Li[0][0], Li[1][0]], [Li[0][1], Li[1][1]]])
                Rt.append((Mt[0][0], (Mt[0][1] + Mt[1][0]) / 2, Mt[1][1]))
            self.pre.append(dict(Jmax=R.amax(Jc), Gl=Gl, RR=RR, dev=dev, Rt=Rt))

    def eps_es(self, k, pqc, pqb, ra):
        """Efron-Stein majorant of eps_k: centre jet and box jet of  main  and  dev  (Jmax folded in),
        eps_k(tau) = (1+tau) main + (1+1/tau) dev,  main >= lam_max(L^{-1} (sum_i pq_i R_i R_i') L^{-T}),  Gl = L L'.
        With Mt = sum_i pq_i Rt_i = [[a, b], [b, d]] (Rt_i = L^{-1} R_i R_i' L^{-T}, constant):
           lam_max = (a+d)/2 + sqrt(q),  q = ((a-d)/2)^2 + b^2,
        and sqrt(q) <= (q + s^2) / (2 s) for any constant s > 0 (AM-GM), so
           main = (a+d)/2 + (q + s^2)/(2 s)   is a POLYNOMIAL majorant in pq (no branch, smooth everywhere),
        exact where q = s^2.  s is fixed per box from the centre: s = max(sqrt(q(centre)), 0.02 (a+d)(centre)).
        (v1/v2 used the closed form with a square root; near-isotropic Mt (lam_1 ~ lam_2, which happens at id71) makes
        the radicand ~ 0, so the box had to fall back to the trace = up to 2x loose.)"""
        P = self.pre[k]
        def pieces(pqs):
            a = Jet(arb(0)); b = Jet(arb(0)); d = Jet(arb(0)); dv = Jet(arb(0))
            for i in range(N):
                ai, bi, di = P['Rt'][i]
                a = a + pqs[i] * ai; b = b + pqs[i] * bi; d = d + pqs[i] * di
                dv = dv + pqs[i] * P['dev'][i]
            hd = (a - d) * arb(0.5)
            return a + d, hd * hd + b * b, dv
        trc, qc, dvc = pieces(pqc)
        trb, qb, dvb = pieces(pqb)
        J = P['Jmax']
        mains = [(trc * J, trb * J)]                      # trace majorant (a, d >= 0)
        t0 = float(trc.v.upper()); q0 = float(qc.v.upper()) ** 0.5
        for sv in sorted({max(q0, 0.02 * t0), max(2 * q0, 0.1 * t0), 0.25 * t0, 0.4 * t0}):
            if not sv > 0:
                continue
            s_ = arb(sv); inv2s = 1 / (2 * s_); s2 = s_ * s_
            mains.append(((trc * arb(0.5) + (qc + s2) * inv2s) * J, (trb * arb(0.5) + (qb + s2) * inv2s) * J))
        return mains, (dvc * J, dvb * J)

    def jets(self, cen, rad):
        """center jets (exact derivatives at u=cen) and box jets (enclosures over the box) of Q_k and eps parts.
        Multiaffine pieces (J_k, A, N_k) use their exact corner tensors: on a sub-box a multiaffine f, each df/du_i and
        each d2f/du_i du_j (i != j) is multiaffine, so its range is the hull of its values at the sub-box corners
        (obtained by finite differences of the sub-box corner values, exact for multiaffine f); d2f/du_i^2 = 0."""
        ca = [arb(float(c)) for c in cen]; ra = [arb(float(r)) for r in rad]
        Jc_, Jb_ = zip(*[ma_jets(self.Jt[k], ca, ra) for k in range(5)])
        Nc_, Nb_ = zip(*[ma_jets(self.Nt[k], ca, ra) for k in range(5)])
        Ac_, Ab_ = ma_jets(self.At, ca, ra)
        Jic = [Jc_[k].recip() for k in range(5)]; Jib = [Jb_[k].recip() for k in range(5)]
        Qc = [Ac_ * Nc_[k] * Jic[k] for k in range(5)]
        Qb = [Ab_ * Nb_[k] * Jib[k] for k in range(5)]
        Xc = Xb = None
        if self.emode in ('exact', 'both'):
            Ec = []; Eb = []
            for k in range(5):
                Sc, Sb = zip(*[ma_jets(t, ca, ra) for t in self.St[k]])
                Dc, Db = zip(*[ma_jets(t, ca, ra) for t in self.Adt[k]])
                trs = []; parts = []
                for (S_, A_, Ji) in ((Sc, Dc, Jic[k]), (Sb, Db, Jib[k])):
                    g11 = (A_[0] * A_[0] + A_[1] * A_[1]) * Ji
                    g12 = (A_[0] * A_[2] + A_[1] * A_[3]) * Ji
                    g22 = (A_[2] * A_[2] + A_[3] * A_[3]) * Ji
                    h11 = S_[0] - g11; h12 = S_[1] - g12; h22 = S_[2] - g22
                    dSi = (S_[0] * S_[2] - S_[1] * S_[1]).recip()
                    tr = (S_[2] * h11 + S_[0] * h22 - S_[1] * h12 * 2) * dSi
                    dt = (h11 * h22 - h12 * h12) * dSi
                    parts.append((tr, dt))
                # e = lam_max(S^-1 H) = tr/2 + sqrt(q), q = tr^2/4 - det >= 0;  sqrt(q) <= (q + s^2)/(2s)  (AM-GM, s const)
                (trc, dtc), (trb, dtb) = parts
                qc = trc * trc * arb(0.25) - dtc; qb = trb * trb * arb(0.25) - dtb
                cands = [(trc, trb)]                                   # trace majorant (both eigenvalues >= 0)
                t0 = float(trc.v.upper()); q0 = max(float(qc.v.upper()), 0.0) ** 0.5
                for sv in sorted({max(q0, 0.02 * t0), max(2 * q0, 0.1 * t0), 0.25 * t0, 0.4 * t0}):
                    if not sv > 0:
                        continue
                    s_ = arb(sv); i2s = 1 / (2 * s_); s2 = s_ * s_
                    cands.append((trc * arb(0.5) + (qc + s2) * i2s, trb * arb(0.5) + (qb + s2) * i2s))
                Ec.append([c for c, _ in cands]); Eb.append([b for _, b in cands])
            Xc, Xb = Ec, Eb
        Ec = Eb = None
        if self.emode in ('es', 'both'):
            ubx = [c.union(c - r).union(c + r) for c, r in zip(ca, ra)]
            pqb = [pq_jet(ubx, i) for i in range(N)]; pqc = [pq_jet(ca, i) for i in range(N)]
            Ec = [self.eps_es(k, pqc, pqb, ra) for k in range(5)]        # per k: (mains [(c, b)], (dev c, dev b))
            Eb = True
        return (Qc, Ec, Xc), (Qb, Eb, Xb)


def hull(xs):
    """ball containing all balls xs, using EXACT arb endpoints (v6: was via float(), which can round inward)."""
    lo = xs[0].lower(); hi = xs[0].upper()
    for x in xs[1:]:
        l = x.lower(); h = x.upper()
        if l < lo:
            lo = l
        if h > hi:
            hi = h
    return lo.union(hi)


def ma_jets(T, ca, ra):
    """T: object array (2,)*6 of corner values (index 0 <-> u_i=-1, 1 <-> u_i=+1).  returns (center Jet, box Jet)."""
    S = T
    for i in range(N):
        a = ca[i] - ra[i]; b = ca[i] + ra[i]
        S0 = S.take(0, axis=i); S1 = S.take(1, axis=i)
        n0 = S0 * ((1 - a) / 2) + S1 * ((1 + a) / 2)
        n1 = S0 * ((1 - b) / 2) + S1 * ((1 + b) / 2)
        S = np.stack([n0, n1], axis=i)
    flat = list(S.ravel())
    vc = sum(flat[1:], flat[0]) / 64
    vb = hull(flat)
    gc = []; gb = []
    for i in range(N):
        d = list(((S.take(1, axis=i) - S.take(0, axis=i)) / (2 * ra[i])).ravel())
        gc.append(sum(d[1:], d[0]) / 32); gb.append(hull(d))
    hc = [[arb(0)] * N for _ in range(N)]; hb = [[arb(0)] * N for _ in range(N)]
    for i in range(N):
        for j in range(i + 1, N):
            Si1 = S.take(1, axis=i); Si0 = S.take(0, axis=i)
            jj = j - 1                                   # axis index of j after removing axis i (j > i)
            d2 = (Si1.take(1, axis=jj) - Si1.take(0, axis=jj) - Si0.take(1, axis=jj) + Si0.take(0, axis=jj)) / (4 * ra[i] * ra[j])
            d2 = list(d2.ravel())
            hc[i][j] = hc[j][i] = sum(d2[1:], d2[0]) / 16
            hb[i][j] = hb[j][i] = hull(d2)
    return Jet(vc, gc, hc), Jet(vb, gb, hb)


def _taylor(jc, jb, ra, sign):
    """sign=-1: lower bound over the box; sign=+1: upper bound.  returns (bound, per-coordinate contributions)."""
    b = jc.v; con = np.zeros(N)
    for i in range(N):
        gi = abs(jc.g[i]) * ra[i]
        b = b + sign * gi
        con[i] += float(gi.upper())
        for j in range(N):
            hij = abs(jb.h[i][j]) * ra[i] * ra[j] / 2
            b = b + sign * hij
            con[i] += float(hij.upper()) / 2; con[j] += float(hij.upper()) / 2
    return b, con


def box_lower(cd, cen, rad, taus):
    """rigorous lower bound of B = (1 - max_k' eps_k') * min_k Q_k over the box, and per-coordinate contributions.
    v3 (2026-10-31): B = min over the 25 pairs (k, k') of P_kk' = Q_k * (1 - eps_k'), each pair bounded by its own
    second-order Taylor form, so the first-order correlation between Q (smallest at the box corner) and eps
    (zero at the cell corners, largest inside) is kept.  The decoupled form (min Q_k^lower) * (1 - max eps^upper)
    is valid but multiplies extremes taken at different points and stalls the branch-and-bound."""
    ra = [arb(float(r)) for r in rad]
    (Qc, Ec, Xc), (Qb, Eb, Xb) = cd.jets(cen, rad)
    best = None
    for tau in (taus if Ec is not None else (None,)):
        ecands = []                                   # per k': list of (centre jet, box jet) majorants of eps_k'
        for k in range(5):
            cs = []
            if tau is not None:
                t = arb(tau)
                mains, (dvc, dvb) = Ec[k]
                for mc, mb in mains:
                    cs.append((mc * (1 + t) + dvc * (1 + 1 / t), mb * (1 + t) + dvb * (1 + 1 / t)))
            if Xc is not None:
                cs.extend(zip(Xc[k], Xb[k]))
            # keep the majorant with the smallest Taylor upper bound on this box
            ups = [float(_taylor(jc, jb, ra, +1)[0].upper()) for jc, jb in cs]
            ecands.append(cs[int(np.argmin(ups))])
        lows = []; cons = []
        for k in range(5):
            for kk in range(5):
                jc = Qc[k] * (1 - ecands[kk][0]); jb = Qb[k] * (1 - ecands[kk][1])
                b, con = _taylor(jc, jb, ra, -1)
                lows.append(b); cons.append(con)
        i = int(np.argmin([float(x.lower()) for x in lows]))
        lb = arb(lows[i].lower())
        if best is None or lb.lower() > best[0].lower():
            best = (lb, cons[i], None, tau)
    return best


def certify_T1(cd, Fstar_up, max_boxes=20000, min_rad=1e-3, taus=(0.03, 0.1, 0.3), log_every=200):
    """branch-and-bound over u in [-1,1]^6 (bisect the coordinate with the largest Taylor-loss contribution)."""
    t0 = time.time()
    heap = [(0.0, 0, tuple([0.0] * N), tuple([1.0] * N))]
    cnt = 0; nacc = 0; worst = None; seq = 1; vol = 0.0; leaves = []
    while heap:
        _, _, cen, rad = heapq.heappop(heap)
        cnt += 1
        if cnt > max_boxes:
            return dict(status='budget', boxes=cnt, accepted=nacc, vol_frac=vol, sec=round(time.time() - t0, 1))
        lb, con, emax, tau = box_lower(cd, cen, rad, taus)
        if lb > Fstar_up:
            nacc += 1; vol += float(np.prod(rad)); leaves.append([list(cen), list(rad), tau, R.f_down(lb)])
            worst = lb if worst is None else arb(min(worst.lower(), lb.lower()))
            continue
        i = int(np.argmax(con * (np.array(rad) > min_rad)))
        if rad[i] <= min_rad:
            return dict(status='fail', boxes=cnt, accepted=nacc, at=list(cen), rad=list(rad), lower=float(lb.lower()),
                        emax=emax, vol_frac=vol, sec=round(time.time() - t0, 1))
        r2 = list(rad); r2[i] = rad[i] / 2
        for sgn in (-1, 1):
            c2 = list(cen); c2[i] = cen[i] + sgn * rad[i] / 2
            heapq.heappush(heap, (float(lb.lower()), seq, tuple(c2), tuple(r2))); seq += 1
        if cnt % log_every == 0:
            print('  ...boxes', cnt, 'accepted', nacc, 'vol', round(vol, 4), 'queue', len(heap), 'last lb', float(lb.lower()),
                  round(time.time() - t0), 's', file=sys.stderr, flush=True)
    return dict(status='certified', boxes=cnt, accepted=nacc, vol_frac=vol, min_lower=float(worst.lower()),
                sec=round(time.time() - t0, 1), leaves=leaves,
                cover_note='leaves are the accepted nodes of a halving tree rooted at [-1,1]^6: an exact partition (vol_frac=1)')


def certify_cell_T1(p, rho, K=8, nseg=1000, workers=1, corners_from=None, emode='es'):
    """rigorous corners (rig_box via rig_paqc._cf) + min-corner bound (for comparison) + T1 branch-and-bound."""
    t0 = time.time()
    lo, hi = R.cell_endpoints(p, rho)
    corners = R.cell_corners(lo, hi)
    res = dict(rho=rho, p=[float(v) for v in p], cell_lo=[repr(float(v)) for v in lo], cell_hi=[repr(float(v)) for v in hi],
               Fstar_up=float(R.FSTAR_UP.upper()), method='T1 v3: min_{k,k2} Q_k (1-eps_k2), 2nd-order Taylor B&B over u, Efron-Stein eps with AM-GM polynomial lam_max majorant')
    if corners_from:
        certs = json.load(open(corners_from))['corner_certificates']
        for c, x in zip(certs, corners):
            assert [float(v) for v in c['vertices_x']] == [float(v) for v in x], 'corner mismatch'
    else:
        if workers > 1:
            from multiprocessing import Pool
            with Pool(workers) as pool:
                outs = pool.starmap(R._cf, [(list(x), K, nseg) for x in corners])
        else:
            outs = [R._cf(list(x), K, nseg) for x in corners]
        if any(o is None for o in outs):
            res.update(status='fail_corner_index'); return res
        certs = [o[2] for o in outs]
    res['corner_certificates'] = certs
    Fl = [arb(R.cert_F_lo(c)) for c in certs]          # v6: recomputed from the saved balls, rounded down
    res['corner_sec'] = round(time.time() - t0, 1)
    g = R.geometry_bound(p, rho)
    if g['ok']:
        res['mincorner_LB'] = float((g['factor'] * R.amin(Fl)).lower())
        res['mincorner_certified'] = bool(g['factor'] * R.amin(Fl) > R.FSTAR_UP)
    cd = CellData(p, rho, Fl, emode=emode)
    res['emode'] = emode
    t = certify_T1(cd, R.FSTAR_UP)
    res['T1'] = t
    res['certified'] = t['status'] == 'certified'
    res['sec'] = round(time.time() - t0, 1)
    return res


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'cell':                 # rig_t1.py cell P rho workers [corners_json]
        p = np.array(json.loads(sys.argv[2])); rho = float(sys.argv[3]); w = int(sys.argv[4]) if len(sys.argv) > 4 else 1
        cf = sys.argv[5] if len(sys.argv) > 5 and sys.argv[5] != '-' else None
        em = sys.argv[6] if len(sys.argv) > 6 else 'es'
        print(json.dumps(certify_cell_T1(p, rho, workers=w, corners_from=cf, emode=em)), flush=True)
    if mode == 'floatprobe':           # float corner values (allowance) -> B&B test, for choosing rho before rigorous corners
        p = np.array(json.loads(sys.argv[2])); rho = float(sys.argv[3]); allow = float(sys.argv[4]) if len(sys.argv) > 4 else 3e-4
        sys.path.insert(0, R.C456)
        from mps import PolyMPS, poly_area
        lo, hi = R.cell_endpoints(p, rho)
        Fc = []
        cache = sys.argv[5] if len(sys.argv) > 5 else None
        if cache:
            Fc = [arb(v - allow) for v in json.load(open(cache))['%.6f' % rho]]
        for x in ([] if cache else R.cell_corners(lo, hi)):
            V = np.array([[0, 0], [1, 0], [x[0], x[1]], [x[2], x[3]], [x[4], x[5]]], float); a = poly_area(V)
            Fc.append(arb(PolyMPS(V, K=11).solve(18.95 / a, width=0.04, ngrid=13) * a - allow))
        cd = CellData(p, rho, Fc)
        print(json.dumps(certify_T1(cd, R.FSTAR_UP)), flush=True)
