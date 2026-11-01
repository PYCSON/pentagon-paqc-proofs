"""Rigorous single-box certificate (M3 pilot for C4), own code, MEAS-01.

Claim certified:  for EVERY pentagon L with v0=(0,0), v1=(1,0) and
(v2,v3,v4) in the sup-norm box  p +- d  (6-D),   lambda_1(L)*Area(L) > F*_up,
where F*_up is a rigorous upper bound for lambda_1*Area of the regular pentagon.

Chain (all inequalities enclosed with arb ball arithmetic, python-flint):
 1. Outer polygon Qt (float vertices) built from the centre pentagon's normals;
    containment L subset Qt checked exactly for all box corners  -> domain
    monotonicity  lambda_1(L) >= lambda_1(Qt).
 2. Area(L) >= A_min = min over the 64 box corners (area is multi-affine).
 3. lambda_1(Qt): float MPS gives (lam~, u); u solves Delta u + lam~ u = 0 in Qt
    exactly (Bessel-sine particular solutions at the 5 corners).  Fox-Henrici-Moler:
    some eigenvalue lam_k has |lam~ - lam_k|/lam_k <= mu,
    mu = sqrt|Qt| sup_bdry|u| / ||u||_2.   sup_bdry|u|: Taylor-2 on edge segments
    (value+slope at segment midpoint, interval bound on 2nd derivative);
    ||u||_2: lower bound by interval enclosures on small interior triangles.
 4. Index: Hong-Krahn-Szego  lambda_2 |Om| >= 2 pi j01^2;  if lam~/(1-mu)|Qt| is
    below that, lam_k = lambda_1.  Hence lambda_1(Qt) >= lam~/(1+mu).
 5. Same FHM+HKS at the regular pentagon gives F*_up = lam~_R5/(1-mu_R5) * A_R5.
"""
import itertools, sys, time
import numpy as np
from flint import arb, ctx
from mps import PolyMPS, regular, poly_area

ctx.prec = 80
PI = arb.pi()
J01 = arb("2.404825557695772768621631879326454643124244909145")  # first zero of J0
HKS = 2 * PI * J01 ** 2


def A(x):
    return arb(float(x))


class ArbU:
    """u = sum_j sum_m c_jm J_{m pi/a_j}(k r_j) sin(m pi th_j/a_j) in ball arithmetic."""

    def __init__(self, m):
        self.V = m.V
        self.n = len(m.V)
        self.K = m.K
        self.k = arb(float(m.lam)).sqrt()
        self.c = [float(x) for x in m.c]
        self.fr = []
        for j in range(self.n):
            p = [A(self.V[j][0]), A(self.V[j][1])]
            a = [A(self.V[(j + 1) % self.n][0]) - p[0], A(self.V[(j + 1) % self.n][1]) - p[1]]
            b = [A(self.V[(j - 1) % self.n][0]) - p[0], A(self.V[(j - 1) % self.n][1]) - p[1]]
            phi0 = arb.atan2(a[1], a[0])
            phib = arb.atan2(b[1], b[0])
            al = phib - phi0
            if al < 0:
                al += 2 * PI
            nus = [m_ * PI / al for m_ in range(1, self.K + 1)]
            self.fr.append((p, phi0.cos(), phi0.sin(), al, nus))

    def polar(self, j, x, y):
        p, c0, s0, al, nus = self.fr[j]
        dx, dy = x - p[0], y - p[1]
        wx = c0 * dx + s0 * dy           # rotate by -phi0
        wy = -s0 * dx + c0 * dy
        r = (wx * wx + wy * wy).sqrt()
        th = arb.atan2(wy, wx)
        return r, th, wx, wy

    def value(self, x, y, skip=()):
        tot = arb(0)
        for j in range(self.n):
            if j in skip:
                continue
            r, th, _, _ = self.polar(j, x, y)
            z = self.k * r
            for mm, nu in enumerate(self.fr[j][4]):
                cc = self.c[j * self.K + mm]
                tot += cc * z.bessel_j(nu) * (nu * th).sin()
        return tot

    def d1d2(self, x, y, ex, ey, skip=()):
        """(f', f'') along direction e at (x,y): f'=grad u.e, f''=e^T Hess u e."""
        f1 = arb(0); f2 = arb(0)
        for j in range(self.n):
            if j in skip:
                continue
            p, c0, s0, al, nus = self.fr[j]
            r, th, wx, wy = self.polar(j, x, y)
            # direction in local frame
            lex = c0 * ex + s0 * ey
            ley = -s0 * ex + c0 * ey
            cr, sr = wx / r, wy / r                 # e_r components
            a = lex * cr + ley * sr                 # e . e_r
            b = -lex * sr + ley * cr                # e . e_theta
            z = self.k * r
            k = self.k
            for mm, nu in enumerate(nus):
                cc = self.c[j * self.K + mm]
                J0 = z.bessel_j(nu); Jm1 = z.bessel_j(nu - 1); Jp1 = z.bessel_j(nu + 1)
                Jm2 = z.bessel_j(nu - 2); Jp2 = z.bessel_j(nu + 2)
                Jd = (Jm1 - Jp1) / 2
                Jdd = (Jm2 - 2 * J0 + Jp2) / 4
                S = (nu * th).sin(); C = (nu * th).cos()
                ur = k * Jd * S
                ut_r = nu * J0 * C / r
                Hrr = k * k * Jdd * S
                Hrt = k * nu * Jd * C / r - nu * J0 * C / (r * r)
                Htt = k * Jd * S / r - nu * nu * J0 * S / (r * r)
                f1 += cc * (a * ur + b * ut_r)
                f2 += cc * (a * a * Hrr + 2 * a * b * Hrt + b * b * Htt)
        return f1, f2


def sup_boundary(U, nseg):
    best = arb(0)
    V = U.V; n = U.n
    for i in range(n):
        q0 = [A(V[i][0]), A(V[i][1])]
        q1 = [A(V[(i + 1) % n][0]), A(V[(i + 1) % n][1])]
        ex, ey = q1[0] - q0[0], q1[1] - q0[1]
        skip = (i, (i + 1) % n)
        h = arb(1) / nseg
        for s in range(nseg):
            tm = (arb(s) + arb("0.5")) * h
            xm, ym = q0[0] + tm * ex, q0[1] + tm * ey
            f0 = U.value(xm, ym, skip)
            f1, _ = U.d1d2(xm, ym, ex, ey, skip)
            tI = tm + arb(0, 0.5 / nseg + 1e-18)          # whole segment
            xI, yI = q0[0] + tI * ex, q0[1] + tI * ey
            _, f2 = U.d1d2(xI, yI, ex, ey, skip)
            bound = abs(f0) + abs(f1) * h / 2 + abs(f2) * h * h / 8
            best = arb(max(best.upper(), bound.upper()))
    return best


def l2_lower(U, ngrid=None):
    """Rigorous lower bound of ||u||_2^2 via the Helmholtz mean-value identity:
    for a disk B(x0,R) inside the polygon, int_B u = u(x0) 2 pi R J1(kR)/k exactly,
    hence ||u||^2 >= int_B u^2 >= (int_B u)^2/(pi R^2).  x0 = vertex centroid,
    R = 0.999 * (distance from x0 to the nearest edge line)."""
    V = U.V; n = U.n
    c = V.mean(axis=0)
    dmin = None
    for i in range(n):
        q0 = [A(V[i][0]), A(V[i][1])]; q1 = [A(V[(i + 1) % n][0]), A(V[(i + 1) % n][1])]
        ex, ey = q1[0] - q0[0], q1[1] - q0[1]
        dist = abs(ex * (A(c[1]) - q0[1]) - ey * (A(c[0]) - q0[0])) / (ex * ex + ey * ey).sqrt()
        dmin = dist if dmin is None else arb(min(dmin.lower(), dist.lower()))
    R = arb(dmin.lower()) * arb("0.999")
    u0 = U.value(A(c[0]), A(c[1]))
    kR = U.k * R
    intB = u0 * 2 * PI * R * kR.bessel_j(1) / U.k
    lb = intB * intB / (PI * R * R)
    return arb(lb.lower())


def certify_polygon(V, K, guess, nseg=500, ngrid=None, width=0.03):
    m = PolyMPS(V, K=K)
    lam = m.solve(guess, width=width, ngrid=11)
    U = ArbU(m)
    sup = sup_boundary(U, nseg)
    l2 = l2_lower(U, ngrid)
    xs = [A(v[0]) for v in V]; ys = [A(v[1]) for v in V]
    area = abs(sum(xs[i] * ys[(i + 1) % 5] - xs[(i + 1) % 5] * ys[i] for i in range(5))) / 2
    mu = area.sqrt() * sup / l2.sqrt()
    lamt = arb(float(lam))
    ok_index = (lamt / (1 - mu) * area) < HKS
    lam_lo = lamt / (1 + mu)
    lam_hi = lamt / (1 - mu)
    return dict(lam=lam, mu=mu, ok_index=bool(ok_index), lam_lo=lam_lo, lam_hi=lam_hi, area=area)


CORN = np.array(list(itertools.product([-1, 1], repeat=6)), float)


def pent(p):
    return np.array([[0, 0], [1, 0], [p[0], p[1]], [p[2], p[3]], [p[4], p[5]]], float)


def certify_box(p, d, K=8, nseg=500, ngrid=None):
    V = pent(p)
    # normals of the centre pentagon
    nr = []
    for i in range(5):
        e = V[(i + 1) % 5] - V[i]
        nr.append(np.array([e[1], -e[0]]) / np.hypot(*e))
    nr = np.array(nr)
    # support numbers over the box (vertex boxes for j=2,3,4)
    H = np.zeros(5)
    for i in range(5):
        vals = [V[0] @ nr[i], V[1] @ nr[i]]
        for j in (2, 3, 4):
            vals.append(V[j] @ nr[i] + d * (abs(nr[i][0]) + abs(nr[i][1])))
        H[i] = max(vals) + 1e-12
    Q = np.array([np.linalg.solve(np.array([nr[i - 1], nr[i]]), [H[i - 1], H[i]]) for i in range(5)])
    # rigorous containment: every box-corner vertex lies left of every edge of Q
    for i in range(5):
        q0, q1 = Q[i], Q[(i + 1) % 5]
        for s in CORN:
            L = pent(p + d * s)
            for v in L:
                cr = (A(q1[0]) - A(q0[0])) * (A(v[1]) - A(q0[1])) - (A(q1[1]) - A(q0[1])) * (A(v[0]) - A(q0[0]))
                if not cr >= 0:
                    return dict(fail="containment")
    # rigorous min area over the box (multi-affine -> corners)
    Amin = None
    for s in CORN:
        L = pent(p + d * s)
        xs = [A(v[0]) for v in L]; ys = [A(v[1]) for v in L]
        ar = sum(xs[i] * ys[(i + 1) % 5] - xs[(i + 1) % 5] * ys[i] for i in range(5)) / 2
        Amin = ar if Amin is None else arb(min(Amin.lower(), ar.lower()))
    guess = 18.92 / poly_area(Q)
    cQ = certify_polygon(Q, K, guess, nseg, ngrid, width=0.06)
    return dict(Q=cQ, Amin=Amin, F_lo=cQ["lam_lo"] * Amin)


if __name__ == "__main__":
    t0 = time.time()
    # (a) rigorous upper bound for F* at the regular pentagon (unit circumradius)
    R5 = regular(5)
    cR = certify_polygon(R5, 8, 18.9191042862 / poly_area(R5), nseg=int(sys.argv[1]) if len(sys.argv) > 1 else 500)
    Fstar_up = cR["lam_hi"] * cR["area"]
    print("R5: mu=%s index_ok=%s  F* in [%s, %s]  (%.0fs)" % (cR["mu"].str(3), cR["ok_index"],
          (cR["lam_lo"] * cR["area"]).str(15), Fstar_up.str(15), time.time() - t0), flush=True)
    # (b) one box: the r=0.1 test centre direction of paqc (seed 7 first draw) is not
    #     reproducible here, so use a fixed off-centre pentagon at distance 0.1
    p0 = np.array([1.309017, 0.951057, 0.5, 1.538842, -0.309017, 0.951057])
    rng = np.random.default_rng(11)
    dirn = rng.normal(size=6); dirn /= np.linalg.norm(dirn)
    p = p0 + 0.1 * dirn
    for d in [float(x) for x in sys.argv[2:]] or [1e-4, 3e-4]:
        t1 = time.time()
        res = certify_box(p, d, nseg=int(sys.argv[1]) if len(sys.argv) > 1 else 500)
        if "fail" in res:
            print("box d=%.1e: %s" % (d, res["fail"]))
            continue
        q = res["Q"]
        print("box d=%.1e: lam~(Q)=%.12f mu=%s index_ok=%s  A_min=%s  F_lo=%s  F_lo>F*_up: %s  (%.0fs)" % (
            d, q["lam"], q["mu"].str(3), q["ok_index"], res["Amin"].str(12), res["F_lo"].str(12),
            bool(res["F_lo"] > Fstar_up), time.time() - t1), flush=True)
    print("total %.0fs" % (time.time() - t0))
