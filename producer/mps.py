"""Method of particular solutions (Betcke-Trefethen subspace-angle variant) for the
first Dirichlet eigenvalue of a convex polygon, plus edge integrals of (d_n u)^2.

Own code, MEAS-01 / C4.  Float64 (scipy.special.jv).  Not rigorous.

Basis: at every vertex j with interior angle a_j, functions
    J_{k pi/a_j}(sqrt(lam) r_j) sin(k pi th_j / a_j),  k = 1..K
which vanish on the two edges adjacent to vertex j.  For convex polygons every
point of the closed polygon lies in every vertex wedge (th_j in [0, a_j]).
"""
import numpy as np
from scipy.special import jv
from scipy.optimize import minimize_scalar


def poly_area(V):
    x, y = V[:, 0], V[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)


def make_ccw(V):
    V = np.asarray(V, float)
    if poly_area(V) < 0:
        V = V[::-1].copy()
    return V


def vertex_frames(V):
    n = len(V)
    frames = []
    for j in range(n):
        p, a, b = V[j], V[(j + 1) % n], V[(j - 1) % n]
        e1 = a - p
        e2 = b - p
        phi0 = np.arctan2(e1[1], e1[0])
        ang = np.arctan2(e2[1], e2[0]) - phi0
        ang = np.mod(ang, 2 * np.pi)          # interior angle (CCW from e1 to e2)
        frames.append((p, phi0, ang))
    return frames


def local_coords(P, frame):
    p, phi0, ang = frame
    d = P - p
    r = np.hypot(d[:, 0], d[:, 1])
    th = np.mod(np.arctan2(d[:, 1], d[:, 0]) - phi0, 2 * np.pi)
    th = np.where(th > ang + 1e-9 + np.pi, th - 2 * np.pi, th)  # tiny negative angles
    return r, th


class PolyMPS:
    def __init__(self, V, K=14, nb=None, ni=None, seed=0):
        self.V = make_ccw(V)
        self.n = len(self.V)
        self.frames = vertex_frames(self.V)
        self.K = K
        self.rtol = 1e-14
        n = self.n
        nb = nb or 3 * K
        # boundary points: Chebyshev-like clustering toward corners on each edge
        t = 0.5 * (1 - np.cos(np.pi * (np.arange(nb) + 0.5) / nb))
        B = []
        for j in range(n):
            a, b = self.V[j], self.V[(j + 1) % n]
            B.append(a[None, :] + t[:, None] * (b - a)[None, :])
        self.B = np.vstack(B)
        # interior points: random points in fan triangles (away from boundary)
        rng = np.random.default_rng(seed)
        c = self.V.mean(axis=0)
        ni = ni or 4 * K * n
        I = []
        for k in range(ni):
            j = rng.integers(n)
            u, v = rng.random(2)
            if u + v > 1:
                u, v = 1 - u, 1 - v
            s = 0.85
            pt = c + s * (u * (self.V[j] - c) + v * (self.V[(j + 1) % n] - c))
            I.append(pt)
        self.I = np.array(I)
        self.coordsB = [local_coords(self.B, f) for f in self.frames]
        self.coordsI = [local_coords(self.I, f) for f in self.frames]

    def basis(self, k, coords):
        cols = []
        for (r, th), f in zip(coords, self.frames):
            ang = f[2]
            for m in range(1, self.K + 1):
                nu = m * np.pi / ang
                cols.append(jv(nu, k * r) * np.sin(nu * th))
        return np.array(cols).T

    def _Q(self, lam):
        k = np.sqrt(lam)
        AB = self.basis(k, self.coordsB)
        AI = self.basis(k, self.coordsI)
        A = np.vstack([AB, AI])
        sc = np.linalg.norm(A, axis=0)
        sc[sc == 0] = 1
        A = A / sc
        U, S, Vt = np.linalg.svd(A, full_matrices=False)
        keep = S > self.rtol * S[0]
        return U[:, keep], S[keep], Vt[keep], sc, len(AB)

    def sigma(self, lam):
        U, S, Vt, sc, nB = self._Q(lam)
        s = np.linalg.svd(U[:nB], compute_uv=False)
        return s[-1]

    def solve(self, lam_guess, width=0.05, tol=1e-15, ngrid=41):
        lo, hi = lam_guess * (1 - width), lam_guess * (1 + width)
        grid = np.linspace(lo, hi, ngrid)
        vals = [self.sigma(g) for g in grid]
        i = int(np.argmin(vals))
        a, b = grid[max(i - 1, 0)], grid[min(i + 1, len(grid) - 1)]
        # golden section on the V-shaped sigma
        gr = (np.sqrt(5) - 1) / 2
        c, d = b - gr * (b - a), a + gr * (b - a)
        fc, fd = self.sigma(c), self.sigma(d)
        while (b - a) > 1e-9 * lam_guess:
            if fc < fd:
                b, d, fd = d, c, fc
                c = b - gr * (b - a)
                fc = self.sigma(c)
            else:
                a, c, fc = c, d, fd
                d = a + gr * (b - a)
                fd = self.sigma(d)
        # V-fit refinement: sigma ~ s*|lam-lam*| (floor ignored)
        x0 = 0.5 * (a + b)
        h = 1e-9 * lam_guess
        for _ in range(3):
            xl, xr = x0 - h, x0 + h
            fl, fr = self.sigma(xl), self.sigma(xr)
            f0 = self.sigma(x0)
            # slopes from outer points; assume x* between xl and xr
            # use symmetric V: f = s|x-x*|; with fl, fr: x* = x0 + h*(fl-fr)/(fl+fr)
            x0 = x0 + h * (fl - fr) / (fl + fr) if (fl + fr) > 0 else x0
            h = h * 1e-2
            if h < tol * lam_guess:
                break
        self.lam = x0
        self.sig = self.sigma(x0)
        self._coeffs()
        return self.lam

    def _coeffs(self):
        U, S, Vt, sc, nB = self._Q(self.lam)
        _, s, wt = np.linalg.svd(U[:nB])
        y = wt[-1]
        # A_scaled = U S Vt ; want A_scaled x = U y  ->  x = Vt^T S^-1 y
        x = Vt.T @ (y / S)
        self.c = x / sc

    # ---- evaluation of u and grad u ----
    def eval(self, P):
        P = np.atleast_2d(P)
        k = np.sqrt(self.lam)
        u = np.zeros(len(P))
        gx = np.zeros(len(P))
        gy = np.zeros(len(P))
        idx = 0
        for f in self.frames:
            r, th = local_coords(P, f)
            p, phi0, ang = f
            # unit vectors
            phi = th + phi0
            er = np.stack([np.cos(phi), np.sin(phi)], 1)
            et = np.stack([-np.sin(phi), np.cos(phi)], 1)
            rs = np.where(r == 0, 1e-300, r)
            for m in range(1, self.K + 1):
                nu = m * np.pi / ang
                cm = self.c[idx]
                idx += 1
                J = jv(nu, k * r)
                dJ = 0.5 * (jv(nu - 1, k * r) - jv(nu + 1, k * r))
                s, co = np.sin(nu * th), np.cos(nu * th)
                u += cm * J * s
                dr = cm * k * dJ * s
                dt = cm * nu * J * co / rs
                gx += dr * er[:, 0] + dt * et[:, 0]
                gy += dr * er[:, 1] + dt * et[:, 1]
        return u, gx, gy

    def l2norm2(self, nq=40):
        """int_Omega u^2 by fan triangulation from centroid, collapsed (Duffy) at
        the centroid; Gauss-Legendre tensor rule."""
        c = self.V.mean(axis=0)
        xg, wg = np.polynomial.legendre.leggauss(nq)
        xg = 0.5 * (xg + 1)
        wg = 0.5 * wg
        tot = 0.0
        for j in range(self.n):
            a, b = self.V[j], self.V[(j + 1) % self.n]
            # point = c + s*((1-t) a + t b - c), jacobian = s * |det(a-c, b-c)|
            det = abs((a - c)[0] * (b - c)[1] - (a - c)[1] * (b - c)[0])
            S, T = np.meshgrid(xg, xg, indexing="ij")
            W = np.outer(wg, wg)
            pts = c + S[..., None] * ((1 - T)[..., None] * a + T[..., None] * b - c)
            u, _, _ = self.eval(pts.reshape(-1, 2))
            tot += np.sum(W.reshape(-1) * S.reshape(-1) * u ** 2) * det
        return tot

    def edge_w(self, nq=60):
        """w_i = int_{edge i} (d_n u)^2 ds for u normalised in L2.  Edge i runs
        from V[i] to V[i+1].  Graded substitution s = 3t^2-2t^3 near both ends."""
        nrm2 = self.l2norm2()
        xg, wg = np.polynomial.legendre.leggauss(nq)
        t = 0.5 * (xg + 1)
        wt = 0.5 * wg
        # cubic grading
        s = t ** 3 / (t ** 3 + (1 - t) ** 3)
        ds = 3 * t ** 2 * (1 - t) ** 2 / (t ** 3 + (1 - t) ** 3) ** 2
        w = np.zeros(self.n)
        normals = np.zeros((self.n, 2))
        hs = np.zeros(self.n)
        for i in range(self.n):
            a, b = self.V[i], self.V[(i + 1) % self.n]
            e = b - a
            L = np.hypot(*e)
            nvec = np.array([e[1], -e[0]]) / L     # outward for CCW
            pts = a + s[:, None] * e
            _, gx, gy = self.eval(pts)
            dn = gx * nvec[0] + gy * nvec[1]
            w[i] = np.sum(wt * ds * dn ** 2) * L / nrm2
            normals[i] = nvec
            hs[i] = np.dot(a, nvec)
        self.nrm2 = nrm2
        return w, normals, hs

    def bdry_residual(self, m=400):
        """FHM quantity mu = sqrt|A| max_bdry|u| / ||u||_2 (sampled, not rigorous)."""
        t = np.linspace(0, 1, m)
        mx = 0
        for j in range(self.n):
            a, b = self.V[j], self.V[(j + 1) % self.n]
            pts = a + t[:, None] * (b - a)
            u, _, _ = self.eval(pts)
            mx = max(mx, np.max(np.abs(u)))
        nrm2 = self.l2norm2()
        return np.sqrt(abs(poly_area(self.V))) * mx / np.sqrt(nrm2)


def regular(n, R=1.0, rot=0.0):
    th = rot + 2 * np.pi * np.arange(n) / n
    return np.stack([R * np.cos(th), R * np.sin(th)], 1)


if __name__ == "__main__":
    import time
    t0 = time.time()
    sq = np.array([[0, 0], [1, 0], [1, 1], [0, 1.]])
    m = PolyMPS(sq, K=10)
    print("square", m.solve(2 * np.pi ** 2), 2 * np.pi ** 2, time.time() - t0)
    tri = regular(3)
    side = np.sqrt(3)
    ex = 16 * np.pi ** 2 / (3 * side ** 2)
    m = PolyMPS(tri, K=10)
    print("triangle", m.solve(ex), ex)
    P5 = regular(5)
    s5 = 2 * np.sin(np.pi / 5)
    for K in (8, 12, 16, 20):
        t0 = time.time()
        m = PolyMPS(P5, K=K)
        lam = m.solve(10.9964270845598 / s5 ** 2)
        print("pentagon K=%d" % K, "lam*side^2=%.16f" % (lam * s5 ** 2), "sigma=%.2e" % m.sig,
              "t=%.2fs" % (time.time() - t0))
    print("Jones 2017 unit-edge: 10.996427084559806648")
