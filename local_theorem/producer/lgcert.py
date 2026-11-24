"""Polynomial lower bound for the first Dirichlet eigenvalue of a convex polygon (Claude, 2026-11-16).  DRAFT for review.

Lemma LG (a one-vector Lehmann-Goerisch bound).  Let Omega be a bounded domain, 0 < lambda_1 <= lambda_2 <= ... its
Dirichlet eigenvalues, and nu a number with nu <= lambda_2.  Let u in H^1_0(Omega), u != 0, and w in H(div; Omega) with
        -div w = u.
Put  A0 = |grad u|^2,  A1 = |u|^2,  A2 = |w|^2  (L^2 norms squared),  D = A0 - nu A1,  G = A0 - 2 nu A1 + nu^2 A2.
If D < 0 then  G - D > 0  and
        lambda_1 >= nu (-D) / (G - D).
Proof.  Let psi in H^1_0 solve -Lap psi = u.  For f in H^1_0:  int grad f . w = int f u = int grad f . grad psi, so
w - grad psi is L^2-orthogonal to gradients of H^1_0 functions and A2 = |grad psi|^2 + |w - grad psi|^2 >= B2 := (u, A^-1 u).
Expand u = sum c_k phi_k in the eigenfunctions, p_k = c_k^2:  A1 = sum p_k,  A0 = sum lambda_k p_k,  B2 = sum p_k / lambda_k.
With G0 := A0 - 2 nu A1 + nu^2 B2 = sum p_k (lambda_k - nu)^2 / lambda_k >= 0 we have G >= G0, and t -> nu (-D) / (t - D)
is decreasing for t > D, so it suffices to show lambda_1 (G0 - D) >= nu (-D), i.e. lambda_1 G0 >= (-D)(nu - lambda_1).
D < 0 gives lambda_1 <= A0 / A1 < nu.  Since lambda_k >= nu for k >= 2,
        -D = sum p_k (nu - lambda_k) <= p_1 (nu - lambda_1),       lambda_1 G0 >= p_1 (nu - lambda_1)^2,
hence lambda_1 G0 >= p_1 (nu - lambda_1)^2 >= (-D)(nu - lambda_1).  (G - D >= G0 - D > 0 because -D > 0.)   qed

Use for a convex polygon P (vertices exact binary64, counter-clockwise):
  1. frame: p -> (p - c) / s with dyadic c, s (lambda scales by 1 / s^2);
  2. P+ : a polygon with vertices on the grid 2^-GRID containing the image of P (exact containment test); by domain
     monotonicity lambda_1(P) >= lambda_1(P+).  All polynomial algebra is done on P+ (small rational numbers);
  3. u = b q,  b = product of the affine functions of the edges of P+ (so u = 0 on the boundary), q a polynomial;
  4. w = polynomial vector field with -div w = u EXACTLY: a floating-point field wt is corrected by (R, 0),
     R = the x-antiderivative of the polynomial  -u - div wt;  the identity is then verified exactly;
  5. nu = the larger of
        PI2_LO * |d|^2 * min(4 / e1^2 + 1 / e2^2, 1 / e1^2 + 4 / e2^2):  the second eigenvalue of the rectangle spanned by
        the extents e1, e2 of P+ along a direction d and d-perp (lambda_2 is monotone under inclusion; best of a few d),
        KS_LO / |P+|:  the Krahn-Szego inequality lambda_2(Omega) >= 2 pi j_{0,1}^2 / |Omega| (the second eigenfunction
        has two nodal domains, one of area <= |Omega| / 2; Faber-Krahn on it);
  6. A0, A1, A2 are integrals of polynomials over P+, computed exactly (rational arithmetic; Green's theorem on the
     edges), and the bound is evaluated exactly.
The floating-point part (Ritz vector for q, least-squares flux for wt) only chooses q and wt; nothing in the proof
depends on its accuracy.  Returned: F_lo <= lambda_1(P) |P| (binary64, rounded down) and the data to re-check it."""
import math
from fractions import Fraction
import numpy as np
from numpy.polynomial import legendre as Lg
from scipy.linalg import eigh, null_space
from flint import fmpq, fmpz, fmpq_poly, fmpq_mpoly_ctx

GRID = 16                    # P+ has vertices on the grid 2^-GRID (normalised frame)
EPS0 = 2.0 ** -13            # first inflation factor 1 + EPS0 of P about its centroid
CBITS = 40                   # coefficients of q and wt are rounded to multiples of 2^-CBITS (relative to the largest)
PI2_LO = fmpq(98696044, 10 ** 7)                                   # 9.8696044 < pi^2 = 9.8696044010...
KS_LO = fmpq(3633682907, 10 ** 8)                                # 36.33682907 < 2 pi j_{0,1}^2 = 36.336829071...
CTX = fmpq_mpoly_ctx.get(('x', 'y'))
X, Y = CTX.gens()
_cache = {}


def _fq(v):
    a, b = float(v).as_integer_ratio(); return fmpq(a, b)


def _leg_exact(n):
    """L_0..L_n as exact polynomials in x and in y"""
    if ('leg', n) not in _cache:
        out = []
        for g in (X, Y):
            L = [CTX.from_dict({(0, 0): fmpq(1)}), g]
            for k in range(1, n):
                L.append((fmpq(2 * k + 1, k + 1) * g * L[k]) - fmpq(k, k + 1) * L[k - 1])
            out.append(L[:n + 1])
        _cache[('leg', n)] = out
    return _cache[('leg', n)]


def _legvals(n, t):
    V = np.empty((n + 1, len(t))); D = np.zeros((n + 1, len(t))); V[0] = 1.0
    if n >= 1:
        V[1] = t; D[1] = 1.0
    for k in range(1, n):
        V[k + 1] = ((2 * k + 1) * t * V[k] - k * V[k - 1]) / (k + 1); D[k + 1] = D[k - 1] + (2 * k + 1) * V[k]
    return V, D


def _structure(n):
    """polygon-independent linear algebra for the flux of tensor degree n: divergence matrix and its null space"""
    if ('str', n) not in _cache:
        Dm = np.zeros((n + 1, n + 1))
        for k in range(1, n + 1):
            dc = Lg.legder([0] * k + [1]); Dm[:len(dc), k] = dc
        I = np.eye(n + 1); N1 = (n + 1) ** 2; C = np.zeros((N1, 2 * N1)); C[:, :N1] = np.kron(Dm, I); C[:, N1:] = np.kron(I, Dm)
        _cache[('str', n)] = (C, np.linalg.pinv(C), null_space(C, rcond=1e-10))
    return _cache[('str', n)]


def _quad(U, nq):
    if ('gl', nq) not in _cache:
        g, w = Lg.leggauss(nq); _cache[('gl', nq)] = ((g + 1) / 2, w / 2)
    g, w = _cache[('gl', nq)]; Xs = []; Ys = []; Ws = []
    for k in range(1, len(U) - 1):
        a, b, c = U[0], U[k], U[k + 1]; J = abs((b - a)[0] * (c - a)[1] - (b - a)[1] * (c - a)[0])
        S, T = np.meshgrid(g, g, indexing='ij')
        Xs.append((a[0] + S * (b - a)[0] + (1 - S) * T * (c - a)[0]).ravel()); Ys.append((a[1] + S * (b - a)[1] + (1 - S) * T * (c - a)[1]).ravel()); Ws.append((np.outer(w, w) * (1 - S) * J).ravel())
    return np.concatenate(Xs), np.concatenate(Ys), np.concatenate(Ws)


def float_stage(U, d, extra):
    """U: float vertices of P+ (ccw).  Returns (q coefficients (d+1, d+1), flux coefficients a, b (n+1, n+1)) in the
    tensor Legendre basis, or None."""
    U = np.asarray(U, float); nv = len(U); m = nv + d; n = m + extra
    x, y, wq = _quad(U, n + 14)
    B = np.ones_like(x); Bx = np.zeros_like(x); By = np.zeros_like(x)
    for k in range(nv):
        p, q = U[k], U[(k + 1) % nv]; nx, ny = -(q - p)[1], (q - p)[0]; l = nx * (x - p[0]) + ny * (y - p[1])
        Bx, By = Bx * l + B * nx, By * l + B * ny; B = B * l
    LX, DX = _legvals(n, x); LY, DY = _legvals(n, y)
    P0 = np.einsum('ip,jp->ijp', LX[:d + 1], LY[:d + 1]).reshape((d + 1) ** 2, -1)
    Px = np.einsum('ip,jp->ijp', DX[:d + 1], LY[:d + 1]).reshape((d + 1) ** 2, -1); Py = np.einsum('ip,jp->ijp', LX[:d + 1], DY[:d + 1]).reshape((d + 1) ** 2, -1)
    F = B * P0; Fx = Bx * P0 + B * Px; Fy = By * P0 + B * Py
    M = (F * wq) @ F.T; K = (Fx * wq) @ Fx.T + (Fy * wq) @ Fy.T
    # Ritz vector in a well-conditioned subspace: eigen-decomposition of the mass matrix, tiny directions dropped
    try:
        em, Qm = np.linalg.eigh(M); keep = em > 1e-11 * em.max(); T0 = Qm[:, keep] / np.sqrt(em[keep])
        ev, vec = np.linalg.eigh(T0.T @ K @ T0)
    except Exception:
        return None
    cu = T0 @ vec[:, 0]
    if not np.all(np.isfinite(cu)) or not (np.abs(cu).max() > 0):
        return None
    cu = cu / np.abs(cu).max()
    # Legendre coefficients of u = b q on [-1, 1]^2 (exact for polynomials: tensor Gauss rule)
    if ('gl1', m + 2) not in _cache:
        _cache[('gl1', m + 2)] = Lg.leggauss(m + 2)
    g, wg = _cache[('gl1', m + 2)]; GX, GY = np.meshgrid(g, g, indexing='ij'); gx = GX.ravel(); gy = GY.ravel(); Bg = np.ones_like(gx)
    for k in range(nv):
        p, q = U[k], U[(k + 1) % nv]; Bg = Bg * (-(q - p)[1] * (gx - p[0]) + (q - p)[0] * (gy - p[1]))
    lg, _ = _legvals(n, g); lgx, _ = _legvals(d, gx); lgy, _ = _legvals(d, gy)
    ug = (Bg * (cu @ np.einsum('ip,jp->ijp', lgx, lgy).reshape((d + 1) ** 2, -1))).reshape(len(g), len(g))
    T = (lg * wg) * ((2 * np.arange(n + 1) + 1) / 2)[:, None]; uc = T @ ug @ T.T; uc[m + 1:, :] = 0; uc[:, m + 1:] = 0
    C, Cp, Nn = _structure(n); N1 = (n + 1) ** 2
    Phi = np.einsum('ip,jp->ijp', LX, LY).reshape(N1, -1); Gm = (Phi * wq) @ Phi.T
    zp = Cp @ (-uc.ravel()); N_a = Nn[:N1]; N_b = Nn[N1:]
    H = N_a.T @ Gm @ N_a + N_b.T @ Gm @ N_b; rhs = N_a.T @ (Gm @ zp[:N1]) + N_b.T @ (Gm @ zp[N1:])
    yv = -np.linalg.lstsq(H, rhs, rcond=None)[0]; z = zp + Nn @ yv
    return cu.reshape(d + 1, d + 1), z[:N1].reshape(n + 1, n + 1), z[N1:].reshape(n + 1, n + 1)


def inflate(V):
    """frame and P+.  V: exact binary64 vertices (ccw, convex).  Returns (c (2 Fractions), s (Fraction), Uplus (list of
    pairs of Fractions on the grid), area of P (Fraction)) or None.  Containment of the image of P in P+ is checked exactly."""
    Vf = [(Fraction(float(a)), Fraction(float(b))) for a, b in V]; nv = len(Vf)
    for i in range(nv):                                            # P itself must be strictly convex and counter-clockwise
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
    cg = (sum(float(a) for a, b in U) / nv, sum(float(b) for a, b in U) / nv); g = 2 ** GRID; eps = EPS0
    for _ in range(4):
        Up = [(Fraction(round((cg[0] + (1 + eps) * (float(a) - cg[0])) * g), g), Fraction(round((cg[1] + (1 + eps) * (float(b) - cg[1])) * g), g)) for a, b in U]
        ok = True
        for i in range(nv):
            p, q, r = Up[i], Up[(i + 1) % nv], Up[(i + 2) % nv]
            if (q[0] - p[0]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[0] - q[0]) <= 0:
                ok = False; break                                 # P+ not strictly convex and counter-clockwise
            if any((q[0] - p[0]) * (v[1] - p[1]) - (q[1] - p[1]) * (v[0] - p[0]) < 0 for v in U):
                ok = False; break                                 # a vertex of P outside the half-plane of this edge of P+
        if ok:
            return c, s, Up, area
        eps *= 2
    return None


def _moments(Up, Dx, Dy):
    """exact m[i][j] = int_{P+} x^i y^j for i <= Dx, j <= Dy (Green: (1 / (i+1)) * boundary integral of x^(i+1) y^j dy)"""
    nv = len(Up); xp = []; yp = []; dys = []
    for k in range(nv):
        p, q = Up[k], Up[(k + 1) % nv]
        px = fmpq_poly([fmpq(p[0].numerator, p[0].denominator), fmpq((q[0] - p[0]).numerator, (q[0] - p[0]).denominator)])
        py = fmpq_poly([fmpq(p[1].numerator, p[1].denominator), fmpq((q[1] - p[1]).numerator, (q[1] - p[1]).denominator)])
        a = [fmpq_poly([1])]; b = [fmpq_poly([1])]
        for _ in range(Dx + 1):
            a.append(a[-1] * px)
        for _ in range(Dy):
            b.append(b[-1] * py)
        xp.append(a); yp.append(b); dys.append(fmpq((q[1] - p[1]).numerator, (q[1] - p[1]).denominator))
    one = fmpq(1); m = [[None] * (Dy + 1) for _ in range(Dx + 1)]
    for i in range(Dx + 1):
        for j in range(Dy + 1):
            tot = fmpq(0)
            for k in range(nv):
                if dys[k] != 0:
                    tot += dys[k] * (xp[k][i + 1] * yp[k][j]).integral()(one)
            m[i][j] = tot / (i + 1)
    return m


CT = fmpq_mpoly_ctx.get(('t',))
T_ = CT.gens()[0]


def integrate(f, Up):
    """exact integral of the polynomial f over P+ (ccw): Green's theorem with F = x-antiderivative of f,
    int_{P+} f = boundary integral of F dy; on each edge x, y are affine in t in [0, 1] and F(x(t), y(t)) is a polynomial
    in t (composition in exact arithmetic).  Equal to the moment formula (_moments / _integrate, kept for the tests)."""
    F = f.integral(0); tot = fmpq(0); nv = len(Up)
    for k in range(nv):
        p, q = Up[k], Up[(k + 1) % nv]; dy = q[1] - p[1]
        if dy == 0:
            continue
        xt = fmpq(p[0].numerator, p[0].denominator) + fmpq((q[0] - p[0]).numerator, (q[0] - p[0]).denominator) * T_
        yt = fmpq(p[1].numerator, p[1].denominator) + fmpq(dy.numerator, dy.denominator) * T_
        g = F.compose(xt, yt, ctx=CT); a = fmpq(0)
        for e, cf in g.to_dict().items():
            a += cf / (e[0] + 1)
        tot += fmpq(dy.numerator, dy.denominator) * a
    return tot


def _integrate(f, m):
    tot = fmpq(0)
    for (i, j), cf in f.to_dict().items():
        tot += cf * m[i][j]
    return tot


def _round_coeffs(A):
    """the table A rounded to multiples of 2^(e - CBITS), 2^(e-1) <= max |A| < 2^e, as exact rationals"""
    A = np.asarray(A, float); sc = np.abs(A).max()
    if not (sc > 0) or not np.isfinite(sc):
        return None
    den = Fraction(2) ** (CBITS - math.frexp(sc)[1]); out = []
    for row in A:
        r = []
        for v in row:
            f = Fraction(round(Fraction(float(v)) * den)) / den; r.append(fmpq(f.numerator, f.denominator))
        out.append(r)
    return out


def exact_stage(Up, qc, ac, bc):
    """Up: vertices of P+ (Fractions).  qc, ac, bc: coefficient tables (floats) in the tensor Legendre basis.
    Returns the exact rational lower bound of lambda_1(P+) (fmpq) and nu, or None if the hypotheses fail."""
    nv = len(Up); d = len(qc) - 1; n = len(ac) - 1; LXe, LYe = _leg_exact(max(d, n))
    b = CTX.from_dict({(0, 0): fmpq(1)})
    for k in range(nv):
        p, q = Up[k], Up[(k + 1) % nv]; nx = -(q[1] - p[1]); ny = q[0] - p[0]
        b = b * (fmpq(nx.numerator, nx.denominator) * (X - fmpq(p[0].numerator, p[0].denominator)) + fmpq(ny.numerator, ny.denominator) * (Y - fmpq(p[1].numerator, p[1].denominator)))
    qr = _round_coeffs(qc); ar = _round_coeffs(ac); br = _round_coeffs(bc)
    if qr is None or ar is None or br is None:
        return None
    zero = CTX.from_dict({})

    def build(tab, deg):
        f = zero
        for i in range(deg + 1):
            row = zero
            for j in range(deg + 1):
                if tab[i][j] != 0:
                    row = row + tab[i][j] * LYe[j]
            f = f + LXe[i] * row
        return f
    u = b * build(qr, d)
    w1 = build(ar, n); w2 = build(br, n)
    r = -u - w1.derivative(0) - w2.derivative(1); w1 = w1 + r.integral(0)
    if not (-(w1.derivative(0) + w2.derivative(1)) - u).is_zero():
        return None                                                # -div w = u must hold exactly
    ux = u.derivative(0); uy = u.derivative(1); f0 = ux * ux + uy * uy; f1 = u * u; f2 = w1 * w1 + w2 * w2
    A0 = integrate(f0, Up); A1 = integrate(f1, Up); A2 = integrate(f2, Up)
    if not (A1 > 0):
        return None
    # nu: second eigenvalue of a bounding rectangle, the best of the edge directions and the axes
    best = fmpq(0)
    dirs = [(Fraction(1), Fraction(0))] + [(Up[(k + 1) % nv][0] - Up[k][0], Up[(k + 1) % nv][1] - Up[k][1]) for k in range(nv)]
    for dx, dy in dirs:
        a1 = [dx * p[0] + dy * p[1] for p in Up]; a2 = [-dy * p[0] + dx * p[1] for p in Up]; e1 = max(a1) - min(a1); e2 = max(a2) - min(a2); dd = dx * dx + dy * dy
        v = dd * min(4 / (e1 * e1) + 1 / (e2 * e2), 1 / (e1 * e1) + 4 / (e2 * e2)); v = fmpq(v.numerator, v.denominator)
        if v > best:
            best = v
    nu = PI2_LO * best
    ar = sum(Up[i][0] * Up[(i + 1) % nv][1] - Up[(i + 1) % nv][0] * Up[i][1] for i in range(nv)) / 2          # |P+|, exact
    ks = KS_LO / fmpq(ar.numerator, ar.denominator)                # Krahn-Szego: lambda_2(Omega) >= 2 pi j_{0,1}^2 / |Omega|
    if ks > nu:
        nu = ks
    D = A0 - nu * A1
    if not (D < 0):
        return None
    G = A0 - 2 * nu * A1 + nu * nu * A2
    return nu * (-D) / (G - D), nu, A0 / A1


SLACK = 1e-6                 # the stated F_lo is the proved bound times (1 - SLACK): a re-check with its own flux passes


def F_lower(V, d=4, extra=1):
    """Rigorous lower bound of F(P) = lambda_1(P) |P| for a convex ccw polygon with exact binary64 vertices.
    Returns a dict (F_lo float rounded down, and the certificate data) or None."""
    fr = inflate(V)
    if fr is None:
        return None
    c, s, Up, area = fr
    fs = float_stage([(float(a), float(b)) for a, b in Up], d, extra)
    if fs is None:
        return None
    qc, ac, bc = fs
    ex = exact_stage(Up, qc, ac, bc)
    if ex is None:
        return None
    lam, nu, rho = ex
    Flo = lam * fmpq(area.numerator, area.denominator) / fmpq((s * s).numerator, (s * s).denominator)
    v = float(Fraction(int(Flo.p), int(Flo.q))) * (1 - SLACK); v = math.nextafter(v, -math.inf)
    return dict(F_lo=v, d=d, extra=extra, grid=GRID, cbits=CBITS, c=[str(c[0]), str(c[1])], s=str(s), Uplus=[[str(a), str(b)] for a, b in Up],
                q=np.asarray(qc).tolist(), nu=float(Fraction(int(nu.p), int(nu.q))), rho_plus=float(Fraction(int(rho.p), int(rho.q))))


def recheck(V, rec):
    """re-proves a stored record: same P+ and q, the flux is recomputed.  True if the proved bound is >= rec['F_lo']."""
    fr = inflate(V)
    if fr is None:
        return False
    c, s, Up0, area = fr; Up = [(Fraction(a), Fraction(b)) for a, b in rec['Uplus']]; nv = len(Up)
    Vf = [((Fraction(float(a)) - Fraction(rec['c'][0])) / Fraction(rec['s']), (Fraction(float(b)) - Fraction(rec['c'][1])) / Fraction(rec['s'])) for a, b in V]
    for i in range(nv):                                            # P+ strictly convex, ccw, and contains the image of P
        p, q, r = Up[i], Up[(i + 1) % nv], Up[(i + 2) % nv]
        if (q[0] - p[0]) * (r[1] - q[1]) - (q[1] - p[1]) * (r[0] - q[0]) <= 0 or any((q[0] - p[0]) * (v[1] - p[1]) - (q[1] - p[1]) * (v[0] - p[0]) < 0 for v in Vf):
            return False
    fs = float_stage([(float(a), float(b)) for a, b in Up], rec['d'], rec['extra'])
    if fs is None:
        return False
    ex = exact_stage(Up, rec['q'], fs[1], fs[2])
    if ex is None:
        return False
    s2 = Fraction(rec['s']) ** 2; Flo = ex[0] * fmpq(area.numerator, area.denominator) / fmpq(s2.numerator, s2.denominator)
    return Fraction(int(Flo.p), int(Flo.q)) >= Fraction(float(rec['F_lo']))
