#!/usr/bin/env python3
"""Claude's pointwise check of GPT's PENT-24 LG witnesses (written from LG_LEMMA.md / FAR_THEOREM.md, NOT from GPT's
code; GPT's programs are never imported or run).  A witness is read as data only: box, chart, cut index j, degree N,
target L, and the integer coefficient lists of the affine trial family u = u0 + sum x_i u_i, chi = chi0 + sum x_i chi_i
(x in [-1, 1]^d the normalised box coordinate).  The coefficient order is GPT's documented layout (triangles t, Bernstein
triples (i, j, l) with i = exponent at O ascending, then j ascending; spoke / edge / interior / centre keys in order of
first appearance; u has zeros on the boundary edge, chi has chi(O) = 0).

At a rational point x of the box this program, independently:
  - builds the polygon from the chart formulas (FAR_THEOREM.md sections 2-3; near: (0,0), (1,0), (x0,x1), (x2,x3), (x4,x5)),
  - takes O = the vertex average of the polygon at the box centre,
  - builds u and chi on every fan triangle as polynomials in the physical coordinates (X, Y) = p - O,
  - checks: u = 0 on the boundary edges, u and chi continuous across the spokes, -div w = u for
    w = -(X, Y) h + curl chi with h = sum_d u_d / (d + 2) (homogeneous parts about O), w.n continuous across the spokes,
  - integrates |grad u|^2, u^2, |w|^2 exactly over every triangle,
  - and evaluates the Lehmann-Goerisch bound with nu = KS/|P|:  F_lb = |P| Dc / (KS A2 - |P| A1),
    Dc = KS A1 - |P| A0;  required: all fan determinants > 0, Dc > 0, KS A2 - |P| A1 > 0, F_lb > L.

    python3 -I lg_point_check.py WITNESS.jsonl [--points K] [--limit M] [--seed S]
"""
import sys, json, random
from fractions import Fraction as Q
from math import factorial

KS = Q('36.33682907')


# ---- bivariate polynomials {(p, q): Fraction}
def padd(*ps):
    r = {}
    for p in ps:
        for k, v in p.items():
            r[k] = r.get(k, 0) + v
    return {k: v for k, v in r.items() if v}


def pscale(p, c):
    return {k: v * c for k, v in p.items() if v * c}


def pmul(p, q):
    r = {}
    for (a, b), v in p.items():
        for (c, d), w in q.items():
            r[a + c, b + d] = r.get((a + c, b + d), 0) + v * w
    return {k: v for k, v in r.items() if v}


def ppow(p, n, cache):
    key = (id(p), n)
    if key not in cache:
        r = {(0, 0): Q(1)}
        for _ in range(n):
            r = pmul(r, p)
        cache[key] = r
    return cache[key]


def pdx(p):
    return {(a - 1, b): v * a for (a, b), v in p.items() if a}


def pdy(p):
    return {(a, b - 1): v * b for (a, b), v in p.items() if b}


def peval(p, x, y):
    return sum(v * x ** a * y ** b for (a, b), v in p.items())


def integrate(p, e1, e2):
    """integral over the triangle (0, e1, e2) of p(X, Y): substitute X = a e1x + b e2x, Y = a e1y + b e2y and use
    int_ref a^m b^n = m! n! / (m + n + 2)!  times |det|."""
    det = e1[0] * e2[1] - e1[1] * e2[0]
    X = {(1, 0): e1[0], (0, 1): e2[0]}; Y = {(1, 0): e1[1], (0, 1): e2[1]}
    cache = {}; tot = Q(0)
    for (a, b), v in p.items():
        s = pmul(ppow(X, a, cache), ppow(Y, b, cache))
        for (m, n), w in s.items():
            tot += v * w * Q(factorial(m) * factorial(n), factorial(m + n + 2))
    return tot * abs(det)


# ---- polygon from the chart
def polygon(chart, j, z):
    if chart == 'near':
        return [(Q(0), Q(0)), (Q(1), Q(0)), (z[0], z[1]), (z[2], z[3]), (z[4], z[5])]
    if chart in ('std', 'qstd'):
        q2, q3 = (z[0], z[1]), (z[2], z[3])
    elif chart in ('E2', 'qE2'):
        q3 = (z[0], z[1]); a, h = z[2], z[3]; q2 = (1 + a * (q3[0] - 1), a * q3[1] + h)
    elif chart in ('E3', 'qE3'):
        q2 = (z[0], z[1]); b, h = z[2], z[3]; q3 = (b * q2[0], b * q2[1] + h)
    else:
        raise ValueError(chart)
    q0, q1 = (Q(0), Q(0)), (Q(1), Q(0))
    if chart.startswith('q'):
        return [q0, q1, q2, q3]
    t, u = z[4], z[5]
    if j == 3:
        return [q0, q1, q2, (q3[0] + t * (q2[0] - q3[0]), q3[1] + t * (q2[1] - q3[1])), ((1 - u) * q3[0], (1 - u) * q3[1])]
    if j == 0:
        return [(t, Q(0)), q1, q2, q3, (u * q3[0], u * q3[1])]
    raise ValueError(j)


# ---- GPT's coefficient layout (data format)
def layout(n, k, bdy):
    ids = {}; rows = []
    for t in range(k):
        row = []
        for i in range(n + 1):
            for jj in range(n - i + 1):
                l = n - i - jj
                if bdy and i == 0:
                    row.append(None); continue
                if i == n:
                    if not bdy:
                        row.append(None); continue
                    key = ('c',)
                elif jj == 0:
                    key = ('s', (t + 1) % k, i)
                elif l == 0:
                    key = ('s', t, i)
                elif i == 0:
                    key = ('e', t, jj)
                else:
                    key = ('i', t, i, jj)
                ids.setdefault(key, len(ids)); row.append((i, jj, l, ids[key]))
        rows.append(row)
    return len(ids), rows


def bernstein_poly(n, row, coef, e1, e2):
    """sum over the nodes of the triangle of coef * n!/(i! j! l!) lam^i a^j b^l in (X, Y); a, b barycentric of e1, e2."""
    det = e1[0] * e2[1] - e1[1] * e2[0]
    a = {(1, 0): e2[1] / det, (0, 1): -e2[0] / det}
    b = {(1, 0): -e1[1] / det, (0, 1): e1[0] / det}
    lam = padd({(0, 0): Q(1)}, pscale(a, -1), pscale(b, -1))
    cache = {}; r = {}
    for node in row:
        if node is None:
            continue
        i, jj, l, g = node
        c = coef[g]
        if not c:
            continue
        m = Q(factorial(n), factorial(i) * factorial(jj) * factorial(l)) * c
        r = padd(r, pscale(pmul(pmul(ppow(lam, i, cache), ppow(a, jj, cache)), ppow(b, l, cache)), m))
    return r


def check_point(w, x):
    lo = [Q(v) for v in w['lo']]; hi = [Q(v) for v in w['hi']]; d = len(lo)
    mid = [(a + b) / 2 for a, b in zip(lo, hi)]; half = [(b - a) / 2 for a, b in zip(lo, hi)]
    z = [m + h * xi for m, h, xi in zip(mid, half, x)]
    V = polygon(w['chart'], w['j'], z); V0 = polygon(w['chart'], w['j'], mid); k = len(V)
    O = (sum(v[0] for v in V0) / k, sum(v[1] for v in V0) / k)
    E = [(v[0] - O[0], v[1] - O[1]) for v in V]
    dets = [E[t][0] * E[(t + 1) % k][1] - E[t][1] * E[(t + 1) % k][0] for t in range(k)]
    if min(dets) <= 0:
        return 'DET', None
    n = w['N']
    nu_, urows = layout(n, k, True); nc_, crows = layout(n + 1, k, False)
    fam_u = w['u']; fam_c = w['chi']
    if any(len(r) != nu_ for r in fam_u) or any(len(r) != nc_ for r in fam_c):
        return 'LAYOUT', None
    xs = [Q(1)] + list(x[:len(fam_u) - 1])
    ucoef = [sum(xi * r[g] for xi, r in zip(xs, fam_u)) for g in range(nu_)]
    ccoef = [sum(xi * r[g] for xi, r in zip(xs, fam_c)) for g in range(nc_)]
    A0 = A1 = A2 = Q(0); U = []; W = []; CH = []
    for t in range(k):
        e1, e2 = E[t], E[(t + 1) % k]
        u = bernstein_poly(n, urows[t], ucoef, e1, e2); ch = bernstein_poly(n + 1, crows[t], ccoef, e1, e2)
        h = {kk: v / (kk[0] + kk[1] + 2) for kk, v in u.items()}
        wx = padd(pscale(pmul({(1, 0): Q(1)}, h), -1), pdy(ch))
        wy = padd(pscale(pmul({(0, 1): Q(1)}, h), -1), pscale(pdx(ch), -1))
        if padd(pdx(wx), pdy(wy), u):          # -div w = u  <=>  div w + u = 0
            return 'DIV', None
        ux, uy = pdx(u), pdy(u)
        A0 += integrate(padd(pmul(ux, ux), pmul(uy, uy)), e1, e2)
        A1 += integrate(pmul(u, u), e1, e2)
        A2 += integrate(padd(pmul(wx, wx), pmul(wy, wy)), e1, e2)
        U.append(u); W.append((wx, wy)); CH.append(ch)
    # boundary and continuity, by evaluation at n + 3 points of each segment (polynomials of degree <= n + 2 on a line)
    for t in range(k):
        e1, e2 = E[t], E[(t + 1) % k]
        for s in range(n + 3):
            r = Q(s + 1, n + 4)
            pb = (e1[0] + r * (e2[0] - e1[0]), e1[1] + r * (e2[1] - e1[1]))       # boundary edge
            if peval(U[t], *pb):
                return 'BOUNDARY', None
            ps = (r * e2[0], r * e2[1])                                         # spoke t+1, shared with t+1
            tn = (t + 1) % k
            if peval(U[t], *ps) != peval(U[tn], *ps) or peval(CH[t], *ps) != peval(CH[tn], *ps):
                return 'CONTINUITY', None
            nx, ny = e2[1], -e2[0]
            if (peval(W[t][0], *ps) - peval(W[tn][0], *ps)) * nx + (peval(W[t][1], *ps) - peval(W[tn][1], *ps)) * ny:
                return 'FLUX', None
    area = sum(dets) / 2
    Dc = KS * A1 - area * A0; den = KS * A2 - area * A1
    if Dc <= 0 or den <= 0:
        return 'LG_SIGN', None
    F = area * Dc / den
    return ('OK' if F > Q(w['L']) else 'BELOW'), F


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
    opt = dict(a[2:].split('=') for a in sys.argv[1:] if a.startswith('--'))
    K = int(opt.get('points', 3)); M = int(opt.get('limit', 0)); rng = random.Random(int(opt.get('seed', 20261122)))
    tot = {}; worst = None; n = 0
    for line in open(args[0]):
        r = json.loads(line); w = r.get('witness')
        if not w:
            continue
        n += 1
        if M and n > M:
            break
        d = len(w['lo'])
        pts = [[Q(rng.choice((-1, 1))) for _ in range(d)]]                       # one random corner
        pts += [[Q(rng.randint(-2 ** 20, 2 ** 20), 2 ** 20) for _ in range(d)] for _ in range(K - 1)]
        for x in pts:
            st, F = check_point(w, x)
            tot[st] = tot.get(st, 0) + 1
            if F is not None:
                m = (F - Q(w['L'])) / Q(w['L'])
                worst = m if worst is None or m < worst else worst
            if st != 'OK':
                print('FAIL', r.get('id'), st, [str(v) for v in x], flush=True)
        if n % 25 == 0:
            print('progress', n, tot, 'smallest relative margin %.3e' % float(worst), flush=True)
    print('RESULT witnesses', n if not M else min(n, M), 'points', tot, 'smallest relative margin %.3e' % float(worst))


if __name__ == '__main__':
    main()
