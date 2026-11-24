#!/usr/bin/env python3
"""Standalone checker of the final positivity proof of Theorem L (far/LOCAL_LEMMA.md section 7, far/LOCAL_CERT.md).

    python3 lglocal_cert_check.py CERT_DIR [--leaves K] [--seed S]

Uses only the Python standard library (exact integers / fractions); imports nothing of the producer.  Input (plain or .xz):
  CERT_DIR/witness.json   Tm, Em (maps "a0,...,a5" -> "p/q"), rho_in
  CERT_DIR/leaves.jsonl   first line {"t": "root", R, W, rho_in, keep}; then one leaf per line (t = I, O or P, m, r, witness)
Claim checked ("level A"):  IF the true polynomial Q(d) = sum_{|alpha| >= 2} s_alpha d^alpha has s_alpha in
Tm_alpha +- Em_alpha for every alpha (Tm_alpha = 0 if absent), THEN Q(d) > 0 for every d != 0 with |d|_2 < R.
  1. inner cube |d_i| <= rho_in: tau and the 64 vertex matrices H/2 + sum_i (+-rho_in) T_i - tau I are positive definite
     (exact LDL^T);
  2. the leaves tile the half cube {|d_i| <= W, d_2 >= 0} exactly (inside, pairwise interior-disjoint, volumes add up);
  3. every I leaf lies in the inner cube, every O leaf has sum_i max(|m_i| - r_i, 0)^2 >= R^2, and on every P leaf
       p(x) = sum_{2 <= |alpha| <= keep} Tm_alpha (m + x)^alpha  satisfies  min_{|x_i| <= r_i} p(x) > err,
       err = sum_alpha E'_alpha prod_i (|m_i| + r_i)^alpha_i,  E' = Em + |Tm_alpha| for |alpha| > keep,
     by the bound  p(x) >= c0 + g.x + x^T H x / 2 - tail,  tail = sum_{|beta| >= 3} |p_beta| r^beta,  and for the quadratic
     part the Lagrangian bound with the stored multipliers mu >= 0 (A = H + 2 diag(mu) positive definite by exact LDL^T,
     value c0 - g^T A^-1 g / 2 - sum mu_i r_i^2) or the crude bound c0 - sum |g_i| r_i - sum |H_ij| r_i r_j / 2;
  4. the half d_2 < 0 is the mirror image (not checked here: far/LOCAL_LEMMA.md section 7).
--leaves K: check only K random P leaves (items 1, 2 and the I / O leaves are always checked in full)."""
import sys, os, json, itertools, random, time, lzma
from fractions import Fraction as F

ND = 6


def ldl_pd(A):
    """pivots of the LDL^T factorisation of the symmetric matrix A (Fractions), or None if one is not positive"""
    n = len(A); L = [[F(0)] * n for _ in range(n)]; D = []
    for j in range(n):
        dj = A[j][j] - sum(L[j][k] * L[j][k] * D[k] for k in range(j))
        if dj <= 0:
            return None
        D.append(dj)
        for i in range(j + 1, n):
            L[i][j] = (A[i][j] - sum(L[i][k] * L[j][k] * D[k] for k in range(j))) / dj
    return D, L


def inner(Tm, Em, rho):
    H = [[F(0)] * ND for _ in range(ND)]; T = [[[F(0)] * ND for _ in range(ND)] for _ in range(ND)]
    for a, v in Tm.items():
        ix = [i for i in range(ND) for _ in range(a[i])]
        if len(ix) == 2:
            if ix[0] == ix[1]:
                H[ix[0]][ix[0]] += 2 * v
            else:
                H[ix[0]][ix[1]] += v; H[ix[1]][ix[0]] += v
        elif len(ix) == 3:
            ps = set(itertools.permutations(ix))
            for p in ps:
                T[p[0]][p[1]][p[2]] += v / len(ps)
    w = lambda a: F(1) if max(a) == sum(a) else F(1, 2)
    e2 = lambda i, j: Em.get(tuple((2 if i == j else 1) if k in (i, j) else 0 for k in range(ND)), F(0))
    tau2 = max(e2(i, i) + sum(e2(i, j) for j in range(ND) if j != i) / 2 for i in range(ND))
    tau3 = rho * sum(v * w(a) for a, v in Em.items() if sum(a) == 3)
    tau4 = sum(rho ** (sum(a) - 2) * (abs(Tm.get(a, F(0))) + Em.get(a, F(0))) * w(a) for a in set(Tm) | set(Em) if sum(a) >= 4)
    tau = tau2 + tau3 + tau4; mn = None
    for sg in itertools.product((-1, 1), repeat=ND):
        M = [[H[i][j] / 2 + rho * sum(sg[k] * T[k][i][j] for k in range(ND)) - (tau if i == j else 0) for j in range(ND)] for i in range(ND)]
        r = ldl_pd(M)
        if r is None:
            return False, dict(tau=tau, vertex=sg)
        mn = min(r[0]) if mn is None else min(mn, min(r[0]))
    return True, dict(tau=tau, tau2=tau2, tau3=tau3, tau4=tau4, min_pivot=mn, H=H)


def shift(P, m):
    """the coefficients of x -> P(m + x)  (P: dict exponent tuple -> Fraction)"""
    for i in range(ND):
        if m[i] == 0:
            continue
        Q = {}; pw = [F(1)]
        for a, v in P.items():
            ai = a[i]
            while len(pw) <= ai:
                pw.append(pw[-1] * m[i])
            c = 1
            for k in range(ai, -1, -1):                              # binomial(ai, k) m^(ai - k) x^k
                b = a[:i] + (k,) + a[i + 1:]; Q[b] = Q.get(b, 0) + v * c * pw[ai - k]
                c = c * k // (ai - k + 1)
        P = Q
    return P


def leaf_positive(P0, Eterms, m, r, wit, mu):
    P = shift(P0, m); c = [abs(m[i]) + r[i] for i in range(ND)]
    pw = [[F(1)] for _ in range(ND)]
    for i in range(ND):
        for _ in range(Eterms['maxdeg']):
            pw[i].append(pw[i][-1] * c[i])
    err = F(0)
    for a, v in Eterms['E'].items():
        t = v
        for i in range(ND):
            if a[i]:
                t *= pw[i][a[i]]
        err += t
    c0 = P.get((0,) * ND, F(0)) - err; g = [F(0)] * ND; H = [[F(0)] * ND for _ in range(ND)]; tail = F(0)
    rp = [[F(1)] for _ in range(ND)]
    for a, v in P.items():
        d = sum(a)
        if d == 1:
            g[a.index(1)] = v
        elif d == 2:
            ix = [i for i in range(ND) if a[i]]
            if len(ix) == 1:
                H[ix[0]][ix[0]] = 2 * v
            else:
                H[ix[0]][ix[1]] = H[ix[1]][ix[0]] = v
        elif d >= 3:
            t = abs(v)
            for i in range(ND):
                while len(rp[i]) <= a[i]:
                    rp[i].append(rp[i][-1] * r[i])
                if a[i]:
                    t *= rp[i][a[i]]
            tail += t
    crude = c0 - sum(abs(g[i]) * r[i] for i in range(ND)) - sum(abs(H[i][j]) * r[i] * r[j] for i in range(ND) for j in range(ND)) / 2
    best = crude
    if wit == 'dual' and mu is not None and all(x >= 0 for x in mu):
        A = [[H[i][j] + (2 * mu[i] if i == j else 0) for j in range(ND)] for i in range(ND)]; f = ldl_pd(A)
        if f is not None:
            D, L = f; y = [None] * ND
            for i in range(ND):
                y[i] = g[i] - sum(L[i][k] * y[k] for k in range(i))
            best = max(best, c0 - sum(y[i] * y[i] / D[i] for i in range(ND)) / 2 - sum(mu[i] * r[i] * r[i] for i in range(ND)))
    return best - tail > 0, float(best - tail), float(err)


def main():
    cert = sys.argv[1]; K = int(sys.argv[sys.argv.index('--leaves') + 1]) if '--leaves' in sys.argv else -1
    seed = int(sys.argv[sys.argv.index('--seed') + 1]) if '--seed' in sys.argv else 1; t0 = time.time(); bad = []
    op = lambda name: open(os.path.join(cert, name)) if os.path.exists(os.path.join(cert, name)) else lzma.open(os.path.join(cert, name + '.xz'), 'rt')   # plain or .xz
    wj = json.load(op('witness.json')); pm = lambda d: {tuple(int(x) for x in k.split(',')): F(v) for k, v in d.items()}
    Tm = pm(wj['Tm']); Em = pm(wj['Em'])
    if not (all(len(a) == ND and sum(a) >= 2 for a in set(Tm) | set(Em)) and all(v >= 0 for v in Em.values())):
        bad.append('witness: exponents of degree < 2 or negative error bounds')
    lines = [json.loads(l) for l in op('leaves.jsonl')]; root = lines[0]; leaves = lines[1:]
    assert root['t'] == 'root'; R2 = F(root['R']) ** 2; W = F(root['W']); rho = F(root['rho_in']); keep = int(root['keep'])
    if F(wj['rho_in']) != rho:
        bad.append('rho_in differs between witness.json and leaves.jsonl')
    print('witness: %d coefficients, %d error bounds; R = %s (%.9f), W = %s, rho_in = %s, keep = %d; %d leaves' % (len(Tm), len(Em), root['R'], float(F(root['R'])), W, rho, keep, len(leaves)), flush=True)
    if not (W * W > R2 and rho < W):
        bad.append('the cube of half-width W does not contain the ball of radius R')
    # 1. inner cube
    ok, inf = inner(Tm, Em, rho)
    if ok:
        sc = inf['H']; print('1. inner cube: positive definite at the 64 vertices; tau = %.6e (tau2 %.6e, tau3 %.6e, tau>=4 %.6e), smallest pivot %.6e'
                             % (float(inf['tau']), float(inf['tau2']), float(inf['tau3']), float(inf['tau4']), float(inf['min_pivot'])), flush=True)
    else:
        bad.append('inner cube not positive definite at vertex %s' % (inf['vertex'],)); print('1. inner cube: FAILED', flush=True)
    # 2. tiling of the half cube
    boxes = []
    for d in leaves:
        m = [F(x) for x in d['m']]; r = [F(x) for x in d['r']]
        if d['t'] not in ('I', 'O', 'P') or any(x <= 0 for x in r):
            bad.append('leaf with an unknown label or an empty box'); continue
        boxes.append(([m[i] - r[i] for i in range(ND)], [m[i] + r[i] for i in range(ND)], d, m, r))
    rlo = [-W, -W, F(0), -W, -W, -W]; rhi = [W] * ND; vol = F(0)
    for lo, hi, d, m, r in boxes:
        if not all(rlo[i] <= lo[i] and hi[i] <= rhi[i] for i in range(ND)):
            bad.append('leaf outside the half cube'); break
        v = F(1)
        for i in range(ND):
            v *= hi[i] - lo[i]
        vol += v
    rv = F(1)
    for i in range(ND):
        rv *= rhi[i] - rlo[i]
    if vol != rv:
        bad.append('the volumes of the leaves do not add up to the half cube (ratio %.9f)' % float(vol / rv))
    order = sorted(range(len(boxes)), key=lambda i: boxes[i][0][0]); act = []; overlap = False
    for i in order:
        lo, hi = boxes[i][0], boxes[i][1]; act = [k for k in act if boxes[k][1][0] > lo[0]]
        for k in act:
            if all(lo[t] < boxes[k][1][t] and boxes[k][0][t] < hi[t] for t in range(ND)):
                overlap = True; break
        if overlap:
            bad.append('overlapping leaves'); break
        act.append(i)
    cnt = {t: sum(1 for b in boxes if b[2]['t'] == t) for t in 'IOP'}
    print('2. tiling: %d leaves (I %d, O %d, P %d); inside the half cube, pairwise interior-disjoint, volumes add up: %s  (%.0fs)' % (len(boxes), cnt['I'], cnt['O'], cnt['P'], 'ok' if not bad else 'see below', time.time() - t0), flush=True)
    # 3. the leaves
    P0 = {a: v for a, v in Tm.items() if sum(a) <= keep}; E = dict(Em)
    for a, v in Tm.items():
        if sum(a) > keep:
            E[a] = E.get(a, F(0)) + abs(v)
    Eterms = dict(E={a: v for a, v in E.items() if v != 0}, maxdeg=max(sum(a) for a in E))
    pl = [b for b in boxes if b[2]['t'] == 'P']
    if K >= 0:
        random.Random(seed).shuffle(pl); pl = pl[:K]
    nI = nO = 0
    for lo, hi, d, m, r in boxes:
        if d['t'] == 'I':
            nI += 1
            if not all(abs(m[i]) + r[i] <= rho for i in range(ND)):
                bad.append('I leaf not inside the inner cube')
        elif d['t'] == 'O':
            nO += 1
            if not sum(max(abs(m[i]) - r[i], 0) ** 2 for i in range(ND)) >= R2:
                bad.append('O leaf not outside the ball')
    worst = None; nP = 0; ncrude = 0
    for lo, hi, d, m, r in pl:
        mu = [F(x) for x in d['mu']] if d.get('w') == 'dual' else None
        ok, margin, err = leaf_positive(P0, Eterms, m, r, d.get('w'), mu); nP += 1; ncrude += d.get('w') != 'dual'
        if not ok:
            bad.append('P leaf not proved: m = %s r = %s' % (d['m'], d['r']))
            if len(bad) > 20:
                break
        worst = margin if worst is None else min(worst, margin)
        if nP % 500 == 0:
            print('   %d / %d P leaves, %.0fs' % (nP, len(pl), time.time() - t0), flush=True)
    print('3. leaves: I %d and O %d checked; P leaves proved %d of %d checked (%s; crude bound %d); smallest margin %.3e  (%.0fs)'
          % (nI, nO, nP - sum(1 for b in bad if b.startswith('P leaf')), nP, 'all' if K < 0 else 'random sample', ncrude, worst if worst is not None else float('nan'), time.time() - t0), flush=True)
    for b in bad[:10]:
        print('  PROBLEM:', b)
    full = K < 0
    print(('CERTIFICATE OK (level A, all leaves)' if full else 'CERTIFICATE OK (level A, sample of %d P leaves)' % nP) if not bad else 'CERTIFICATE FAILED: %d problems' % len(bad))
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
