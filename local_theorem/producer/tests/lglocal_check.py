#!/usr/bin/env python3
"""tests of lglocal.py (the local theorem at the regular pentagon; far/LOCAL_LEMMA.md).  Run from v2_run:
    python3 tests/lglocal_check.py [quick]
quick: degree 8 only (the proof itself needs degree 16; about 2 minutes more)."""
import os, sys, random, itertools
from fractions import Fraction
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
import lglocal as LL, lgpoly as LP          # noqa: E402
from flint import fmpq, arb                 # noqa: E402

FAIL = []


def check(name, cond, info=''):
    print('%-100s %s %s' % (name, 'ok' if cond else 'FAIL', info), flush=True)
    if not cond:
        FAIL.append(name)


# 1. the algebra of the residual form (float, arbitrary polynomial error)
w = LL.selftest(6, verbose=False); check('residual form: exact identity with an arbitrary error function (degree 6)', w < 1e-9, '%.1e' % w)
w = LL.selftest(5, seed=7, verbose=False); check('residual form: the same, degree 5, another error function', w < 1e-9, '%.1e' % w)

# 2. polynomial balls: products contain the products of perturbed polynomials
rnd = random.Random(3)


def rpoly(nt, pos=False):
    return LL.CX.from_dict({tuple(rnd.randint(0, 2) for _ in range(LL.ND)): fmpq(rnd.randint(0 if pos else -9, 9), rnd.randint(1, 5)) for _ in range(nt)})


ok = True
for trial in range(20):
    a = LL.PB(rpoly(6), rpoly(3, True)); b = LL.PB(rpoly(5), rpoly(3, True)); c = LL.PB(rpoly(4))
    pa = LL.CX.from_dict({k: v + fmpq(rnd.choice((-1, 1))) * LL.pdict(a.r).get(k, fmpq(0)) for k, v in {**{k: fmpq(0) for k in LL.pdict(a.r)}, **LL.pdict(a.m)}.items()})
    pb = LL.CX.from_dict({k: v + fmpq(rnd.choice((-1, 1))) * LL.pdict(b.r).get(k, fmpq(0)) for k, v in {**{k: fmpq(0) for k in LL.pdict(b.r)}, **LL.pdict(b.m)}.items()})
    for X, px in ((a * b + c, pa * pb + c.m), (a * b - a * c, pa * pb - pa * c.m), ((a - b) * (a + b), (pa - pb) * (pa + pb))):
        d = LL.pdict(px - X.m); r = LL.pdict(X.r) if X.r is not None else {}
        ok = ok and all(abs(v) <= r.get(k, fmpq(0)) for k, v in d.items())
check('polynomial balls: sums and products contain the exact results for extreme members', ok)
x = arb(2).sqrt() / 3; m, r = LL.ball(x); check('ball(): rational centre and radius enclose the ball', abs(x - LL.A(m)) < LL.A(r) * 1.0000001 + arb(10) ** -100 and LL.up(x) >= m and LL.lo(x) <= m)

# 3. the rigorous data at degree 8 (too coarse for the theorem, exercises every part)
fam, ce, rng = LL._test_family(8, 1); X = LL.Local(8, test_list=ce, verbose=False)
GR = [np.array([[LL._f(v) for v in row] for row in g]) for g in X.GR]
check('Gram matrices symmetric', all(np.allclose(g, g.T, rtol=1e-12, atol=1e-30) for g in GR))
def psd(g):
    i = np.nonzero(np.diag(g) > 0)[0]; h = g[np.ix_(i, i)]; s_ = np.sqrt(np.diag(h)); return np.linalg.eigvalsh(h / np.outer(s_, s_))[0] > -1e-8


check('Gram matrices positive semidefinite (float, after scaling)', all(psd(g) for g in GR))
d0 = X.debug(seed=5, size=0.0); d1 = X.debug(seed=5, size=1e-2)
check('computed polynomial K D equals the pointwise formula (centre and |d| ~ 0.01), relative to the scale', max(d0, d1) < 1e3, 'see the lines above (float noise of the pointwise formula about 1e-12)')
check('symmetry: relabelling maps fix c, mean differential has norm < 1', X.symmetry_check(verbose=False))
check('degree 8: F0 enclosure contains the degree-16 enclosure [18.9191028, 18.9191048]', X.Flo < Fraction(189191028, 10 ** 7) and X.Fup > Fraction(189191048, 10 ** 7))
check('degree 8: inner box of half-width 0.01 NOT proved (error of u0 too large) - the test can fail', not X.certify_inner(Fraction(1, 100), verbose=False))
check('constant and linear coefficients are not used', all(sum(a) >= 2 for a in X.Tm) and all(sum(a) >= 2 for a in X.Em))
check('error bounds non-negative, tight bounds for degrees 0 to 3', all(v >= 0 for v in X.Em.values()) and all(sum(a) <= 3 for a in X.Elin))
check('degree 8: the computed coefficients of degree 0 and 1 lie within their error bounds (the true ones are zero)', X.low_check(verbose=False))
tb = X.tb; q8 = X.q
ray = all(LL._f(tb[k]['qaa'][b][b] + 2 * tb[k]['qab'][b][b] + tb[k]['qbb'][b][b]) <= 0.25 * LL._f(tb[k]['m'][b][b]) * (1 + 1e-12) for k in range(LL.NV) for b in range(q8 + 1))
check('ray bound: int (a + b)^2 h^2 <= (1/4) int u^2 on every triangle, for the 7 lists and the test function', ray)
check('tight bounds of the linear part do not exceed the crude majorant', all(X.Elin[a] <= X.parts['crude'].get(a, 0) * (1 + Fraction(1, 10 ** 6)) for a in X.Elin))
H, T3 = X.quadratic(); ev = np.linalg.eigvalsh(np.array([[float(v) for v in row] for row in H]) / X.scale)
check('Hessian of the computed part: eigenvalues about 0.53, 0.56, 5.03, 5.69, 12.8, 15.3 (units of F)', np.allclose(ev, [0.5265, 0.5563, 5.0293, 5.6929, 12.7553, 15.3028], atol=5e-3), str(np.round(ev, 4)))
Y = LL.Local.__new__(LL.Local); Y.__dict__.update(X.__dict__); Y.Em = {a: Fraction(0) for a in X.Em}
check('degree 8 with all error bounds set to zero: the inner box 0.005 passes (so the failure above is the error budget)', Y.certify_inner(Fraction(1, 200), verbose=False))

if len(sys.argv) > 1 and sys.argv[1] == 'quick':
    print('PASS (quick)' if not FAIL else 'FAILED: %s' % FAIL); sys.exit(1 if FAIL else 0)

# 4. degree 16: the theorem on a small region
X = LL.Local(16, verbose=False)
check('degree 16: F0 in [18.919102, 18.919105]', X.Flo > Fraction(18919102, 10 ** 6) and X.Fup < Fraction(18919105, 10 ** 6), '[%.10f, %.10f]' % (float(X.Flo), float(X.Fup)))
check('degree 16: the computed coefficients of degree 0 and 1 lie within their error bounds', X.low_check())
check('degree 16: inner box of half-width 1/128 proved', X.certify_inner(Fraction(1, 128), verbose=False))
check('degree 16: side conditions on the cube of half-width 3/128', X.side_conditions(Fraction(3, 128), verbose=False))
ok, st = X.certify_shell(Fraction(1, 128), Fraction(1, 50), Fraction(3, 128), verbose=False)
check('degree 16: shell 1/128 <= |e|_inf, |e|_2 < 1/50 covered', ok, str(st))
Y = LL.Local.__new__(LL.Local); Y.__dict__.update(X.__dict__); Y.Em = {a: 30 * v for a, v in X.Em.items()}
check('degree 16 with the error bounds multiplied by 30: the inner box fails', not Y.certify_inner(Fraction(1, 128), verbose=False))
Y.Tm = {a: -v if sum(a) == 2 else v for a, v in X.Tm.items()}; Y.Em = X.Em
ok2, st2 = Y.certify_shell(Fraction(1, 128), Fraction(1, 64), Fraction(1, 64), verbose=False, max_sec=20)
check('degree 16 with the quadratic form negated: the shell is not proved', not ok2, str(st2))
print('PASS' if not FAIL else 'FAILED: %s' % FAIL); sys.exit(1 if FAIL else 0)
