import json, lzma, sys
from fractions import Fraction as F
arch = json.load(open(sys.argv[1])); g = json.load(lzma.open(sys.argv[2], 'rt'))
print('gpt keys', list(g), 'gbits', g['gbits'], 'rho_in', g['rho_in'], 'arch rho', arch['rho_in'], arch['gbits'])
Ta = {k: F(v) for k, v in arch['Tm'].items()}; Ea = {k: F(v) for k, v in arch['Em'].items()}
Tg = {k: F(v) for k, v in g['Tm'].items()}; Eg = {k: F(v) for k, v in g['Em'].items()}
print('counts arch Tm %d Em %d ; gpt Tm %d Em %d' % (len(Ta), len(Ea), len(Tg), len(Eg)))
keys = set(Ta) | set(Ea) | set(Tg) | set(Eg)
print('keys union', len(keys), 'gpt-only', len(keys - set(Ta) - set(Ea)), 'arch-only', len(keys - set(Tg) - set(Eg)))
disj = 0; incl = 0; worst = 0; rmax = 0; rmin = 10; negE = sum(1 for v in Eg.values() if v < 0)
for k in keys:
    ta, ea, tg, eg = Ta.get(k, F(0)), Ea.get(k, F(0)), Tg.get(k, F(0)), Eg.get(k, F(0))
    d = abs(ta - tg)
    if d > ea + eg: disj += 1
    if d + eg <= ea: incl += 1
    if ea > 0:
        r = float(eg / ea); rmax = max(rmax, r); rmin = min(rmin, r)
    if ea + eg > 0: worst = max(worst, float(d / (ea + eg)))
print('negative Em in gpt', negE)
print('DISJOINT (both cannot be valid):', disj, '; archived interval contains gpt interval:', incl, 'of', len(keys))
print('max |Tm-Tm\'|/(Em+Em\')', worst, '; Em\'/Em range', rmin, rmax)
w = dict(format='pent-local-cert/1', rho_in=g['rho_in'], gbits=g['gbits'], Tm=g['Tm'], Em=g['Em'])
json.dump(w, open(sys.argv[3], 'w'))
