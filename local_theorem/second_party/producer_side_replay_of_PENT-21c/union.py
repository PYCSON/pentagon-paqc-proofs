import json, lzma, sys
from fractions import Fraction as F
arch = json.load(open(sys.argv[1])); g = json.load(lzma.open(sys.argv[2], 'rt'))
Ta = {k: F(v) for k, v in arch['Tm'].items()}; Ea = {k: F(v) for k, v in arch['Em'].items()}
Tg = {k: F(v) for k, v in g['Tm'].items()}; Eg = {k: F(v) for k, v in g['Em'].items()}
Em = {}; wider = 0
for k in set(Ea) | set(Eg):
    ta, tg = Ta.get(k, F(0)), Tg.get(k, F(0))
    e = max(Ea.get(k, F(0)), abs(ta - tg) + Eg.get(k, F(0)))   # Ta +- e contains both intervals
    wider += e > Ea.get(k, F(0)); Em[k] = e
print('union enclosure: centre = archived Tm, radius = max(Em, |Tm - Tm\'| + Em\'); wider than archived for', wider, 'of', len(Em))
json.dump(dict(format='pent-local-cert/1', rho_in=arch['rho_in'], gbits=arch['gbits'], Tm=arch['Tm'], Em={k: str(v) for k, v in Em.items()}), open(sys.argv[3], 'w'))
