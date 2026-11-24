import json, lzma, sys, collections
from fractions import Fraction as F
arch = json.load(open(sys.argv[1])); g = json.load(lzma.open(sys.argv[2], 'rt'))
parts = json.load(lzma.open(sys.argv[3], 'rt'))
Ea = {k: F(v) for k, v in arch['Em'].items()}; Eg = {k: F(v) for k, v in g['Em'].items()}
b = collections.Counter(); lo = []
for k in Ea:
    r = float(Eg[k] / Ea[k]); deg = sum(map(int, k.split(',')))
    b[(deg, 'lt0.9' if r < 0.9 else ('le1' if r <= 1 else 'gt1'))] += 1
    if r < 0.6: lo.append((r, k))
for d in range(2, 14): print(d, {t: b[(d, t)] for t in ('lt0.9', 'le1', 'gt1')})
lo.sort(); print('lowest', lo[:5])
print(type(parts), list(parts)[:3] if isinstance(parts, dict) else parts[:1])
for r, k in lo[:3]:
    print(k, parts[k] if isinstance(parts, dict) and k in parts else '-', float(Ea[k]), float(Eg[k]))
