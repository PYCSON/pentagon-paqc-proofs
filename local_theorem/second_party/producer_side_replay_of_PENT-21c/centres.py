# |Tm - Tm'| <= K_rad + 2 grid steps for every alpha (archived witness.json, GPT own_levelB_enclosure_refined.json.xz, archived parts.json.xz)
import json, lzma, sys
from fractions import Fraction as F
a = json.load(open(sys.argv[1])); g = json.load(lzma.open(sys.argv[2], 'rt')); p = json.load(lzma.open(sys.argv[3], 'rt'))['parts']
U = F(1, 2**160); worst = 0; bad = 0
for k, v in p.items():
    d = abs(F(a['Tm'].get(k, '0')) - F(g['Tm'][k]))
    lim = F(v[0]) * U + 2 * U
    bad += d > lim; worst = max(worst, float(d / lim))
print('centres: |Tm - Tm\'| <= K_rad + 2 grid steps for', len(p) - bad, 'of', len(p), '; worst ratio', worst)
