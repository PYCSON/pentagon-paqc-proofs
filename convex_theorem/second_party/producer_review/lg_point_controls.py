#!/usr/bin/env python3
"""Negative controls for lg_point_check.py:  python3 -I lg_point_controls.py WITNESS.jsonl"""
import sys, json, os, copy
from fractions import Fraction as Q
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lg_point_check as L
w = next(json.loads(l)['witness'] for l in open(sys.argv[1]) if json.loads(l).get('witness'))
d = len(w['lo']); x = [Q(1, 3)] * d
print('unchanged      ', L.check_point(w, x)[0])
a = copy.deepcopy(w); a['L'] = '19'; print('target 19      ', L.check_point(a, x)[0])
b = copy.deepcopy(w); b['u'][0] = b['u'][0][1:] + b['u'][0][:1]; print('u list rotated ', L.check_point(b, x)[0])
c = copy.deepcopy(w); c['chi'][0] = [0] * len(c['chi'][0]); c['chi'] = [c['chi'][0]] * len(c['chi']); print('chi = 0        ', L.check_point(c, x)[0])
e = copy.deepcopy(w); e['hi'] = [str(Q(lo) + 40 * (Q(hi) - Q(lo))) for lo, hi in zip(e['lo'], e['hi'])]; print('box 40x wider  ', L.check_point(e, [Q(1)] * d)[0])
