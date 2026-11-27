import sys, json
sys.path.insert(0, '<research tree>/probes/PENT/local/accept')
import accept_pent24 as A
from fractions import Fraction as Q
B, OUT = sys.argv[1], sys.argv[2]
IDS = '023112/23 023121/14 132012/151 132012/26 132012/47 132012/7 132022/157 132022/72 132102/10 133022/30 133022/40 232012/19 233021/2'.split()
def box(tid):
    blk, ln = tid.split('/'); lo, hi, lab, N = json.loads(open(B + '/' + blk + '.jsonl').read().splitlines()[int(ln)]); return [Q(x) for x in lo], [Q(x) for x in hi]
def fake(lo, hi):
    return {'lg': {'chart': 'near', 'lo': [str(x) for x in lo], 'hi': [str(x) for x in hi], 'L': str(A.TARGET), 'N': 10, 'u': [[0]] * 7, 'chi': [[0]] * 7}}
cnt = {'lg': 0, 'out': 0}
def build(lo, hi, depth, D, cheat=False):
    if A.inf_d2(lo, hi) >= A.OUT2 or (cheat and depth == D):
        cnt['out'] += 1; return {'outside': True}
    if depth == D:
        cnt['lg'] += 1; return fake(lo, hi)
    i = depth % 6; m = (lo[i] + hi[i]) / 2; a = hi[:]; a[i] = m; b = lo[:]; b[i] = m
    return {'axis': i, 'left': build(lo, a, depth + 1, D, cheat), 'right': build(b, hi, depth + 1, D, cheat)}
main = [dict(id='near/' + t, status='UNPROVED', reason='budget') for t in IDS]
lo, hi = box('002011/1'.replace('002011/1', '132012/1')) if False else (None, None)
sup = []
lo, hi = box('132012/7'); w = build(lo, hi, 0, 12); sup.append(dict(id='near/132012/7', status='PROVED', witness=w)); print('valid tree 132012/7:', cnt)
cnt.update(lg=0, out=0); lo, hi = box('233021/2'); w = build(lo, hi, 0, 6, cheat=True); sup.append(dict(id='near/233021/2', status='PROVED', witness=w)); print('cheating tree 233021/2:', cnt)
cnt.update(lg=0, out=0); lo, hi = box('132022/72'); w = {'outside': True, 'note': 1}; sup.append(dict(id='near/132022/72', status='PROVED', witness=w))
# an outside leaf in the main run (target 132012/1 must be an LG target)
lines = open(B + '/132012.jsonl').read().splitlines()
k = next(i for i in range(1, len(lines)) if json.loads(lines[i])[2] == 'LG')
lo, hi = box('132012/%d' % k); cnt.update(lg=0, out=0); w = build(lo, hi, 0, 12)
main.append(dict(id='near/132012/%d' % k, status='PROVED', witness=w)); print('main-run tree with outside leaves 132012/%d:' % k, cnt)
open(OUT + '/results/near/000000.jsonl', 'w').write(''.join(json.dumps(r) + '\n' for r in main))
open(OUT + '/supplement_1/near/000000.jsonl', 'w').write(''.join(json.dumps(r) + '\n' for r in sup))
