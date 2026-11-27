#!/usr/bin/env python3
"""Claude's acceptance of the LG results of GPT's independent long run (package PENT24b r3).  Standard library only.
GPT's programs are not imported or run: the frozen inputs and the result files are read as data.

    python3 -I accept_pent24b_r3.py --near NEAR_DIR --far FAR_OUT_DIR --results RESULTS_DIR
                                    [--stages near,q_lg,far_std,far_E2,far_E3] [--sample M] [--points K] [--seed S]
                                    [--partial] [--report FILE] [--shards-out DIR]

NEAR_DIR     frozen near cover (contains blocks/*.jsonl)            = work/data          of the long run
FAR_OUT_DIR  frozen far cover (std|E2|E3/blocks/*/cert.jsonl, qtable/cells.jsonl) = work/data/out
RESULTS_DIR  work/results of the long run (one directory per stage with *.jsonl shards)

For every LG stage (near, q_lg, far_std, far_E2, far_E3):
 1. the targets are enumerated again from the frozen inputs by this program (identifier, chart, cut index, box, target
    value F*_up or the value stored in the quadrilateral table);
 2. every result row must belong to a target of its stage (by identifier), at most once;
 3. a row counts as proved only if status = PROVED and its witness is a tree of this form: a node {"axis", "left",
    "right"} halves the box in coordinate "axis" (exact dyadic midpoint), a leaf {"lg": w} carries a single-box
    witness w whose chart and cut index are those of the target, whose box is exactly the box of that node, whose
    bound L is not below the target value, which has no "vertices" field and whose family has 1 or d + 1 members;
    at most 64 leaves.  So the leaves tile the target box and each leaf claims the full statement on its part;
 4. coverage: every target has a proved row; all others are listed (UNPROVED / ERROR / MISSING / REJECTED);
 5. spot check: for M leaves drawn uniformly from all proved leaves (reservoir sampling), the pointwise bound of
    lg_point_check.py at K points each (one corner and K - 1 random points): an exact computation by another route
    (polygon from the chart formulas, u and chi as polynomials in physical coordinates, exact integrals, Lemma LG).
The positivity of the box polynomials themselves is GPT's replay and is not repeated here.
Exit status 0 only if every target of every requested stage is proved, nothing is rejected and the spot check has no
failure (with --partial: missing rows are allowed and reported, and an incomplete last line of a shard that is still
being written is skipped).
--shards-out DIR writes DIR/<stage>.shards.sha256 with one line "<sha256 of the shard file>  main/<stage>/<file>" per
shard, the format of GPT's replay report, so that both checks can be seen to have read the same bytes; the SHA-256 of
that listing is printed."""
import sys, os, json, glob, random, hashlib, time
from fractions import Fraction as Q

TARGET = Q('18.919134195733914')
STAGES = ('near', 'q_lg', 'far_std', 'far_E2', 'far_E3')


def targets(stage, near, far):
    """identifier -> (chart, j or None, lo, hi, L); lo, hi as read (binary64), L a Fraction"""
    T = {}
    if stage == 'near':
        for p in sorted(glob.glob(os.path.join(near, 'blocks', '*.jsonl'))):
            with open(p) as f:
                next(f)
                for i, line in enumerate(f, 1):
                    lo, hi, lab, N = json.loads(line)
                    if lab == 'LG':
                        T['near/%s/%d' % (os.path.basename(p)[:-6], i)] = ('near', None, lo, hi, TARGET)
    elif stage == 'q_lg':
        for i, line in enumerate(open(os.path.join(far, 'qtable', 'cells.jsonl')), 1):
            r = json.loads(line)
            for jj, (lo, hi, L, met) in enumerate(r.get('lg', {}).get('parts', [])):
                if met == 'lg':
                    T['qstd/%d/%d' % (i, jj)] = ('qstd', None, lo, hi, Q(L))
    else:
        chart = stage[4:]
        for p in sorted(glob.glob(os.path.join(far, chart, 'blocks', '*', 'cert.jsonl'))):
            for i, line in enumerate(open(p), 1):
                q = json.loads(line)
                T['%s/%s/%d' % (chart, os.path.basename(os.path.dirname(p)), i)] = (chart, q['j'], q['lo'], q['hi'], TARGET)
    return T


def walk(w, lo, hi, tgt, leaves, depth=0):
    """returns None if the tree is a valid tiling of [lo, hi] by leaves bound to the target, else a reason"""
    if depth > 12 or not isinstance(w, dict):
        return 'depth or type'
    if 'lg' in w:
        if set(w) != {'lg'}:
            return 'extra keys at a leaf'
        p = w['lg']; chart, j, _, _, L = tgt; d = len(lo)
        if 'vertices' in p:
            return 'vertices field'
        if p.get('chart') != chart or (j is not None and p.get('j') != j):
            return 'chart or cut index'
        if len(p['lo']) != d or len(p['hi']) != d or [Q(x) for x in p['lo']] != lo or [Q(x) for x in p['hi']] != hi:
            return 'leaf box is not the box of its node'
        if Q(p['L']) < L:
            return 'bound below the target value'
        if len(p['u']) not in (1, d + 1) or len(p['chi']) != len(p['u']):
            return 'family size'
        leaves.append(p)
        return None
    if set(w) != {'axis', 'left', 'right'} or type(w['axis']) is not int or not 0 <= w['axis'] < len(lo):
        return 'bad inner node'
    i = w['axis']; m = (lo[i] + hi[i]) / 2
    a = hi[:]; a[i] = m; b = lo[:]; b[i] = m
    return walk(w['left'], lo, a, tgt, leaves, depth + 1) or walk(w['right'], b, hi, tgt, leaves, depth + 1)


def main():
    a = sys.argv[1:]; opt = {}
    while a:
        x = a.pop(0)
        if x in ('--partial',):
            opt[x[2:]] = True
        else:
            opt[x[2:]] = a.pop(0)
    stages = opt.get('stages', ','.join(STAGES)).split(',')
    M = int(opt.get('sample', 0)); K = int(opt.get('points', 2)); rng = random.Random(int(opt.get('seed', 20261122)))
    out = []
    def say(*s):
        line = ' '.join(str(x) for x in s); print(line, flush=True); out.append(line)
    say('accept_pent24b_r3', time.strftime('%Y-%m-%d %H:%M:%S'), 'stages', stages)
    plan = os.path.join(opt['results'], 'plan.json')
    if os.path.exists(plan):
        say('plan.json', json.dumps(json.load(open(plan))['config'], sort_keys=True))
    ok = True; reservoir = []; seen_leaves = 0; total_leaf = 0
    for st in stages:
        T = targets(st, opt.get('near'), opt.get('far')); state = {}; rej = {}; nleaf = 0; maxleaf = 0; files = sorted(glob.glob(os.path.join(opt['results'], st, '*.jsonl')))
        h = hashlib.sha256(); listing = []
        for fn in files:
            fh = hashlib.sha256()
            with open(fn, 'rb') as f:
                for raw in f:
                    h.update(raw); fh.update(raw)
                    if not raw.endswith(b'\n'):
                        if opt.get('partial'):
                            say('         incomplete last line skipped in', os.path.basename(fn)); continue
                        rej['incomplete line'] = rej.get('incomplete line', 0) + 1; continue
                    r = json.loads(raw); i = r.get('id')
                    if i not in T or i in state:
                        rej['foreign or duplicate row'] = rej.get('foreign or duplicate row', 0) + 1; ok = False; continue
                    s = r.get('status')
                    if s != 'PROVED':
                        state[i] = s if s in ('UNPROVED', 'ERROR') else 'REJECTED'; continue
                    chart, j, lo, hi, L = T[i]; leaves = []
                    try:
                        why = walk(r.get('witness'), [Q(x) for x in lo], [Q(x) for x in hi], T[i], leaves)
                        if why is None and len(leaves) > 64:
                            why = 'more than 64 leaves'
                    except Exception as e:
                        why = 'malformed witness %r' % (e,)
                    if why:
                        state[i] = 'REJECTED'; rej[why] = rej.get(why, 0) + 1
                        if sum(rej.values()) <= 20:
                            say('  REJECTED', i, why)
                        continue
                    state[i] = 'PROVED'; nleaf += len(leaves); maxleaf = max(maxleaf, len(leaves))
                    if M:
                        for p in leaves:
                            seen_leaves += 1
                            if len(reservoir) < M:
                                reservoir.append((i, L, p))
                            else:
                                k = rng.randrange(seen_leaves)
                                if k < M:
                                    reservoir[k] = (i, L, p)
            listing.append('%s  main/%s/%s' % (fh.hexdigest(), st, os.path.basename(fn)))
        ltxt = '\n'.join(listing) + '\n'
        if opt.get('shards-out'):
            os.makedirs(opt['shards-out'], exist_ok=True)
            open(os.path.join(opt['shards-out'], st + '.shards.sha256'), 'w').write(ltxt)
        ct = {}
        for i in T:
            s = state.get(i, 'MISSING'); ct[s] = ct.get(s, 0) + 1
        total_leaf += nleaf
        say('%-8s targets %d  %s  leaves %d (largest tree %d)  shard files %d  sha256 of the rows %s  sha256 of the shard listing %s'
            % (st, len(T), json.dumps(ct, sort_keys=True), nleaf, maxleaf, len(files), h.hexdigest(), hashlib.sha256(ltxt.encode()).hexdigest()))
        if rej:
            say('         rejected:', json.dumps(rej, sort_keys=True))
        notp = sorted(i for i in T if state.get(i, 'MISSING') not in ('PROVED', 'MISSING'))
        for i in notp[:40]:
            say('         not proved:', i, state[i])
        if len(notp) > 40:
            say('         ... and %d more' % (len(notp) - 40))
        if ct.get('PROVED', 0) != len(T) and not (opt.get('partial') and set(ct) <= {'PROVED', 'MISSING'}):
            ok = False
        if rej:
            ok = False
    if M and reservoir:
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import lg_point_check as P
        fails = 0; worst = None; t0 = time.time()
        for n, (i, L, p) in enumerate(reservoir, 1):
            d = len(p['lo'])
            pts = [[Q(rng.choice((-1, 1))) for _ in range(d)]] + [[Q(rng.randint(-2 ** 20, 2 ** 20), 2 ** 20) for _ in range(d)] for _ in range(K - 1)]
            for x in pts:
                stt, F = P.check_point(p, x)
                if stt != 'OK':
                    fails += 1; say('  SPOT CHECK FAIL', i, stt)
                else:
                    m = (F - L) / L
                    worst = m if worst is None or m < worst else worst
            if n % 50 == 0:
                print('  spot check %d / %d, %.0f s' % (n, len(reservoir), time.time() - t0), flush=True)
        say('spot check: %d leaves (uniform over %d proved leaves) x %d points: failures %d, smallest relative margin over the target %s'
            % (len(reservoir), seen_leaves, K, fails, 'n/a' if worst is None else '%.3e' % float(worst)))
        if fails:
            ok = False
    say('ACCEPT: coverage and binding complete for the requested stages' if ok else 'NOT COMPLETE (see above)')
    if opt.get('report'):
        open(opt['report'], 'w').write('\n'.join(out) + '\n')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
