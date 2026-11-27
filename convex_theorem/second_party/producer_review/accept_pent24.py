#!/usr/bin/env python3
"""Claude's acceptance of the results of GPT's independent long run (package PENT24b r3, supplements of PENT24c).
Version 2 (2026-11-23): all nine stages, supplements, reasons of the unproved items.  Standard library only.
GPT's programs are not imported or run: the frozen inputs and the result files are read as data.

    python3 -I accept_pent24.py --near NEAR_DIR --far FAR_OUT_DIR --results RESULTS_DIR [--supplements DIR1,DIR2]
                                [--stages all | near,q_lg,...] [--sample M] [--points K] [--seed S]
                                [--partial] [--report FILE] [--shards-out DIR]

NEAR_DIR     frozen near cover (contains blocks/*.jsonl)                       = work/data      of the long run
FAR_OUT_DIR  frozen far cover (std|E2|E3/blocks/*, qtable/cells.jsonl, roots)  = work/data/out
RESULTS_DIR  work/results (one directory per stage with *.jsonl shards); supplements: work/supplement_<n>

For every stage the targets are enumerated again from the frozen inputs by this program, with the identifiers of the
long run:  near/<block>/<line> (LG leaves of the near cover), qstd/<cell>/<part> (table parts proved by LG: stage q_lg;
all other parts: q_other; a cell without parts: qstd/<cell>, q_other), <chart>/<block>/<line> (certified boxes,
far_std|far_E2|far_E3), B/<chart>/<block>/<line> (every leaf with a closed-form label, B_std|B_E2|B_E3).
Every result row must belong to a target of its stage, at most once per directory.

LG stages (near, q_lg, far_*): a row counts as proved only if status = PROVED and its witness is a tree whose inner
nodes {"axis", "left", "right"} halve the box at the exact dyadic midpoint and whose leaves {"lg": w} carry a
single-box witness with the chart and cut index of the target, exactly the box of that node, a bound L not below the
target value (F*_up, or the value stored in the table), no "vertices" field, a family of 1 or d + 1 members; at most 64
leaves in the main run, at most 256 in a supplement (where the degree N may differ from the frozen one).  So the
leaves tile the target box.  Version 2.1 (2026-11-23): in a supplement of the near stage a leaf may also be
{"outside": true}; it is accepted only if inf |z - C|^2 >= (3/5 + 2^-40)^2 on its box (exact), because Theorem N claims
only the points with |z - C| < 3/5; such leaves are counted apart, are not LG leaves and are not spot-checked, and a
tree needs at least one LG leaf.  The positivity of the box polynomials is GPT's replay and is not repeated here; instead
M leaves drawn uniformly from all proved leaves are checked at K points each by lg_point_check.py (polygon from the
chart formulas, u and chi rebuilt in physical coordinates, exact integrals, Lemma LG): another route.
q_other and B stages: coverage only (status PROVED with a witness; for B_std also CONDITIONAL_TABLE, which counts only
if q_lg and q_other are complete).  Their criteria were reviewed in the code and against GPT's LABEL_ROUTES.md; the
closed-form labels of the frozen cover were also re-proved in full by the producer's second implementation.

A target that is not proved in the main run is looked up in the supplements (only there).  Every target that remains
unproved is listed with the reason and the diagnostic fields its row carries.
Exit status 0 only if every target of every requested stage is proved, nothing is rejected and the spot check has no
failure (with --partial: missing rows are allowed and reported, and an incomplete last line is skipped).
--shards-out DIR writes DIR/<stage>.shards.sha256 ("<sha256>  main/<stage>/<file>" per shard, the format of GPT's
replay report) for the main run and DIR/<stage>.<supplement>.shards.sha256 for each supplement."""
import sys, os, json, glob, random, hashlib, time
from fractions import Fraction as Q

TARGET = Q('18.919134195733914')
LG = ('near', 'q_lg', 'far_std', 'far_E2', 'far_E3')
ALL = ('near', 'q_lg', 'q_other', 'far_std', 'far_E2', 'far_E3', 'B_std', 'B_E2', 'B_E3')
DIAG = ('reason', 'lower_center', 'ritz_center', 'parts', 'attempts', 'seconds', 'terms_PD', 'terms_PC', 'route')


def targets(stage, near, far):
    """identifier -> (chart, j or None, lo, hi, L) for LG stages, identifier -> None otherwise"""
    T = {}
    if stage == 'near':
        for p in sorted(glob.glob(os.path.join(near, 'blocks', '*.jsonl'))):
            with open(p) as f:
                next(f)
                for i, line in enumerate(f, 1):
                    lo, hi, lab, N = json.loads(line)
                    if lab == 'LG':
                        T['near/%s/%d' % (os.path.basename(p)[:-6], i)] = ('near', None, lo, hi, TARGET)
    elif stage in ('q_lg', 'q_other'):
        for i, line in enumerate(open(os.path.join(far, 'qtable', 'cells.jsonl')), 1):
            r = json.loads(line); parts = r.get('lg', {}).get('parts', [])
            if not parts and stage == 'q_other':
                T['qstd/%d' % i] = None
            for jj, (lo, hi, L, met) in enumerate(parts):
                if (met == 'lg') == (stage == 'q_lg'):
                    T['qstd/%d/%d' % (i, jj)] = ('qstd', None, lo, hi, Q(L)) if stage == 'q_lg' else None
    elif stage.startswith('far_'):
        chart = stage[4:]
        for p in sorted(glob.glob(os.path.join(far, chart, 'blocks', '*', 'cert.jsonl'))):
            for i, line in enumerate(open(p), 1):
                q = json.loads(line)
                T['%s/%s/%d' % (chart, os.path.basename(os.path.dirname(p)), i)] = (chart, q['j'], q['lo'], q['hi'], TARGET)
    elif stage.startswith('B_'):
        chart = stage[2:]
        blocks = sorted(json.load(open(os.path.join(far, chart, 'roots.json')))['blocks'])
        for b in blocks:
            with open(os.path.join(far, chart, 'blocks', b, 'leaves.jsonl')) as f:
                for i, _ in enumerate(f, 1):
                    T['B/%s/%s/%d' % (chart, b, i)] = None
    else:
        raise ValueError(stage)
    return T


# near chart: the binary64 point C of Theorem N (v2_run/near_lg.py, CF) and the squared radius of the "outside" test.
# GPT's programs use decimal literals that differ from C by one unit in the last place in five coordinates; the margin
# 2^-40 covers both points and the regular pentagon itself.
C_NEAR = [Q(float.fromhex(v)) for v in ('0x1.4f1bbcdcbfa53p+0', '0x1.e6f0e134454ffp-1', '0x1.0000000000000p-1',
                                        '0x1.89f188bdcd7aep+0', '-0x1.3c6ef372fe94cp-2', '0x1.e6f0e134454ffp-1')]
OUT2 = (Q(3, 5) + Q(1, 2 ** 40)) ** 2


def inf_d2(lo, hi):
    """exact inf over the box of |z - C_NEAR|^2"""
    s = Q(0)
    for l, h, c in zip(lo, hi, C_NEAR):
        if l > c:
            s += (l - c) ** 2
        elif h < c:
            s += (h - c) ** 2
    return s


def walk(w, lo, hi, tgt, leaves, depth=0, outs=None):
    """None if the tree tiles [lo, hi] by leaves bound to the target, else a reason.  outs is None: only LG leaves;
    otherwise (near stage, supplements only) also leaves {"outside": true} with inf |z - C|^2 >= (3/5 + 2^-40)^2 on
    their box, which are appended to outs: Theorem N claims only the points with |z - C| < 3/5 (near_lg.py)."""
    if depth > 16 or not isinstance(w, dict):
        return 'depth or type'
    if 'outside' in w:
        if outs is None:
            return 'outside leaf not allowed here'
        if set(w) != {'outside'} or w['outside'] is not True:
            return 'extra keys at a leaf'
        if len(lo) != 6 or inf_d2(lo, hi) < OUT2:
            return 'outside leaf meets the near region'
        outs.append((lo, hi))
        return None
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
        if not 1 <= int(p['N']) <= 14:
            return 'degree'
        leaves.append(p)
        return None
    if set(w) != {'axis', 'left', 'right'} or type(w['axis']) is not int or not 0 <= w['axis'] < len(lo):
        return 'bad inner node'
    i = w['axis']; m = (lo[i] + hi[i]) / 2
    a = hi[:]; a[i] = m; b = lo[:]; b[i] = m
    return walk(w['left'], lo, a, tgt, leaves, depth + 1, outs) or walk(w['right'], b, hi, tgt, leaves, depth + 1, outs)


class Run:
    def __init__(self, opt):
        self.opt = opt; self.out = []; self.ok = True; self.nout = 0
        self.M = int(opt.get('sample', 0)); self.K = int(opt.get('points', 2))
        self.rng = random.Random(int(opt.get('seed', 20261122))); self.reservoir = []; self.seen = 0

    def say(self, *s):
        line = ' '.join(str(x) for x in s); print(line, flush=True); self.out.append(line)

    def scan(self, st, T, root, tag, eligible, state, info, rej, maxleaves):
        """read the shards of one directory; returns (files, sha of the rows, listing text, leaves, largest tree)"""
        files = sorted(glob.glob(os.path.join(root, st, '*.jsonl'))); h = hashlib.sha256(); listing = []
        nleaf = 0; big = 0; seen_here = set()
        for fn in files:
            fh = hashlib.sha256()
            with open(fn, 'rb') as f:
                for raw in f:
                    h.update(raw); fh.update(raw)
                    if not raw.endswith(b'\n'):
                        if self.opt.get('partial'):
                            self.say('         incomplete last line skipped in', tag, os.path.basename(fn)); continue
                        rej['incomplete line'] = rej.get('incomplete line', 0) + 1; continue
                    r = json.loads(raw); i = r.get('id')
                    if i not in T or i in seen_here or (eligible is not None and i not in eligible):
                        rej['foreign, duplicate or not eligible (%s)' % tag] = rej.get('foreign, duplicate or not eligible (%s)' % tag, 0) + 1
                        continue
                    seen_here.add(i)
                    s = r.get('status')
                    good = s == 'PROVED' or (st == 'B_std' and s == 'CONDITIONAL_TABLE')
                    if not good or r.get('witness') is None:
                        if state.get(i) not in ('PROVED', 'CONDITIONAL_TABLE'):
                            state[i] = s if s in ('UNPROVED', 'ERROR') else 'REJECTED'
                            info[i] = {k: r[k] for k in DIAG if k in r}; info[i]['in'] = tag
                        continue
                    if st in LG:
                        leaves = []; outs = [] if (st == 'near' and tag != 'main') else None
                        try:
                            why = walk(r['witness'], [Q(x) for x in T[i][2]], [Q(x) for x in T[i][3]], T[i], leaves, 0, outs)
                            if why is None and len(leaves) > maxleaves:
                                why = 'more than %d leaves' % maxleaves
                            if why is None and not leaves:
                                why = 'no LG leaf'
                        except Exception as e:
                            why = 'malformed witness %r' % (e,)
                        if why:
                            state[i] = 'REJECTED'; rej[why] = rej.get(why, 0) + 1
                            if sum(rej.values()) <= 20:
                                self.say('  REJECTED', i, tag, why)
                            continue
                        nleaf += len(leaves); big = max(big, len(leaves)); self.nout += len(outs or ())
                        for p in leaves if self.M else ():
                            self.seen += 1
                            if len(self.reservoir) < self.M:
                                self.reservoir.append((i, T[i][4], p))
                            else:
                                k = self.rng.randrange(self.seen)
                                if k < self.M:
                                    self.reservoir[k] = (i, T[i][4], p)
                    state[i] = s; info.pop(i, None)
            listing.append('%s  %s/%s/%s' % (fh.hexdigest(), tag, st, os.path.basename(fn)))
        return files, h.hexdigest(), '\n'.join(listing) + '\n', nleaf, big

    def stage(self, st):
        opt = self.opt; T = targets(st, opt.get('near'), opt.get('far')); state = {}; info = {}; rej = {}
        files, hrows, ltxt, nleaf, big = self.scan(st, T, opt['results'], 'main', None, state, info, rej, 64)
        self.listing(st, 'main', ltxt)
        self.say('%-8s targets %d  main: shard files %d, leaves %d (largest tree %d), sha256 of the rows %s, of the listing %s'
                 % (st, len(T), len(files), nleaf, big, hrows, hashlib.sha256(ltxt.encode()).hexdigest()))
        for sd in [x for x in opt.get('supplements', '').split(',') if x]:
            tag = os.path.basename(os.path.normpath(sd))
            eligible = {i for i in T if state.get(i) in ('UNPROVED', 'ERROR')}
            if not os.path.isdir(os.path.join(sd, st)):
                continue
            n0 = self.nout
            files, hrows, ltxt, nleaf, big = self.scan(st, T, sd, tag, eligible, state, info, rej, 256)
            self.listing(st, tag, ltxt)
            self.say('         %s: shard files %d, leaves %d (largest tree %d)%s, sha256 of the rows %s, of the listing %s'
                     % (tag, len(files), nleaf, big, ', outside leaves %d' % (self.nout - n0) if st == 'near' else '',
                        hrows, hashlib.sha256(ltxt.encode()).hexdigest()))
        ct = {}
        for i in T:
            s = state.get(i, 'MISSING'); ct[s] = ct.get(s, 0) + 1
        self.say('         result %s' % json.dumps(ct, sort_keys=True))
        if rej:
            self.say('         rejected:', json.dumps(rej, sort_keys=True))
        notp = sorted(i for i in T if state.get(i, 'MISSING') not in ('PROVED', 'CONDITIONAL_TABLE', 'MISSING'))
        for i in notp[:60]:
            self.say('         not proved:', i, state[i], json.dumps(info.get(i, {}), sort_keys=True)[:300])
        if len(notp) > 60:
            self.say('         ... and %d more' % (len(notp) - 60))
        allowed = {'PROVED', 'CONDITIONAL_TABLE'} | ({'MISSING'} if opt.get('partial') else set())
        if not set(ct) <= allowed or rej:
            self.ok = False
        return ct

    def listing(self, st, tag, ltxt):
        d = self.opt.get('shards-out')
        if d:
            os.makedirs(d, exist_ok=True)
            open(os.path.join(d, '%s.shards.sha256' % st if tag == 'main' else '%s.%s.shards.sha256' % (st, tag)), 'w').write(ltxt)

    def spot(self):
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        import lg_point_check as P
        fails = 0; worst = None; t0 = time.time()
        for n, (i, L, p) in enumerate(self.reservoir, 1):
            d = len(p['lo'])
            pts = [[Q(self.rng.choice((-1, 1))) for _ in range(d)]]
            pts += [[Q(self.rng.randint(-2 ** 20, 2 ** 20), 2 ** 20) for _ in range(d)] for _ in range(self.K - 1)]
            for x in pts:
                stt, F = P.check_point(p, x)
                if stt != 'OK':
                    fails += 1; self.say('  SPOT CHECK FAIL', i, stt)
                else:
                    m = (F - L) / L; worst = m if worst is None or m < worst else worst
            if n % 50 == 0:
                print('  spot check %d / %d, %.0f s' % (n, len(self.reservoir), time.time() - t0), flush=True)
        self.say('spot check: %d leaves (uniform over %d proved leaves) x %d points: failures %d, smallest relative margin %s'
                 % (len(self.reservoir), self.seen, self.K, fails, 'n/a' if worst is None else '%.3e' % float(worst)))
        if fails:
            self.ok = False


def main():
    a = sys.argv[1:]; opt = {}
    while a:
        x = a.pop(0)
        if x == '--partial':
            opt['partial'] = True
        else:
            opt[x[2:]] = a.pop(0)
    st = opt.get('stages', 'all'); stages = ALL if st == 'all' else tuple(st.split(','))
    R = Run(opt)
    R.say('accept_pent24 v2', time.strftime('%Y-%m-%d %H:%M:%S'), 'stages', ','.join(stages), 'supplements', opt.get('supplements', '-'))
    plan = os.path.join(opt['results'], 'plan.json')
    if os.path.exists(plan):
        R.say('plan.json', json.dumps(json.load(open(plan))['config'], sort_keys=True))
    counts = {s: R.stage(s) for s in stages}
    if 'B_std' in counts and counts['B_std'].get('CONDITIONAL_TABLE'):
        closed = all(s in counts and set(counts[s]) == {'PROVED'} for s in ('q_lg', 'q_other'))
        R.say('B_std: %d labels depend on the table; the table is %s in this check'
              % (counts['B_std']['CONDITIONAL_TABLE'], 'complete (q_lg and q_other all proved)' if closed else 'NOT complete'))
        if not closed:
            R.ok = False
    if R.M and R.reservoir:
        R.spot()
    R.say('ACCEPT: coverage and binding complete for the requested stages' if R.ok else 'NOT COMPLETE (see above)')
    if opt.get('report'):
        open(opt['report'], 'w').write('\n'.join(R.out) + '\n')
    sys.exit(0 if R.ok else 1)


if __name__ == '__main__':
    main()
