#!/usr/bin/env python3
"""Synthetic controls for compare_reports.py v2 (Claude, 2026-11-24).  Builds full-size fake report archives in a
temporary directory (15,223 main shard lines per side, fixed targets) and checks that the comparator accepts the two
legal formats (replay24c combined listings, replay26a per-round listings) and rejects every listed defect, including
the eight defects of GPT's PENT-26b controls.  Synthetic data only; no real report is involved.

    python3 -I test_compare_reports_v2.py"""
import os, sys, io, json, tarfile, tempfile, hashlib, subprocess, copy

HERE = os.path.dirname(os.path.abspath(__file__)); CMP = os.path.join(os.path.dirname(HERE), 'compare_reports.py')
sys.path.insert(0, os.path.dirname(HERE))
import compare_reports as K

def h(s): return hashlib.sha256(s.encode()).hexdigest()

def base():
    lines = {st: ['%s  main/%s/%06d.jsonl' % (h(st + str(i)), st, i) for i in range(K.FIXED[st][1])] for st in K.STAGES}
    supp = {'near': {'supplement_1': ['%s  supplement_1/near/000000.jsonl' % h('s1')], 'supplement_2': ['%s  supplement_2/near/000000.jsonl' % h('s2')]}}
    stats = {st: {'PROVED': K.FIXED[st][0]} for st in K.STAGES}
    stats['B_std'] = {'PROVED': K.FIXED['B_std'][0] - 1000, 'CONDITIONAL_TABLE': 1000}
    summ = dict(N_computational_gate=True, F_computational_gate=True, q_table_closed=True, structure=True, all_shards_enumerated=True,
                attempt_reject_stages=[], sample=0, stats=stats)
    plan = dict(config=dict(code=K.CODE_ID, zip_sha256=K.ZIP_SHA, sample=0), full_totals={st: K.FIXED[st][0] for st in K.STAGES},
                totals={st: K.FIXED[st][0] for st in K.STAGES})
    rep = ['accept_pent24 v2 2026-11-26 10:00:00 stages %s supplements work/supplement_1,work/supplement_2' % ','.join(K.STAGES),
           'plan.json {}']
    for st in K.STAGES:
        rep += ['%-8s targets %d  main: shard files %d ...' % (st, K.FIXED[st][0], K.FIXED[st][1]),
                '         result %s' % json.dumps(stats[st])]
    rep += ['B_std: 1000 labels depend on the table; the table is complete (q_lg and q_other all proved) in this check',
            'spot check: 400 leaves (uniform over 200000 proved leaves) x 2 points: failures 0, smallest relative margin 9.251e-04',
            'ACCEPT: coverage and binding complete for the requested stages']
    return dict(lines=lines, supp=supp, summ=summ, plan=plan, rep=rep, gpt_format='24c', extra_gpt=[], extra_claude=[])

def tgz(path, files):
    with tarfile.open(path, 'w:gz') as t:
        for name, data in files:
            ti = tarfile.TarInfo(name); b = data.encode() if isinstance(data, str) else data; ti.size = len(b); t.addfile(ti, io.BytesIO(b))

def build(d, F):
    g = [('replay/replay_summary.json', json.dumps(F['summ'])), ('results/plan.json', json.dumps(F['plan']))]
    c = [('claude_accept_all.txt', '\n'.join(F['rep']) + '\n')]
    for st in K.STAGES:
        main = F['lines'][st]; sp = F['supp'].get(st, {})
        if F['gpt_format'] == '24c':
            g.append(('replay/%s.shards.sha256' % st, '\n'.join(main + [l for r in sorted(sp) for l in sp[r]]) + '\n'))
        else:
            g.append(('replay/%s.shards.sha256' % st, '\n'.join(main) + '\n'))
            for r in sorted(sp):
                g.append(('replay/%s.%s.shards.sha256' % (st, r), '\n'.join(sp[r]) + '\n'))
        c.append(('claude_shards/%s.shards.sha256' % st, '\n'.join(main) + '\n'))
        for r in sorted(sp):
            c.append(('claude_shards/%s.%s.shards.sha256' % (st, r), '\n'.join(sp[r]) + '\n'))
    g += F['extra_gpt']; c += F['extra_claude']
    tgz(os.path.join(d, 'g.tgz'), g); tgz(os.path.join(d, 'c.tgz'), c)
    r = subprocess.run([sys.executable, '-I', CMP, '--gpt', os.path.join(d, 'g.tgz'), '--claude', os.path.join(d, 'c.tgz')], capture_output=True, text=True)
    return r.returncode, r.stdout.strip().splitlines()

def case(name, mutate, expect_ok):
    F = base(); mutate(F)
    with tempfile.TemporaryDirectory() as d:
        rc, out = build(d, F)
    ok = (rc == 0) == expect_ok
    why = [l for l in out if l.startswith('  - ')][:2]
    print('%-4s %-46s expected %-7s got %-15s %s' % ('PASS' if ok else 'FAIL', name, 'ACCEPT' if expect_ok else 'REJECT',
          'CONSISTENT' if rc == 0 else 'NOT CONSISTENT', ' | '.join(x.strip() for x in why)))
    return ok

def m(fn): return fn
C = [
 ('legal: replay24c combined listings', lambda F: None, True),
 ('legal: replay26a per-round listings', lambda F: F.update(gpt_format='26a'), True),
 ('both sides claim target 1 (wrong totals)', lambda F: (F['summ']['stats'].update(near={'PROVED': 1}), F['rep'].__setitem__(2, 'near     targets 1  main: ...'), F['rep'].__setitem__(3, '         result {"PROVED": 1}')), False),
 ('extra MISSING and ERROR', lambda F: (F['summ']['stats']['far_E2'].update(MISSING=3, ERROR=1)), False),
 ('CONDITIONAL_TABLE in near', lambda F: F['summ']['stats'].update(near={'PROVED': 38279, 'CONDITIONAL_TABLE': 1}), False),
 ('sample run', lambda F: (F['plan']['config'].update(sample=5), F['summ'].update(sample=5)), False),
 ('attempt_reject_stages absent', lambda F: F['summ'].pop('attempt_reject_stages'), False),
 ('MISSING pseudo-hash in main listing', lambda F: F['lines']['q_lg'].__setitem__(0, 'MISSING  main/q_lg/000000.jsonl'), False),
 ('duplicate listing line', lambda F: F['lines']['far_std'].append(F['lines']['far_std'][5]), False),
 ('ACCEPT after a rejected line', lambda F: F['rep'].insert(-3, '         rejected: {"bound below the target value": 1}'), False),
 ('spot check with too few leaves', lambda F: F['rep'].__setitem__(-2, 'spot check: 10 leaves (uniform over 200000 proved leaves) x 2 points: failures 0, smallest relative margin 1e-3'), False),
 ('Claude report not for all stages', lambda F: F['rep'].__setitem__(0, 'accept_pent24 v2 2026-11-26 10:00:00 stages near,q_lg supplements -'), False),
 ('one main shard missing on both sides', lambda F: F['lines']['B_E3'].pop(), False),
 ('one shard hash differs on one side', lambda F: F['extra_claude'].append(('claude_shards/near.supplement_3.shards.sha256', '%s  supplement_3/near/000000.jsonl\n' % h('x'))), False),
 ('a gate is the string "true"', lambda F: F['summ'].update(structure='true'), False),
 ('wrong code identity', lambda F: F['plan']['config'].update(code='0' * 64), False),
 ('unsafe member name in an archive', lambda F: F['extra_gpt'].append(('../evil.txt', 'x')), False),
 ('table not closed but B_std conditional counted', lambda F: F['summ'].update(q_table_closed=False), False),
 ('line in the wrong stage file', lambda F: F['lines']['near'].append('%s  main/q_lg/999999.jsonl' % h('w')), False),
]
res = [case(*c) for c in C]
print('%d / %d controls behave as expected' % (sum(res), len(res)))
sys.exit(0 if all(res) else 1)
