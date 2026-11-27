#!/usr/bin/env python3
"""Compare GPT's replay report with Claude's acceptance report -- version 2 (Claude, 2026-11-24; rewritten after GPT's
PENT-26b review, which found 7 false accepts and 1 false reject in version 1).  Standard library only; nothing is
extracted to disk.

    python3 -I compare_reports.py --gpt REPLAY_REPORT.tgz --claude claude_accept_report.tgz [--format-test]

REPLAY_REPORT.tgz = work/PENT24b_replay_report.tgz (replay24c) or work/PENT26a_replay_report.tgz (replay26a);
claude_accept_report.tgz = claude_accept_all.txt + claude_shards/ (accept_pent24.py v2.1 with --shards-out).

Archive safety: only regular files and directories; no absolute names, no '..', no links, no duplicate names.
Files are read at fixed names (never "the first file of that name somewhere").

Contract (all must hold for VERDICT: CONSISTENT):
 1. identity: the replay's results plan has code 7cf79bee..., zip_sha256 1aba5f0a..., sample 0; its totals and the
    summary's sample agree with the fixed full targets below (1,947,824 tasks, 15,223 main shards);
 2. shard listings: every line of either side is "<64 hex>  <round>/<stage>/<6 digits>.jsonl" with round main or
    supplement_<k>, in the file of its stage (and of its round, for per-round files); the multiset of
    (stage, round, path, sha256) is the same on both sides; no tuple occurs twice; every stage has exactly its fixed
    number of distinct main shards; lines "MISSING ..." are tolerated only for supplement rounds and are reported;
 3. counts: on both sides, for every stage, PROVED (+ CONDITIONAL_TABLE for B_std only, and only when the table is
    closed -- q_table_closed on GPT's side, the table-complete line on Claude's side) equals the fixed target, and every
    other status (UNPROVED, ERROR, MISSING, REJECT, REJECTED, CONDITIONAL_TABLE elsewhere, ...) is 0;
 4. GPT's gates N_computational_gate, F_computational_gate, q_table_closed, structure, all_shards_enumerated are the
    boolean true, and attempt_reject_stages is present and empty;
 5. Claude's report: all nine stages requested; no rejected rows, no incomplete lines, no "NOT COMPLETE", no spot-check
    failure; spot check with at least 400 leaves x 2 points, 0 failures, positive smallest margin; the last line is the
    single ACCEPT line.
--format-test replaces the fixed targets by the replay plan's own totals and allows a sample run, to test formats on
small samples; its verdict is printed as FORMAT TEST and is never a sign-off.
This comparison shows that both programs read the same bytes and report the same counts; it is not a mathematical
replay, and matching hashes cannot exclude two forged reports: the raw shards stay on the user's machine."""
import sys, re, json, tarfile, collections

STAGES = ('near', 'q_lg', 'q_other', 'far_std', 'far_E2', 'far_E3', 'B_std', 'B_E2', 'B_E3')
FIXED = {'near': (38280, 300), 'q_lg': (66138, 517), 'q_other': (118962, 930), 'far_std': (293846, 2296),
         'far_E2': (47804, 374), 'far_E3': (35865, 281), 'B_std': (934801, 7304), 'B_E2': (185823, 1452),
         'B_E3': (226305, 1769)}
CODE_ID = '7cf79beefbc81c796969279b7c4d8bf35c7ecb0525e41f2fe3d88d2c629e2d3b'
ZIP_SHA = '1aba5f0a693313b515d9b3b0952747cd61ced2b4b5fb2d3c9c5783333ff2f958'
GATES = ('N_computational_gate', 'F_computational_gate', 'q_table_closed', 'structure', 'all_shards_enumerated')
LINE = re.compile(r'^([0-9a-f]{64})  ((main|supplement_[0-9]+)/([A-Za-z_0-9]+)/([0-9]{6}\.jsonl))$')
LFILE = re.compile(r'^(?:replay|claude_shards)/([A-Za-z_0-9]+?)(?:\.(supplement_[0-9]+))?\.shards\.sha256$')


def read_tar(path, issues, who):
    files = {}
    with tarfile.open(path) as t:
        for m in t.getmembers():
            n = m.name
            while n.startswith('./'):
                n = n[2:]
            if n.startswith('/') or '..' in n.split('/'):
                issues.append('%s archive: unsafe member name %r' % (who, m.name)); continue
            if m.isdir():
                continue
            if not m.isreg():
                issues.append('%s archive: member %r is not a regular file' % (who, m.name)); continue
            if n in files:
                issues.append('%s archive: duplicate member %r' % (who, n)); continue
            files[n] = t.extractfile(m).read()
    return files


def listings(files, who, issues, notes):
    """Counter of (stage, round, path, sha) from all listing files of one side"""
    C = collections.Counter(); seen = 0
    for name, data in sorted(files.items()):
        m = LFILE.match(name)
        if not m:
            continue
        st, rnd = m.group(1), m.group(2)
        if st not in STAGES:
            issues.append('%s: listing for unknown stage %s' % (who, name)); continue
        seen += 1
        for raw in data.decode('utf-8', 'replace').splitlines():
            if not raw.strip():
                continue
            x = LINE.match(raw)
            if x:
                sha, path, r, st2 = x.group(1), x.group(2), x.group(3), x.group(4)
                if st2 != st or (rnd and r != rnd):
                    issues.append('%s: line in the wrong listing file %s: %s' % (who, name, raw[:100])); continue
                C[(st, r, path, sha)] += 1
            elif raw.startswith('MISSING') and 'supplement_' in raw and 'main/' not in raw:
                notes.append('%s: sentinel in %s: %s' % (who, name, raw[:100]))
            else:
                issues.append('%s: malformed listing line in %s: %r' % (who, name, raw[:100]))
    if not seen:
        issues.append('%s: no shard listings found' % who)
    for k, n in C.items():
        if n > 1:
            issues.append('%s: duplicate listing entry %s x%d' % (who, k[2], n))
    return C


def check_counts(who, st, counts, target, table_closed, issues):
    good = counts.get('PROVED', 0)
    if st == 'B_std' and table_closed:
        good += counts.get('CONDITIONAL_TABLE', 0)
    other = {k: v for k, v in counts.items() if v and k != 'PROVED' and not (k == 'CONDITIONAL_TABLE' and st == 'B_std' and table_closed)}
    if good != target:
        issues.append('%s %s: proved %d, target %d' % (who, st, good, target))
    if other:
        issues.append('%s %s: other statuses %s' % (who, st, json.dumps(other, sort_keys=True)))


def main():
    a = sys.argv[1:]; opt = {}
    while a:
        k = a.pop(0)
        if k == '--format-test':
            opt['format-test'] = True
        else:
            opt[k[2:]] = a.pop(0)
    ft = opt.get('format-test', False)
    issues = []; notes = []; out = []
    G = read_tar(opt['gpt'], issues, 'GPT'); C = read_tar(opt['claude'], issues, 'Claude')
    # 1. identity and targets
    summ = json.loads(G['replay/replay_summary.json']) if 'replay/replay_summary.json' in G else None
    plan = json.loads(G.get('results/plan.json') or G.get('replay/results_plan.json') or b'null')
    if summ is None or plan is None:
        issues.append('GPT archive: replay/replay_summary.json or the results plan is missing')
        summ = summ or {}; plan = plan or {}
    cfg = plan.get('config', {})
    if cfg.get('code') != CODE_ID or cfg.get('zip_sha256') != ZIP_SHA:
        issues.append('identity: code %s / zip %s' % (cfg.get('code'), cfg.get('zip_sha256')))
    target = {st: FIXED[st][0] for st in STAGES}; nmain = {st: FIXED[st][1] for st in STAGES}
    if ft:
        target = {st: plan.get('totals', {}).get(st, -1) for st in STAGES}; nmain = None
        out.append('FORMAT TEST: targets taken from the replay plan %s' % json.dumps(target, sort_keys=True))
    else:
        if cfg.get('sample') != 0 or summ.get('sample') != 0:
            issues.append('sample run (plan sample %s, summary sample %s)' % (cfg.get('sample'), summ.get('sample')))
        for key in ('full_totals', 'totals'):
            if plan.get(key) is not None and plan[key] != {st: FIXED[st][0] for st in STAGES}:
                issues.append('plan %s differ from the fixed targets' % key)
    # 4. gates
    for g in GATES:
        if summ.get(g) is not True:
            (notes if ft else issues).append('GPT gate %s is %r' % (g, summ.get(g)))
    if 'attempt_reject_stages' not in summ or summ['attempt_reject_stages'] != []:
        issues.append('GPT attempt_reject_stages is %r' % summ.get('attempt_reject_stages', '<absent>'))
    # 5. Claude's report
    rep = C.get('claude_accept_all.txt')
    if rep is None:
        issues.append('Claude archive: claude_accept_all.txt missing'); text = []
    else:
        text = rep.decode('utf-8', 'replace').splitlines()
    head = text[0] if text else ''
    m = re.match(r'^accept_pent24 v2 \S+ \S+ stages (\S+) supplements (\S+)', head)
    if not m or m.group(1).split(',') != list(STAGES):
        issues.append('Claude report: header does not request all nine stages: %r' % head[:120])
    cres = {}; ctarg = {}; cur = None
    for l in text:
        x = re.match(r'^(\w+)\s+targets (\d+)', l)
        if x and x.group(1) in STAGES:
            cur = x.group(1); ctarg[cur] = int(x.group(2)); continue
        x = re.match(r'^\s+result (\{.*\})\s*$', l)
        if x and cur:
            cres[cur] = json.loads(x.group(1))
        for bad in ('REJECTED', 'rejected:', 'incomplete last line', 'NOT COMPLETE', 'SPOT CHECK FAIL', 'malformed'):
            if bad in l:
                issues.append('Claude report: %r in line %r' % (bad, l[:120])); break
    tclosed_c = any('complete (q_lg and q_other all proved)' in l for l in text)
    sp = [l for l in text if l.startswith('spot check:')]
    if len(sp) != 1:
        issues.append('Claude report: %d spot-check summary lines' % len(sp))
    else:
        x = re.match(r'^spot check: (\d+) leaves .*? x (\d+) points: failures (\d+), smallest relative margin(?: over the target)? (\S+)', sp[0])
        if not x:
            issues.append('Claude report: spot-check line not understood')
        else:
            nl, npnt, nf, mg = int(x.group(1)), int(x.group(2)), int(x.group(3)), x.group(4)
            if nf != 0 or (not ft and (nl < 400 or npnt < 2)) or mg == 'n/a' or float(mg) <= 0:
                issues.append('Claude report: spot check %d leaves x %d points, %d failures, margin %s' % (nl, npnt, nf, mg))
    last = [l for l in text if l.strip()][-1:] if text else []
    if last != ['ACCEPT: coverage and binding complete for the requested stages'] or sum(l.startswith('ACCEPT') for l in text) != 1:
        issues.append('Claude report: the last line is not the single ACCEPT line')
    # 2. listings
    LG = listings(G, 'GPT', issues, notes); LC = listings(C, 'Claude', issues, notes)
    # 3. counts and the per-stage table
    for st in STAGES:
        g = collections.Counter({k: n for k, n in LG.items() if k[0] == st}); c = collections.Counter({k: n for k, n in LC.items() if k[0] == st})
        rounds = sorted({k[1] for k in list(g) + list(c)}, key=lambda r: (r != 'main', r))
        mains = len({k[2] for k in g if k[1] == 'main'})
        same = g == c and len(g) > 0
        gs = summ.get('stats', {}).get(st, {}); cs = cres.get(st, {})
        out.append('%-8s target %d | listings GPT %d, Claude %d, %s (rounds %s, main shards %d) | GPT %s | Claude %s' % (
            st, target[st], sum(g.values()), sum(c.values()), 'IDENTICAL' if same else 'DIFFERENT', ','.join(rounds) or '-',
            mains, json.dumps({k: v for k, v in gs.items() if v}, sort_keys=True), json.dumps(cs, sort_keys=True)))
        if not same:
            onlyg = sorted((g - c).elements()); onlyc = sorted((c - g).elements())
            issues.append('%s: shard listings differ (only GPT %d, only Claude %d)' % (st, len(onlyg), len(onlyc)))
            for k in onlyg[:3]:
                out.append('           only GPT:    %s  %s' % (k[3], k[2]))
            for k in onlyc[:3]:
                out.append('           only Claude: %s  %s' % (k[3], k[2]))
        if nmain is not None and mains != nmain[st]:
            issues.append('%s: %d main shards, fixed %d' % (st, mains, nmain[st]))
        if ctarg.get(st) != target[st]:
            issues.append('Claude %s: targets %s, fixed %d' % (st, ctarg.get(st), target[st]))
        check_counts('GPT', st, gs, target[st], summ.get('q_table_closed') is True, issues)
        check_counts('Claude', st, cs, target[st], tclosed_c, issues)
    out.append('GPT gates: %s; attempt_reject_stages %s' % (json.dumps({g: summ.get(g) for g in GATES}, sort_keys=True), summ.get('attempt_reject_stages', '<absent>')))
    out += ['Claude: ' + l for l in text if l.startswith('ACCEPT') or l.startswith('spot check:')]
    out += ['note: ' + n for n in notes]
    print('\n'.join(out))
    tag = 'FORMAT TEST (never a sign-off): ' if ft else ''
    print(tag + ('VERDICT: CONSISTENT' if not issues else 'VERDICT: NOT CONSISTENT\n  - ' + '\n  - '.join(issues)))
    sys.exit(0 if not issues else 1)


if __name__ == '__main__':
    main()
