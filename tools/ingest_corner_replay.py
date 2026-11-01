"""Ingest the second-party corner replay (summary_N300.json written by second_party/q1_batch.py) for cells A and B.

Checks, with our own code (the second-party checker itself is not executed here):
  * exactly the 2 x 64 corners of cells A and B, each once;
  * every record has status PASS, HKS_pass true, mu_upper < 1 and F_lo >= the producer value;
  * the producer value the checker compared against equals OUR v6 corner bound, recomputed from the certificate
    (rig_paqc.cert_F_lo, i.e. lam_tilde/(1+mu) * area rounded down) -- so the comparison was against the right numbers;
  * cell lower bounds from the independent corner values: factor * min_c F_lo(second party) > F*_up.
On success (--apply) it copies the summary to second_party/results/, updates replay_status.json and prints what to do next.
Usage (from the release root):  python3 tools/ingest_corner_replay.py PATH/summary_N300.json [--apply]"""
import os, sys, json, hashlib, shutil
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'producer'))
import rig_paqc as R                      # noqa: E402  (producer code: used only to recompute OUR corner bounds)
from flint import arb                     # noqa: E402
import numpy as np                        # noqa: E402

CELLS = {'mps_cell_liu_point_rho0.025_v6.json': ('A', 'cell_liu_point_rho0.025_v6.json'),
         'mps_cell_r0.1_id71_rho0.012_v6.json': ('B', 'cell_r0.1_id71_rho0.012_v6.json')}


def main(path, apply):
    rows = json.load(open(path))
    sha = hashlib.sha256(open(path, 'rb').read()).hexdigest()
    print('summary', path, 'sha256', sha, 'records', len(rows))
    problems = []; seen = {}
    for r in rows:
        if r.get('file') not in CELLS:
            problems.append('unexpected file %r' % r.get('file')); continue
        key = (r['file'], r['corner'])
        if key in seen:
            problems.append('duplicate %s corner %d' % key)
        seen[key] = r
    for f in CELLS:
        for c in range(64):
            if (f, c) not in seen:
                problems.append('missing %s corner %d' % (f, c))
    out = {}
    for f, (name, certfile) in CELLS.items():
        d = json.load(open(os.path.join(ROOT, 'certificates', certfile)))
        ours = [R.cert_F_lo(c) for c in d['corner_certificates']]
        recs = [seen.get((f, c)) for c in range(64)]
        npass = 0; margins = []; theirs = []
        for c, rr in enumerate(recs):
            if rr is None:
                continue
            ok = (rr.get('status') == 'PASS' and rr.get('HKS_pass') is True and float(rr.get('mu_upper', 9)) < 1
                  and float(rr['F_lo']) >= float(rr['producer_F_lo_v6']))
            if float(rr['producer_F_lo_v6']) != ours[c]:
                problems.append('%s corner %d: checker compared against %r, our v6 bound is %r' % (name, c, rr['producer_F_lo_v6'], ours[c]))
                ok = False
            if not ok:
                problems.append('%s corner %d: not PASS (%s)' % (name, c, {k: rr.get(k) for k in ('status', 'HKS_pass', 'mu_upper', 'margin')}))
            npass += ok; margins.append(float(rr['F_lo']) - ours[c]); theirs.append(float(rr['F_lo']))
        res = dict(cell=name, corners_pass=npass, corners_total=64)
        if len(theirs) == 64:
            g = R.geometry_bound(np.array(d['p']), d['rho'])
            LB = g['factor'] * R.amin([arb(v) for v in theirs])
            res.update(min_margin=min(margins), LB_second_party_corners=R.f_down(LB),
                       certified_with_second_party_corners=bool(LB > R.FSTAR_UP))
        out[name] = res
        print(json.dumps(res))
    if problems:
        print('PROBLEMS (%d):' % len(problems)); [print('  ' + p) for p in problems[:40]]
    allok = not problems and all(v['corners_pass'] == 64 and v.get('certified_with_second_party_corners') for v in out.values())
    print('ALL 128 CORNERS PASS' if allok else 'NOT COMPLETE')
    if apply and allok:
        dst = os.path.join(ROOT, 'second_party', 'results'); os.makedirs(dst, exist_ok=True)
        shutil.copy(path, os.path.join(dst, 'summary_N300_cells_AB.json'))
        st_path = os.path.join(ROOT, 'replay_status.json'); st = json.load(open(st_path))
        for it in st['items']:
            if it['item'].startswith('Corner eigenvalue lower bounds, cells A and B') or it['item'].startswith('Corner eigenvalue lower bounds, certificates A and B'):
                it['status'] = 'PASS 128/128'
                it['method'] += ('; results in \\texttt{second\\_party/results/}; smallest margin over the producer bound '
                                 '%.1e' % min(v['min_margin'] for v in out.values()))
        st['conclusion'] = ('Boxes 1 and 3 (certificates A and B) are fully replayed by the second party: geometry factors '
                            '(exact fractions), all 128 corner eigenvalue lower bounds and the regular-pentagon upper bound. '
                            'Certificates C--F are producer certificates whose whole-leaf bounds have not yet been replayed '
                            'independently, and should be regarded accordingly.')
        st['abstract'] = ('for boxes 1 and 3 (certificates A and B) every load-bearing component -- geometry factors, all 128 '
                          'corner eigenvalue bounds and the regular-pentagon bound -- has been recomputed by independent code.')
        json.dump(st, open(st_path, 'w'), indent=1, ensure_ascii=False)
        print('applied: results copied, replay_status.json updated. Next: python3 tools/make_paper.py VERSION DATE, then the MANIFEST.')
    return allok


if __name__ == '__main__':
    sys.exit(0 if main(sys.argv[1], '--apply' in sys.argv) else 1)
