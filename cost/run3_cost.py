"""Cost of a rigorous PAQC cover of the near region from the run3 measurements (Claude, 2026-11-01, revised after
GPT PENT-07 NEW-01 and PENT-08 NEW-03/04).

Input: the run3 output folder (results.jsonl, summary.json, SHA256SUMS from the user's machine) and the run2 centre
file (directions d; the centre is x = x_* + r d).
Model (an ESTIMATE, not a proof: run3 uses float corner eigenvalues minus a 2e-5 relative margin):
  * centre i at chart distance r_i has certified T1 half-widths rho tested on a ~9% grid; q = rho / r;
  * the centres are log-uniform in r on [0.02, 0.6] with uniform directions, so for any region W of the chart and the
    size field rho(x) defined by the rule below, (pi^3/64) ln(b/a) mean_band(1_W q^-6), averaged over ALL sampled
    centres of the band, is an unbiased Monte Carlo estimate of the continuum proxy  int_{W cap shell} (2 rho(x))^-6 dx
    (shell measure pi^3 r^5 dr).  It is not the number of boxes of a constructed cover (boundary crossings, grid
    interfaces, overlaps and certification overhead are not counted).  A centre outside W contributes 0 whether or not
    it passed; a centre inside W that did not pass makes the estimate undefined, and the program stops (GPT PENT-08
    NEW-03: rc2 dropped the three failed centres, all outside D, from the denominator, +2.2% on [0.02, 0.05));
  * a proof only needs one labelling of each pentagon.  We take W = D, the Dirichlet domain of the D5 relabelling
    action: the points whose identity labelling is (weakly) the closest of their 10 labellings to x_*.  The action is
    NOT an isometry of the chart (at x_* the Jacobian of a cyclic shift has singular values 2.66, 2.66, 1, 1, 0.375,
    0.375), so D does not take 1/10 of a shell; the earlier '/10' discount (cost_estimate rc1) was wrong.  About 71%
    of the sampled centres lie in D;
  * cost per cell = corners amortised on a grid (C corners per cell, t_c CPU seconds per strict corner) + the
    branch-and-bound time measured at that size; per centre the tested passing size with the smallest cost per
    volume is used (sizes below the tested ones were not measured and may be cheaper);
  * uncertainty: 95% bootstrap interval over centres (resampled within each band); a uniform-shrink model
    s * rho_max with bb(s) = bb(rho_max) s^beta (beta from GPT PENT-07: quartiles 11.26 / 12.88 pooled / 14.20).
Usage: python3 run3_cost.py OUT_DIR [--centres FILE (default: centres_run2.jsonl next to this script or in measurement/)] [--corner-sec 19.5] [--corners-per-cell 1.5] [--boot 4000]
                           [--json cost.json]"""
import sys, os, json, hashlib, math
import numpy as np

XS = np.array([1.3090169943749472, 0.9510565162951535, 0.49999999999999994, 1.5388417685876266,
               -0.30901699437494745, 0.9510565162951536])
BANDS = [(0.02, 0.05), (0.05, 0.2), (0.2, 0.6)]


def orbit(x):
    """The 10 chart points of the pentagon x under relabelling (V0 = 0, V1 = 1, CCW); orbit(x)[0] is x itself."""
    z = np.array([0, 1, x[0] + 1j * x[1], x[2] + 1j * x[3], x[4] + 1j * x[5]])
    out = []
    for k in range(5):
        for rev in (False, True):
            idx = [(k - j) % 5 for j in range(5)] if rev else [(k + j) % 5 for j in range(5)]
            w = z[idx]
            u = (w - w[0]) / (w[1] - w[0])
            if rev:
                u = np.conj(u)
            out.append(np.array([u[2].real, u[2].imag, u[3].real, u[3].imag, u[4].real, u[4].imag]))
    return out


def in_D(x):
    d = [float(np.linalg.norm(y - XS)) for y in orbit(x)]
    return d[0] <= min(d) + 1e-15, min(d)


def shift_jacobian_sv():
    J = np.zeros((6, 6)); h = 1e-6
    for a in range(6):
        e = np.zeros(6); e[a] = h
        J[:, a] = (orbit(XS + e)[2] - orbit(XS - e)[2]) / (2 * h)
    return np.linalg.svd(J)[1], float(np.linalg.det(J))


def main():
    a = sys.argv[1:]
    out = a[0]
    here = os.path.dirname(os.path.abspath(__file__))
    cfile = a[a.index('--centres') + 1] if '--centres' in a else next(
        f for f in (os.path.join(here, 'centres_run2.jsonl'), os.path.join(here, 'measurement', 'centres_run2.jsonl'))
        if os.path.exists(f))
    tc = float(a[a.index('--corner-sec') + 1]) if '--corner-sec' in a else 19.5
    cpc = float(a[a.index('--corners-per-cell') + 1]) if '--corners-per-cell' in a else 1.5
    nboot = int(a[a.index('--boot') + 1]) if '--boot' in a else 4000
    jpath = a[a.index('--json') + 1] if '--json' in a else None
    shas = {}
    if os.path.exists(os.path.join(out, 'SHA256SUMS')):
        for line in open(os.path.join(out, 'SHA256SUMS')):
            h, name = line.split()
            shas[name] = h
    for name in ('results.jsonl', 'summary.json'):
        p = os.path.join(out, name)
        if os.path.exists(p) and name in shas:
            h = hashlib.sha256(open(p, 'rb').read()).hexdigest()
            print('%-14s sha256 %s  %s' % (name, h[:16], 'OK' if h == shas[name] else 'MISMATCH'))
    recs = {}
    for line in open(os.path.join(out, 'results.jsonl')):
        x = json.loads(line); recs[x['id']] = x          # last record wins
    cent = {}
    for line in open(cfile):
        c = json.loads(line); cent[c['id']] = c
    sv, det = shift_jacobian_sv()
    print('cyclic shift Jacobian at x_*: singular values', np.round(sv, 4), 'det %.4f' % det)
    base = cpc * tc
    rows = []                                             # one per centre
    for i, x in sorted(recs.items()):
        c = cent[i]
        assert abs(c['r'] - x['r']) < 1e-12
        p = XS + c['r'] * np.array(c['d'])
        inD, dmin = in_D(p)
        z = np.array([0, 1, p[0] + 1j * p[1], p[2] + 1j * p[3], p[4] + 1j * p[5]]); e = np.abs(np.roll(z, -1) - z)
        best = None; rmax = None
        for t in x['tests']:
            if not t.get('ok'):
                continue
            qq = t['rho'] / x['r']; tt = base + (t.get('bb_sec') or 0.0)
            if best is None or tt / qq ** 6 < best[0]:
                best = (tt / qq ** 6, qq, t.get('bb_sec') or 0.0)
            if rmax is None or t['rho'] > rmax[0]:
                rmax = (t['rho'], t.get('bb_sec') or 0.0)
        rows.append(dict(id=i, r=x['r'], edge_ratio=float(e.max() / e.min()), ok=x.get('status') == 'ok' and best is not None, inD=bool(inD), dmin=dmin,
                         dens=best[0] if best else None, q=best[1] if best else None, bb=best[2] if best else None,
                         qmax=rmax[0] / x['r'] if rmax else None, bbmax=rmax[1] if rmax else None))
    nD = sum(w['inD'] for w in rows)
    print('centres %d, ok %d, no_pass %d; in D %d (%.1f%%)' % (len(rows), sum(w['ok'] for w in rows),
                                                                sum(not w['ok'] for w in rows), nD, 100 * nD / len(rows)))
    K = math.pi ** 3 / 64
    rng = np.random.default_rng(20261101)
    bands = []; boot_tot = {}
    def integrand(sel):
        """per sampled centre: (1_D, 1_D q^-6, 1_D * cost density, bb); outside D -> zeros; in D and failed -> stop"""
        bad_in = [w['id'] for w in sel if w['inD'] and not w['ok']]
        if bad_in:
            raise SystemExit('estimate undefined: centres %s lie in D but passed at no tested size' % bad_in)
        ind = np.array([1.0 if w['inD'] else 0.0 for w in sel])
        qi = np.array([w['q'] ** -6.0 if w['inD'] else 0.0 for w in sel])
        de = np.array([w['dens'] if w['inD'] else 0.0 for w in sel])
        bb = np.array([w['bb'] if w['inD'] else 0.0 for w in sel])
        return ind, qi, de, bb
    for lo_, hi_ in BANDS:
        sel = [w for w in rows if lo_ <= w['r'] < hi_]
        n = len(sel); L = math.log(hi_ / lo_)
        ind, qi6, de, bb = integrand(sel)
        cells = K * L * float(np.mean(qi6))
        days = K * L * float(np.mean(de)) / 86400
        okw = [w for w in sel if w['ok']]
        full = None if len(okw) < n else (K * L * float(np.mean([w['q'] ** -6.0 for w in okw])),
                                          K * L * float(np.mean([w['dens'] for w in okw])) / 86400)
        bbw = float(np.sum(qi6 * bb) / np.sum(qi6))
        bs = np.array([np.mean(de[rng.integers(0, n, n)]) for _ in range(nboot)]) * K * L / 86400
        boot_tot[(lo_, hi_)] = bs
        order = np.argsort(-de)
        top = [dict(id=sel[j]['id'], r=sel[j]['r'], edge_ratio=sel[j]['edge_ratio'], share=float(de[j] / np.sum(de)))
               for j in order[:3]]
        qD = [w['q'] for w in sel if w['inD']]
        bands.append(dict(band=[lo_, hi_], centres=n, passed=len(okw), in_D=int(ind.sum()),
                          q_median_D=float(np.median(qD)), cells=cells,
                          bb_sec_weighted=bbw, sec_per_cell=base + bbw, cpu_days=days,
                          cpu_days_ci95=[float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                          cells_full_shell=full[0] if full else None, cpu_days_full_shell=full[1] if full else None,
                          top=top, no_pass_outside_D=n - len(okw)))
    tot = {}
    for key in ('cells', 'cpu_days'):
        tot['all_' + key] = sum(x[key] for x in bands)
        tot['L4_' + key] = sum(x[key] for x in bands if x['band'][0] >= 0.05)
    l4b = boot_tot[(0.05, 0.2)] + boot_tot[(0.2, 0.6)]
    tot['L4_cpu_days_ci95'] = [float(np.percentile(l4b, 2.5)), float(np.percentile(l4b, 97.5))]
    # sensitivity: leave the single largest contributor of [0.2,0.6) out (biased low; shows the dependence on one shape)
    b3 = bands[2]; worst = b3['top'][0]
    sel = [w for w in rows if 0.2 <= w['r'] < 0.6 and w['id'] != worst['id']]
    d3 = K * math.log(3) * float(np.mean(integrand(sel)[2])) / 86400
    tot['L4_cpu_days_without_top'] = dict(id=worst['id'], cpu_days=bands[1]['cpu_days'] + d3)
    # uniform-shrink model (GPT PENT-07): cell s * rho_max, bb(s) = bb(rho_max) s^beta
    model = []
    for s in (0.7, 0.8, 0.9, 1.0):
        for beta in (11.26, 12.88, 14.20):
            dd = 0.0
            for lo_, hi_ in BANDS[1:]:
                sel = [w for w in rows if lo_ <= w['r'] < hi_]
                integrand(sel)                      # same validity check
                dd += K * math.log(hi_ / lo_) * np.mean([(base + w['bbmax'] * s ** beta) / (s * w['qmax']) ** 6
                                                         if w['inD'] else 0.0 for w in sel]) / 86400
            model.append(dict(s=s, beta=beta, L4_cpu_days=float(dd)))
    bad = [w for w in rows if not w['ok']]
    res = dict(model=__doc__.split('Usage')[0].strip(), corner_sec=tc, corners_per_cell=cpc, n_centres=len(rows),
               n_ok=len(rows) - len(bad), n_not_ok=len(bad), n_in_D=nD, shift_jacobian_sv=[float(v) for v in sv],
               bands=bands, totals=tot, shrink_model=model,
               not_ok=[dict(id=w['id'], r=w['r'], in_D=w['inD'], closest_labelling_distance=w['dmin'],
                            whys=sorted({t.get('why') for t in recs[w['id']]['tests']}),
                            smallest_rho_tested=min(t['rho'] for t in recs[w['id']]['tests'])) for w in bad])
    for x in bands:
        print('r in [%.2f, %.2f): %3d centres (%3d in D), q median(D) %.3f, cells(D) %.2e, %.0f s/cell, %.0f CPU-days '
              '[95%% %.0f, %.0f]; full shell %s' % (x['band'][0], x['band'][1], x['centres'], x['in_D'], x['q_median_D'],
              x['cells'], x['sec_per_cell'], x['cpu_days'], x['cpu_days_ci95'][0], x['cpu_days_ci95'][1],
              '%.0f' % x['cpu_days_full_shell'] if x['cpu_days_full_shell'] is not None else 'undefined (no_pass)'), ' top', [(t['id'], round(t['share'], 3)) for t in x['top']])
    print('L4 band r in [0.05, 0.6), domain D: cells %.2e, %.0f CPU-days [95%% %.0f, %.0f]; without id %d: %.0f' % (
        tot['L4_cells'], tot['L4_cpu_days'], tot['L4_cpu_days_ci95'][0], tot['L4_cpu_days_ci95'][1], worst['id'],
        tot['L4_cpu_days_without_top']['cpu_days']))
    print('r in [0.02, 0.6): cells %.2e, %.0f CPU-days' % (tot['all_cells'], tot['all_cpu_days']))
    for m in model:
        print('  shrink s=%.1f beta=%.2f: %.0f CPU-days' % (m['s'], m['beta'], m['L4_cpu_days']))
    for w in res['not_ok']:
        print('  no_pass id %d r %.4f in D %s, closest labelling at %.4f' % (w['id'], w['r'], w['in_D'],
                                                                           w['closest_labelling_distance']))
    if jpath:
        json.dump(res, open(jpath, 'w'), indent=1)
        print('wrote', jpath)


if __name__ == '__main__':
    main()
