"""Fill the result / replay tables of paper/main.tex from the certificate JSON files and replay_status.json, then build the
PDF.  Usage (from the release root): python3 tools/make_paper.py VERSION DATE"""
import json, math, os, sys, subprocess
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CERT = os.path.join(ROOT, 'certificates')
XS = np.array([1.3090169943749472, 0.9510565162951535, 0.49999999999999994, 1.5388417685876266, -0.30901699437494745, 0.9510565162951536])
FSTAR_UP = 18.919134195733918
rows = [('A', 'cell_liu_point_rho0.025_v6.json', 'min-corner', '1'), ('B', 'cell_r0.1_id71_rho0.012_v6.json', 'min-corner', '3'),
        ('C', 't1_liu_point_rho0.03_rigorous_v6.json', 'T1', '2'), ('D', 't1_r0.1_id71_rho0.015_rigorous_v6.json', 'T1', '4'),
        ('E', 't1_id301_rho0.0239_rigorous_v6.json', 'T1', '5'), ('F', 't1_r0.1_id71_rho0.012_rigorous_v6.json', 'T1', '3')]
boxes = {}
for name, f, form, bx in rows:          # R01: verify the box relations claimed in the text from the certificate endpoints
    d = json.load(open(os.path.join(CERT, f)))
    boxes.setdefault(bx, []).append((name, tuple(d['cell_lo']), tuple(d['cell_hi'])))
assert len(boxes) == 5 and boxes['3'][0][1:] == boxes['3'][1][1:], 'B and F must be the same box'
def inside(a, b):
    return all(float(x) >= float(y) for x, y in zip(a[1], b[1])) and all(float(x) <= float(y) for x, y in zip(a[2], b[2]))
assert inside(boxes['1'][0], boxes['2'][0]) and inside(boxes['3'][0], boxes['4'][0]), 'nesting A in C, B=F in D'
tab = [r'\begin{tabular}{lllrrrll}', r'\toprule',
       r'cert. & box & form & $\|x_C-x_*\|_2$ & $\rho$ & certified lower bound & margin & B\&B leaves\\', r'\midrule']
for name, f, form, bx in rows:
    d = json.load(open(os.path.join(CERT, f)))
    dist = float(np.linalg.norm(np.array(d['p']) - XS))
    if form == 'min-corner':
        lb = d['LB']; leaves = '--'
    else:
        lb = d['T1']['min_lower']; leaves = '%d' % d['T1']['accepted']
    tab.append(r'%s & %s & %s & %.4f & %.4g & %.9f & %.2e & %s\\' % (name, bx, form, dist, d['rho'], lb, lb - FSTAR_UP, leaves))
tab += [r'\bottomrule', r'\end{tabular}']
res = ('There are six certificates (A--F) for five distinct boxes (1--5) at three centres. Each box is '
       '$C=\\prod_i[\\ell_i,h_i]$ with the binary64 endpoints listed in its certificate file (centre $\\pm\\rho$ in every '
       'coordinate, sup norm); $x_*$ is the regular pentagon in the chart. ``Margin'' is the certified lower bound of $F$ on '
       'the whole box minus $F_*^{\\rm up}=18.919134195733918$. Boxes 1 and 2 share a centre, the cell identified as hardest '
       'in the outward-continuation prototype of \\cite{liu}, and $1\\subset2$; boxes 3 and 4 share a centre at chart '
       'distance $0.1$ and $3\\subset4$; certificates B and F certify the same box 3 by the two methods; box 5 is centred '
       'at the hardest of 480 random centres of our own cost sweep. The nested and repeated certificates are kept as method '
       'comparisons.\n\n'
       '\\begin{center}\\small\n' + '\n'.join(tab) + '\n\\end{center}\n\n'
       'For boxes 4 and 5 the min-corner form fails (its bound is $18.9099$ resp.\\ $18.6492$) and only the T1 form certifies '
       'the box; box 5 needs the exact covariance majorant $e_k$ in addition to Lemma~\\ref{lem:es}.')
st = json.load(open(os.path.join(ROOT, 'replay_status.json')))
rt = [r'\begin{tabular}{p{6.2cm}p{3.2cm}p{5.4cm}}', r'\toprule', r'item & status & method of the second party\\', r'\midrule']
for it in st['items']:
    rt.append(r'%s & %s & %s\\' % (it['item'], it['status'], it['method']))
rt += [r'\bottomrule', r'\end{tabular}']
rep = (st['intro'] + '\n\n\\begin{center}\\small\n' + '\n'.join(rt) + '\n\\end{center}\n\n' + st['conclusion'])
cpath = os.path.join(ROOT, 'cost_estimate.json')
if os.path.exists(cpath):
    c = json.load(open(cpath))
    assert 'n_in_D' in c and all('no_pass_outside_D' in b for b in c['bands']), \
        'cost_estimate.json must come from the PENT-08 version of run3_cost.py (all sampled centres in the mean)'
    def sci(v):                                      # 3.0e6 -> 3.0\cdot10^{6}
        m, e = ('%.1e' % v).split('e')
        return r'%s\cdot10^{%d}' % (m, int(e))
    ct = [r'\begin{tabular}{lrrrrrr}', r'\toprule',
          r'band of $r$ & centres & in $D$ & median $\rho/r$ in $D$ & cells (proxy) & s/cell & CPU-days (95\%)\\', r'\midrule']
    for b in c['bands']:
        ct.append(r'$[%.2f,%.2f)$ & %d & %d & %.3f & $%s$ & %.0f & %.0f\ (%.0f--%.0f)\\' % (
            b['band'][0], b['band'][1], b['centres'], b['in_D'], b['q_median_D'], sci(b['cells']), b['sec_per_cell'],
            b['cpu_days'], b['cpu_days_ci95'][0], b['cpu_days_ci95'][1]))
    ct += [r'\bottomrule', r'\end{tabular}']
    sv = c['shift_jacobian_sv']
    tot = c['totals']; b1, b2, b3 = c['bands']
    top3 = b3['top'][0]; top1 = b1['top'][0]
    no = c.get('not_ok', [])
    assert all(not x['in_D'] for x in no), 'text below says the failed centres lie outside D'
    sh = [m for m in c['shrink_model'] if m['s'] == 0.8]
    cost = ('To estimate the cost of covering a whole region we measured, at %d centres log-uniformly distributed in chart '
            'distance $r\\in[0.02,0.6]$ (directions uniform), for which half-widths $\\rho$ the rigorous T1 branch and bound '
            'succeeds (a search on a grid of ratio $0.85$ up to the largest passing size). In this measurement the corner '
            'eigenvalues are floating-point values reduced by a relative margin $2\\cdot10^{-5}$ (a strict corner certificate '
            'loses about that much), so the numbers below are an \\emph{estimate}, not a certificate.\n\n'
            'A proof needs only one labelling of each pentagon. We use the Dirichlet domain $D$ of the relabelling action of '
            '$D_5$: the points of the chart whose own labelling is (weakly) the closest of their ten labellings to $x_*$. '
            'The action is not an isometry of the chart (at $x_*$ the Jacobian of a cyclic relabelling has singular values '
            '$%.2f,%.2f,1,1,%.3f,%.3f$), so $D$ does not occupy a tenth of a shell around $x_*$; %d of the %d centres lie in '
            '$D$. The shell $a\\le r<b$ has measure $\\pi^3r^5\\,dr$. Let $\\rho(x)=q(x)\\,r$ be the half-width chosen at $x$ '
            'by the rule below. For $r$ log-uniform on $[a,b)$ and uniform directions\n'
            '\\[ N_D=\\frac{\\pi^3}{64}\\ln\\frac ba\\;\\mathbb E\\bigl[\\mathbf 1_D\\,q^{-6}\\bigr]'
            '=\\int_{D\\cap\\{a\\le r<b\\}}(2\\rho(x))^{-6}\\,dx, \\]\n'
            'and the mean over all sampled centres of the band, with those outside $D$ counted as zero, is an unbiased Monte '
            'Carlo estimate of this integral. The integral is a continuum proxy for the number of cells, not the number of '
            'boxes of a constructed rigorous cover: cells that cross the boundary of $D$, the interfaces of an adaptive grid, '
            'overlaps and the certification overhead are not counted. The cost is estimated in the same way, with the '
            'integrand multiplied by the cost per cell: %.1f amortised strict '
            'corners (%.1f CPU seconds each, shared between neighbouring grid cells) plus the branch-and-bound time '
            'measured at that size. For each centre we use the tested passing size with the smallest cost per volume; '
            'smaller sizes were not measured and may be cheaper. Times are CPU seconds on the author\'s laptop (Intel '
            'i9-14900HX, WSL2); the intervals are 95\\%% bootstrap intervals over the centres of each band.\n\n'
            '\\begin{center}\\small\n' % (c['n_centres'], sv[0], sv[1], sv[4], sv[5], c['n_in_D'], c['n_centres'],
                                          c['corners_per_cell'], c['corner_sec'])
            + '\n'.join(ct) + '\n\\end{center}\n\n'
            'For $r\\in[0.05,0.6)$ this gives a proxy of about $%s$ cells and %.0f CPU-days (bootstrap interval %.0f--%.0f). The '
            'estimate rests on few shapes: one centre ($r=%.3f$, one edge about %.1f times shorter than the longest) '
            'carries %.0f\\%% of the band $[0.2,0.6)$, and without it the total would be about %.0f CPU-days. Shrinking '
            'every cell to $0.8$ times its largest passing size, with a branch-and-bound time proportional to '
            '$\\rho^\\beta$ ($\\beta=11.3$--$14.2$, the interquartile range fitted on pairs of passing sizes by the second party), gives %.0f--%.0f '
            'CPU-days, so a uniform change of the size does not help. Including $[0.02,0.05)$ the estimate is $%s$ cells '
            'and %.0f CPU-days, but that band is dominated by one centre ($r=%.3f$, %.0f\\%% of the band) and its interval '
            '(%.0f--%.0f) is too wide to be useful. The region $r\\ge0.6$ and a neighbourhood $r<0.02$ of the regular '
            'pentagon are not included, and a different fundamental domain, or a choice of labelling per region, could be '
            'cheaper; this was not measured.\n\n'
            '%d of the %d centres did not pass at any tested size. All of them lie at $r<%.3f$ and outside $D$ (their '
            'closest labellings are at $r=%.3f$--$%.3f$), and every failed test failed because the corner values, reduced '
            'by the $2\\cdot10^{-5}$ margin, were already below $F_*^{\\rm up}$. This concerns the tested sizes and the '
            'chosen margin, not every possible refinement. Lying outside $D$, these centres contribute zero to the estimate; '
            'a centre in $D$ that passed at no tested size would leave the estimate undefined, and the estimator stops in '
            'that case. '
            'Close to the regular pentagon the accuracy of the corner eigenvalues, not the transplantation bound, limits '
            'the method.' % (
                sci(tot['L4_cells']), tot['L4_cpu_days'], tot['L4_cpu_days_ci95'][0], tot['L4_cpu_days_ci95'][1],
                top3['r'], top3['edge_ratio'], 100 * top3['share'], tot['L4_cpu_days_without_top']['cpu_days'],
                min(m['L4_cpu_days'] for m in sh), max(m['L4_cpu_days'] for m in sh),
                sci(tot['all_cells']), tot['all_cpu_days'], top1['r'], 100 * top1['share'],
                b1['cpu_days_ci95'][0], b1['cpu_days_ci95'][1],
                c['n_not_ok'], c['n_centres'], math.ceil(1000 * max([x['r'] for x in no] or [0])) / 1000,
                min([x['closest_labelling_distance'] for x in no] or [0]),
                max([x['closest_labelling_distance'] for x in no] or [0])))
else:
    cost = 'The measurement of the cost of a cover is in progress and will be reported in this section.'
src = open(os.path.join(ROOT, 'paper', 'main.tex')).read()
src = src.replace('__COST_SECTION__', cost)
src = src.replace('__RESULTS_TABLE__', res).replace('__REPLAY_TABLE__', rep).replace('__REPLAY_ABSTRACT__', st['abstract'])
src = src.replace('__VERSION__', sys.argv[1]).replace('__DATE__', sys.argv[2])
open(os.path.join(ROOT, 'paper', 'paper.tex'), 'w').write(src)
import re                                            # the README carries the same version line
rp = os.path.join(ROOT, 'README.md'); rd = open(rp).read()
rd, nsub = re.subn(r'^Version \S+ \([^)]*\)\. Author', 'Version %s (%s). Author' % (sys.argv[1], sys.argv[2]), rd, flags=re.M)
assert nsub == 1, 'README version line not found'
open(rp, 'w').write(rd)
for _ in range(2):
    subprocess.run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error', 'paper.tex'], cwd=os.path.join(ROOT, 'paper'),
                   stdout=subprocess.DEVNULL, check=True)
print('built', os.path.join(ROOT, 'paper', 'paper.pdf'))
