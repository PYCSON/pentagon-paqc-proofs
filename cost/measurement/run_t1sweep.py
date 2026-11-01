"""run3: rigorous-geometry T1 cell-size sweep over the 480 run2 centres (Claude, 2026-10-31).  Own code.
For each run2 centre p, find (geometric search, ~8% resolution) the largest cube half-width rho such that the
RIGOROUS-GEOMETRY T1 test of rig_t1.py passes: ball-arithmetic convexity gate, sector Jacobians, matrix Efron-Stein eps(u)
and second-order Taylor branch-and-bound over the whole cell, but with FLOAT corner values
F_c = MPS(K=11) value * (1 - allow_rel)  (allow_rel models the loss of the rigorous corner certificates, ~2e-5).
So the result is NOT a proof; it measures the cell sizes the strict pipeline will reach, for the cost estimate.
Also reports the min-corner (old) bound at each tested rho, for comparison.
Resumable (results appended to OUT/results.jsonl, finished ids skipped).  Guards on free disk / available memory.
Usage: python3 run_t1sweep.py [params.json] [--test N]
"""
import os
for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[k] = "1"
import sys, json, time, shutil, hashlib, itertools
import multiprocessing as mp
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mps import PolyMPS, regular, poly_area
from flint import arb
import rig_paqc as R
import rig_t1 as T

FSTAR_UP = R.FSTAR_UP
P0 = regular(5); A0 = poly_area(P0)
CORN = R.CORN
PRM = None


def r5_params():
    Rg = regular(5, R=1 / (2 * np.sin(np.pi / 5)), rot=-np.pi / 2 - np.pi / 5)
    Rg = Rg - Rg[0]
    e = Rg[1] - Rg[0]; ang = np.arctan2(e[1], e[0])
    Rm = np.array([[np.cos(-ang), -np.sin(-ang)], [np.sin(-ang), np.cos(-ang)]])
    return (Rg @ Rm.T)[2:].reshape(-1)


def Ffloat(x):
    """lambda_1 * area (scale invariant) of the raw-chart pentagon x, float MPS with retries (as run2)."""
    V = np.array([[0, 0], [1, 0], [x[0], x[1]], [x[2], x[3]], [x[4], x[5]]], float)
    V = V * np.sqrt(A0 / poly_area(V))
    last = None
    for K, sd in ((11, 0), (12, 1), (13, 2)):
        try:
            m = PolyMPS(V, K=K, seed=sd)
            g = 18.95 / A0; w = 0.04
            lam = m.solve(g, width=w, ngrid=13)
            if abs(lam - g * (1 - w)) < 0.02 * w * g or abs(lam - g * (1 + w)) < 0.02 * w * g:
                lam = m.solve(27.0 / A0, width=1.0 / 3.0, ngrid=181)
            return lam * A0
        except np.linalg.LinAlgError as ex:
            last = ex
    raise last


def precheck(cd_float, nsamp, rng):
    """cheap float minimum of the T1 bound (correct global-eps form): samples (incl. corner-biased ones) + local
    refinement (L-BFGS-B) from the 4 best.  Below F*_up => this rho cannot pass, skip the B&B."""
    from scipy.optimize import minimize
    J, lam, pre, a0, Gc, Dc = cd_float

    def B(u):
        u = np.clip(u, -1, 1)
        a = np.prod((1 + CORN * u) / 2, axis=1); pq = (1 - u * u) / 4
        JP = a @ J; AP = JP @ a0
        Q = min(AP * (a @ (lam * J[:, k])) / JP[k] for k in range(5))
        e = 0.0
        for k in range(5):
            Jmax, RR, dev = pre[k]
            main = np.linalg.eigvalsh(sum(pq[i] * RR[i] for i in range(6))).max(); dv = pq @ dev
            ees = Jmax * min((1 + t) * main + (1 + 1 / t) * dv for t in (0.03, 0.1, 0.3))
            if PRM.get("emode", "es") == "both":     # exact sector e = lam_max(S^-1 (S - G_P)), as in rig_t1 emode 'both'
                Sk = np.tensordot(a, Gc[:, k], axes=1); DP = np.tensordot(a, Dc[:, k], axes=1)
                adj = np.array([[DP[1, 1], -DP[0, 1]], [-DP[1, 0], DP[0, 0]]]); GP = adj @ adj.T / np.linalg.det(DP)
                ees = min(ees, np.linalg.eigvals(np.linalg.solve(Sk, Sk - GP)).real.max())
            e = max(e, ees)
        return (1 - e) * Q
    us = [np.zeros(6)] + [rng.uniform(-1, 1, 6) for _ in range(nsamp)] + \
         [rng.choice([-1, 1], 6) * rng.uniform(0.5, 1, 6) for _ in range(nsamp // 3)]
    vals = sorted(((B(u), i) for i, u in enumerate(us)))
    best = vals[0][0]
    for _, i in vals[:4]:
        r = minimize(B, us[i], method='L-BFGS-B', bounds=[(-1, 1)] * 6, options=dict(maxiter=60))
        best = min(best, float(r.fun))
    return best


def float_cell(lo, hi, Fc):
    c0 = P0.mean(0)
    E = [np.column_stack([P0[k] - c0, P0[(k + 1) % 5] - c0]) for k in range(5)]
    Ei = [np.linalg.inv(e) for e in E]; a0 = np.array([abs(np.linalg.det(e)) / 2 for e in E])
    xs = [np.where(s > 0, hi, lo) for s in CORN]
    Ds = []
    for x in xs:
        V = np.array([[0, 0], [1, 0], [x[0], x[1]], [x[2], x[3]], [x[4], x[5]]]); c = V.mean(0)
        Ds.append([np.column_stack([V[k] - c, V[(k + 1) % 5] - c]) @ Ei[k] for k in range(5)])
    J = np.array([[np.linalg.det(D[k]) for k in range(5)] for D in Ds])
    lam = np.array(Fc) / (J @ a0)
    pre = []
    for k in range(5):
        X = [np.linalg.inv(D[k]) for D in Ds]; G = [J[c, k] * X[c] @ X[c].T for c in range(64)]
        Gcen = sum(G) / 64; eta = max(np.abs(np.linalg.eigvalsh(g - Gcen)).max() for g in G)
        lmc = np.linalg.eigvalsh(Gcen).min()
        mu = min(np.linalg.eigvals(np.linalg.solve(Gcen, g)).real.min() for g in G)
        Gl = (Gcen - eta * np.eye(2)) if (1 - eta / lmc) >= mu else mu * Gcen       # same rule as rig_paqc.lower_metric
        Li = np.linalg.inv(np.linalg.cholesky(Gl))
        RR = []; dev = []
        for i in range(6):
            bit = 1 << (5 - i); d = [X[c] - X[c ^ bit] for c in range(64) if not c & bit]; Rm = sum(d) / 32
            RR.append(Li @ Rm @ Rm.T @ Li.T); dev.append(max(np.linalg.eigvalsh(Li @ (q - Rm) @ (q - Rm).T @ Li.T).max() for q in d))
        pre.append((J[:, k].max(), RR, np.array(dev)))
    Gc = np.array([[J[c, k] * np.linalg.inv(Ds[c][k]) @ np.linalg.inv(Ds[c][k]).T for k in range(5)] for c in range(64)])
    Dc = np.array([[Ds[c][k] for k in range(5)] for c in range(64)])
    return J, lam, pre, a0, Gc, Dc


def test_rho(p, rho, rng):
    t0 = time.time()
    lo, hi = R.cell_endpoints(p, rho)
    out = dict(rho=rho)
    try:
        g = R.convexity_gate(lo, hi)
        if not (g > 0):
            out.update(ok=False, why='convexity'); return out
        Fc = [Ffloat(x) * (1 - PRM["allow_rel"]) for x in R.cell_corners(lo, hi)]
        out['minF'] = float(min(Fc))
        if min(Fc) <= float(FSTAR_UP):
            out.update(ok=False, why='corner<=F*', sec=round(time.time() - t0, 1)); return out
        gb = R.geometry_bound(p, rho)
        out['mincorner_ok'] = bool(gb['ok'] and float(gb['factor'].lower()) * min(Fc) > float(FSTAR_UP))
        if gb['ok']:
            out['mincorner_LB'] = float(gb['factor'].lower()) * min(Fc)
        pc = precheck(float_cell(lo, hi, Fc), PRM["nsamp_precheck"], rng)
        out['precheck'] = float(pc)
        if not (pc > float(FSTAR_UP)):
            out.update(ok=False, why='precheck', sec=round(time.time() - t0, 1)); return out
        cd = T.CellData(p, rho, [arb(v) for v in Fc], emode=PRM.get("emode", "es"))
        r = T.certify_T1(cd, FSTAR_UP, max_boxes=PRM["max_boxes"], min_rad=PRM["min_rad"], log_every=10 ** 9)
        out.update(ok=r['status'] == 'certified', why=r['status'], boxes=r['boxes'], bb_sec=r['sec'])
    except (AssertionError, ZeroDivisionError, ValueError, np.linalg.LinAlgError) as ex:
        out.update(ok=False, why='exception:' + repr(ex)[:120])
    out['sec'] = round(time.time() - t0, 1)
    return out


def task(rec2):
    t0 = time.time()
    i = rec2['id']
    rng = np.random.default_rng(PRM["seed"] * 1000003 + i)
    p = r5_params() + rec2['r'] * np.array(rec2['d'])
    rec = dict(id=i, r=rec2['r'], rho_run2=rec2.get('rho_pass'))
    rp = rec2.get('rho_pass') or 0.1 * rec2['r']
    tests = []; best = None; fail = None
    m = PRM["start_mult"]
    o = test_rho(p, m * rp, rng); tests.append(o)
    if o['ok']:
        best = o['rho']
        for _ in range(PRM["max_up"]):
            o = test_rho(p, best / PRM["step"], rng); tests.append(o)
            if not o['ok']:
                fail = o['rho']; break
            best = o['rho']
    else:
        fail = o['rho']
        for _ in range(PRM["max_down"]):
            o = test_rho(p, fail * PRM["step"], rng); tests.append(o)
            if o['ok']:
                best = o['rho']; break
            fail = o['rho']
    if best and fail:
        o = test_rho(p, float(np.sqrt(best * fail)), rng); tests.append(o)
        if o['ok']:
            best = o['rho']
        else:
            fail = o['rho']
    rec.update(status='ok' if best else 'no_pass', rho_T1=best, rho_fail=fail, tests=tests, sec=round(time.time() - t0, 1))
    mc = [t['rho'] for t in tests if t.get('mincorner_ok')]
    rec['rho_mincorner_best_tested'] = max(mc) if mc else None
    return rec


def init(prm):
    global PRM
    PRM = prm
    try:
        os.nice(19)
    except Exception:
        pass


def mem_available_gb():
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemAvailable:"):
                    return int(line.split()[1]) / 1048576
    except Exception:
        pass
    return 99.0


def main():
    args = sys.argv[1:]
    pfile = args[0] if args and not args[0].startswith("--") else os.path.join(HERE, "params.json")
    prm = json.load(open(pfile))
    centres = [json.loads(l) for l in open(os.path.join(HERE, prm["centres"]))]
    centres = [c for c in centres if c.get('status') == 'ok']
    if "--test" in args:
        n = int(args[args.index("--test") + 1])
        centres = centres[:n]; prm = dict(prm, out_dir="out_test", workers=min(2, prm["workers"]))
    if "--ids" in args:                      # e.g. --ids 301,433  (test specific centres)
        want = {int(x) for x in args[args.index("--ids") + 1].split(",")}
        centres = [c for c in centres if c["id"] in want]; prm = dict(prm, out_dir="out_ids", workers=min(2, prm["workers"]))
    out = os.path.join(HERE, prm["out_dir"]); os.makedirs(out, exist_ok=True)
    res = os.path.join(out, "results.jsonl")
    done = set()
    if os.path.exists(res):
        for line in open(res):
            try:
                done.add(json.loads(line)["id"])
            except Exception:
                pass
    todo = [c for c in centres if c['id'] not in done]
    json.dump(dict(prm, started=time.strftime("%F %T"), resumed_with_done=len(done), n_centres=len(centres),
                   fstar_up=float(FSTAR_UP)), open(os.path.join(out, "run_info.json"), "w"), indent=1)
    print(f"total {len(centres)}  already done {len(done)}  to do {len(todo)}  workers {prm['workers']}", flush=True)
    if todo:
        ctx = mp.get_context("spawn")
        pool = ctx.Pool(prm["workers"], initializer=init, initargs=(prm,))
        fh = open(res, "a")
        try:
            for rec in pool.imap_unordered(task, todo):
                fh.write(json.dumps(rec) + "\n"); fh.flush(); os.fsync(fh.fileno())
                print(time.strftime("%T"), rec["id"], rec["status"], rec.get("rho_T1"), rec.get("rho_run2"), rec.get("sec"), flush=True)
                free_gb = shutil.disk_usage(out).free / 2 ** 30
                mem_gb = mem_available_gb()
                if free_gb < prm["min_free_disk_gb"] or mem_gb < prm["min_mem_available_gb"]:
                    open(os.path.join(out, "STOPPED_guard.txt"), "w").write(
                        f"{time.strftime('%F %T')} stopped safely: free disk {free_gb:.2f} GB, MemAvailable {mem_gb:.2f} GB\n")
                    print("GUARD STOP", free_gb, mem_gb, flush=True)
                    pool.terminate(); fh.close(); sys.exit(3)
        finally:
            fh.close()
        pool.close(); pool.join()
    recs = [json.loads(l) for l in open(res)]
    ids = sorted({r["id"] for r in recs})
    from collections import Counter
    summ = dict(finished=time.strftime("%F %T"), n_records=len(recs), n_unique=len(ids), complete=len(ids) == len(centres),
                status=dict(Counter(r["status"] for r in recs)))
    ok = [r for r in recs if r["status"] == "ok"]
    if ok:
        q = np.array([r["rho_T1"] / r["r"] for r in ok])
        summ["q_T1_quantiles_5_25_50_75_95"] = [float(x) for x in np.percentile(q, [5, 25, 50, 75, 95])]
        summ["q_T1_eff"] = float(np.mean(q ** -6.0) ** (-1 / 6))
        ratio = np.array([r["rho_T1"] / r["rho_run2"] for r in ok if r.get("rho_run2")])
        if len(ratio):
            summ["rho_T1_over_run2_median"] = float(np.median(ratio))
    json.dump(summ, open(os.path.join(out, "summary.json"), "w"), indent=1)
    with open(os.path.join(out, "SHA256SUMS"), "w") as f:
        for name in ("results.jsonl", "summary.json", "run_info.json"):
            h = hashlib.sha256(open(os.path.join(out, name), "rb").read()).hexdigest()
            f.write(f"{h}  {name}\n")
    open(os.path.join(out, "DONE.txt"), "w").write(time.strftime("%F %T") + "\n")
    print("DONE", json.dumps(summ), flush=True)


if __name__ == "__main__":
    main()
