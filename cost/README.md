# cost — data and code of the cost estimate (paper, Section 6)

This is an **estimate, not a certificate**: the measurement uses floating-point corner eigenvalues reduced by a
relative margin of 2·10⁻⁵.

| path | content |
|---|---|
| `run3_cost.py` | the estimator. `python3 cost/run3_cost.py cost/results_run3 --json cost_estimate.json` (run from the release root) rewrites `../cost_estimate.json`, from which `tools/make_paper.py` builds Section 6 |
| `results_run3/` | the measurement as returned from the author's laptop: `results.jsonl` (one line per centre: the tested half-widths, pass or fail and the reason, branch-and-bound boxes and CPU seconds), `summary.json`, `SHA256SUMS` (it also lists `run_info.json`, a machine description that was not copied) |
| `measurement/` | the package that produced it, unchanged (`sha256sum -c SHA256SUMS` inside it): `run_t1sweep.py` (driver), `params.json`, the copies of `rig_t1.py`, `rig_paqc.py`, `rig_box.py`, `mps.py` it ran with, `fstar_exact.json`, and `centres_run2.jsonl` (the 480 centres x = x\* + r·d) |

The modules in `measurement/` are the versions the sweep ran with. They are not the certificate producer: the
certificates of this release are replayed with `producer/`.

The estimator needs NumPy only. It checks the hashes of `results.jsonl` and `summary.json`, decides for each centre
whether it lies in the Dirichlet domain D of the D5 relabelling action (its own labelling is the closest of its ten
labellings to x\*), and computes the per-band cell counts, CPU-days, 95% bootstrap intervals (4,000 resamples, fixed
seed), the cost without the largest contributor, and a uniform-shrink model.
