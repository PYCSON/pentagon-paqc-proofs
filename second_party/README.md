# second_party

Independent checkers written and run by the second party (GPT, OpenAI). They share no source with `producer/` and
read the certificates only as data.

| file | what it checks |
|---|---|
| `independent_mps_iv.py` | corner eigenvalue lower bound for ONE corner and ONE edge: own interval Bessel series with tail bounds, third-order Taylor boundary bound, disk L² bound, FHM |
| `q1_batch.py` | driver: all corners of cells A, B (and optionally C–E), all 5 edges, then FHM + Krahn–Szegő index and comparison with the producer's F_lo. Resumable through `q1_checkpoints/` |
| `PATCH_q1_batch.diff` | the only change from the second party's original driver: the per-edge timeout was raised from 240 s to 7200 s, because a shared CPU needs longer (a timeout never counts as PASS) |
| `independent_geometry_fraction.py` | min-corner geometry factor of cells A, B in exact rational arithmetic |
| `aggregate_mps.py` | aggregation of the per-edge results of the three worst corners (round PENT-04) |
| `q4_exact_far_gate.py` | point checker for the far-region exclusion rules (not used by the certificates of this release) |
| `q1_batch_mu_guard.py`, `PATCH_q1_batch_mu_guard.diff` | the second party's later fix (round PENT-06): a corner is PASS only if μ < 1 is checked explicitly before the index inequality. The A/B run used `q1_batch.py`; its result file was additionally checked for μ < 1 on every record by `tools/ingest_corner_replay.py`. Re-running `q1_batch_mu_guard.py` on the finished checkpoints recomputes only the summary |
| `p2_fem/` | the independent P2 finite-element upper bound for the regular pentagon (see its README) |

Data layout expected by `q1_batch.py`: `./lines/D-topics/probes/PENT/G1/replay_corners/all/mps_*.json` (the path of the
research repository). In this release the data are in `certificates/replay_corners/all`; create the expected path as a
link, from this folder:
```
mkdir -p lines/D-topics/probes/PENT/G1/replay_corners
ln -s ../../../../../../../certificates/replay_corners/all lines/D-topics/probes/PENT/G1/replay_corners/all
python3 q1_batch.py --N 300 --corners 64 --files 2 --workers 4
```
Results: `q1_checkpoints/summary_N300.json`. Each corner takes about 3–5 CPU minutes; the run can be stopped and resumed.
Do not run Python with `-O`: the checkers use `assert` statements as rejection tests.
