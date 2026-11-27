# The complete far cover `pent_far5` (user's machine, 2026-11-18/19) and its checks

Run: package `pent_far5` (T1-free), COMMIT f30df709be7ad03480fca3c6543760ca373693e9, eight processes on the author's
computer, started 11-18 10:44 CST, `ALL STAGES DONE` 11-19 01:57:46 CST; 1,260 / 1,260 blocks, no unresolved leaf.
Acceptance: `../local/accept/FAR5_USER_ACCEPT_20261119.md`.

| file | content | SHA256 |
|---|---|---|
| `pent_far5_cert.tar.xz` (13,749,176 bytes; 262 MB unpacked) | `out/`: `config.json`, `<chart>/roots.json`, every block's `leaves.jsonl`, `cert.jsonl`, `done.json`, `block.log`; `qtable/cells.jsonl`, `table_info.json` | 01cec8f9ddacf8fb39b6b188d75cb848c7d18ba2b37ba9b354d8577dd5fed174 |
| `far5_results.tgz` (1,860,051 bytes) | `run.log`, `status.txt`, `COMMIT`, `verify.txt` (empty: the self-check crashed, see below), per-block summaries | 4760db3ee7cf018700d2e177c9dc01ca5a4cc189f52bdc5b4bf617cf0ee3af4c |

**Unpacking.** Extract into an empty directory (`tar xJf pent_far5_cert.tar.xz --no-same-owner --no-same-permissions`).
The archive lacks the per-block `corners.jsonl` because the packing command omitted it (fixed in `../local/accept/README.md`).
Every `done.json` records the SHA256 of the empty file for it (eigen corners: 0 everywhere), so recreate them with
`for d in out/*/blocks/*/; do : > "$d/corners.jsonl"; done` before running `far_verify.py`.

## Checks (`sandbox_checks/`, all on these files, in the sandbox)
| log | program | result |
|---|---|---|
| `verify_edge.txt` | production `far_verify.py --replay 3 --relabel 300 --table 20 --strict-table --complete --chart E2,E3` | VERIFY OK |
| `verify_std.txt` | the same for std (before the `lg_check` fix) | everything OK except 10 table LG parts "not re-proved": the non-monotone test, see the acceptance note |
| `table_A.txt`, `table_B.txt` | the table check for std after the fix (two halves of the blocks) | VERIFY OK; 43,103 LG entries re-proved |
| `root_spot.txt` | `v2_run/tests/far_root_spot.py` (second implementation) | no missing root |
| `tiling_spot.txt` | `v2_run/tests/far_tiling_spot.py` (second implementation) | no failure |
| `label_*of8.txt` | `v2_run/tests/far_label_spot.py --sample 0 --part i/8` (second implementation) | 1,346,916 re-proved in one piece, 13 not |
| `../local/accept/tools/label_subdivide_far5_log.txt` | `local/accept/tools/label_subdivide.py` | the 13: every one re-proved after one halving into 64 pieces |
| `replay/replay_<block>.txt` | production `far_verify.py --charts std --blocks B --replay -1` | the nine std blocks around the regular pentagon (16_4 … 18_6): **all 29,805 certified boxes re-proved from scratch, VERIFY OK** (about 0.45 s per box); std blocks bordering the near ball: **all 59 blocks, 137,384 boxes, every one VERIFY OK** (complete 11-20) (13_10, 13_11, 13_12, 13_6, 13_7, 13_8, 13_9, 14_10, 14_11, 14_12, 14_5, 14_6, 14_7, 14_8, 14_9, 15_10, 15_11, 15_12, 15_13, 15_4, 15_5, 15_6, 15_7, 15_8, 15_9, 16_10, 16_11, 16_12, 16_13, 16_7, 16_8, 16_9, 17_10, 17_11, 17_12, 17_7, 17_8, 17_9, 18_10, 18_11, 18_7, 18_8, 18_9, 19_10, 19_11, 19_3, 19_4, 19_5, 19_6, 19_7, 19_8, 19_9, 20_10, 20_4, 20_5, 20_6, 20_7, 20_8, 20_9) |
| `replay_E2/E2_<block>.txt` | production `far_verify.py --charts E2 --blocks=B --replay -1` | 31 E2 blocks replayed in full (1,871 boxes), VERIFY OK; E2 runs at about 1.5–4 s per box, so a full replay of E2/E3 does not fit the sandbox; the remaining blocks of all three charts get a sampled replay (25 per block, below) |
| `replay_std_more/replay_<block>.txt` | production `far_verify.py --charts std --blocks B --replay -1` | further std blocks replayed in full, nearest to the regular pentagon first: 235 blocks, 125,119 boxes, all VERIFY OK (queue finished 11-22). With `replay/` this is 422 of the 428 std blocks (119 of them have no C box) and 292,308 of the 293,846 std C boxes; the 6 blocks left (21_6, 21_4, 12_3, 14_2, 20_2, 18_1; 1,538 boxes) hit the 2-hour job limit and only have the sampled replay below |
| `replay_sample/sample_<job>_<chart>.txt` | production `far_verify.py --charts CH --blocks=LIST --replay 25 --seed 20261120` (block lists in `sample_queue.txt`) | sampled replay of every block not replayed in full: 11 of 11 jobs done, blocks E2 353, E3 448, std 360; 16,905 certified boxes re-proved, all VERIFY OK |
| `e2e_50000_final.txt`, `e2e_hard_30000_final.txt`, `e2e_shell_30000_final.txt` | `v2_run/tests/far_e2e_check.py` (floating point, no proof): random strictly convex pentagon → kept representation → leaf of this cover (through the hand-over) → the leaf's assertion from the vertices | 110,000 pentagons (random; adversarial: tiny / full cuts, flat vertices, ties, elongated; chart distance 0.5–0.8 around the hand-over sphere): **0 violations, 0 inconclusive** |

Not covered by a second implementation: the 377,515 certified boxes (C; 2,532 of them replayed from scratch by the
production self-check), the LG parts of the quadrilateral table
(re-proved by the production code only), the closed forms of the label code (re-derived by hand and by
`tests/far_formula_check.py`). These are the second party's items 3.4 and 3.5 in `handoffs/GPT-PENT-22-20261117.md`.
