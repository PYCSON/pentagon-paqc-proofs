# Acceptance of the user's complete pent_far5 run (Claude, 2026-11-19)

Run: `pent_far5` (T1-free) on the user's machine, COMMIT f30df709be7ad03480fca3c6543760ca373693e9, started 11-18 10:44 CST, `ALL STAGES DONE` 11-19 01:57:46 CST.

| file | SHA256 | size |
|---|---|---|
| far5_results.tgz (summary) | 4760db3ee7cf018700d2e177c9dc01ca5a4cc189f52bdc5b4bf617cf0ee3af4c | 1,860,051 |
| pent_far5_cert.tar.xz (full certificate: leaves.jsonl, cert.jsonl, done.json of every block, roots, qtable) | 01cec8f9ddacf8fb39b6b188d75cb848c7d18ba2b37ba9b354d8577dd5fed174 | 13,749,176 (262 MB unpacked) |

## The self-check on the user's machine did not run
`far_verify.py` ended with a segmentation fault after 5 minutes (run.log line `./far5_run.sh: line 47: 4800 Segmentation fault`);
`verify.txt` is empty.  `far5_run.sh` nevertheless printed `ALL STAGES DONE` (it does not test the exit code of far_verify: to be fixed
in the next package).  The self-check is therefore re-run in the sandbox on the user's data (below).

## Results
1. Summary (`accept_far_summary.py`): E2 384/384, E3 448/448, std 428/428 done, complete root lists; CPU 32.9 + 25.1 + 58.1 h;
   labels C 377,515, X 578,965, T 461,545, N 85,147, S 64,950, TX 60,758, O 58,040, E3 18,266, E2 14,798, V 2,336, R 2,124;
   **unresolved U: 0**.  Compared with the 204 sandbox reference blocks: 156 identical files, 48 with different files but the
   same label counts (float pre-filters), 0 with different counts.
2. Second tiling check (`tests/far_tiling_spot.py`, no production module) on the full certificate: 1,260 blocks, 196,798 columns,
   1,724,444 boxes (= the label total), **no failure**.
3. Second root-range check (`tests/far_root_spot.py`) on the user's roots.json: no cell of a strictly convex Q missing
   (std 0, E2 0, E3 0); std has 1,592 unnecessary roots (harmless); no root with other bounds.
4. Production self-check (`far_verify.py --replay 3 --relabel 300 --table 20 --strict-table --complete`) in the sandbox on the
   user's data: **VERIFY OK on all three charts** (std after the `lg_check` fix; details below).
5. Second label check (`tests/far_label_spot.py --sample 0`) on all 1,346,929 non-C leaves: **all re-proved**
   (1,346,916 in one piece, 13 after one halving of all six coordinates into 64 pieces; 0 not re-proved, 0 not re-checked).

6. End-to-end test (`tests/far_e2e_check.py`, floating point): 110,000 strictly convex pentagons (50,000 random, 30,000
   adversarial, 30,000 with chart distance 0.5–0.8) followed through the kept representation into a leaf of this cover
   (and through the hand-over into E2 / E3); every final leaf's assertion holds at the point (pointwise lower bound of F
   above F*_up, or chart distance < 3/5 for N); no point in X, V or outside the cover. Logs in
   `../../far5_run_20261119/sandbox_checks/e2e_*_final.txt`.

## Verdict
The user's complete far5 run is accepted: summary, tiling, root range, production self-check and a second implementation of
every leaf label all pass.  Remaining for Theorem F: GPT's review and independent checker (PENT-22 + addendum).

## 11-19, sandbox self-check on the user's data (production `far_verify.py`, same options as far5_run.sh, per chart)
The certificate archive did not contain the per-block `corners.jsonl` (the packing command omitted them; fixed in README.md).
Every done.json records for it the SHA256 of the empty file (eigen corners 0 everywhere), so empty files were recreated
and the hash check passes.
- E2: 384 blocks, 31,292 columns, 233,627 boxes, 753 certificates replayed, 92,515 labels recomputed, root range complete: no problem.
- std: 428 blocks, 114,186 columns, 1,228,647 boxes, 904 certificates replayed, 117,456 labels recomputed, root range complete,
  8,288 table LG entries re-proved: **10 problems, all of one kind**: "LG bound L of the part ... not re-proved"
  (roots of std 11_10, 14_3, 16_4, 2_6, 5_10, 5_9, 6_7, 6_9 (two), 8_9).
  Diagnosis: the table bounds are reproduced bit for bit by the producer here (`qtable.lg_lower`), and the 13,969 table cells
  computed in the sandbox batches equal the user's exactly, so it is not a machine difference.  The producer proves
  F > L_found and records L = (1 - 2^-10) L_found; the verifier `qtable.lg_check` tests the exact sufficient condition
  P1 - L' P2 > 0 only at L' = L, and that test is not monotone in L'.  At every one of the 10 parts the same exact test passes
  for L' = L (1 + k 2^-14) with the smallest k between 1 and 16 (1, 3, 4, 5, 5, 6, 11, 12, 15, 16), which proves F > L' >= L: **all 10 bounds hold**.
  `lg_check` now also tries L' = L (1 + k 2^-14), k = 1..63 (`v2_run/qtable.py`).  `far5_run.sh` now stops with
  "SELF-CHECK NOT PASSED" when far_verify fails (the user's run printed ALL STAGES DONE after a segmentation fault).
- Full table check (all std blocks in two halves, patched verifier, `--table -1 --strict-table`): 57,093 roots with a bound,
  **all 43,103 LG table entries re-proved**, no tangent-bound entry, VERIFY OK in both halves (tiling of all 1,228,647 std
  boxes checked again).
- E3 (job together with E2): 448 blocks, 51,320 columns, 262,170 boxes, 875 certificates replayed, 113,358 labels recomputed,
  root range complete: no problem; **E2,E3 VERIFY OK**.

So the production self-check of far5_run.sh holds on the user's complete data, with the table check widened from 20
entries per block to all of them and `lg_check` patched as above (the only change).
- Second label check (`tests/far_label_spot.py --sample 0`, 8 parts):
  - part 1/8 (105 min): 190,992 leaves re-proved, 1 not reproduced: E3 block 7_3, w-box key [7,3,-1,0], j = 0, label T,
    w in [7/16,1/2]x[3/16,1/4]x[-1/16,0]x[1/32,1/16], t in [0,1/2], u in [0,1]. The checker's own quadrilateral bound there is
    22.027 (fibre) against a need of 22.399 (sup tu = 1/2, tau_bar = 0.3107); the producer uses a sharper bound.
    Re-proved independently by halving all six coordinates once: all 64 pieces pass the same T test
    (`tools/label_subdivide.py`, exact rationals). So 190,993 / 190,993 re-proved.
  - part 0/8 (110 min): 210,525 leaves re-proved, 0 not reproduced.
  - part 3/8 (103 min): 192,671 re-proved in one piece, 3 not reproduced (all E3, label T): block 16_5 key [16,5,8,0] j = 3
    (two leaves, t-halves [1/2,1] x u [0,1/2] and t [0,1/2] x u [0,1]) and block 21_11 key [21,11,10,1] j = 0.
    Each re-proved by one halving of all six coordinates (64 pieces, all pass): 192,674 / 192,674.
  - part 2/8 (112 min): 196,662 re-proved in one piece, 2 not reproduced (chart E2): block 15_8 key [15,8,11,1] j = 3 label TX,
    and block 1_4 key [1,4,9,0] j = 0 label T. Each re-proved by one halving of all six coordinates (64 pieces, all pass):
    196,664 / 196,664.
  - part 5/8 (85 min): 153,910 re-proved in one piece, 4 not reproduced (label T): E2 block -9_11 key [-9,11,7,0] j = 0 (two leaves)
    and E3 block 15_7 key [15,7,5,0] j = 3 (two leaves). Each re-proved by one halving of all six coordinates (64 pieces):
    153,914 / 153,914.
  - part 4/8 (95 min): 178,283 re-proved in one piece, 1 not reproduced: E3 block 16_6 key [16,6,7,1] j = 0 label T;
    re-proved by one halving (64 pieces): 178,284 / 178,284.
  - part 6/8 (110 min): 97,760 re-proved in one piece, 1 not reproduced: E2 block 8_4 key [8,4,1,1] j = 0 label T;
    re-proved by one halving (64 pieces): 97,761 / 97,761.
  - part 7/8 (114 min): 126,101 re-proved in one piece, 1 not reproduced: E2 block -4_5 key [-4,5,8,1] j = 3 label T;
    re-proved by one halving (64 pieces): 126,102 / 126,102.
  - **Total: 1,346,929 / 1,346,929 re-proved** (per chart and label: std E2 14,798, E3 18,266, N 85,147, O 17,522, R 531,
    S 19,832, T 324,858, TX 15,120, V 83, X 438,644; E2 O 21,002, R 314, S 11,568, T 59,179, TX 21,935, X 71,825;
    E3 O 19,516, R 1,279, S 33,550, T 77,508, TX 23,703, V 2,253, X 68,496).  The 13 subdivided leaves: 12 T, 1 TX, all in E2/E3.
  - Command: `python3 -I tools/label_subdivide.py TESTS OUT E3 BLOCK KEY J T 2` (log of all 13: `tools/label_subdivide_far5_log.txt`).
