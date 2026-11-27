# far_theorem — the pentagons far from the regular one (part of version 2)

Status of the statement: computer-assisted proof candidate, not peer reviewed; re-certified in full by an independently
written program of a second-party AI reviewer and signed by it (PENT-26, 2026-11-28) (see "Second party" below).
Author: Puyi Sun (ORCID [0009-0004-1644-5300](https://orcid.org/0009-0004-1644-5300)). Roles of AI systems: see below.

## Result
F(Ω) = λ₁(Ω)·|Ω|, F\* = 18.919134195733914 (exact terminating decimal; the common comparison threshold, not the value F(P_reg) ≈ 18.919104); C = the binary64 chart point of the regular pentagon (as in `near_theorem/`).

> **Theorem F.** Every strictly convex pentagon P has F(P) > F\*, or a labelling of its vertices whose chart point z
> (P congruent to the polygon (0,0), (1,0), (z₁,z₂), (z₃,z₄), (z₅,z₆) up to scaling) satisfies |z − C| < 3/5.

Together with Theorem N (the annulus 1/20 ≤ |z − C| < 3/5, `../near/`) and Theorem L (the ball around the regular
pentagon, `../../local_theorem/`) this gives F(P) ≥ F(P_reg) for every convex pentagon, with equality only for pentagons
similar to the regular one (the assembly is Section 1 of the paper).

## Files
| path | content |
|---|---|
| `proof/FAR_THEOREM.md` | the cut chart, the kept representation, the three charts, the labels and their criteria, the cover argument, file formats |
| `proof/FAR_PROOFS.md` | proofs of the criteria; `proof/FIBRE_LEMMA.md` (one-dimensional fibre bound); `proof/LG_LEMMA.md` (the box certificate) |
| `certificate/pent_far5_cert.tar.xz` | the cover (262 MB unpacked): root lists of the three charts; for each of the 1,260 blocks the leaves (box + label) and the certified boxes; the quadrilateral table |
| `certificate/far5_results.tgz` | run log, status and per-block summaries of the run |
| `producer/pent_far5.zip` | the package that produced the cover (commit f30df709), with its run script and README; a public copy in which local paths in 14 text members (the run README and probe and test scripts) were replaced, the production programs unchanged; `producer/ABOUT_pent_far5.md` gives the SHA-256 of the original and of every changed member |
| `producer/qtable.py` | `qtable.py` with the corrected table check of 19 November 2026 (see below); replaces the file of the package for the check |
| `checkers/` | second implementations (no producer module): `far_root_spot.py`, `far_tiling_spot.py`, `far_label_spot.py`, `label_subdivide.py`; the symbolic check of the closed forms `far_formula_check.py` (needs SymPy) |
| `logs/` | logs of every check listed below |
| `second_party/replay_status.json` | what has been checked, by which program, on which part of the data |
| `MANIFEST-SHA256.txt` | SHA-256 of every file of this directory |

## The run
`pent_far5` on the author's computer, 18–19 November 2026, eight processes: 1,260 / 1,260 blocks (std 428, E2 384,
E3 448), 116 CPU hours, **no unresolved leaf**. 1,724,444 boxes: 377,515 certified by the box certificate (label C) and
1,346,929 leaves with closed-form labels (X 578,965, T 461,545, N 85,147, S 64,950, TX 60,758, O 58,040, E3 18,266,
E2 14,798, V 2,336, R 2,124). The quadrilateral table has 57,093 cells (43,103 with Lehmann–Goerisch parts).

## Checks done
| check | program | data | result |
|---|---|---|---|
| summaries of 204 blocks against a second computation | the same package in a second environment | 204 blocks | 156 byte-identical, 48 with the same label counts (floating-point heuristics choose different boxes), 0 different |
| root ranges complete | `checkers/far_root_spot.py` (second implementation) | all | no missing root |
| tiling | `checkers/far_tiling_spot.py` (second implementation) | all | no failure |
| every closed-form label | `checkers/far_label_spot.py --sample 0` (second implementation) | all 1,346,929 | 1,346,916 re-proved directly, 13 after one bisection of all six coordinates (`label_subdivide.py`) |
| producer self-check | `far_verify.py --replay 3 --relabel 300 --table 20 --strict-table --complete` | all | VERIFY OK (std with `producer/qtable.py`) |
| table | the same, `--table` over every root | all 43,103 entries | re-proved |
| end to end (floating point) | `far_e2e_check.py` (uses the producer modules; `logs/e2e_*`): random convex pentagon → leaf → assertion | 110,000 pentagons (random, adversarial, around the hand-over sphere) | 0 violations |

The table check in the package tests positivity of P₁ − L·P₂ only at the recorded L (= (1 − 2⁻¹⁰) × the value found).
That positivity test is sufficient but not monotone in L, and 10 entries failed it at the recorded value while passing at
L(1 + k·2⁻¹⁴), k ≤ 16, which still proves F > L. `producer/qtable.py` tries k = 0…63.

Not yet re-proved by a second implementation: the 377,515 certified boxes (a sample of 2,532 was replayed by the
producer) and the Lehmann–Goerisch parts of the table.

## How to verify
Tested with Python 3.13 / 3.14, python-flint 0.9.0, NumPy, SciPy (Linux x86_64, WSL2 Ubuntu). Do not run Python with `-O`.
```
sha256sum -c MANIFEST-SHA256.txt
mkdir pkg && cd pkg && unzip ../producer/pent_far5.zip && cp ../producer/qtable.py pent_far5/PENT/v2_run/
tar xJf ../certificate/pent_far5_cert.tar.xz -C pent_far5 --no-same-owner --no-same-permissions    # pent_far5/out
for d in pent_far5/out/*/blocks/*/; do : > "$d/corners.jsonl"; done   # empty corner files, not in the archive (logs/README_RUN_DATA.md)
# second implementations, standard library only
python3 -I ../checkers/far_root_spot.py pent_far5/out                   # about 2 minutes
python3 -I ../checkers/far_tiling_spot.py pent_far5/out                 # about 5 minutes
python3 -I ../checkers/far_label_spot.py pent_far5/out --sample 0 --part 0/8   # one eighth of the leaves, 1.5-2 hours; parts 0/8 ... 7/8
python3 -I ../checkers/far_label_spot.py pent_far5/out --chart E2 --sample 20   # 20 random leaves per E2 block: 20 minutes
# producer self-check (python-flint 0.9.0, NumPy, SciPy)
cd pent_far5/PENT/v2_run && python3 far_verify.py ../../out --replay 3 --relabel 300 --table 20 --strict-table --complete
# a quick try on one block: python3 far_verify.py ../../out --charts E2 --blocks=-10_10 --replay 2 --relabel 50   (13 seconds)
```
Replayed from a copy of this directory on 2026-11-19: manifest OK; root ranges, tiling, the E2 label sample (7,161 leaves) and
the one-block producer check pass with the times given.

## Second party
Lemma LG and the box certificate were reviewed and re-implemented by a second party (review PENT-21); the cut chart,
the kept representation, the corner principle, the ear comparisons, the fibre lemmas and the cover argument in
PENT-22, the label and table criteria with a ledger of statements and proofs in PENT-24b/24c. The second party's
independently written exact program re-proved on the author's computer all 377,515 certified boxes, all 185,100 table
entries and all 1,346,929 labels (22–24 November 2026); every stored witness was replayed (finished 27 November 2026),
and the producer's acceptance program agrees with the replay shard by shard. Theorem F is signed by the second party
(PENT-26, 2026-11-28). `second_party/replay_status.json` and `../second_party/final_reports/` record the details. This is a review by
an AI system, not peer review.

## Roles of AI systems
Claude (Anthropic) found the method, wrote the proofs, the programs and the checks listed above (including the second
implementations), and ran the sample computations. GPT (OpenAI) reviews the mathematics and writes independent
programs. The author initiated and coordinated the work, ran the computation on his own computer and decides on the
release; he did not establish the mathematical correctness by hand. Responsibility for the publication lies with the
author.

## License
Text: CC BY 4.0. Code: MIT.
