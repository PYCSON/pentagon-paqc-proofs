# near_theorem — the annulus around the regular pentagon (part of version 2)

Status of the statement: computer-assisted proof candidate, not peer reviewed; re-certified in full by an independently
written program of a second-party AI reviewer and signed by it (PENT-26, 2026-11-28) (see "Second party" below).
Author: Puyi Sun (ORCID [0009-0004-1644-5300](https://orcid.org/0009-0004-1644-5300)). Roles of AI systems: see below.

## Result
Chart: z ∈ ℝ⁶ ↦ P(z) = the polygon with vertices (0,0), (1,0), (z₁,z₂), (z₃,z₄), (z₅,z₆). C = the binary64 chart point
of the regular pentagon used by all programs (`near_lg.CF`; |C − c|_∞ ≤ 2⁻⁵⁰ for the exact point c; C₃ = 1/2,
C₅ = 1 − C₁, C₆ = C₂ exactly). F(Ω) = λ₁(Ω)·|Ω|, F\* = 18.919134195733914 (exact terminating decimal; the common comparison threshold, not the value F(P_reg) ≈ 18.919104).

> **Theorem N.** For every z ∈ ℝ⁶ with 1/20 ≤ |z − C|₂ < 3/5, P(z) is a simple pentagon, star-shaped with respect to
> an interior point, and F(P(z)) > F\*.

F\* is larger than F of the regular pentagon (Theorem L of `local_theorem/`: F(P(c)) < 18.9191047951), and Theorem L
covers the ball |z − c| < 52429/2²⁰ ⊃ {|z − C| < 1/20}. Together: **F(P(z)) > F(P(c)) for all z ≠ c with |z − C| < 3/5.**
By itself this does not prove the Pólya–Szegő conjecture for pentagons; the pentagons outside this ball are the
subject of Theorem F.

## Files
| path | content |
|---|---|
| `proof/NEAR_THEOREM.md` | statement, the cover and its labels, the proof from a complete cover, file format |
| `proof/LG_LEMMA.md` | the lower bound (one-vector Lehmann–Goerisch), the fan flux, the certificate for a box of polygons |
| `certificate/blocks.tgz` | the cover: 1,344 block files with 51,434 leaves (38,280 certified boxes, 13,118 outside the annulus, 36 inside the inner ball); a leaf is only a box and a label, a checker computes its own trial functions |
| `certificate/run.log`, `verify_a.txt`, `verify_b.txt`, `reproof_full_20261118.txt` | logs of the run of 16 November 2026 and of the producer-side checks |
| `producer/` | `near_lg.py` (driver), `near_lg_verify.py` (check), `lgbox.py`, `lgpoly.py`, `lgfan.py`, `lgcert.py` (used by the tests), and two tests |
| `second_party/replay_status.json` | what has been checked independently, by whom, on which part of the data |
| `MANIFEST-SHA256.txt` | SHA-256 of every file of this directory |

## What is proved by which program
1. `near_lg.py` cuts a cube of half-width 5/8 around C into boxes with dyadic ends. Each leaf has one label:
   `mirror` (the box lies in z₃ ≤ 1/2; the mirror image of P(z) has the same F and the same distance to C), `far`
   (the box does not meet the ball |z − C| < 3/5), `inner` (the box lies in |z − C| < 1/20), or `LG`.
   The first three tests are exact rational comparisons.
2. `LG`: for all z of the box, the five fan determinants are positive and two polynomials P_D, P_C(F\*) in z with
   rational coefficients are positive (Lemma on positivity on a box: exact Lagrangian bound of the quadratic part, sum
   of absolute values of the higher terms). By Lemma LG this gives the two claims of the theorem on the box.
3. `near_lg_verify.py` recomputes the list of blocks, checks that the leaves of every block tile it (exact), repeats
   the three exact tests and re-proves LG leaves from scratch with new trial functions. It uses the modules of the
   producer: it is a re-run, not an independent check.

## How to verify
Tested with Python 3.13.17 and 3.14, python-flint 0.9.0, NumPy 2.5.3, SciPy 1.18.1 (Linux x86_64, WSL2 Ubuntu).
Do not run Python with `-O`.
```
sha256sum -c MANIFEST-SHA256.txt
pip install python-flint==0.9.0 numpy scipy
mkdir run && tar xzf certificate/blocks.tgz -C run            # run/blocks/<id>.jsonl
cd producer
python3 tests/lgpoly_check.py                                 # test of the box certificate against pointwise values
python3 near_lg_verify.py ../run --sample 1                   # tiling, skip labels, one LG leaf per block: 27 minutes on one core
python3 near_lg_verify.py ../run                              # every LG leaf re-proved: about as long as the run
python3 near_lg.py ../new run --jobs 8                        # the whole cover from scratch (7-10 CPU hours in our run)
python3 near_lg.py ../new status
```
The leaf lists of two runs need not be identical (the choice of boxes depends on floating-point heuristics); each run
is a proof of the same statement if it ends without unresolved leaves and passes the check.

## Runs so far
- 16 November 2026, sandbox (two processes): complete, no unresolved leaf; check with one leaf per block: VERIFY OK;
  re-proof of every leaf of the 1,201 smallest blocks (18,650 of the 38,280 certified boxes): all passed.
- 17 November 2026, the author's computer (8 processes): complete, no unresolved leaf; check with three leaves per
  block: VERIFY OK.
- 18 November 2026: replay of this directory from a copy (manifest, unpack, test, check with one leaf per block):
  `blocks checked 1344 of 1344 … LG leaves re-proved 1344`, VERIFY OK.

## Second party
The lemma and the formulas of the box certificate were reviewed and re-implemented by a second party (review
PENT-21); the simplicity lemma, the mirror step, the skip tests and the cover argument of `proof/NEAR_THEOREM.md`
were reviewed in PENT-23. The second party's independently written exact program re-proved all 38,280 near targets on
the author's computer (22–26 November 2026): 38,267 in the main run, 1 in a first supplementary round with the trial
degree raised by two, and 12 in a second supplementary round in which the box is split exactly and subboxes lying
entirely outside the near region are excluded (86 LG and 109 outside subboxes; see the paper, Section 15). Every stored
witness was replayed, and the producer's acceptance program agrees with the replay shard by shard. Theorem N is
signed by the second party (PENT-26, 2026-11-28). `second_party/replay_status.json` and `../second_party/final_reports/` record the
details. This is a review by an AI system, not peer review.

## Roles of AI systems
Claude (Anthropic) found the method, wrote the proofs, the programs and the producer-side checks, and ran the first
computation. GPT (OpenAI) reviews the mathematics and writes independent programs. The author initiated and
coordinated the work, repeated the computation on his own computer and decides on the release; he did not establish
the mathematical correctness by hand. Responsibility for the publication lies with the author.

## License
Text: CC BY 4.0. Code: MIT.
