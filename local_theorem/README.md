# local_theorem — an explicit ball of strict minimality of the regular pentagon (version 1.1)

**Status: computer-assisted proof candidate, not peer reviewed.**
Author: Puyi Sun (ORCID [0009-0004-1644-5300](https://orcid.org/0009-0004-1644-5300)). Roles of AI systems: see below.
This directory is new in version 1.1; the material of version 1.0 (the PAQC cell certificates) is unchanged.

## Problem
For a planar domain Ω let F(Ω) = λ₁(Ω)·|Ω|, λ₁ the first Dirichlet eigenvalue of the Laplacian. The Pólya–Szegő
conjecture (1951) says that among all pentagons the regular one minimises F.
A preprint by H. Wang, "The Pólya–Szegő inequalities for polygons" (arXiv:2611.22305, dated 17 November 2026), claims a
proof that for every n ≥ 5 the regular n-gon minimises λ₁ among simple polygons with at most n sides and the same
area, which contains the pentagon case. It has not been peer reviewed, and we have not checked it; we learned of it
on 23 November 2026. Theorem L below does not depend on it.
Known before this work: B. Bogosel and D. Bucur reduced the case of a fixed number of sides to a finite computation
(J. Éc. polytech. Math. 11 (2024) 19–105, doi:10.5802/jep.250) and proved with validated (interval) computing that the
regular pentagon is a local minimiser (arXiv:2406.11575; we did not find an explicit radius there). A public repository
maintained under the GitHub account xfliu (`github.com/xfliu/polya-szego-conjecture`, read at commit `bf7c9c2`,
September 2026) gives a computer-assisted proof of strict local minimality along every transversal direction for
parameters 0 < t ≤ 1/6250 in its own parametrisation, using the one-vector Lehmann–Goerisch bound with Raviart–Thomas
fluxes; that is an explicit neighbourhood in other coordinates.

## Result
Chart: z ∈ ℝ⁶ ↦ P(z) = the pentagon with vertices (0,0), (1,0), (z₁,z₂), (z₃,z₄), (z₅,z₆); c = the point of the
regular pentagon.

> **Theorem L.** For every z with 0 < |z − c|₂ < R = 52429/2²⁰ (> 1/20), P(z) is a simple pentagon and
> F(P(z)) > F(P(c)). Moreover 18.9191028942 < F(P(c)) < 18.9191047951.

This does **not** prove the conjecture: nothing is claimed outside the ball.
The radius 52429/2²⁰ refers to the Euclidean norm of the fixed-first-edge chart above; the explicit neighbourhood of
the xfliu repository is expressed in different coordinates, and we do not compare the maximal neighbourhoods of the two
methods.
Our contribution: the explicit radius in this chart; the particular residual identity (Proposition R) for a trial family
that contains the *true* eigenfunction of the regular pentagon, known only up to an error with a proved H¹ bound — the
first-order error terms are paired with residual fields, and the quadratic error term is bounded separately; and a
rational coefficient-enclosure certificate whose final positivity proof can be checked independently. We make no
priority claim.

## Files
| path | content |
|---|---|
| `paper/paper.pdf`, `paper/paper.tex` | the note (8 pages) |
| `proof/LOCAL_LEMMA.md` | the full proof with all formulas as implemented; `proof/LOCAL_CERT.md`: the certificate format and what to check; `proof/LG_LEMMA.md`: the lower-bound lemma and the fan flux |
| `certificate/` | the certificate of the run of 17 November 2026: `witness.json.xz` (enclosures Tm ± Em of the coefficients of the polynomial Q), `leaves.jsonl.xz` (7,939 boxes with multipliers), `inner.json`, and the data behind the enclosures (`inputs.json`, `gram.json.xz`, `constants.json`, `lowdeg.json`, `rho_c.json`; `parts.json.xz`: the four parts of every error radius, diagnostic only); `local.json` (summary with the SHA-256 of each file in the encoding named by its entry: `gram.json.xz` and
`parts.json.xz` as compressed bytes, the others uncompressed), `run.log`, `cert_check_standalone.log` |
| `checker/lglocal_cert_check.py` | standalone checker of the final positivity proof ("level A"); Python standard library only, shares no code with the producer |
| `producer/` | the program that computes everything (`lglocal.py`, with `lgfan.py`, `lgbox.py`, `lgpoly.py`, `near_lg.py`, `lgcert.py`) and its tests |
| `second_party/` | report and programs of the independent review; `replay_status.json` says exactly what has been checked independently and what has not |
| `MANIFEST-SHA256.txt` | SHA-256 of every file of this directory except itself (the root manifest of the repository lists this file too) |

## What is proved by which program
1. *Level A* (`checker/`): **if** the coefficients of the polynomial Q(d) = (Π det_k)·Ψ(c + d) lie in the stored
   intervals Tm ± Em, **then** Q > 0 on the punctured ball. Exact rational arithmetic, no floating point in any decision.
2. *Level B* (`producer/lglocal.py`): the coefficients do lie in these intervals (Proposition R of the proof, exact
   polynomial arithmetic with FLINT, coefficientwise error radii), and the side conditions.
   Recomputed by the reviewer with independently written programs (`second_party/`, round PENT-21c): its own
   enclosure of all 13,024 coefficients (of degree 2 to 13; 12,968 of our centres are nonzero) intersects ours
   everywhere, and level A was re-proved with its enclosure (route (ii)); inclusion in our intervals holds for 7,086
   coefficients only, so the inclusion route (i) does not pass.
3. The analysis (Lemma LG, the fan flux, Lemma S on the vanishing gradient, Proposition R) is in the paper and in
   `proof/LOCAL_LEMMA.md`.

## How to verify
Tested with Python 3.13.17 and 3.14, python-flint 0.9.0, NumPy 2.5.3, SciPy 1.18.1 (Linux x86_64, WSL2 Ubuntu).
Do not run Python with `-O` (checks are `assert` statements).
```
sha256sum -c MANIFEST-SHA256.txt

# 1. the final positivity proof from the stored certificate; standard library only
#    (936 s = about 16 minutes on the machine of the run; timings vary)
python3 checker/lglocal_cert_check.py certificate
#    (add  --leaves 300  for a random sample of the boxes: about 2 minutes)

# 2. everything from scratch (needs network access for pip): the producer run took 383 s (about 6 minutes) on one
#    core of the machine of the run; checking its output (last line) takes another ~16 minutes
pip install python-flint==0.9.0 numpy scipy
cd producer
python3 tests/lglocal_check.py            # tests of the identities and bounds; last line PASS
python3 lglocal.py run ../out 16          # last line: LOCAL THEOREM: PROVED (computation complete)
python3 ../checker/lglocal_cert_check.py ../out/cert
cd ..
```
The trial functions come from floating-point linear algebra, so a recomputation may differ from the stored
certificate in the last digits of the lists; the proof does not depend on them. On the two machines used so far the
list of boxes was identical.

## Roles of AI systems
Claude (Anthropic) found the method, wrote the proof, the producer code and the level-A checker, and ran the
computation. GPT (OpenAI) reviewed the mathematics and wrote independent programs for the second-party checks
(`second_party/`); in the review rounds PENT-21, 21b and 21c it did not run or import any producer program, and in
the release review PENT-25 it also ran the published checker as a reader would. The author initiated and coordinated the work, ran the computation on his own computer, and decides
on the release. The author did not establish the mathematical correctness by hand; responsibility for the publication
lies with the author.

## Licence
Text and data: CC BY 4.0. Code: MIT. (See the licence files in the repository root.)
