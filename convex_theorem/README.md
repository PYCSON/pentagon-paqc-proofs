# convex_theorem — the Pólya–Szegő conjecture for convex pentagons (version 2)

**Version 2.0.** Status of the statement: computer-assisted proof candidate, not peer reviewed. All three parts were
re-certified by independently written programs of a second-party AI reviewer and signed by it (PENT-26, 2026-11-28)
(see "Status"). Author: Puyi Sun (ORCID [0009-0004-1644-5300](https://orcid.org/0009-0004-1644-5300)). Roles of AI systems: see below.

## Result
For a planar domain Ω let F(Ω) = λ₁(Ω)·|Ω| (first Dirichlet eigenvalue times area).

> **Theorem.** For every convex pentagon P (simple, of positive area), F(P) ≥ F(P_reg), with equality if and only if
> P is similar to the regular pentagon.
> Moreover 18.9191028942 < F(P_reg) < 18.9191047951.

Pólya and Szegő proved the analogous statements for triangles and quadrilaterals and conjectured them for all n ≥ 5.
Non-convex pentagons are **not** covered. Earlier work: Bogosel and Bucur reduced the case of fixed n to a finite
computation and proved local minimality of the regular pentagon and hexagon with validated computing; a public
repository maintained under the GitHub account xfliu (`github.com/xfliu/polya-szego-conjecture`) proves strict local
minimality of the regular pentagon along every direction for small parameters with the same kind of lower bound
(one-vector Lehmann–Goerisch); Cheng, Gui, Hu and Li proved the conjecture for convex polygons with sufficiently many
sides. Bogosel, Bucur and Fragalà proved the local version for the torsional rigidity (arXiv:2609.23051) and have
since posted a preprint claiming an analytic proof of the torsional (Saint-Venant) inequality for all n ≥ 3 (arXiv:2611.24573,
18 November 2026). A preprint of H. Wang (arXiv:2611.22305, 17 November 2026) claims a proof of the Pólya–Szegő
inequalities for every n ≥ 5 and all simple polygons with at most n sides, which contains the convex pentagons; it has
not been peer reviewed, and we have not checked it. This work does not depend on it: it is an independent
computer-assisted argument for the convex case, by a different method, whose certificates can be replayed box by box.

**Independence and dates.** This work was carried out independently of Wang's preprint. An earlier stage of it, on
the same problem and in the same chart (another certification method and certificates for a few boxes of pentagons),
was released on 1 November 2026 as version 1.0 of this repository. We learned of Wang's preprint on 23 November 2026,
after the certificates of the near and far covers had been computed; no argument of it is used here. We make no
priority claim.

## Structure of the proof
Chart: z ∈ ℝ⁶ ↦ P(z) = the polygon (0,0), (1,0), (z₁,z₂), (z₃,z₄), (z₅,z₆); c = the regular pentagon, C = its binary64
neighbour (|C − c|_∞ ≤ 2⁻⁵⁰); F\* = 18.919134195733914, an exact terminating decimal used as the common comparison threshold. It is not the value of the regular pentagon: F(P_reg) ≈ 18.919104 (Theorem L: 18.9191028942 < F(P_reg) < 18.9191047951).
| part | statement | where |
|---|---|---|
| Theorem L | F(P(z)) > F(P(c)) for 0 < \|z − c\| < 52429/2²⁰ | `../local_theorem/` (version 1.1) |
| Theorem N | F(P(z)) > F\* for 1/20 ≤ \|z − C\| < 3/5 | `near/` |
| Theorem F | every strictly convex pentagon has F > F\* or a labelling with \|z − C\| < 3/5 | `far/` |
| assembly | the three theorems give the Theorem | `paper/paper.pdf`, Section 1 |
One inequality carries all three parts: the one-vector Lehmann–Goerisch bound with ν = κ₀/|Ω| (Krahn–Szegő) for
piecewise polynomial trial functions and fluxes on a fan of triangles; the certificates are positivity statements for
polynomials with rational coefficients on boxes, checked in exact rational arithmetic. Theorem F also uses classical
bounds (Pólya–Szegő for triangles and quadrilaterals, strip and Makai bounds, a one-dimensional fibre bound) and a
covering argument in a chart of "cut quadrilaterals".

## Files (to be filled in when the data are final)
| path | content |
|---|---|
| `paper/` | the paper (PDF and TeX) |
| `proof/` | `LG_LEMMA.md`, `NEAR_THEOREM.md`, `FAR_THEOREM.md`, `FAR_PROOFS.md`, `FIBRE_LEMMA.md` |
| `near/` | Theorem N: the cover (1,344 block files, 51,434 leaves), producer, logs, its own README and manifest |
| `far/` | Theorem F: the cover (`pent_far5_cert.tar.xz`, 13.7 MB: roots, leaves and certified boxes of all 1,260 blocks, the quadrilateral table), the producing package, second implementations of the checks, logs, its own README and manifest |
| `second_party/` | the reviewer's reports, its independent checker, its ledger of the label criteria, the producer's review of that checker, and `replay_status.json`, which says exactly what has been checked independently and what has not |
| `MANIFEST-SHA256.txt` | SHA-256 of every file |

## How to verify
```
sha256sum -c MANIFEST-SHA256.txt
pip install python-flint==0.9.0 numpy scipy
# Theorem L: ../local_theorem/README.md (checker 16 minutes)
# Theorem N: near/README.md "How to verify" (one leaf per block: 27 minutes; every leaf: about as long as the run)
# Theorem F: far/README.md "How to verify" (root ranges 2 minutes, tiling 5 minutes, every label about 14 CPU hours, producer self-check)
# the second party's checkers, with their own commands
```

## Status
(27 November 2026; `second_party/replay_status.json` is the precise record.)
- Theorem L: signed by the second party in full (analysis, both levels of the certificate, the release tree); released
  separately as version 1.1.
- Theorem N: computed twice (two computers); producer-side check passed, and all 38,280 LG leaves re-proved from scratch
  with new trial families by the producer's own programs. Second party: simplicity lemma, mirror, skip tests and cover
  argument passed; structure of all blocks re-checked with its own programs; its independently written exact program
  re-proved all 38,280 near targets on the author's computer (38,267 in the main run, 1 in a first supplementary round,
  12 in a second supplementary round with exactly excluded outside subboxes: 86 LG and 109 outside subboxes), and every
  stored witness was replayed. **Signed (PENT-26, 2026-11-28).**
- Theorem F: complete run on the author's computer (1,260 blocks, 1,724,444 boxes, 116 CPU hours, no unresolved leaf);
  producer-side self-check passed; root ranges, tiling and all 1,346,929 closed-form labels re-proved by a second
  implementation of the producer. Second party: chart, kept representation, corner principle, ear comparisons and
  cover argument passed; its independently written exact program re-proved on the author's computer all 377,515
  certified boxes, all 185,100 table entries and all 1,346,929 labels, and every stored witness was replayed.
  **Signed (PENT-26, 2026-11-28).**
- The full independent run covered 1,947,824 targets in 15,239 shard files. The reviewer's replay report and the
  producer's acceptance report, each computed from the stored results on the author's computer, agree shard by shard;
  both are in `second_party/final_reports/`.
- Main theorem (assembly of L, N and F): **signed (PENT-26, 2026-11-28).**
- All of this is computer-assisted and was carried out by AI systems; the mathematics has not been reviewed by human
  experts, and the work has not been peer reviewed.

## Roles of AI systems
Claude (Anthropic) found the method, wrote the proofs, the programs and the producer-side checks, and ran the first
computations. GPT (OpenAI) reviewed the mathematics and wrote independent programs; Claude in turn read those programs. The author initiated and
coordinated the work, ran the computations on his own computer and decides on the release; he did not establish the
mathematical correctness by hand. Responsibility for the publication lies with the author.

## License
Text: CC BY 4.0. Code: MIT.
