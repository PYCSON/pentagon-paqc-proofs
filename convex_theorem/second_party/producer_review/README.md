# The producer's review of the second party's independent checker (Claude, 22 November 2026)

The second party (GPT, an AI reviewer) wrote its own programs for the certificates of Theorems N and F from the written
specification, without running or importing any producer program (`independent_checker/`). Because the signatures on
Theorems N and F rest on these programs, the producer in turn reviewed them. The programs were only read, never run by
the producer; the checks below use separate code.

1. **Exact kernel of the box certificate** (`lg_exact.py`, 208 lines, standard library). Read line by line against
   `proof/LG_LEMMA.md` and the chart formulas of `proof/FAR_THEOREM.md`: the coefficient layout (u in H¹₀, chi
   continuous), the monomial moments, h = (2 + E)⁻¹u, the fourteen bilinear forms, the four chart maps and the cut
   j = 3, 0, the assembly of A₀, A₁, A₂ with all scalings, P_D and P_C, the positivity test (crude bound; quadratic
   dual bound with exact LDLᵀ, multipliers ≥ 0 and the absolute sum of the higher terms), the exact binomial split.
   No error found. Three hardening requests (reject witnesses that carry fixed vertices; bind every witness to the
   frozen inputs; guard the packed exponents) were implemented by the reviewer in its replay program.
2. **Pointwise cross-check** (`producer_review/lg_point_check.py`): a separate exact computation that never forms the
   box polynomials. A witness is read as data; at rational points of its box the polygon is built from the chart
   formulas, u and chi are rebuilt as polynomials in physical coordinates, the conditions u = 0 on the boundary,
   continuity, −div w = u and continuity of the normal flux are tested, the three integrals are computed exactly and
   Lemma LG is evaluated. 1,222 witnesses of the reviewer (all 582 near leaves that touch the inner sphere, the
   largest leaf of 13 blocks, 16 leaves with non-convex members, 11 far boxes of the three charts, 600 table entries),
   2,484 points: no failure; smallest relative margin over the target 1.8·10⁻⁵. Five negative controls behave as they
   should.
3. **Label and table criteria** (`label_exact.py` and the corresponding part of `proofs24b.py`): every criterion was
   compared with section 5 and 6 of `proof/FAR_THEOREM.md`, `proof/FIBRE_LEMMA.md` and the formulas listed in the
   producer's `qchart.py` and `cutchart.py` (empty boxes, hand-over E2/E3, T, the four untouched-ear cases and the four
   stub polynomials of X, the homogeneous face of E3, TX, the strips S(a) and S(b), R, the fibre bound O, the near
   label N, the quadrilateral bounds and the closure of the table). The reviewer's two new routes (a joint Bernstein
   inequality for the enclosing triangle; the fibre bound with the breakpoints of the hull) were re-derived. No error
   found. The reviewer's own ledger of statements and proofs (`LABEL_ROUTES.md`) agrees with this reading.
4. **Structure checkers** (`structure.py`, `near_robust.py`, the geometry part of the replay): exact dyadic tiling,
   completeness of the root cells of the three charts, the near blocks and skip labels with the uncertainty of the
   centre: read, no error found.
5. **Acceptance tool** (`producer_review/accept_pent24b_r3.py`): re-enumerates the targets from the frozen inputs,
   checks that every target has a witness tree whose leaves tile exactly its box and are bound to its chart, cut index
   and target value, spot-checks uniformly drawn leaves with the pointwise computation, and writes the list of SHA-256
   values of the result shards in the format of the reviewer's replay report, so that both checks can be seen to have
   read the same bytes.

Not read: the programs that only propose trial coefficients (`lg_search.py`, `lg_stdlib.py`); they do not enter the
acceptance of a certificate.

This is a review by the producer, an AI system, of programs written by another AI system. It is not peer review.
