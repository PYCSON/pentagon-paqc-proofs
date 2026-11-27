# Theorem N (near region): statement, proof structure, files (Claude, 2026-11-17)

**Status: proof candidate.** The certificate used (Lemma LG, computation C, exact polynomial positivity) was reviewed
by a second party (GPT, PENT-21: formulas passed); the cover argument of this note and an independent checker of the
leaves have not been. Programs: `v2_run/near_lg.py` (driver), `near_lg_verify.py` (producer-side check), `lgbox.py`,
`lgpoly.py`, `lgfan.py`. History: `v2_run/NEAR_SPEC.md` section 19.

## 0. Statement
Chart: x ∈ ℝ⁶ ↦ P(x) = the polygon (0,0), (1,0), (x₀,x₁), (x₂,x₃), (x₄,x₅). C = the binary64 point `near_lg.CF`
(= `plan_cells.C`), |C − c|_∞ ≤ 2⁻⁵⁰ for the chart point c of the regular pentagon; C₂ = 1/2 and C₄ = 1 − C₀, C₅ = C₁
exactly. F*_up = 18.919134195733914 (the exact terminating decimal; note of 2026-11-23 after PENT-26a M1: producer certificates that used the slightly larger binary64 number imply the same inequality).

> **Theorem N.** For every x ∈ ℝ⁶ with 1/20 ≤ |x − C|₂ < 3/5, the polygon P(x) is a simple pentagon, star-shaped with
> respect to an interior point, and F(P(x)) > F*_up.

Together with Theorem L (`LOCAL_LEMMA.md`: F(P(x)) > F(P(c)) for 0 < |x − c|₂ < R = 52429/2²⁰; and
|x − C| < 1/20 implies |x − c| < 1/20 + √6·2⁻⁵⁰ < R; F(P(c)) < 18.9191047951 < F*_up) this gives: **for every x with
|x − C| < 3/5 and x ≠ c, F(P(x)) > F(P(c)).** Theorem F (`FAR_THEOREM.md`) reduces every strictly convex pentagon to
this ball or proves F > F*_up directly.

## 1. Inputs
Lemma LG with ν = κ₀/|Ω| (Krahn–Szegő), computation C of `LG_LEMMA.md` (one trial family for a box of polygons whose
vertices are affine in the parameters; here the vertices are the parameters); exact positivity of polynomials on boxes.
Reflection invariance of F.

## 2. The cover
- Mirror: m(x) = (1 − x₄, x₅, 1 − x₂, x₃, 1 − x₀, x₁). P(m(x)) is the mirror image of P(x) under (X, Y) ↦ (1 − X, Y)
  with the vertices relabelled, so F(P(m(x))) = F(P(x)); m is an isometry of ℝ⁶ with m(C) = C exactly. Hence it
  suffices to prove the statement for x₂ ≥ 1/2.
- Root cube: centre ROOT_C = C rounded to multiples of 2⁻¹⁰, half-width 5/8; it contains the ball |x − C| < 3/5.
- The root cube is bisected; a box [lo, hi] (binary64 ends, exact dyadic rationals) becomes a leaf with the label
  | label | test (exact rational arithmetic) | meaning |
  |---|---|---|
  | mirror | hi₂ ≤ 1/2 | every point is the mirror image of a point with x₂ ≥ 1/2 at the same distance from C |
  | far | inf_box |x − C|² ≥ 9/25 | no point of the annulus in the box |
  | inner | sup_box |x − C|² < 1/400 | inside the ball of Theorem L; nothing claimed here |
  | LG | box certificate | for all x of the box the five fan determinants are positive (P(x) is a simple polygon, star-shaped about O) and F(P(x)) > F*_up |
  The first three are tested in this order before any certificate; inf and sup of |x − C|² over a box are exact
  (coordinatewise).
- Blocks: the 4⁶ cubes of side 5/16 of the root grid; a block (or a larger cube of the grid) that is closed by one of
  the three skip tests is not stored; the others are stored with all their leaves: 1,344 blocks.
**Lemma S (simplicity; 2026-11-19).** For every x with |x − C| < 3/5, P(x) is a simple closed polygon: its five
vertices are distinct, two non-adjacent edges are disjoint, and two adjacent edges meet only in their common vertex.
*Proof.* In the regular pentagon P(c) (side 1, angles 108° > 90°) the distance between two non-adjacent edges and the
distance from a vertex to an edge not containing it are both at least 1 (equal to 1 for the nearest pairs, whose closest
points are two vertices joined by an edge; e.g. a vertex and the opposite edge are at distance cos 18° + sin 36° > 1;
correction of the wording after review PENT-23), and two vertices are at distance ≥ 1. Moving the endpoints of a segment by at most δ moves the segment by at
most δ in the Hausdorff distance, and the distance between two sets is 1-Lipschitz in each set for the Hausdorff
distance. In each of the three statements the two objects involve disjoint sets of vertices, say with largest
displacements a and b; V₀, V₁ do not move and a² + b² ≤ |x − c|² < (3/5 + √6·2⁻⁵⁰)², so a + b < √2·0.6001 < 0.85 < 1.
Hence the distances stay positive. Two adjacent edges V_{k−1}V_k, V_kV_{k+1} that met in a further point would be
collinear and overlapping, so V_{k+1} would lie on V_{k−1}V_k or V_{k−1} on V_kV_{k+1}: excluded by the vertex–edge
case. ∎
*Why it is needed.* The box certificate proves that the five fan triangles about O are positively oriented on the whole
box; that makes P(x) star-shaped with respect to O only when the winding number about O is 1, and it is constant on the
box (`LG_LEMMA.md`, computation C). A pentagram has positive fan determinants and winding number 2. Every LG leaf
contains a point x* with |x* − C| < 3/5 (it is not `far`); there P(x*) is simple by Lemma S, so the winding number is 1
on the whole box. (The data are not affected; the argument was missing in the note of 2026-11-17.)

**Proof of Theorem N from a complete run without unresolved leaves.** Let 1/20 ≤ |x − C| < 3/5; replacing x by m(x)
if necessary, x₂ ≥ 1/2. The second coordinate of the root cube is [1/2 − 5/8, 1/2 + 5/8] and the grid of side 5/16 has
the plane x₂ = 1/2 as a cell boundary, so x lies in a closed grid cube B with lo₂ ≥ 1/2. Neither B nor a larger grid
cube containing it is closed by a skip test: not `mirror` (it contains points with x₂ > 1/2), not `far` (it contains
x with |x − C|² < 9/25), not `inner` (|x − C|² ≥ 1/400). So B is a stored block and x lies in one of its leaves. That
leaf has lo₂ ≥ 1/2 and positive width, so it is not `mirror`; it contains x, so it is not `far` or `inner`. Hence it
is `LG`, and the certificate with Lemma S (winding number 1 on the leaf) gives both claims at x. ∎

## 3. Files of a run
`OUT/blocks/<id>.jsonl`: first line `{"block": id, "lo", "hi", "stats"}`; then one line per leaf
`[lo, hi, label, N]` (label ∈ mirror, far, inner, LG; N = the degree of the producer's family for LG leaves, a checker
may use any). The id is the path in the grid (six base-4 digits). `near_lg.blocks()` lists the stored blocks.
Labels `open` (time limit) and `U` (minimal size reached) mean "not proved"; a complete run has none.

## 4. What the producer-side check does (`near_lg_verify.py`)
The stored blocks are exactly those not closed by a skip test (recomputed); the leaves of each block tile it (exact);
the skip labels are re-tested; LG leaves are re-proved from scratch (a sample per block, or all: a full re-proof of
the sandbox run of 2026-11-16, smallest blocks first, covered 1,201 of the 1,344 blocks with 18,650 of the 38,280 LG
leaves, all passing, before it was stopped); the root cube contains the ball; m(C) = C.
It uses the production modules; an independent checker has to implement the box certificate itself.

## 5. What a reviewer should attack
0. Lemma S and its use: the winding number of the fan on an LG leaf (2026-11-19).
1. The mirror step and the role of the exact symmetry of C; the boundary case x₂ = 1/2.
2. That "far" uses the infimum and "inner" the supremum of |x − C|², and the order of the tests.
3. The skipped cubes of the grid (not stored): re-enumerate them.
4. The box certificate on a sample of leaves, with an own implementation; in particular leaves close to the inner
   sphere (smallest margins: F − F*_up of order 10⁻⁴ relative) and leaves that contain non-convex polygons (allowed:
   the certificate needs only positive fan determinants).
