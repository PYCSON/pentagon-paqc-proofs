# Theorem F (far region): statement, proof structure, and what a complete run proves (Claude, 2026-11-17)

**Status: proof candidate, unreviewed as a whole.** This note is the clean statement of the far cover in its T1-free
form ("far5": `PENT_LGBOX=1 PENT_QLG=1`). `FAR_SPEC.md` is the history of how it was found; where the two differ, this
note is the intended statement. Reviewed so far by a second party: Lemma LG and the fan certificate (B), the formulas of
the box certificate (C) (GPT, PENT-21). Not yet reviewed: everything specific to the cut chart (sections 2–7).
**Proofs of every criterion and derivations of every closed form: `FAR_PROOFS.md`** (2026-11-18; its identities are
checked symbolically by `v2_run/tests/far_formula_check.py`).
Programs (all in `v2_run/`): `cutchart.py`, `qchart.py`, `cut_run.py` (labels), `qtable.py`, `fibre.py`, `cutcell.py`,
`lgbox.py`, `lgpoly.py`, `far_block.py`, `far_runner.py`, `far_verify.py`.

## 0. Statement
Pentagon chart: x ∈ ℝ⁶ ↦ the pentagon (0,0), (1,0), (x₀,x₁), (x₂,x₃), (x₄,x₅). C = the binary64 point `plan_cells.C`
(= `near_lg.CF`; |C − c|_∞ ≤ 2⁻⁵⁰, c the regular pentagon; C is exactly invariant under the mirror of the chart).
For a convex pentagon P, d(P) = min over its ten labellings h (five rotations, with and without reflection) of
|x_h(P) − C|₂, x_h(P) the chart point of P in the labelling h. F*_up = 18.919134195733914 (the exact terminating decimal; note of 2026-11-23 after PENT-26a M1: producer certificates that used the slightly larger binary64 number imply the same inequality) (> F(regular pentagon), by
Theorem L: F₀ < 18.9191047951).

> **Theorem F.** Every strictly convex pentagon P with d(P) ≥ 3/5 has F(P) = λ₁(P)·|P| > F*_up.

What is actually proved is the disjunction: **every strictly convex pentagon P has F(P) > F*_up, or has a labelling h
with |x_h(P) − C|₂ < 3/5.** The second case is the domain of the near cover (`near_lg.py`: all chart points with
1/20 ≤ |x − C| < 3/5) and of Theorem L (|x − c| < 52429/2²⁰).

## 1. Inputs from analysis
1. Domain monotonicity of λ₁; similarity and reflection invariance of F.
2. Pólya: among triangles the equilateral one minimises F (F ≥ 4π²/√3); among quadrilaterals the square (F ≥ 2π²).
3. Strip bound: a domain inside a strip of width W has λ₁ ≥ π²/W².
4. Lemma OD / PS (fibre bound, `FIBRE_LEMMA.md`).
5. Lemma LG with ν = κ₀/|Ω| ≤ λ₂ (Krahn–Szegő), in the forms B (one polygon) and C (a box of polygons) of `LG_LEMMA.md`.
6. Makai's inequality for convex domains, λ₁ ≥ π²·Per²/(16|Ω|²) (used only for quadrilaterals, in the std chart).

## 2. The cut chart
Let P be a strictly convex pentagon, counter-clockwise, e_i = V_{i+1} − V_i, and for an edge k
A_k = e_{k−1} × e_k > 0, B_k = e_k × e_{k+1} > 0, C_k = e_{k−1} × e_{k+1}.
- If C_k > 0, the lines of the edges k − 1 and k + 1 meet beyond the edge k in a point W, the triangle (V_k, W, V_{k+1})
  has area ear_k = A_k B_k/(2C_k), and Q_k = P ∪ ear is a strictly convex quadrilateral. Some edge has C_k > 0 (the
  exterior angles sum to 2π, so two consecutive ones sum to less than π).
- **(F1)** P ⊂ Q_k, so F(P) ≥ F(Q_k)·|P|/|Q_k|.
- **Kept representation.** Call k *kept* for P if C_k > 0 and ear_k ≤ ear_k′ for every k′ with C_k′ > 0. Every P has a
  kept edge. It suffices to treat every P in the representation of one kept edge.
- **Coordinates.** Normalise Q = Q_k by a similarity, and a reflection if necessary, to Q(q) = (q₀ = 0, q₁ = (1,0), q₂,
  q₃), q = (q₂ₓ, q₂ᵧ, q₃ₓ, q₃ᵧ), with the base q₀q₁ a longest edge of Q and W = q_j, j ∈ {3, 0} (W = q₂ or q₁ are the
  mirror images of these two cases). With t, u ∈ (0, 1) the fractions of the two edges at W that are cut away,
  j = 3: P = [q₀, q₁, q₂, q₃ + t(q₂ − q₃), (1 − u)q₃];  j = 0: P = [(t, 0), q₁, q₂, q₃, u·q₃];
  |P| = |Q|(1 − t·u·τ_j(q)), τ₃ = D/(q₂ᵧ + D), τ₀ = q₃ᵧ/(q₂ᵧ + D).
  **Admissible** points: q with Q(q) strictly convex and |q₂ − q₁|, |q₃ − q₂|, |q₃| ≤ 1; then
  q ∈ 𝔅 := [0,2] × [0,1] × [−1,1] × [0,1]. So: every strictly convex pentagon is, up to similarity and reflection,
  P(q, t, u; j) for an admissible q, t, u ∈ (0,1), j ∈ {3, 0}, with the cut edge kept.
- **Invariants of Q** (cross products of its edge vectors e₀ = (1,0), e₁ = q₂ − q₁, e₂ = q₃ − q₂, e₃ = −q₃):
  q₃ᵧ = e₃×e₀, q₂ᵧ = e₀×e₁, K = e₁×e₂, D = e₂×e₃, M = e₃×e₁ = D − q₃ᵧ, L = e₀×e₂ = q₃ᵧ − q₂ᵧ, 2|Q| = q₂ᵧ + D.
  Q strictly convex ⟺ q₂ᵧ, K, D, q₃ᵧ > 0.

## 3. Three charts for q, and the corner principle
std: w = q.  E2: w = (q₃ₓ, q₃ᵧ, a, h), q₂ = q₁ + a(q₃ − q₁) + h(0,1) (then K = h(1 − q₃ₓ)).
E3: w = (q₂ₓ, q₂ᵧ, b, h), q₃ = b·q₂ + h(0,1) (then D = h·q₂ₓ). A box is a product of intervals in z = (w, t, u).
**Corner principle.** In each chart every vertex of P (and of Q) is multi-affine in z, and each coordinate moves all
vertices along one common direction (depending on the other coordinates). Hence every cross product of two vertex
differences is multi-affine in z (its second derivative in z_m is 2(∂_m u)×(∂_m v) = 0), and so are the area and every
linear functional of the vertices. A multi-affine function on a box takes its extreme values at the corners. All box
criteria below are of this kind: exact corner values (the corner vertices are exact dyadic rationals, asserted), then
monotone combinations. No interval enclosure of a 6-dimensional map is used anywhere.

## 4. Domains
- std roots: all cells of side 1/16 of 𝔅 except those that contain no strictly convex Q (some invariant q₂ᵧ, q₃ᵧ, K, D
  is ≤ 0 at all 16 corners) or no admissible q (an edge 1, 2, 3 is longer than 1 on the whole cell). 57,093 roots.
- Hand-over: a q-box of the std chart is labelled `E2` if q₃ᵧ ≥ 1/4 on it and K < H₂(1 − q₃ₓ) at its 16 corners
  (H₂ = 1/8), `E3` if q₂ₓ ≥ 1/4 and D < H₃·q₂ₓ (H₃ = 1/8). For an admissible strictly convex q of such a box: in E2,
  1 − q₃ₓ > 0, h = K/(1 − q₃ₓ) ∈ (0, H₂), and a(q₃ − q₁) = e₁ − h(0,1) with |q₃ − q₁| ≥ q₃ᵧ ≥ 1/4 gives |a| ≤ 4.5; in
  E3, h = D/q₂ₓ ∈ (0, H₃) and b = q₃ₓ/q₂ₓ, |b| ≤ 4.
- E2 roots: cells of side 1/16 of [−1,1] × [1/4,1] × [−4.5,4.5] × [0,1/8] that may contain a strictly convex Q
  (15,646); E3 roots: of [1/4,2] × [0,1] × [−4,4] × [0,1/8] (25,660).
- For every root c and j ∈ {3, 0} the *column* c × [0,1]² is cut by bisection into leaves; every leaf carries a label.

## 5. Labels: what each one asserts about the box, and why
Every assertion is about the members z = (w, t, u) of the box with Q(w) strictly convex and t, u ∈ (0,1).
LB denotes a number with F(Q(w)) ≥ LB for every strictly convex Q(w) of the w-box (section 6).
| label | assertion | criterion |
|---|---|---|
| V | no strictly convex Q in the w-box | an invariant q₂ᵧ, q₃ᵧ, K or D is ≤ 0 at the 16 corners |
| E2, E3 | (std) covered in the edge chart | section 4 |
| T | F(P) > F*_up | the certified quadrilateral bound LB > F*_up, sup τ_j > 0, and sup(t)·sup(u) < θ_T := (1 − F*_up/LB)/sup τ_j; then by (F1) F(P) ≥ LB(1 − tuτ_j) > F*_up |
| X | the cut edge is not kept | another edge k′ has C_k′ > 0 and ear_k′ < ear_cut on the whole box: (i) k′ untouched by the cut: t·u > g(w), e.g. j = 3, k′ = 0: M > 0, g = q₂ᵧq₃ᵧ/(DM); the four cases are listed in `qchart.py`; evaluated with inf(t)·inf(u) and the corner ranges of the invariants; (ii) k′ one of the two shortened edges: two multi-affine polynomials (`cutchart.STUB`) have a sign at all 64 corners; (iii) chart E3 at the face q₂ᵧ = 0 = h: the comparisons (i) are homogeneous in (q₂ᵧ, h) and are verified for all ratios h/q₂ᵧ ∈ (0, ∞) (`qchart.face_not_kept`) |
| TX | every member: F > F*_up or not kept | whole column, strict separation θ_T > sup g, with the positivity assumptions of T and of X: t·u < θ_T is T, t·u ≥ θ_T > g is X(i) |
| S | F(P) > F*_up | strip bounds: (a) for a non-zero vector n, W_n = max_{i,k}|⟨V_i − V_k, n⟩|: choose a rational τ ≥ 0 with |P| − τ|n|⁻²⟨V_i − V_k, n⟩² ≥ 0 at every corner and for every vertex pair (no division by a corner width, which may be 0 at a degenerate corner); this function is concave in each coordinate, so it is ≥ 0 on the box; a strictly convex member has W_n > 0, lies in a strip of width W_n/|n| and has F ≥ π²τ; the label requires π²τ > F*_up. (b) for an edge of length l > 0 with max_k c_k > 0: F ≥ π²l²/(2 max_k c_k); the label requires the verified lower bound > F*_up. (Text revised 2026-11-24 after GPT PENT-26a M5 and PENT-26b; the computation is unchanged.) |
| R | F(P) > F*_up | two non-adjacent edges m, m+2 with C_m > 0 and C_{m+2} > 0: P lies in the triangle T of the other three edge lines, |T| = |P| + ear_m + ear_{m+2}, F ≥ (4π²/√3)|P|/|T|; the label requires this verified lower bound > F*_up |
| O | F(P) > F*_up | fibre bound of the convex hull H of the 320 corner vertices (P(z) ⊂ H), times min corner area |
| N | d(P) < 3/5 | for a labelling h: Σ_p (|r_p − C_p| + ρ_p)² < 9/25, where u_p = (W_p − W₀)/(W₁ − W₀) is the chart coordinate, r_p its value at the centre, ρ_p = max_corners |N_p − r_pD| / min_corners Re(D ē) (N_p, D complex multi-affine) |
| H | — | **not used** in the far5 run (`PENT_NO_H=1`). In older runs: closed-form tests of `far_exclude.classify(rep=False)` on the bounding box in the pentagon chart; it did not occur in about 60,000 sampled leaves |
| C | F(P) > F*_up | LG box certificate (computation C of `LG_LEMMA.md`): the fan determinants, PD and PC are positive on the box, with one trial family; the fan determinants are positive on the whole box (part of the certificate); proves it for every member that is a simple polygon, in particular for every strictly convex member (the fan has winding number 1 there; `LG_LEMMA.md`, computation C) |

## 6. The quadrilateral bound LB of a w-box
LB = max of: 2π² (Pólya); the triangle bound u3, (4π²/√3)/(1 + ear_k/|Q|) for an edge k of Q with C_k > 0 on the box
(the four cases in `qchart.py`); the strip bounds π²|Q|/max(q₂ᵧ, q₃ᵧ)² and π²/(2 max(q₂ᵧ, q₃ᵧ)); the fibre bound of
the hull of the corner quadrilaterals; in the std chart also π²·Per²/(16|Q|) (Makai; perimeter from below, area from above,
by interval arithmetic on the box); and, for the root cell containing the w-box (std chart), the table value:
**LG table** (`qtable.py`, `PENT_QLG=1`): for the root cell, or for each of its 16 sub-cells that may contain a strictly
convex Q, either a closed-form bound or a number L with "F(Q) > L for every strictly convex Q of the (sub-)cell" proved by the box certificate with target L (the fan determinants positive on the whole (sub-)cell is part
of the certificate); the table value is the minimum over the parts.

## 7. Proof of Theorem F from a complete run
Let P be strictly convex with kept edge k, and (q, t, u; j) its coordinates (section 2), q admissible.
1. q lies in a cell of 𝔅; that cell is a std root (it contains the admissible, strictly convex q). So (q, t, u) lies in
   a leaf of the column of that root for the index j.
2. If the leaf is `E2` (`E3`): by section 4 the E2 (E3) coordinates of q lie in the root range of that chart, in a
   root (the cell contains a strictly convex Q), and (w, t, u) lies in a leaf of that column. Continue with that leaf.
3. The leaf is not V (Q is strictly convex) and not X (k is kept; X asserts a strictly smaller ear). If it is TX, the
   X alternative is excluded for the same reason. So the label is T, TX, S, R, O or C, and F(P) > F*_up; or it is
   N, and P has a labelling in the ball |x − C| < 3/5. ∎
Points on the boundaries of leaves belong to several closed leaves; each closed leaf's assertion holds for them.

## 8. What the run and the self-check establish (`far_verify.py`, producer side)
For each chart: (a) `--complete`: the roots cover the root range (every other cell re-classified); every block
finished. (b) Tiling: the leaves of every column are interior-disjoint boxes inside the column with total volume equal
to the column's (exact rationals). (c) only known labels, no U. (d) `--relabel`: closed-form labels recomputed with the
production label code. (e) `--replay`: C boxes re-proved from scratch. (f) `--table`: the root bounds against the
table; closed-form entries recomputed, LG entries re-proved part by part, parts complete.
**Not established by the self-check**: the correctness of the label code itself and of the closed forms (it re-runs
the same code); that needs an independent implementation from this note.

**Second implementations written from this note (Claude, 2026-11-19; no production module), run on the complete
cover of `pent_far5`** (`../far5_run_20261119/`): `tests/far_root_spot.py` (section 4: no missing root),
`tests/far_tiling_spot.py` (section 7 step 1: 1,724,444 boxes, no failure), `tests/far_label_spot.py` (section 5: all
1,346,929 non-C leaves re-proved, 13 of them after one bisection, `local/accept/tools/label_subdivide.py`). Not covered:
the certified boxes (C) and the LG parts of the table. A floating-point end-to-end test of section 7
(`tests/far_e2e_check.py`: random strictly convex pentagon → kept representation → leaf of the stored cover, through
the hand-over → the leaf's assertion checked from the vertices with a pointwise lower bound of F) complements them.

## 9. What a reviewer should attack first
1. Section 2: the kept-representation argument (in particular ties, the reflection, and that "base = a longest edge"
   is only used to bound the domain); (F1).
2. Section 3: the corner principle for the three charts and both j (list the direction along which each coordinate
   moves each vertex).
3. Section 4: the hand-over to E2 / E3 and the bounds for a, b; the std root range.
4. Section 5, X: the closed forms of the ear comparisons (four untouched-edge cases, four stub polynomials, the
   homogeneous form at the E3 corner face) — re-derive them independently. (Claude re-derived all eight by hand on
   2026-11-17 from the edge vectors of P; they agree with `cutchart.py` / `qchart.py`; the earlier derivation was by
   computer algebra plus a numerical test on 8,000 points.)
5. Section 5, S(a): the concavity argument; O: P(z) ⊂ hull of the corner vertices; N: the two hull arguments.
6. Section 6: the table (the box certificate gives "F(Q) > L for every simple Q of the cell"; a strictly convex Q is
   simple and the determinants are proved positive on the whole cell, so nothing more is needed; see the winding-number
   remark of `LG_LEMMA.md`, computation C).
7. An independent checker of leaves: tiling, and every label from its criterion above.

## Appendix A. Files of a run (for an independent checker)
`RUN/config.json`: the driver settings (heuristics; H₂, H₃). `RUN/<chart>/roots.json`: `{"chart", "blocks": {name:
[[key, lo4, hi4, LB], …]}, "complete", …}`; key = (0, i₀, i₁, i₂, i₃) in the std chart (cell [i/16, (i+1)/16]), four
indices in E2 / E3; LB = the table value of the root (0 if none).
`RUN/<chart>/blocks/<name>/leaves.jsonl`: one line `[key, j, lo, hi, label]` per closed-form leaf (lo, hi: six binary64
numbers, dyadic; the coordinates are those of the chart, then t, u). `cert.jsonl`: one line per certified box,
`{"key", "j", "lo", "hi", "status": "certified", "how": "LG", "N", "sub"}` (N: degree of the family the producer used;
sub: number of sub-boxes its positivity test used; a checker may use any family and any subdivision).
`done.json`: counts and SHA256 of the two files. `RUN/qtable/cells.jsonl`: one line per std root,
`{"key", "LB", "how", "closed", "lg": {"N", "parts": [[lo4, hi4, bound, "lg" | closed-form name], …]}}`.
The leaves and certified boxes of a root and an index j tile root × [0,1]².
Vertex maps: `cutchart.to_std` (chart → q) and `cutchart.verts` (q, t, u, j → five vertices), reproduced in section 2–3.
A checker needs, per label: V, E2/E3, X, N — exact corner evaluations and the formulas of sections 4–5; T, TX — a bound
LB of its own for the w-box (section 6; the table part by re-proving the stored parts); S, R — corner evaluations;
O — an implementation of the fibre bound (`FIBRE_LEMMA.md`) or any other valid lower bound; C — a box certificate.
A checker is free to prove a leaf's assertion by a different valid criterion than the producer used.
