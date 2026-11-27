# A polynomial certificate for λ₁ of a convex polygon (Claude, 2026-11-16; draft for review, not yet checked by a second party)

**Current status (note added 2026-11-22; the text below is the working note of 2026-11-16 and its status lines are
historical).** Lemma LG, the fan flux (computation B) and the box certificate formulas (computation C) were reviewed by
the second party (GPT, an AI reviewer; PENT-21, 2026-11-17), which re-proved 240 of 240 single-polygon certificates in
exact arithmetic with its own code. In release 1.1 only Lemma LG and the fan construction are used (by Theorem L,
`LOCAL_LEMMA.md`); the near and far covers mentioned below are not part of that release. Not peer reviewed.

Purpose: a rigorous lower bound of the first Dirichlet eigenvalue of a given convex polygon from exact arithmetic with
piecewise polynomials: no eigenvalue solver in the proof, no Bessel functions, no boundary supremum, no restriction
λ₁|P| < 2πj₀₁² in the form of an index argument. Three computations:

- **B (fan; production since the afternoon of 2026-11-16)**: `v2_run/lgfan.py`, tests `v2_run/tests/lgfan_check.py`.
  Piecewise polynomials on the five triangles (centre, edge). About 2 ms per polygon at degree 5, within 10⁻⁴–10⁻³ of
  the true value for F < 24; degree 12 reaches the accuracy of the eigenvalue certificates (1.5·10⁻⁶) in 14 ms.
  Integer arithmetic only. It replaces the eigenvalue certificates (0.6 s each) in the far region.
- **C (box certificate; prototype of the afternoon of 2026-11-16)**: `v2_run/lgbox.py`, `lgpoly.py`, `near_lg.py`,
  tests `v2_run/tests/lgpoly_check.py`. One trial family for a whole box of pentagons; the conditions of Lemma LG become
  polynomial inequalities with exact rational coefficients. Intended for the near region (5–20 CPU-hours in total).
- **A (one piece; the first version, morning of 2026-11-16)**: `v2_run/lgcert.py`, tests `v2_run/tests/lg_check.py`.
  One polynomial on the whole polygon. Kept as an independent cross-check; less accurate (0.97–0.999).

## Lemma LG

Let Ω ⊂ ℝ² be a bounded domain with Dirichlet eigenvalues 0 < λ₁ ≤ λ₂ ≤ …, and let ν ≤ λ₂. Let u ∈ H¹₀(Ω), u ≠ 0, and
w ∈ H(div; Ω) with −div w = u. Put

    A₀ = ‖∇u‖²,  A₁ = ‖u‖²,  A₂ = ‖w‖²   (L² norms),      D = A₀ − νA₁,      G = A₀ − 2νA₁ + ν²A₂.

If D < 0, then G − D > 0 and

    λ₁ ≥ ν(−D)/(G − D).

Proof. Let ψ ∈ H¹₀ be the solution of −Δψ = u. For f ∈ H¹₀: ∫∇f·w = −∫f div w = ∫fu = ∫∇f·∇ψ. Hence w − ∇ψ is
L²-orthogonal to ∇H¹₀, in particular to ∇ψ, and A₂ = ‖∇ψ‖² + ‖w − ∇ψ‖² ≥ ‖∇ψ‖² = (u, A⁻¹u) =: B₂.
Expand u = Σ c_k φ_k in an orthonormal eigenbasis and put p_k = c_k². Then A₁ = Σp_k, A₀ = Σλ_k p_k (u is in the form
domain), B₂ = Σp_k/λ_k, and

    G₀ := A₀ − 2νA₁ + ν²B₂ = Σ p_k (λ_k − ν)²/λ_k ≥ 0,        D = Σ p_k (λ_k − ν).

D < 0 gives λ₁ ≤ A₀/A₁ < ν. Because λ_k ≥ ν for k ≥ 2,

    −D = p₁(ν − λ₁) − Σ_{k≥2} p_k(λ_k − ν) ≤ p₁(ν − λ₁),          λ₁G₀ ≥ p₁(ν − λ₁)².

Therefore λ₁G₀ ≥ p₁(ν − λ₁)² ≥ (−D)(ν − λ₁), i.e. λ₁(G₀ − D) ≥ ν(−D), and G₀ − D ≥ −D > 0. Finally G ≥ G₀ and
t ↦ ν(−D)/(t − D) is decreasing for t > D. ∎

Remarks. (i) This is the one-vector case of the Lehmann–Goerisch method (X = L²(Ω)², T = ∇); the proof above is
self-contained. (ii) With u = φ₁ and w = ∇φ₁/λ₁ the bound equals λ₁. If A₂ = (1 + δ)B₂ and u = φ₁, the bound is
λ₁/(1 + δν/(ν − λ₁)): the loss is of second order in the L² error of the flux, times ν/(ν − λ₁). (iii) No H²
regularity of u is used, so corner singularities of the eigenfunction cost little — unlike Temple's inequality with a
polynomial trial function, which needs Δu ∈ L² to be approximated.

## The bound ν ≤ λ₂

Two classical facts, the larger value is used:
1. *Monotonicity.* If Ω ⊂ R then λ₂(Ω) ≥ λ₂(R). For a rectangle R with sides a, b: λ₂(R) = π² min(4/a² + 1/b², 1/a² + 4/b²).
2. *Krahn–Szegő.* λ₂(Ω) ≥ 2πj₀₁²/|Ω| for every bounded open Ω ⊂ ℝ² (a second eigenfunction has a nodal domain of area
   ≤ |Ω|/2 whose first eigenvalue is λ₂(Ω); Faber–Krahn). Constant used: 36.33682907 < 2πj₀₁².

With the second fact ν/λ₁ ≥ 36.34/F(Ω), so the certificate works whenever F(Ω) is clearly below 36.34; elongated
polygons with larger F are covered by the fibre bound (FIBRE_LEMMA.md).

## Computation B: piecewise polynomials on the fan (lgfan.py)

1. *Frame.* p ↦ (p − c)/s with dyadic c and s = 2^k, so that the image lies in about [−1, 1]²; λ₁ scales by s⁻².
2. *P⁺ and O.* A convex polygon P⁺ with vertices on the grid 2⁻²⁴ that contains the image of P (P inflated by
   1 + 2⁻²¹ about its centroid and rounded), and a grid point O. Verified exactly: P⁺ strictly convex and
   counter-clockwise, every vertex of P in P⁺, O strictly inside P⁺. Domain monotonicity: λ₁(P) ≥ s⁻² λ₁(P⁺); the
   relative loss is about 10⁻⁶. Fan: the triangles T_k = (O, v_k, v_{k+1}), v_k the vertices of P⁺; the segments Ov_k
   are called spokes.
3. *Trial function u.* On T_k let (λ, a, b) be the barycentric coordinates of (O, v_k, v_{k+1}), so
   p − O = a e₁ + b e₂ with e₁ = v_k − O, e₂ = v_{k+1} − O. u|T_k = Σ c_m B_m, B_m = N!/(i! j! l!) λ^i a^j b^l,
   i + j + l = N (Bernstein form). The coefficients are integers times 2⁻³⁰; one global list: the centre (i = N), the
   coefficients on the spokes (shared by the two neighbouring triangles: (i, 0, l) on T_k and (i, l, 0) on T_{k+1}),
   the interior ones; the coefficients with i = 0 (on the boundary edge) are zero. Hence u is continuous, piecewise
   polynomial, zero on ∂P⁺: u ∈ H¹₀(P⁺).
4. *Flux w.* h := (2 + E)⁻¹u on every T_k, E = (p − O)·∇ the Euler operator about O. In the coordinates (a, b), which
   are linear in p − O: write u = Σ ũ_qr a^q b^r (λ = 1 − a − b), then h = Σ ũ_qr a^q b^r/(q + r + 2). Further χ:
   continuous, piecewise of degree n = N + 1 in the same Bernstein form (no boundary condition), integer coefficients
   times a power of two. Put

       w := −(p − O) h + curl χ,      curl χ = (∂_y χ, −∂_x χ).

   Then on every T_k: div((p − O)h) = 2h + Eh = u and div curl χ = 0, so −div w = u. Across a spoke: h is
   continuous (h(O + ρe) = ∫₀¹ τ u(O + τρe) dτ depends only on u along the ray, and u is continuous), so the first
   part is continuous; the normal component of curl χ is the tangential derivative of χ, continuous because χ is.
   A piecewise smooth field with continuous normal components is in H(div; P⁺), and its divergence is the piecewise
   one. So w ∈ H(div; P⁺) and −div w = u: the hypotheses of Lemma LG, by construction (nothing to verify at run time).
5. *Integrals.* With n₁₁ = |e₁|², n₁₂ = e₁·e₂, n₂₂ = |e₂|², det = e₁ × e₂ = 2|T_k| > 0, and
   ∫_T λ^p a^q b^r = det · p! q! r!/(p + q + r + 2)!:
   - ∇f = D_a f ∇a + D_b f ∇b, D_a = ∂_a − ∂_λ, D_b = ∂_b − ∂_λ, |∇a|² = n₂₂/det², ∇a·∇b = −n₁₂/det², |∇b|² = n₁₁/det²,
     so ∫_T ∇f·∇g = [n₂₂ S_aa(f, g) − n₁₂ (S_ab(f, g) + S_ab(g, f)) + n₁₁ S_bb(f, g)]/det, S_βγ(f, g) = ∫ D_β f D_γ g / det;
   - ∫_T |p − O|² h² = det [n₁₁ Q_aa + 2 n₁₂ Q_ab + n₂₂ Q_bb](h, h), Q_aa(h, h) = ∫ a² h²/det etc.;
   - (p − O) × ∇χ = {D_aχ (−a n₁₂ − b n₂₂) + D_bχ (a n₁₁ + b n₁₂)}/det, and (p − O)·curl χ = (p − O) × ∇χ, so
     ∫_T (−(p − O)h)·curl χ = n₁₂ X_{a,a} + n₂₂ X_{b,a} − n₁₁ X_{a,b} − n₁₂ X_{b,b},  X_{s,β} = ∫ s h D_βχ / det.
   All forms are fixed rational matrices of the reference triangle (they depend only on N and n); A₀ = ‖∇u‖²,
   A₁ = ‖u‖², A₂ = ‖w‖² = Σ_k [∫|p − O|²h² + 2∫(−(p − O)h)·curl χ + ∫|∇χ|²] are exact rationals.
6. *Result.* ν ≤ λ₂(P⁺) as above (bounding rectangles in the directions of the edges and the axes; Krahn–Szegő with
   |P⁺|; rational lower bounds 9.8696044 < π², 36.33682907 < 2πj₀₁²), D = A₀ − νA₁ < 0 is checked, and
   F(P) ≥ s⁻² · ν(−D)/(G − D) · |P|, |P| the exact area of P. The stated F_lo is this number times (1 − 10⁻⁷), rounded
   down, so that a re-check with a recomputed flux still passes.

The coefficients of u are a floating-point Rayleigh–Ritz vector, those of χ a floating-point least-squares solution
(minimise ‖w‖² over χ: one linear system with the stiffness matrix of degree n); nothing in the proof depends on their
quality. Why piecewise: on T_k only ONE edge function vanishes, so u/λ behaves like r^(π/α − 1) at a vertex of angle α
(bounded), where the one-piece trial function b·q of computation A needs q ≈ u/b ~ r^(π/α − 2) (unbounded; worst for
nearly flat vertices).

### The stored record (B)

JSON under the key `lg` of a corner record in `corners.jsonl`; sample with explicit polygons:
`v2_run/tests/lg_samples/fan_samples.jsonl` (one object per line: `V` = the vertices of P as binary64 numbers,
counter-clockwise; `lg` = the record).

| field | meaning |
|---|---|
| `kind` | `"fan"` |
| `c`, `s` | the frame (rationals as strings): (x, y) = (p − c)/s |
| `Uplus`, `O` | the vertices of P⁺ (counter-clockwise) and the point O in the frame (rationals as strings) |
| `N`, `u` | the degree and the integer coefficients of u (u = 2⁻³⁰ Σ `u[g]` B; the factor does not matter) in the global order: [0] the centre; then for the spokes k = 0, 1, …: i = 1, …, N − 1 (i = exponent of λ, so from the boundary end towards O); then for the triangles k = 0, 1, …: the interior coefficients (i, j, l ≥ 1) in the order i ascending, then j ascending |
| `n` | the degree of χ used by the producer (a checker may use any) |
| `F_lo` | the claimed bound: F(P) = λ₁(P)|P| ≥ `F_lo` |
| `x`, `xe` | optional: the flux, χ = 2^`xe` Σ `x[g]` B in the global order of degree n without boundary condition: [0] centre; spokes as above; the vertices v_0, v_1, …; the edges k = 0, 1, … with j = 1, …, n − 1 (j = exponent of a); the interiors |
| `grid`, `cb`, `nu`, `rho_plus` | information (ν and the Rayleigh quotient of u on P⁺ as floats) |

A checker verifies the geometry of step 2 exactly, builds u from `u`, takes the stored flux or finds its own χ (any
continuous piecewise polynomial χ will do; a worse one only gives a smaller bound), computes A₀, A₁, A₂ exactly and
verifies `F_lo` ≤ s⁻² ν(−D)/(G − D)|P|.

### Tests (B): `tests/lgfan_check.py`

Known shapes; 300 random convex pentagons against a Rayleigh–Ritz value from ANOTHER trial space; records accepted /
tampered records refused; and the algebra by another route, exactly: u and w are rebuilt as polynomials in (x, y) on
every triangle (h from the homogeneous components about O), then u = 0 on ∂P⁺, u and w·n continuous across the
spokes, −div w = u, and A₀, A₁, A₂ equal to the integrals computed by Green's theorem on every triangle
(`lgcert.integrate`): identical rationals in all five cases (pentagons, a quadrilateral, a triangle; N = 3…6).

### Accuracy and cost (B; not part of the proof)

Time 1.7 ms (N = 4), 1.9 ms (5), 2.2 ms (6), 3.7 ms (8), 7 ms (10), 14 ms (12); record 1.0 kB at N = 5.
Deficit against the true value: regular pentagon 1.8·10⁻⁴ (N = 4), 8·10⁻⁵ (5), 3·10⁻⁵ (6), 6·10⁻⁶ (8), 2.4·10⁻⁶ (10),
1.5·10⁻⁶ (12; the eigenvalue certificate has 1.6·10⁻⁶); square with a cut corner (the low-margin valley) 8·10⁻⁵ at
N = 5; nearly flat vertex 2·10⁻⁶ at N = 6. On 363 corners with eigenvalue certificates from two valley blocks: no
bound above the certified upper bound; ratio to the certified LOWER bound between 0.9997 and 1.0002 (N = 5).
On random pentagons the accuracy at fixed N falls as F grows (factor ν/(ν − λ₁) and thinner triangles): worst ratio
at N = 5: 0.9992 for F < 22, 0.990 for F < 27, 0.90 for F up to 34; at N = 8: 0.9998 up to F = 30, 0.998 up to 34.
For F ≥ 36.3 there is no certificate (D ≥ 0); those pentagons are elongated and covered by the fibre bound.

## Computation C: a whole box of polygons with one trial family (lgbox.py, lgpoly.py; prototype of 2026-11-16 afternoon)

**Observation.** In computation B nothing in steps 3–5 depends on the polygon except through the four numbers
(n₁₁, n₁₂, n₂₂, det) of each triangle. So let the vertices v_k(z) depend on parameters z in a box X, keep O fixed, and
use the SAME coefficient lists for every z, or lists that are affine in z:

    u(z) = u⁰ + Σ_i (z_i − m_i) uⁱ,      χ(z) = χ⁰ + Σ_i (z_i − m_i) χⁱ      (m = the centre of X; "order 1"; "order 0": u⁰, χ⁰ only).

If det_k(z) > 0 for all k (all fan triangles positively oriented), the fan is a polygon P(z), star-shaped with respect
to O, as soon as ONE member of the box is a simple polygon (convexity is not needed; 2026-11-19: the winding number w
of the closed polygon about O is constant on the connected box, because O never lies on a line V_kV_{k+1} there; the five
angles at O lie in (0, π) and sum to 2πw, so w ∈ {1, 2}; w = 2 occurs, e.g. for a pentagram, and is excluded only by a
simple member, which has |w| ≤ 1; with w = 1 the five triangles have disjoint interiors and their union is P(z). Each user
of the lemma supplies such a member: Theorem L the regular pentagon, Theorem F the strictly convex members — the only
ones its conclusions are used for — and Theorem N Lemma S of `NEAR_THEOREM.md`), u(z) ∈ H¹₀(P(z)),
w(z) = −(p − O)(2 + E)⁻¹u(z) + curl χ(z) ∈ H(div; P(z)) and −div w(z) = u(z), exactly as in B. For each triangle the
14 numbers of `lgfan.tri_forms` become quadratic polynomials in z − m with exact rational coefficients
(`lgfan.tri_bilinear`: the forms evaluated on pairs of lists), and

    A₁(z) = Σ_k det_k m_k,    A₀(z) = Σ_k (n₂₂ aa_k − 2 n₁₂ ab_k + n₁₁ bb_k)/det_k,
    A₂(z) = Σ_k [det_k (n₁₁ qaa_k + 2 n₁₂ qab_k + n₂₂ qbb_k) + 2 (n₁₂ (xa1_k − xb2_k) + n₂₂ xb1_k − n₁₁ xa2_k) + (n₂₂ raa_k − 2 n₁₂ rab_k + n₁₁ rbb_k)/det_k].

**Lemma LG on the box.** Take ν(z) = KS/|P(z)| ≤ λ₂(P(z)) (Krahn–Szegő; KS = 36.33682907 < 2πj₀₁², |P| = Σ det_k/2).
Then D < 0 reads Dc := KS·A₁ − |P|·A₀ > 0, the bound is λ₁ ≥ (νA₁ − A₀)/(νA₂ − A₁), and

    Dc(z) > 0   and   C(z) := |P|·Dc − F*·(KS·A₂ − |P|·A₁) > 0      ⟹      F(P(z)) = λ₁|P| > F*.

**Exact polynomial form (lgpoly.py).** If the vertices are affine in z (near chart: the vertices ARE the parameters),
then with 𝔇 := det_0 ⋯ det_4 the functions PD := Dc·𝔇 and PC := C·𝔇 are polynomials in d = z − mid with rational
coefficients (degree 13 and about 13,000 terms for order 1; expanded exactly with FLINT in 0.04 s). For the box
|d_i| ≤ ρ_i it suffices to prove det_k > 0 (k = 0…4), PD > 0, PC > 0. Positivity of a polynomial p from its
coefficients:

    p(d) ≥ min_box (p₀ + g·d + ½ dᵀHd) − Σ_{|α| ≥ 3} |p_α| ρ^α,

    p₀ + g·d + ½ dᵀHd ≥ p₀ − ½ gᵀA⁻¹g − Σ_i μ_i ρ_i²     for any μ ≥ 0 with A := H + 2 diag(μ) positive definite

(add Σ μ_i (d_i² − ρ_i²) ≤ 0 and minimise over all d; μ comes from a floating-point solution of the box problem; A is
factorised as LDLᵀ in exact rationals and its pivots must be positive; alternatively every term by its absolute value).
All numbers are exact rationals; only the sum of the higher terms is computed in binary64, from upper bounds of the
|p_α| and with an explicit allowance (factor 1 + 10⁻⁹; fewer than 10⁶ terms, so this covers any order of summation) for the rounding of the powers, products and the sum. No interval
arithmetic.
If a test fails on the box, the box is bisected and the failing polynomials are re-expanded at the centres of the two
halves by an exact Taylor shift d ↦ d ± h e_i (same trial family, same polynomials; at most 64 parts per leaf); a
polynomial proved positive on a box is not tested again on its parts.

**Why it is efficient.** With order 1 the certified bound follows F(z) to fourth order in z − mid (the family is
first-order exact), so the box size is limited by the function itself, not by the method: one test (0.12 s, plus
0.05 s for the family) certifies cubes of half-width 0.156·r at distance r = 0.05…0.4 from the regular pentagon, and
0.05·r at r = 0.6. The tangent box test T1 needs boxes proportional to the MARGIN (∝ r²), and 64 corner certificates
per box. Probe cubes at r = 0.05, 0.1, 0.2, 0.4 (all leaves certified) extrapolate to **5–20 CPU-hours for the whole
near region** 1/20 ≤ |z − C| < 3/5 (the T1 route: a week of planning and then 10⁶–10⁷ cells).

**Driver (near_lg.py).** Bisection cover of the cube of half-width 5/8 around C (C = the binary64 chart point of the
regular pentagon used by all programs); leaves: `mirror` (sup x₂ < 1/2: the mirror image z ↦ (1 − x₄, x₅, 1 − x₂, x₃,
1 − x₀, x₁) is a congruent pentagon at the same distance from C), `far` (inf |z − C|² ≥ 9/25), `inner`
(sup |z − C|² < 1/400: the explicit neighbourhood of the regular pentagon that statement L4 excludes; it is not covered by v1.0; since 2026-11-17 it is the subject of Theorem L, `LOCAL_LEMMA.md`), `LG` (the test above with a family computed at the centre of the leaf).
The three skip tests are exact rational comparisons. Result: every z with 1/20 ≤ |z − C| < 3/5 has F(P(z)) > F*.
Since every pentagon shape with shape distance d in [1/20, 3/5) has a representative z with |z − C| = d, no
"cheapest representative" logic, no planner and no ψ-comparison are needed. `near_lg_verify.py`: exact tiling of every
block, skip labels repeated, every LG leaf re-proved from scratch. A leaf record is just the box: a checker computes
its own family (any family that passes is a proof) and may bisect further where its own attempt fails.

**Tests (`tests/lgpoly_check.py`, PASS).** The polynomial identity by another route, exactly: at random rational points
the pair is rebuilt as polynomials in (x, y) on every triangle, A₀, A₁, A₂ are integrated by Green's theorem, and
PC(d) = C·𝔇, PD(d) = Dc·𝔇 hold as rational numbers (orders 0 and 1). `poly_positive`: no false acceptance on 400
random cubics. Certified near boxes: the pointwise bound at 18 points is above F* and 0.99965–0.99999 of a
Rayleigh–Ritz value from another trial space. A box around the regular pentagon is refused.

**Status.** Prototype, one afternoon old, reviewed by nobody. If it survives review it replaces, for the near region:
the planner and its skip certificates, the corner certificates and the tangent box test T1.

## Computation A: one polynomial piece (lgcert.py)

1. *Frame* and 2. *P⁺* as in B, with the grid 2⁻¹⁶ and inflation 1 + 2⁻¹³ (relative loss about 2.5·10⁻⁴).
3. *Trial function.* u = b·q, b = ∏_k ℓ_k with ℓ_k the affine function vanishing on the line of the k-th edge of P⁺,
   q = Σ q_ij L_i(x)L_j(y), 0 ≤ i, j ≤ d (Legendre polynomials), coefficients rounded to multiples of 2⁻⁴⁰ relative to
   the largest. u is a polynomial vanishing on ∂P⁺.
4. *Flux.* w̃ = (Σ a_ij L_iL_j, Σ b_ij L_iL_j), 0 ≤ i, j ≤ n = (number of edges) + d + 1, rounded in the same way;
   r := −u − div w̃ (a polynomial with tiny coefficients); w := w̃ + (R, 0) with R the antiderivative of r in x.
   Then −div w = u identically; the program verifies this identity exactly and stops otherwise.
5. *Integrals.* A₀, A₁, A₂ are integrals over P⁺ of the polynomials |∇u|², u², |w|², computed exactly by Green's
   theorem (each edge is parametrised affinely).
6. *Result* as in B, with the factor (1 − 10⁻⁶).

Record (A): `c`, `s`, `Uplus` as in B; `d`, `q`: q = Σ `q[i][j]` L_i(x) L_j(y) (L_0 = 1, L_1 = t,
(k+1)L_{k+1} = (2k+1)tL_k − kL_{k−1}; binary64 entries, to be read exactly); `F_lo`; the flux is not stored.
Sample: `v2_run/tests/lg_samples/lg_samples.jsonl`.

Accuracy and cost (A): d = 4: regular pentagon 0.99964; on 900 pentagons with eigenvalue certificates median 0.9994
near the regular pentagon, 0.9945 in the low-margin valley, minimum 0.972 (nearly flat vertices); d = 8: minimum 0.988.
14 ms (d = 4), 42 ms (d = 8). With a fine grid the gap on the regular pentagon is 2.3·10⁻⁴ (d = 4), 1.7·10⁻⁵ (8),
3.3·10⁻⁶ (12): algebraic convergence about d⁻³·⁵.

## Questions for the reviewer

0. (Computation C, added in the afternoon; the most valuable thing to check if it is right.) (a) The argument that for
   every z of the box the transplanted pair is admissible (star-shaped fan, shared coefficients, affine lists).
   (b) The bound with ν = KS/|P| and the two polynomial conditions. (c) The positivity test (Lagrangian bound, tail).
   (d) The covering argument of the driver: annulus in the chart, mirror, the role of C as a binary64 point.
   (e) Please re-implement `lgpoly.build` + `poly_positive` independently (any CAS) and re-prove a sample of leaves.

1. Is the proof of Lemma LG complete (in particular: the integration by parts for w ∈ H(div), the use of
   A₀ = Σλ_k p_k for u ∈ H¹₀, and the case of a multiple λ₂)?
2. Krahn–Szegő for λ₂: please confirm the statement and give the reference you would cite.
3. Computation B, step 4: is the argument that w ∈ H(div; P⁺) with −div w = u complete (continuity of h across the
   spokes; normal continuity of curl χ; no singular part of the divergence at O)?
4. Please write an independent checker for the records of B (from this note only, without reading `lgfan.py`): the
   exact geometry tests, u from the coefficient list, an own flux, exact A₀, A₁, A₂, the final inequality. Run it on
   `fan_samples.jsonl` and on polygons of your own; report any record it refuses and any polygon where the bound
   exceeds a trustworthy value of F.
5. Counterexamples: nearly degenerate pentagons (very short edges, nearly flat vertices, O close to an edge), F close
   to 36.3.
6. The near region (margins of 10⁻⁴ and below): at N = 12 the certificate is as sharp as the eigenvalue certificate
   and 40 times cheaper. Is there an obstacle to using it there? Ideas for faster convergence at the vertices of the
   regular pentagon (u ~ r^(5/3)): graded sub-triangles near the vertices, or singular functions with rigorous
   quadrature?
