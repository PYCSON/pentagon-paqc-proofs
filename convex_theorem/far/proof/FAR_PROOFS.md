# Theorem F: the proofs behind `FAR_THEOREM.md` (Claude, 2026-11-18)

**Status: written by the producer, not reviewed.** `FAR_THEOREM.md` states what every label asserts and by which
criterion; this note proves that each criterion implies its assertion, and derives every closed form from the
definitions. Every algebraic identity below is verified symbolically, from the definitions and independently of the
label code, by `v2_run/tests/far_formula_check.py` (878 exact checks; the last part compares the production formulas
with the ones derived here). Notation: u × v = u₁v₂ − u₂v₁; polygons are counter-clockwise; F(Ω) = λ₁(Ω)|Ω|;
F\* = F*_up = 18.919134195733914.

## 1. Ears
Let P be a strictly convex polygon with vertices V_i, edge vectors e_i = V_{i+1} − V_i, and for an edge k
A_k = e_{k−1} × e_k > 0, B_k = e_k × e_{k+1} > 0, C_k = e_{k−1} × e_{k+1}.

**Lemma 1.** Let C_k > 0. (a) The lines of the edges k − 1 and k + 1 meet in W = V_k + s·e_{k−1} = V_{k+1} − r·e_{k+1}
with s = B_k/C_k > 0, r = A_k/C_k > 0. (b) The triangle (V_k, W, V_{k+1}) has area ear_k = A_kB_k/(2C_k) and lies
outside P along the edge k. (c) The polygon P⁺ obtained from P by replacing the vertices V_k, V_{k+1} by W is strictly
convex, P⁺ = P ∪ triangle, and its edge vectors are those of P with e_k removed, e_{k−1} replaced by (1 + s)e_{k−1}
and e_{k+1} by (1 + r)e_{k+1}.

*Proof.* s·e_{k−1} + r·e_{k+1} = e_k has the unique solution given (cross both sides with e_{k+1} and with e_{k−1}).
The area is ½(s·e_{k−1}) × e_k. The consecutive cross products of the edge vectors of P⁺ are positive multiples of
C_k and of cross products that are positive for P, and the edge directions of P⁺ are those of P without one, so the
total turning is still 2π: P⁺ is strictly convex. ∎

**Lemma 2.** A strictly convex pentagon has an edge with C_k > 0.
*Proof.* With exterior angles θ_i ∈ (0, π), Σθ_i = 2π, one has C_k = |e_{k−1}||e_{k+1}| sin(θ_k + θ_{k+1}), and the five
sums θ_k + θ_{k+1} add up to 4π, so one of them is at most 4π/5 < π. ∎

**(F1).** If P ⊂ Q then λ₁(P) ≥ λ₁(Q), so F(P) ≥ F(Q)·|P|/|Q|.

**Kept edge.** k is *kept* for P if C_k > 0 and ear_k ≤ ear_{k′} for all k′ with C_{k′} > 0. By Lemma 2 every strictly
convex pentagon has a kept edge.

## 2. The cut chart
Let P be a strictly convex pentagon and k a kept edge; Q = P⁺ (Lemma 1) is a strictly convex quadrilateral and
W one of its vertices. Apply the similarity that maps a longest edge of Q to the segment from q₀ = (0,0) to
q₁ = (1,0) with Q above it: Q = Q(q) = (q₀, q₁, q₂, q₃), q = (q₂ₓ, q₂ᵧ, q₃ₓ, q₃ᵧ). The reflection (x, y) ↦ (1 − x, y)
maps Q(q) to Q(q′) with q₀′, q₁′, q₂′, q₃′ the images of q₁, q₀, q₃, q₂; so after a reflection, if necessary,
W ∈ {q₃, q₀}: j ∈ {3, 0}. F and the statement of Theorem F (d(P) is a minimum over all ten labellings) are invariant
under similarities and reflections.

The two vertices of P on the edges of Q at W divide them in ratios t, u ∈ (0, 1) (fraction cut away, measured from W;
by Lemma 1 the fractions are s/(1 + s) and r/(1 + r)):
- j = 3: P = [q₀, q₁, q₂, q₃ + t(q₂ − q₃), (1 − u)q₃], cut edge = edge 3 of this list;
- j = 0: P = [(t, 0), q₁, q₂, q₃, u·q₃], cut edge = edge 4.

Edge vectors of Q: e₀ = (1, 0), e₁ = q₂ − q₁, e₂ = q₃ − q₂, e₃ = −q₃, and the invariants
q₃ᵧ = e₃ × e₀, q₂ᵧ = e₀ × e₁, K = e₁ × e₂ = D − q₃ᵧ + q₂ᵧ, D = e₂ × e₃ = q₂ × q₃, M = e₃ × e₁ = D − q₃ᵧ,
L = e₀ × e₂ = q₃ᵧ − q₂ᵧ, 2|Q| = q₂ᵧ + D. A closed quadrilateral is strictly convex and counter-clockwise iff its four
consecutive cross products are positive: q₂ᵧ, K, D, q₃ᵧ > 0 (the four exterior angles then lie in (0, π) and their sum,
a multiple of 2π below 4π, is 2π).

Edge vectors of P:
- j = 3: E₀ = e₀, E₁ = e₁, E₂ = (1 − t)e₂, E₃ = t·e₂ + u·e₃, E₄ = (1 − u)e₃;
- j = 0: E₀ = (1 − t)e₀, E₁ = e₁, E₂ = e₂, E₃ = (1 − u)e₃, E₄ = t·e₀ + u·e₃.

Cut edge: j = 3: A = (1 − t)u·D, B = t(1 − u)D, C = (1 − t)(1 − u)D, ear_cut = t·u·D/2; j = 0: the same with q₃ᵧ in place
of D, ear_cut = t·u·q₃ᵧ/2. Hence |P| = |Q|(1 − t·u·τ_j), τ₃ = D/(q₂ᵧ + D), τ₀ = q₃ᵧ/(q₂ᵧ + D); 0 < τ_j < 1.

**Admissible q:** Q(q) strictly convex and the base a longest edge, |e₁|, |e₂|, |e₃| ≤ 1. Then q₂ᵧ, q₃ᵧ ∈ (0, 1],
q₂ₓ ∈ [0, 2], q₃ₓ ∈ [−1, 1]: q ∈ 𝔅 = [0,2] × [0,1] × [−1,1] × [0,1].

**Proposition 3.** Every strictly convex pentagon is, up to a similarity and a reflection, P(q, t, u; j) with q
admissible, t, u ∈ (0, 1), j ∈ {3, 0}, and the cut edge kept. ∎

## 3. Charts and the corner principle
std: w = q. E2: w = (q₃ₓ, q₃ᵧ, a, h) with q₂ = q₁ + a(q₃ − q₁) + h(0, 1); then K = h(1 − q₃ₓ).
E3: w = (q₂ₓ, q₂ᵧ, b, h) with q₃ = b·q₂ + h(0, 1); then D = h·q₂ₓ. z = (w, t, u).

**Lemma 4 (ruled charts).** In each chart and for j = 3, 0, every vertex of P(z) is affine in each coordinate of z
separately, and for each coordinate z_m the five derivatives ∂_mV_i are parallel to one common vector d_m (which
depends on the other coordinates). Directions:

| | std | E2 | E3 |
|---|---|---|---|
| 1st, 2nd coordinate | q₂ along (1,0), (0,1) | q₃ along (1,0), (0,1); q₂ by a times the same | q₂ along (1,0), (0,1); q₃ by b times the same |
| 3rd | q₃ along (1,0) | a: q₂ along q₃ − q₁ | b: q₃ along q₂ |
| 4th | q₃ along (0,1) | h: q₂ along (0,1) | h: q₃ along (0,1) |
| t | j = 3: V₃ along q₂ − q₃; j = 0: V₀ along (1,0) | same | same |
| u | j = 3: V₄ along −q₃; j = 0: V₄ along q₃ | same | same |

(the vertices of P are combinations of q₀, …, q₃ with coefficients 1, t, 1 − t, u, 1 − u, so a coordinate of w moves them
along the direction along which it moves q₂ and q₃). *Proof:* symbolic check, all 36 cases. ∎

**Corollary 5 (corner principle).** Let f = (V_a − V_b) × (V_c − V_d) for vertices of P(z) (or of Q(w)), or a linear
functional of the vertices, or the signed area. Then f is multi-affine in z, and on a box its minimum and maximum are
attained at corners. Every vertex V_i(z), z in a box, is a convex combination of its 64 corner values, with weights
that do not depend on i.
*Proof.* ∂²_m of a cross product of two differences is 2(∂_m u) × (∂_m v) = 0 by Lemma 4; the area is a sum of such
products. A multi-affine function is the multilinear interpolant of its corner values: f(z) = Σ_c ω_c(z) f(z_c),
ω_c(z) = Π_m (the barycentric weight of z_m), ω_c ≥ 0, Σω_c = 1. ∎

**Lemma 6.** If g is concave in each coordinate separately on a box, its minimum over the box is attained at a corner.
(Move one coordinate at a time to the end at which g is not larger.) ∎

All corner vertices are exact dyadic rationals (the programs assert that the binary64 evaluation is exact), and corner
values of polynomials are enclosed in ball arithmetic.

## 4. Domains
**std roots.** Cells of side 1/16 of 𝔅. A cell is not a root only if (i) one of q₂ᵧ, q₃ᵧ, K, D is ≤ 0 at its 16
corners — then, by Corollary 5, it is ≤ 0 on the cell and the cell contains no strictly convex Q; or (ii) one of
|e₁|², |e₂|², |e₃|² exceeds 1 on the whole cell (interval arithmetic) — then it contains no admissible q. So every
admissible q lies in a root.

**Hand-over E2.** A q-box of the std chart gets the label E2 if q₃ᵧ ≥ 1/4 on it and K − H₂(1 − q₃ₓ) < 0 at its 16
corners (H₂ = 1/8), hence on the box (multi-affine). Let q be admissible in such a box. From 0 < K < H₂(1 − q₃ₓ):
1 − q₃ₓ > 0. Put a = (1 − q₂ₓ)/(1 − q₃ₓ), h = K/(1 − q₃ₓ) ∈ (0, H₂). Then q₂ = q₁ + a(q₃ − q₁) + h(0,1) (the
x-component is the definition of a; the y-component is the identity K = q₃ᵧ(q₂ₓ − 1) + q₂ᵧ(1 − q₃ₓ)). Moreover
a(q₃ − q₁) = e₁ − h(0,1), |e₁| ≤ 1, |q₃ − q₁| ≥ q₃ᵧ ≥ 1/4, so |a| ≤ 4(1 + 1/8) = 4.5; and |q₃| ≤ 1. So the E2
coordinates of q lie in [−1,1] × [1/4,1] × [−4.5,4.5] × [0,1/8], the root range of the E2 run.

**Hand-over E3.** Label E3 if q₂ₓ ≥ 1/4 on the box and D − H₃q₂ₓ < 0 at the corners (H₃ = 1/8). For an admissible q
of the box: b = q₃ₓ/q₂ₓ, h = D/q₂ₓ ∈ (0, H₃), q₃ = b·q₂ + h(0,1), |b| ≤ |q₃ₓ|/q₂ₓ ≤ 4, q₂ₓ ∈ [1/4, 2], q₂ᵧ ∈ (0, 1]:
the E3 coordinates lie in [1/4,2] × [0,1] × [−4,4] × [0,1/8].

**Edge-chart roots.** Cells of side 1/16 of these ranges, except those on which one of q₂ᵧ, q₃ᵧ, K, D is ≤ 0 at the 16
corners. (Admissibility is not used in the edge charts.)

## 5. Labels
Throughout, a *member* of a box is a point z = (w, t, u) of it with Q(w) strictly convex and t, u ∈ (0, 1); then
P(z) is a strictly convex pentagon. Ranges of invariants over the w-box are corner ranges (Corollary 5); for members
q₂ᵧ, q₃ᵧ, K, D are positive, so their ranges are clipped to [0, ∞).

**V.** One of q₂ᵧ, q₃ᵧ, K, D is ≤ 0 at the 16 corners ⟹ no member. ∎

**T.** Let LB ≤ F(Q(w)) for all strictly convex Q(w) of the w-box (section 6), and τ̄ ≥ τ_j on them
(τ₃ = D/(q₂ᵧ + D) is increasing in D and decreasing in q₂ᵧ: take sup D, inf q₂ᵧ; τ₀ ≤ sup q₃ᵧ/inf(q₂ᵧ + D); always
τ_j ≤ 1). Put θ_T = (1 − F\*/LB)/τ̄. For a member with t·u < θ_T, by (F1),
F(P) ≥ F(Q)(1 − t·u·τ_j) ≥ LB(1 − t·u·τ̄) > LB(1 − θ_Tτ̄) = F\*. Criterion: sup t · sup u < θ_T. ∎

**X.** Claim: for every member some edge k′ ≠ cut has C_{k′} > 0 and ear_{k′} < ear_cut; so no member is the kept
representation of its pentagon. With the edge vectors of section 2 (all eight verified symbolically):

| j | k′ | C_{k′} | ear_{k′} < ear_cut ⟺ (given C_{k′} > 0) |
|---|---|---|---|
| 3 | 0 | (1 − u)M | t·u > q₂ᵧq₃ᵧ/(D·M) |
| 3 | 1 | (1 − t)L | t·u > q₂ᵧK/(L·D) |
| 3 | 2 | t·K − u·M | K(1 − 2t) + t·u·M < 0 |
| 3 | 4 | u·q₃ᵧ + t(q₂ᵧ − q₃ᵧ) | q₃ᵧ(1 − 2u) − t·u(q₂ᵧ − q₃ᵧ) < 0 |
| 0 | 1 | (1 − t)L | t·u > q₂ᵧK/(q₃ᵧL) |
| 0 | 2 | −(1 − u)M | t·u > D·K/(q₃ᵧ(−M)) |
| 0 | 0 | u·M + t·q₂ᵧ | q₂ᵧ(1 − 2t) − t·u·M < 0 |
| 0 | 3 | u·D + t(q₂ᵧ − q₃ᵧ) | D(1 − 2u) − t·u(q₂ᵧ − q₃ᵧ) < 0 |

Derivation of a typical line (j = 3, k′ = 2): A = E₁ × E₂ = (1 − t)K, B = E₂ × E₃ = (1 − t)u·D,
C = E₁ × E₃ = t·K − u·M, so ear = (1 − t)²u·K·D/(2C), and ear < t·u·D/2 ⟺ (1 − t)²K < t(t·K − u·M) ⟺
K(1 − 2t) + t·u·M < 0 (u, D, C > 0). For the untouched edges (rows 1, 2, 5, 6) the ear does not depend on t, u:
e.g. j = 3, k′ = 0: A = (1 − u)q₃ᵧ, B = q₂ᵧ, C = (1 − u)M, ear = q₂ᵧq₃ᵧ/(2M).
Criteria on a box:
(i) untouched edges: the sign condition (M > 0, L > 0 or M < 0) holds for the corner range, and
inf t · inf u > ḡ, ḡ = (sup of the numerator's factors)/(inf of the denominator's factors, all positive);
(ii) stub edges (rows 3, 4, 7, 8): C_{k′} and the polynomial on the right are multi-affine in z in all three charts
(they are polynomials in t, u and multi-affine invariants, of degree one in each variable: checked symbolically);
C_{k′} > 0 and the polynomial < 0 at the 64 corners, hence on the box;
(iii) chart E3, boxes with inf q₂ᵧ = 0 = inf h and inf b ≥ 0. In E3 the invariants are homogeneous of degree one in
(q₂ᵧ, h): with r = h/q₂ᵧ ∈ (0, ∞) and x = q₂ₓ,
K = q₂ᵧ(r(x − 1) + 1 − b), L = q₂ᵧ(r − (1 − b)), M = q₂ᵧ(r(x − 1) − b), q₃ᵧ = q₂ᵧ(b + r), D = q₂ᵧ·x·r,
and the conditions of rows 1, 2, 5, 6 become conditions on (x, b, r, t·u) alone. They are verified in ball arithmetic
for all r ∈ (0, 1] and, after division by r², for all s = 1/r ∈ [0, 1], with t·u replaced by inf t · inf u (each
condition has the form t·u·(positive) > …); a piece of the r- or s-range on which K ≤ 0 has no member. ∎

**TX** (a whole column w-box × [0,1]²): θ_T > ḡ (ḡ from (i)). A member has t·u < θ_T, and then F(P) > F\* as in T, or
t·u ≥ θ_T > ḡ, and then it is not kept as in X(i). ∎

**S(a).** (Revised 2026-11-24 after GPT PENT-26b: no division by a corner width.) For a vector n ≠ 0 let
p_i = ⟨V_i, n⟩ and W_n = max p_i − min p_i. Let τ ≥ 0 be a rational number such that at each of the 64 corners z_c and
for every pair i, k: |P(z_c)| − τ|n|⁻²(p_i(z_c) − p_k(z_c))² ≥ 0 (|P(z_c)| the signed area). The producer takes
τ = min over the corners of |n|²|P(z_c)|/W_n(z_c)², which is such a number whenever every corner width is positive; a
degenerate corner with W_n = 0 needs no division. For each pair i, k the function |P(z)| − τ|n|⁻²(p_i(z) − p_k(z))² is
concave in each coordinate (multi-affine minus the square of a multi-affine function), so by Lemma 6 it is ≥ its
smallest corner value, which is ≥ 0. Hence |P(z)| ≥ τ|n|⁻²W_n(z)² on the box. A strictly convex member has W_n > 0, lies
in a strip of width W_n/|n|, so λ₁(P) ≥ π²|n|²/W_n² and F(P(z)) ≥ π²|n|²|P|/W_n² ≥ π²τ. Criterion: π²τ > F\*. ∎

**S(b).** For the edge from V_a to V_b of length l > 0, with max_k c_k > 0, put c_k = (V_b − V_a) × (V_k − V_a) ≥ 0 and H = max_k c_k/l. P lies
in a strip of width H and contains the triangle over the edge with apex at a farthest vertex: |P| ≥ lH/2,
F(P) ≥ π²·(lH/2)/H² = π²l²/(2 max_k c_k). On a box: max_k c_k ≤ the largest corner value of all c_k (multi-affine), and
l² ≥ (inf |Δx|)² + (inf |Δy|)² from the corner ranges of the two components (multi-affine). ∎

**R.** Let m and m + 2 be edges of P with C_m, C_{m+2} > 0. P is the intersection of its five edge half-planes; the
three half-planes of the edges m + 1, m + 3, m + 4 intersect in a set T ⊃ P. Applying Lemma 1 to the edge m of P and
then to the edge m + 2 of P⁺ (whose A, B, C for that edge are those of P up to the common factor 1 + r in A and C, so
its ear is the same) shows that T is a triangle with |T| = |P| + ear_m + ear_{m+2}. By (F1) and Pólya's theorem for
triangles, F(P) ≥ (4π²/√3)/(1 + (ear_m + ear_{m+2})/|P|). On a box: ear ≤ sup A · sup B/(2 inf C) with inf C > 0, and
|P| ≥ the smallest corner area (all multi-affine). (The area identity and P ⊂ T are also checked on 500 exact random
cases.) ∎

**O.** By Corollary 5 every vertex of P(z) lies in the convex hull H of the 320 corner vertices, so P(z) ⊂ H,
λ₁(P(z)) ≥ λ₁(H) ≥ the fibre bound of H (`FIBRE_LEMMA.md`), and |P(z)| ≥ the smallest corner area. ∎

**N.** For a labelling (a cyclic order W₀, …, W₄ of the vertices of P, possibly reversed) the chart coordinates are
u_p = (W_p − W₀)/(W₁ − W₀) ∈ ℂ, p = 2, 3, 4 (conjugated for the reversed order); they are similarity invariant. N_p =
W_p − W₀ and Δ = W₁ − W₀ are multi-affine, complex valued. For constants r_p ∈ ℂ and a unit vector ε:
|N_p(z) − r_pΔ(z)| ≤ the largest corner modulus (a convex combination of the corner values), and
|Δ(z)| ≥ Re(Δ(z)ε̄) ≥ the smallest corner value of Re(Δε̄) =: δ, required positive. So |u_p(z) − r_p| ≤ ρ_p :=
max_c |N_p − r_pΔ|/δ, and |x(P(z)) − C|² ≤ Σ_p (|r_p − C_p| + ρ_p)². Criterion: this sum is < 9/25. ∎

**C.** The box certificate (`LG_LEMMA.md`, computation C) proves that the fan determinants are positive on the whole
box and F(P(z)) > F\* for every member that is a simple polygon; a strictly convex P(z) is one (with O strictly on the
inner side of every edge line, the winding number is 1). The label is only used for strictly convex members. ∎

## 6. The bound LB for a w-box
Each of the following is a lower bound of F(Q(w)) for every strictly convex Q(w) of the w-box; LB is their maximum.
1. 2π² (Pólya–Szegő for quadrilaterals).
2. Triangle bound: for an edge k of Q with C_k^Q > 0 on the box, Q lies in the triangle Q ∪ ear (Lemma 1), so
   F(Q) ≥ (4π²/√3)/(1 + ear_k/|Q|), ear_k/|Q| = A_kB_k/(C_k(q₂ᵧ + D)) with (A, B, C) = (q₃ᵧ, q₂ᵧ, M), (q₂ᵧ, K, L),
   (K, D, −M), (D, q₃ᵧ, −L) for k = 0, 1, 2, 3; evaluated with sup A, sup B, inf C > 0, inf (q₂ᵧ + D) > 0.
3. Strips: Q lies in 0 ≤ y ≤ w̄ := sup max(q₂ᵧ, q₃ᵧ): F(Q) ≥ π²|Q|/w̄² ≥ π² inf(q₂ᵧ + D)/(2w̄²); and Q contains the
   triangle over the base with apex at its highest vertex, |Q| ≥ max(q₂ᵧ, q₃ᵧ)/2, so F(Q) ≥ π²/(2 max(q₂ᵧ, q₃ᵧ)) ≥ π²/(2w̄).
4. Fibre bound of the hull of the corner quadrilaterals (as O), times inf |Q|.
5. std chart: Makai's inequality for convex domains, λ₁ ≥ π²Per²/(16|Q|²), i.e. F(Q) ≥ π²Per²/(16|Q|), with the
   perimeter from below and the area from above (interval arithmetic on the box).
6. std chart: the LG table value of the root cell containing the w-box (`FAR_THEOREM.md` section 6): a minimum over
   parts that cover all strictly convex Q of the cell, each part with a closed-form bound or a box certificate.

## 7. Where each statement is checked
| statement | check |
|---|---|
| Lemma 1; invariants; edge vectors; area; ear of the cut edge | `far_formula_check.py`, sections 1, 2 (symbolic) |
| the eight comparisons of X; multi-affinity of the stub polynomials | section 5X (symbolic) |
| Lemma 4, all charts and j | section 3 (symbolic) |
| hand-over identities | section 4 (symbolic) |
| homogeneous forms at the E3 face | section 5X3 (symbolic) |
| R: area identity, inclusion | section 5R (exact, 500 random cases) |
| production formulas = the formulas of this note | last section (symbolic) |
| assertions of the labels on real leaves, pointwise | `tests/far_leaf_sample.py` (about 21,000 points, floating point) |
| Proposition 3 and the root ranges, forward | `tests/far_forward_check.py` (20,000 random convex pentagons) |

Not covered by any of these: the inequalities quoted from the literature (Pólya, Makai, Krahn–Szegő), Lemma LG,
the fibre lemma (own notes), and the implementation of the criteria in ball arithmetic (re-run by `far_verify.py`,
to be re-implemented by the independent checker).
