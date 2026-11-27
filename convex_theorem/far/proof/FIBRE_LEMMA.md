# The fibre bound: statements and proofs (Claude, 2026-11-16; draft for review, not yet checked by a second party)

Notation. For a bounded open set Ω ⊂ ℝ², λ₁(Ω) = inf { ∫|∇u|² / ∫u² : u ∈ C_c^∞(Ω), u ≠ 0 } is the first Dirichlet
eigenvalue, and F(Ω) = λ₁(Ω)·|Ω|. For an interval I = (a, b) and a measurable V : I → [0, ∞],

    μ(I, V) = inf { ∫_I (f'² + V f²) : f ∈ C_c^∞(I), ∫_I f² = 1 }.

Two monotonicity facts follow from the definition: (M1) if V ≤ W on I then μ(I, V) ≤ μ(I, W); (M2) if I ⊂ J and V is
any extension of V to J, then μ(J, V) ≤ μ(I, V), because C_c^∞(I) ⊂ C_c^∞(J) by extension by zero.

## Lemma OD (one-dimensional reduction)

Let P ⊂ ℝ² be a bounded open convex set, (x, y) orthonormal coordinates, (a, b) the projection of P to the x-axis, and
for x ∈ (a, b) let w(x) > 0 be the length of the section P_x = { y : (x, y) ∈ P } (an interval, by convexity). Then

    λ₁(P) ≥ μ( (a, b), π²/w² ).

Proof. Let u ∈ C_c^∞(P), extended by zero to ℝ². By Fubini,

    ∫_P |∇u|² = ∫_a^b ( ∫_{P_x} u_y² dy ) dx + ∫_ℝ ( ∫_a^b u_x² dx ) dy.

For fixed x, the function y ↦ u(x, y) is smooth with compact support in the interval P_x of length w(x), so the
one-dimensional Dirichlet (Poincaré) inequality gives ∫_{P_x} u_y² dy ≥ (π²/w(x)²) ∫_{P_x} u² dy. Hence, with
V = π²/w²,

    ∫_P |∇u|² ≥ ∫_ℝ ( ∫_a^b (u_x² + V(x) u²) dx ) dy.

For fixed y, the function x ↦ u(x, y) is smooth with compact support in (a, b) (the support of u is a compact subset
of P ⊂ (a, b) × ℝ), so the inner integral is at least μ((a, b), V) ∫_a^b u(x, y)² dx. Integrating over y gives
∫_P |∇u|² ≥ μ((a, b), V) ∫_P u². ∎

Remarks. (i) V is unbounded near a and b; this does not matter, since only compactly supported f enter the definition
of μ. (ii) Equality holds for rectangles with sides parallel to the axes. (iii) Convexity is used only to make the
sections intervals; the direction of the x-axis is arbitrary.

## Lemma PS (positive solution)

Let a = x₀ < x₁ < … < x_N = b, let Vt be constant on each (x_{k−1}, x_k), and let m ∈ ℝ. Suppose φ ∈ C¹[a, b] is C²
on each [x_{k−1}, x_k], satisfies −φ'' + Vt φ = m φ on each open piece, and φ > 0 on (a, b). Then μ((a, b), Vt) ≥ m.

Proof. Let f ∈ C_c^∞(a, b) and v = f/φ; v is C¹, piecewise C², with compact support in (a, b). Then
f' = φ'v + φv' and

    ∫ f'² = ∫ φ'² v² + ∫ φφ' (v²)' + ∫ φ² v'².

The function φφ' is continuous and piecewise C¹, hence absolutely continuous, and v² has compact support, so
∫ φφ' (v²)' = −∫ (φφ')' v² = −∫ (φ'² + φφ'') v². Therefore ∫ f'² = ∫ φ² v'² − ∫ φφ'' v², and with −φ'' = (m − Vt) φ,

    ∫ (f'² + Vt f²) = ∫ φ² v'² + m ∫ f² ≥ m ∫ f². ∎

## The verified computation (`v2_run/fibre.py`)

Input: a convex polygon with vertices p₁, …, p_n in counter-clockwise order (exact binary64 numbers), and a direction
d = (d₁, d₂) ≠ 0 (exact binary64, not normalised). Put D² = d₁² + d₂², x' = ⟨p, d⟩, y' = d × p = d₁p₂ − d₂p₁. The map
p ↦ (x', y') is a rotation followed by the scaling by D, so the image P' satisfies λ₁(P) = D² λ₁(P').

1. *Upper and lower edges.* For a counter-clockwise convex polygon and an edge from pᵢ to pᵢ₊₁, the polygon lies to the
   left of the directed edge. If x'ᵢ₊₁ < x'ᵢ the polygon lies below the line Uᵢ through the edge ("upper edge"); if
   x'ᵢ₊₁ > x'ᵢ it lies above the line Lᵢ ("lower edge"). Hence for every x' in the projection of P' and every pair
   (upper edge i, lower edge j): w'(x') ≤ Uᵢ(x') − Lⱼ(x'). Edges with x'ᵢ₊₁ = x'ᵢ (or of undecided sign in ball
   arithmetic) are not used.
2. *Step minorant.* Choose nodes a' = t₀ < t₁ < … < t_N = b' with [a', b'] containing the projection of P'. On a segment
   [t_{k−1}, t_k], each Uᵢ − Lⱼ is affine, so its maximum is attained at an end point, and

       w'(x') ≤ W_k := min over a set of pairs (i, j) of max( (Uᵢ − Lⱼ)(t_{k−1}), (Uᵢ − Lⱼ)(t_k) )

   for every x' in the segment that lies in the projection of P'. Any non-empty set of pairs gives a valid W_k (the
   program picks, by a floating-point rule, the edges whose x'-range meets the segment). Define Vt = c_k on the k-th
   segment, where c_k is a number ≤ π²/W_k² (c_k arbitrary if W_k ≤ 0, in which case the segment does not meet the
   projection), further reduced to at most (200/(t_k − t_{k−1}))². Then Vt ≤ π²/w'² on the projection of P', and by
   Lemma OD, (M1) and (M2): λ₁(P') ≥ μ((a', b'), Vt).
3. *Positive solution.* For a trial value m, φ is propagated from φ(a') = 0, φ'(a') = 1 through the segments with the
   exact solution formulas (cosh/sinh if c_k > m, cos/sin if c_k < m, affine if c_k = m), in ball arithmetic. Accepted
   if at every node φ > 0 (certainly), and on every trigonometric segment √(m − c_k)·(t_k − t_{k−1}) < π (certainly).
   This implies φ > 0 on (a', b']: a non-trivial solution on a hyperbolic or affine segment has at most one zero, and
   on a trigonometric segment shorter than half a period at most one zero; a zero of a non-trivial solution is simple,
   so φ changes sign there, which is excluded by positive values at both ends. On the first segment φ = sinh(κs)/κ,
   s, or sin(ks)/k with ks < π, positive for s > 0. By Lemma PS, μ((a', b'), Vt) ≥ m.
4. *Result.* λ₁(P) ≥ D² m and F(P) ≥ D² m |P|, with outward rounding.

The value m is found by a floating-point estimate of μ and tried with the factors 0.9985, 0.994, 0.98, 0.95.

## Families of polygons (boxes)

Let z range over a box B ⊂ ℝ⁶ and let P(z) be a convex polygon whose vertices p₁(z), …, p_n(z) are multi-affine
functions of z (affine in each coordinate separately). A multi-affine function on a box is a convex combination of its
values at the corners of the box, with weights that depend only on z. Hence every pᵢ(z) lies in the convex hull H of
the finite set { pᵢ(c) : c a corner of B, i = 1..n }, and P(z) = conv{pᵢ(z)} ⊂ H. By domain monotonicity,

    λ₁(P(z)) ≥ λ₁(H) ≥ fibre bound of H,       F(P(z)) ≥ (fibre bound of H) · min_{z ∈ B} |P(z)|.

In the charts used here the area |P(z)| is multi-affine in z as well (FAR_SPEC 9.9), so its minimum is attained at a
corner. H is computed in exact integer arithmetic (`fibre.hull`). This is the box label O (`cut_run.fibre_label`) and
the quadrilateral bound `qchart.quad_fibre`.

## Accuracy (not part of the proof)

On 650 pentagons with eigenvalue certificates (F between 1.2 F* and 1.9 F*) the bound is 0.924–0.99 of the certified
lower bound, mean 0.955; on the regular pentagon 0.881; on the equilateral triangle 0.760. It is sharp on rectangles.
The loss is the energy of the x-dependence of the transverse profile, which the reduction discards.

## Questions for the reviewer

1. Is every step above correct as stated? In particular the exchange of the order of integration, the use of the
   support of u for fixed y, and the absolute continuity argument in Lemma PS.
2. A reference for Lemma OD in this form (it is the standard lower bound in the spectral theory of thin domains).
3. Step 3: is "positive at both nodes + shorter than half a period" a complete argument for positivity on the segment?
4. The family statement uses only P(z) ⊂ H. Is there a bound for families that does not lose a factor of the order
   (box size)/(width of P)?
