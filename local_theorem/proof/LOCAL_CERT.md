# Certificate of Theorem L: files and what to check (Claude, 2026-11-17; answer to GPT PENT-21)

**Current status (note added 2026-11-22; the text below is the frozen proof of 2026-11-17/18 and its status lines are
historical).** Theorem L (for 0 < |z − c|₂ < R the polygon P(z) is a *simple* pentagon and F(P(z)) > F₀) has passed the
second-party AI-assisted computational review: analysis PENT-21, level A replayed in full PENT-21b, level B recomputed in
full with level A re-proved on the reviewer's own enclosure PENT-21c (route ii). See `second_party/replay_status.json`.
It remains a computer-assisted proof candidate, not peer reviewed. The near and far covers mentioned below (Theorems N
and F) are **not** part of this release and have not been signed off by the second party.


Produced by `python3 lglocal.py run OUT 16` in `OUT/cert/`. All numbers are exact rationals written as `"p/q"` strings
unless marked `approx`. The proof (far/LOCAL_LEMMA.md) splits into two levels; a checker can do them separately.

## Level A: the final positivity proof (sections 7 of LOCAL_LEMMA.md)
**Claim A.** If the true polynomial Q(d) = Σ_{|α|≥2} s_α d^α has s_α ∈ [Tm_α − Em_α, Tm_α + Em_α] for every α
(Tm_α = 0 if absent), then Q(d) > 0 for all d ≠ 0 with |d|₂ < R.

Files:
- `witness.json`: `{"format", "rho_in", "gbits", "Tm": {"a0,a1,a2,a3,a4,a5": "p/q", …}, "Em": {…}}`. Tm, Em are multiples of
  2^−gbits (the centres were rounded and the rounding added to Em). Only exponents of degree ≥ 2 occur; Em ≥ 0.
- `leaves.jsonl`: first line `{"t": "root", "R", "W", "rho_in", "keep"}`; then every leaf of the bisection of the half
  cube {|d_i| ≤ W, d₂ ≥ 0}, in depth-first order: `{"t": "I"|"O"|"P", "m": [6], "r": [6]}` (box m ± r), and for P the
  witness `"w": "dual"` with `"mu"` (six non-negative binary64 numbers, as fractions) or `"w": "crude"`.
- `inner.json` (for comparison only): H, T3, τ and its three parts, the pivots of the 64 vertex matrices.

What to check (all in exact arithmetic; `tests/lglocal_cert_check.py` does it with the standard library only):
1. *Inner cube* |d_i| ≤ rho_in. H_ij from the degree-2 coefficients (Q₂ = ½dᵀHd), T the symmetrised degree-3 tensor
   (Q₃ = Σ T_ijk d_i d_j d_k). τ = τ₂ + τ₃ + τ₄ with
   τ₂ = max_i [Em(2e_i) + ½Σ_{j≠i} Em(e_i + e_j)], τ₃ = rho_in Σ_{|α|=3} Em_α w_α,
   τ₄ = Σ_{|α|≥4} rho_in^{|α|−2}(|Tm_α| + Em_α) w_α, w_α = 1 for pure powers and ½ otherwise.
   The 64 matrices ½H + rho_in Σ_i (±T_i) − τI must be positive definite. (Then Q(d) ≥ dᵀ[½H + Σ d_iT_i − τ]d > 0.)
2. *Tiling.* The leaves lie in the half cube, are pairwise interior-disjoint and their volumes add up to (2W)⁵W.
3. *Leaves.* I: |m_i| + r_i ≤ rho_in for all i. O: Σ_i max(|m_i| − r_i, 0)² ≥ R². P: with
   p(x) = Σ_{2≤|α|≤keep} Tm_α (m + x)^α and err = Σ_α E′_α Π_i (|m_i| + r_i)^{α_i} (E′ = Em, plus |Tm_α| for |α| > keep):
   min over |x_i| ≤ r_i of p(x) is > err. Sufficient: (c₀ − err) + [quadratic bound] − tail > 0, tail = Σ_{|β|≥3}|p_β|r^β,
   quadratic bound = −½gᵀA⁻¹g − Σμ_ir_i² with A = H_p + 2diag(μ) positive definite (c₀, g, H_p: the coefficients of p of
   degree 0, 1, 2), or the crude −Σ|g_i|r_i − ½Σ|H_p,ij|r_ir_j.
4. The half d₂ < 0 is the mirror image under (x₀,…,x₅) ↦ (1 − x₄, x₅, 1 − x₂, x₃, 1 − x₀, x₁) (an isometry of the chart
   fixing c and preserving F). The inner cube is proved for all d, so only the shell uses the mirror.

**No floating point enters any decision of level A** (since 2026-11-17): the error supremum is the exact value of a
polynomial with non-negative coefficients at the far corner, the tail is an exact rational, the multipliers μ are
arbitrary non-negative numbers whose only role is to be good.

## Level B: that the true coefficients lie in Tm ± Em (sections 3–6)
Files:
- `inputs.json`: N, n, xe, CB, O, C_q, and the integer coefficient lists `CI[b]`, `XI[b]` (b = 0: the computed ũ₀, χ̃₀;
  b = 1…6: the first-order lists): u_b = 2^−CB Σ CI[b]·B (degree N Bernstein, numbering of `lgfan.dofs(5, N, True)`),
  χ_b = 2^xe Σ XI[b]·B (degree n, `dofs(5, n, False)`). Everything else is a function of these lists.
- `gram.json.xz`: per triangle the exact 42 × 42 Gram matrix of the basis functions (order: I·7 + b; I = 0: a h_b,
  1: b h_b, 2: ∂_aχ_b, 3: ∂_bχ_b, 4: ∂_aU_b, 5: ∂_bU_b) and the 7 × 7 mass matrix; reference integrals over
  a, b ≥ 0, a + b ≤ 1. Recomputable from `inputs.json` by quadrature of degree ≥ 2n.
- `constants.json`: κ₀ (KS), HB = 2⁻⁵⁰, the enclosures of c, of the fan geometry at c (n₁₁, n₁₂, n₂₂, det per triangle),
  of A₀, A₁, A₂, λ_lo, λ_up, ν, |P(c)|, ε₁², ε₀² (mid, rad), and F_lo, F_up, ε₁, ε₀ as used.
- `rho_c.json`: per triangle the coefficients (enclosures) of ρ_c and 𝔐_cρ_c on the six basis functions of list 0.
- `lowdeg.json`: for every α with |α| ≤ 3: the computed coefficient of 𝔇K (`K_mid`, `K_rad`), the two dual norms
  (`NL` = (Σ_k∫Λ_α·𝔑_cΛ_α)^{1/2}, `NT` = (Σ_k∫|Θ_α|²/det_k(c))^{1/2}), E_lin = ε₁NL + ½ε₀NT, E_quad (majorant of the
  quadratic part), E_shift, and the crude linear bound. Em_α = K_rad + E_lin + E_quad + E_shift (then rounded up);
  for |α| ≥ 4 the crude linear bound replaces E_lin.

- `parts.json.xz` (2026-11-18; diagnostic, the proof uses Em only): for every α of degree ≥ 2 the four parts of the error
  bound, `[K_rad, E_lin_used, E_quad, E_shift, tight]`, as integers in units of 2^−gbits, each rounded up
  (E_lin_used = ε₁NL + ½ε₀NT if tight = 1, i.e. |α| ≤ 3, else the coefficient of the crude majorant). With S the sum of
  the four integers, S − 3 ≤ Em_α·2^gbits ≤ S + 1.

**Full level B (GPT PENT-21b, 2026-11-18: the spot checks pass; the enclosure as a whole still has one implementation.
Update 2026-11-21: done in PENT-21c by route (ii), see the note at the top).**
An independent implementation computes its own enclosure s_α ∈ K′_α ± E′_α of every coefficient (any valid bounds, not
necessarily the majorants of `lglocal.py`). Sufficient for the theorem, either
(i) *inclusion:* |K′_α − Tm_α| + E′_α ≤ Em_α for every α in `witness.json` (and K′_α, E′_α = 0 outside), or
(ii) *re-proof:* level A run again with the checker's own Tm′ = K′, Em′ = E′ (the leaves and multipliers of `leaves.jsonl`
may be reused as hints; a leaf that fails may be subdivided).
The irrational inputs must then be enclosed by the checker itself: the coordinates of the regular pentagon
(V₂ = (1 + cos 72°, sin 72°), V₃ = (½, ½ tan 72°), V₄ = (−cos 72°, sin 72°)) to more digits than the radii of
`constants.json` (about 10⁻⁴³), F₀ ∈ [F_lo, F_up], λ₀ = F₀/|P(c)|, ε₁, ε₀ from its own A₀, A₁, A₂.

What an independent implementation has to reproduce at least (the spot check): the formulas of Proposition R for 𝔇K, 𝔇Λ, 𝔇Θ as polynomials in
d′ = z − C_q with the list weights d′_i − η_i; at least the coefficients of degree 2 (21 numbers: `K_mid` within
`K_rad`) and their bounds. For degrees 0 and 1 `lowdeg.json` shows the computed values and bounds too: the true
coefficients are zero there, so |K_mid| ≤ K_rad + E_lin + E_quad + E_shift is a consistency test.

## The floating-point question raised in PENT-21
- `lglocal.certify_shell`: no floating point in decisions any more (see level A).
- `lgpoly.poly_positive` (used by the near and far box certificates): the sum of the terms of degree ≥ 3 is still
  computed in binary64, but with IEEE multiplications and additions only (no library `pow`), and the comment at that
  place contains the error analysis: n < 10⁶ non-negative terms, at most 64 multiplications per term, no underflow
  (checked), relative error below (n + 66)·2⁻⁵³ < 1.2×10⁻¹⁰ for any order of summation, covered by the factor 1 + 10⁻⁹.
  `exact_tail=True` computes the same sum as an exact rational; the verifiers can use it.
