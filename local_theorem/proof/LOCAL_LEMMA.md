# Theorem L: the regular pentagon is the strict minimiser of F in an explicit ball (Claude, 2026-11-17)

**Current status (note added 2026-11-22; the text below is the frozen proof of 2026-11-17/18 and its status lines are
historical).** Theorem L (for 0 < |z − c|₂ < R the polygon P(z) is a *simple* pentagon and F(P(z)) > F₀) has passed the
second-party AI-assisted computational review: analysis PENT-21, level A replayed in full PENT-21b, level B recomputed in
full with level A re-proved on the reviewer's own enclosure PENT-21c (route ii). See `second_party/replay_status.json`.
It remains a computer-assisted proof candidate, not peer reviewed. The near and far covers mentioned below (Theorems N
and F) are **not** part of this release and have not been signed off by the second party.

**Status (2026-11-17, 14:30): proof candidate. Second-party review PENT-21 (GPT, `reviews/PENT-21/`): the analytic chain
(Lemma S, admissibility, Proposition R, the inequalities of section 5, the logic of sections 7–8) PASSED; the computation
was returned as CHANGES-REQUIRED because the certificate data were not exported and one error sum used floating point.
Both are fixed (certificate: `LOCAL_CERT.md`; all decisions of the shell in exact rational arithmetic); the independent
replay of the certificate by GPT (PENT-21b) is pending. A standalone checker of ours (standard library only) verifies the
final positivity proof from the exported numbers.**
Program: `v2_run/lglocal.py`; tests: `v2_run/tests/lglocal_check.py`; run record: `local_run_20261117/`.
History of the idea and the failed first attempts: `LOCAL_PLAN.md`.

## 0. Statement
Chart: z = (x₀, …, x₅) ∈ ℝ⁶ ↦ P(z) = the pentagon (0,0), (1,0), (x₀,x₁), (x₂,x₃), (x₄,x₅). c = the chart point of the
regular pentagon, F(P) = λ₁(P)·|P|, F₀ = F(P(c)).

> **Theorem L.** For every z with 0 < |z − c|₂ < R, R = 52429/2²⁰ = 1/20 + 1.9×10⁻⁷: P(z) is a simple pentagon
> (section 2) and F(P(z)) > F₀.

By-product: 18.9191028942 < F₀ < 18.9191047951.

Together with the near cover (`near_lg.py`: every z with 1/20 ≤ |z − C_q| < 3/5, C_q the binary64 point with
|c − C_q|_∞ ≤ 2⁻⁵⁰, has F > 18.9191342 > F₀) and the far cover (d ≥ 3/5) this would give F(P) ≥ F₀ for **all convex
pentagons**, with equality only for the regular one. Each of the three parts is unreviewed; non-convex pentagons are not
covered by any of them.

## 1. Inputs from analysis
1. **Lemma LG** (`LG_LEMMA.md`): Ω a bounded domain, ν ≤ λ₂(Ω), u ∈ H¹₀(Ω), w ∈ H(div; Ω) with −div w = u,
   A₀ = ‖∇u‖², A₁ = ‖u‖², A₂ = ‖w‖². If A₀ < νA₁ then λ₁(Ω) ≥ (νA₁ − A₀)/(νA₂ − A₁) (the denominator is positive).
2. **Krahn–Szegő**: λ₂(Ω)|Ω| ≥ 2π j₀₁² > κ₀ := 36.33682907. So ν = κ₀/|Ω| is admissible.
3. Rayleigh's principle, and the invariance of F under similarities and relabellings.

## 2. Reference coordinates
O: a rational point near the centre (the one of the computed family). For any z in the region the fan triangles
T_k(z) = (O, V_k(z), V_{k+1}(z)), k = 0…4, have det_k(z) := e₁ × e₂ > 0 (e₁ = V_k − O, e₂ = V_{k+1} − O; proved on the
whole cube, section 8); since the cube is connected and contains c, the triangles tile P(z), which is star-shaped
about O. On T_k: p = O + a e₁ + b e₂, (a, b) in the reference triangle a, b ≥ 0, a + b ≤ 1. Functions and fields on P(z)
are identified with functions of (k, a, b); ∇_r = (∂_a, ∂_b). With n₁₁ = |e₁|², n₁₂ = e₁·e₂, n₂₂ = |e₂|², det = det_k:

  𝔐_z = det⁻¹ [[n₂₂, −n₁₂], [−n₁₂, n₁₁]],  𝔑_z = 𝔐_z⁻¹ = det⁻¹ [[n₁₁, n₁₂], [n₁₂, n₂₂]]   (per triangle),
  b_z(v, v′) = Σ_k det_k ∫ v v′  (= ∫_{P(z)} v v′),   a_z(v, v′) = Σ_k ∫ ∇_r v·𝔐_z∇_r v′  (= ∫_{P(z)} ∇v·∇v′).

A flux w on P(z) is represented by 𝔰 = det·(its components in the basis e₁, e₂); then ∫_{P(z)} w·w′ = c_z(𝔰, 𝔰′) :=
Σ_k ∫ 𝔰·𝔑_z𝔰′, the physical gradient of v corresponds to 𝔰 = 𝔐_z∇_r v, and ∫_{P(z)} w·∇f = Σ_k ∫ 𝔰·∇_r f.
*Fan flux.* For v continuous, piecewise smooth, and χ continuous: h := (2 + E)⁻¹v, E = a∂_a + b∂_b (i.e.
h(s, θ) = s⁻² ∫₀ˢ t v(t, θ) dt along rays from O), and

  𝔰_z(v, χ) := (−det_k a h + ∂_bχ, −det_k b h − ∂_aχ)   ⟷   w = −(p − O)h + curl χ  on P(z).

w ∈ H(div; P(z)) and −div w = v (LG_LEMMA.md, Computation B). Hence, for f ∈ H¹₀:

  (IP)  Σ_k ∫ 𝔰_z(v, χ)·∇_r f = b_z(v, f).

## 3. The family
Data on the reference: U₁…U₆ continuous piecewise polynomials vanishing on the boundary, X₁…X₆ continuous piecewise
polynomials (computed; arbitrary as far as the proof goes); u₀ the first eigenfunction of P(c), λ₀ = λ₁(P(c)),
h₀ = (2 + E)⁻¹u₀, and χ₀ ∈ H¹ with 𝔰_c(u₀, χ₀) = 𝔐_c∇_r u₀/λ₀, i.e. curl χ₀ = ∇u₀/λ₀ + (p − O)h₀ on P(c) (the right side
is an L² field with zero divergence on a simply connected domain). For e := z − c:

  v_z := u₀ + Σ e_i U_i,  χ_z := χ₀ + Σ e_i X_i,  𝔰_z := 𝔰_z(v_z, χ_z).

(v_z, 𝔰_z) is admissible for Lemma LG on P(z). With |P| = |P(z)|, ν = κ₀/|P|, and

  𝔅_z((v, 𝔰), (v′, 𝔰′)) := |P|(κ₀ + F₀) b_z(v, v′) − |P|² a_z(v, v′) − F₀κ₀ c_z(𝔰, 𝔰′),   Ψ(z) := 𝔅_z((v_z, 𝔰_z), (v_z, 𝔰_z)),

Lemma LG gives: if κ₀ b_z(v_z) > |P| a_z(v_z), then

  (1)  F(P(z)) − F₀ ≥ Ψ(z)/(κ₀ c_z(𝔰_z) − |P| b_z(v_z)),  denominator > 0.

**Exactness.** At z = c: a = λ₀ b, c_c(𝔰_c) = b/λ₀, so Ψ(c) = b[|P|(κ₀ + F₀) − |P|²λ₀ − F₀κ₀/λ₀] = 0.

**Lemma S (∇Ψ(c) = 0).** Let L(z) := F₀ + Ψ(z)/(κ₀c_z(𝔰_z) − |P|b_z(v_z)) and U(z) := |P(z)| a_z(u₀)/b_z(u₀) (the Rayleigh
quotient of the transplanted u₀). Both are smooth near c (rational functions of z whose coefficients are fixed
integrals), L ≤ F ≤ U near c, L(c) = U(c) = F₀. Let φ_g (g ∈ D₅, ten maps) be the relabelling maps of the chart; they
are smooth near c, φ_g(c) = c, F∘φ_g = F. Then L ≤ F = F∘φ_g ≤ U∘φ_g, and averaging over g:
L(z) − F₀ ≤ ∇U(c)·M (z − c) + O(|z − c|²), M := (1/10) Σ_g Dφ_g(c). **M = 0** (computed in ball arithmetic: norm
≤ 5×10⁻¹¹⁹ < 1; M is the projection onto the vectors fixed by the group, so norm < 1 forces M = 0). Hence
∇L(c)·(z − c) ≤ O(|z − c|²) for all z, so ∇L(c) = 0, and ∇Ψ(c) = (denominator at c)·∇L(c) = 0. ∎
(No differentiability of F is used.)

So **Q(e) := (∏_k det_k(z)) Ψ(c + e) is a polynomial in e without constant and linear terms**, and by (1) the theorem
follows from: (i) det_k > 0, (ii) κ₀b_z(v_z) > |P|a_z(v_z), (iii) Q(e) > 0, for all 0 < |e| < R.

## 4. The residual form (Proposition R)
Let ũ₀, χ̃₀ be the computed (piecewise polynomial) approximations, ṽ_z := ũ₀ + Σ e_iU_i, χ̃_z := χ̃₀ + Σ e_iX_i,
𝔰̃_z := 𝔰_z(ṽ_z, χ̃_z), Ψ̃(z) := 𝔅_z of the computed pair with itself. Put e := ũ₀ − u₀ ∈ H¹₀ (the letter e is used for the
function from here on; the parameter is written d), e_h := (2 + E)⁻¹e, λ̂(z) := F₀/|P(z)|, and

  ρ_z := 𝔑_z𝔰̃_z − ∇_r ṽ_z/λ̂(z)   (residual field; ρ_c is the residual of the computed pair at c).

> **Proposition R.** Ψ(z) = K(z) + Σ_k ∫ ∇_r e·Λ_z + Σ_k ∫ e_h (a, b)·Θ_z + q_z(e, e), where
> K = Ψ̃ + 2F₀κ₀ Σ∫ ρ_z·𝔐_cρ_c − F₀κ₀ Σ∫ (𝔐_cρ_c)·𝔑_z(𝔐_cρ_c),
> Λ_z = 2F₀[(κ₀/λ₀)𝔐_cρ_z − |P(z)|𝔐_zρ_z − (κ₀/λ₀)𝔐_c𝔑_z𝔐_cρ_c],
> Θ_z = −2F₀κ₀ (det_k(z) − det_k(c)) [ρ_z − 𝔑_z𝔐_cρ_c],
> q_z(e, e) = |P|(κ₀ + F₀)b_z(e, e) − |P|²a_z(e, e) − F₀κ₀ c_z(𝔱′, 𝔱′),  𝔱′ = 𝔐_c∇_r e/λ₀ − (det_k(z) − det_k(c))(a, b)e_h.

*Proof.* The true pair is the computed pair minus the error pair (e, 𝔱_z), 𝔱_z := 𝔰̃_z − 𝔰_z = (−det_k(z) a e_h + ∂_b e_χ,
−det_k(z) b e_h − ∂_a e_χ), e_χ := χ̃₀ − χ₀ (the terms with U_i, X_i cancel). At z = c: 𝔱_c = 𝔰̃_c − 𝔐_c∇_r u₀/λ₀ =
𝔐_c(ρ_c + ∇_r e/λ₀); for general z, 𝔱_z = 𝔱_c − (det_k(z) − det_k(c))(a, b)e_h. So
Ψ = Ψ̃ − 2𝔅_z((ṽ, 𝔰̃), (e, 𝔱_z)) + 𝔅_z((e, 𝔱_z), (e, 𝔱_z)).
Cross term: by the definition of ρ_z, c_z(𝔰̃, 𝔱_z) = Σ∫ρ_z·𝔱_z + λ̂⁻¹ Σ∫∇_r ṽ·𝔱_z, and Σ∫∇_r ṽ·𝔱_z = b_z(ṽ, e) by (IP) for the
flux 𝔱_z (a difference of two fan fluxes, of ṽ_z and of v_z; ṽ ∈ H¹₀). Also ∇_r ṽ = λ̂(𝔑_z𝔰̃ − ρ_z), so
a_z(ṽ, e) = λ̂ Σ∫𝔰̃·∇_r e − λ̂ Σ∫ρ_z·𝔐_z∇_r e = λ̂ b_z(ṽ, e) − λ̂ Σ∫ρ_z·𝔐_z∇_r e by (IP) again. With |P|²λ̂ = |P|F₀ and
F₀κ₀/λ̂ = κ₀|P| the four terms with b_z(ṽ, e) cancel:

  𝔅_z((ṽ, 𝔰̃), (e, 𝔱_z)) = F₀|P| Σ∫ ρ_z·𝔐_z∇_r e − F₀κ₀ Σ∫ ρ_z·𝔱_z.

Insert 𝔱_z = 𝔐_cρ_c + 𝔱′ and expand c_z(𝔱_z, 𝔱_z) = c_z(𝔐_cρ_c, 𝔐_cρ_c) + 2c_z(𝔐_cρ_c, 𝔱′) + c_z(𝔱′, 𝔱′); collect the terms
without e (K), linear in e (Λ, Θ; 𝔐, 𝔑 are symmetric) and quadratic in e (q). ∎
The algebra is checked numerically with an arbitrary polynomial in the place of e (`lglocal.selftest`: the identity
then holds with the extra term 2κ₀|P|(Σ∫∇_r ṽ·𝔱_z − b_z(ṽ, e)), which is computed; agreement 2×10⁻¹³).

## 5. Bounds for the unknown function e
Normalise u₀ := the L²(P(c)) projection of ũ₀ on the first eigenspace. Then e ⊥ u₀ in L² and in the energy product,
‖∇e‖² − λ₀‖e‖² = A₀ − λ₀A₁ (A_i: the forms of the computed pair at c) and ‖∇e‖² ≥ λ₂‖e‖² ≥ ν‖e‖². Hence, with λ_lo ≤ λ₀ the
bound of Lemma LG at c and λ_up = A₀/A₁ ≥ λ₀:

  ε₁² := ‖∇e‖²_{L²(P(c))} ≤ (A₀ − λ_lo A₁)/(1 − λ_up/ν),   ε₀² := ‖e‖² ≤ ε₁²/ν.

(u₀ ≠ 0 because A₀ < νA₁.) Degree 16: ε₁ ≤ 1.02×10⁻³ (relative 4.6×10⁻⁴), ε₀ ≤ 2.3×10⁻⁴.
*Ray bound.* On every reference triangle ∫ (a + b)² e_h² ≤ ¼ ∫ e² (in the coordinates s = a + b, θ = a/s: da db = s ds dθ,
e_h = s⁻²∫₀ˢ t e dt, Cauchy–Schwarz). In particular ∫ |(a, b)e_h|² ≤ ¼∫e².

## 6. From the identity to bounds on coefficients
Everything is multiplied by 𝔇(z) = ∏ det_k(z) and expanded in d′ = z − C_q (C_q rational, η := c − C_q,
|η_i| ≤ 2⁻⁵⁰, verified in ball arithmetic). All polynomial arithmetic is exact (FLINT `fmpq_mpoly`); numbers that are
not rational or not known exactly are carried as **polynomial balls** (centre polynomial, coefficientwise radius
polynomial ≥ 0; class `PB`): F₀ ∈ [F_lo, F_up] and λ₀ = F₀/|P(c)|; the geometry of P(c) (𝔐_c, det_k(c)); the weights of
the lists, which are d′_i − η_i (centre d′_i, radius 2⁻⁵⁰); ρ_c (the value of polynomial balls at d′ = η).

- **K**: computed as a polynomial ball, 𝔇K ∈ KD.m ± KD.r.
- **Linear part.** 𝔇Λ_z and 𝔇Θ_z are polynomials in d′ whose coefficients are fields in the span of the 6·7 basis
  functions per triangle (a h_b, b h_b, ∂_aχ_b, ∂_bχ_b, ∂_aU_b, ∂_bU_b; b = 0…6 the lists). Their exact Gram matrices
  come from the 14 forms of `lgfan.tri_bilinear` plus the flux–gradient cross terms (`lglocal.tables`). For the
  coefficient of d′^α:
  |Σ∫∇_r e·Λ_α + Σ∫e_h(a, b)·Θ_α| ≤ ε₁ (Σ_k ∫Λ_α·𝔑_cΛ_α)^{1/2} + ½ε₀ (Σ_k ∫|Θ_α|²/det_k(c))^{1/2}
  (Cauchy–Schwarz in the metric of P(c) for ∇e; the ray bound for e_h). Computed exactly for |α| = 2, 3 (`Elin`); for
  all other α a cruder majorant (absolute values of all coefficients times the norms of the basis functions; `Mlin`).
- **Quadratic part** q_z(e, e): majorant `Mq` (every geometric polynomial by its coefficientwise absolute value; the
  integrals of e by ε₀², ε₁² and the extreme eigenvalues of 𝔐_c).
- **Shift from C_q to c.** If B majorises the coefficients t_β of the true polynomial in d′, the coefficients s_α of its
  expansion in d = d′ − η satisfy |s_α − t_α| ≤ (2⁻⁵⁰/2⁻¹⁰)·[B(d + 2⁻¹⁰𝟙) − B(d)]_α.

Result: Q(d) = Σ_{|α|≥2} s_α d^α with **s_α ∈ Tm_α ± Em_α** (`Local.Tm`, `Local.Em`; the computed coefficients of degree
0 and 1 are discarded — the true ones vanish by section 3).

## 7. Positivity
- **Inner cube** |d_i| ≤ ρ_in = 1/128 (`certify_inner`): Q(d) = dᵀ[½H + Σ_i d_iT_i]d + (rest), |rest| ≤ τ|d|², where H, T
  are the quadratic and (symmetrised) cubic coefficients of Tm, and τ collects: the errors of degree 2 (row-sum bound),
  of degree 3 (× ρ), and all terms of degree ≥ 4 including Tm's own (× ρ^{m−2}; |d^α| ≤ ρ^{|α|−2}|d|² for pure powers,
  half of that otherwise). The matrix ½H + Σ d_iT_i − τI is affine in d; it is positive definite at the 64 vertices of
  the cube (exact rational LDLᵀ), hence on the cube, hence Q(d) > 0 for d ≠ 0.
  Numbers (units of F): eigenvalues of H 0.525, 0.557, 5.03, 5.70, 12.76, 15.31; τ = 0.165 = 0.057 (degree-2 errors)
  + 0.013 + 0.095; smallest pivot 0.32.
- **Shell** (`certify_shell`): the cube |d_i| ≤ W = 13/256, d₂ ≥ 0, is bisected; a leaf is `inner` (inside the inner
  cube), `out` (inf |d|² ≥ R²), or proved by the exact polynomial test (`lgpoly.poly_positive`) of
  Σ_{2≤|α|≤6} Tm_α d^α − sup_box Σ E′_α|d^α| > 0 on the box, E′ = Em plus the absolute values of Tm's coefficients of
  degree > 6. 6,942 proved boxes (plus 32 inner and 965 outside leaves: 7,939 leaves tile the half cube), 370 s. Since
  the review PENT-21 every decision is in exact rational arithmetic: the error supremum is the exact value of a polynomial
  with non-negative coefficients at the far corner of the box, and the tail of the polynomial test is an exact rational. The half d₂ < 0 follows from the reflection symmetry
  (x₀,…,x₅) ↦ (1 − x₄, x₅, 1 − x₂, x₃, 1 − x₀, x₁), an isometry of the chart that fixes c and preserves F.

## 8. Side conditions (`side_conditions`, on the cube of half-width W around c)
(i) det_k > 0: exact polynomial test. (ii) Lemma LG for the true family:
κ₀b_z(v) − |P|a_z(v) ≥ (κ₀b_z(ṽ) − |P|a_z(ṽ)) − 2κ₀(b_z(ṽ) b_z(e))^{1/2} − |P|(2(a_z(ṽ) a_z(e))^{1/2} + a_z(e)), with
b_z(e) ≤ 2ε₀² (det_k(z) ≤ 2det_k(c)), a_z(e) ≤ 2ε₁² (𝔐_z ⪯ 2𝔐_c; both checked in ball arithmetic on the cube), upper
bounds of a_z(ṽ), b_z(ṽ), |P| in ball arithmetic, and the exact polynomial test of 𝔇(κ₀b̃ − |P|ã) − rhs·𝔇 > 0 on the
cube. Allowance 0.03 against 7.8.

## 9. The run (2026-11-17; first version 209 s, exact version 383 s, one core)
The numbers below are from the exact version (`local_run_20261117/local.json`, certificate in `local_run_20261117/cert/`,
format in `LOCAL_CERT.md`). The coefficients Tm, Em used by the final proof are rounded to multiples of 2⁻¹⁶⁰ (rounding
added to Em). `tests/lglocal_cert_check.py` (standard library only, no producer code) re-proves the inner cube (64
vertices, smallest pivot 0.0349 = 0.32 in units of F), the exact tiling by the 7,939 leaves and the positivity on all
6,942 P leaves from `witness.json` and `leaves.jsonl`: CERTIFICATE OK (16 minutes; `cert_check_standalone.log`).
The first version (floating-point error sum with a safety factor) is kept in `local_run_20261117/first_run/`; it proved
the same 6,942 boxes.
`python3 lglocal.py run OUT 16`: degree 16 (flux 17), order-1 family at C_q = `near_lg.CF`. symmetry True, inner True,
side True, shell True → `local_run_20261117/local.json`, leaf boxes in `leaves.txt.xz`
(SHA256 of the uncompressed list 313b5767…a7076).
Replicated on a second machine the same day (user's computer, 302 s): PROVED, leaf list byte-identical
(`local_run_20261117/user_machine/`).
Exact tests (`tests/lglocal_exact.py`, PASS): with F₀ fixed to a point, the polynomial 𝔇K and the linear functional
applied to an arbitrary test function equal the formulas of Proposition R evaluated in exact rational arithmetic
(5×10⁻¹⁸ and 3×10⁻¹², the latter limited by c versus C_q); the dual-norm bound holds for the test function on all 84
coefficients of degree ≤ 3 (largest ratio 0.27). Consistency (`Local.low_check`): the computed coefficients of degree 0
and 1 lie within their error bounds (ratio 0.003). End-to-end (`probes_20261117/local_sanity.py`): at 72 points the
certified lower bound of F(z) − F₀ is 0.50–0.99 of the independently computed difference, never above it.
Tests (`tests/lglocal_check.py`, PASS): the identity of Proposition R; containment of polynomial-ball products; Gram
matrices; agreement of the polynomial K with the pointwise formula; the tests fail when they should (degree 8: error of
u₀ too large; error bounds × 30; quadratic form negated).

## 10. What a reviewer should attack first
1. Lemma S (the gradient argument) and the admissibility of the true family (H(div) regularity of the fan flux of u₀;
   existence of χ₀ ∈ H¹).
2. Proposition R: the two uses of (IP), the definition of 𝔱_z, signs.
3. Section 5: the normalisation of u₀ and the three inequalities.
4. Section 6: that the dual-norm bound uses the right metrics (𝔑_c for Λ, 1/det_k(c) for Θ), the polynomial-ball
   bookkeeping (F₀, c, η), the shift.
5. An independent implementation of the final numbers (the 21 + 56 coefficients of degree 2 and 3, their error
   bounds, τ) would be the strongest check; the reference integrals can be recomputed by quadrature.
6. Whether Lemma LG's hypothesis ν ≤ λ₂ is used correctly for P(c) in section 5 and for P(z) in (1).
