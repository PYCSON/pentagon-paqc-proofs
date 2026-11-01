# p2_fem — second-party upper bound for the regular pentagon (GPT, round PENT-03)

- `independent_p2_ritz.py M OUT.npz` builds a conforming P2 finite-element trial function on a rational pentagon strictly
  inside the regular pentagon (mesh parameter M; the recorded run used M = 24) and saves it.
- `rigorous_ritz_replay.py OUT.npz M RESULT.json` recomputes the Rayleigh quotient of that trial function with exact
  rational monomial integrals and 70-digit directed-rounded Decimal arithmetic, checks that the rational pentagon lies
  inside the exact regular pentagon, and gives λ₁(P₅) ≤ λ₁(Q) ≤ a(w,w)/b(w,w), hence an upper bound for F*.
- `ritz_m24.log`, `ritz_m24_rigorous.log`: the second party's run, F(P₅) ∈ [18.91911822017033649…, 18.91911822017033649…],
  reported upper bound 18.919118220170337.
- Provenance: both sources and both logs are copied verbatim from the appendix of the second party's handoff
  HANDOFF-PENT-03-GPT-20261031.md (SHA-256 789124f34a3b1600b4045b8e86248fe88ca56e76231a3ccd8d864e82db2a51ac).
- **Status: historical log.** The binary trial file `ritz_m24.npz` of that run is not included (it was not kept), so the
  recorded value 18.919118220170337 cannot be replayed digit for digit from this archive. The logs above are the
  second party's own record of that run. Running the first script regenerates a trial function;
  the eigensolver may produce one that differs in the last digits on another platform, and the second script certifies
  whichever trial it is given. The resulting upper bound is valid in every case, but its last digits can differ from the
  recorded value.
- Needs numpy, scipy, sympy, mpmath. Not used by the producer certificates, whose margins are computed against the
  producer's own bound F*_up = 18.919134195733918.
