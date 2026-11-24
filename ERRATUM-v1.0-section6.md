# Erratum to version 1.0, Section 6 (cost estimate)

Puyi Sun, version 1.1 (23 November 2026). This file corrects the exploratory cost estimate of Section 6 of
`paper/paper.pdf` (version 1.0, 1 November 2026); the paper itself and the dated version 1.0 release are left unchanged.
The correction was prepared during the second-party review (GPT, rounds PENT-11/12, 6 November 2026) and accepted by
the producer (Claude).

Section 6 of v1.0 reported an estimated **4,174 CPU-days** for the near-region continuum cell-cost proxy, with a bootstrap 95% interval of **1,350–9,305 CPU-days**. Subsequent review identified two systematic accounting effects. First, elapsed-time measurements were influenced by machine load; replacing them with a standardized count of branch-and-bound boxes and a calibrated per-box work coefficient reduces the estimate by a factor of approximately 1.68. Second, the original Dirichlet domain D counts two isometrically reflected representatives of each generic shape, so the corresponding continuum integral is halved.

The corrected **identity-labelling continuum proxy is approximately 1,241 CPU-days**. A constructed and rigorously certified atlas also incurs boundary-cell, reflection, splitting, retry and cache costs. With the measured boundary multiplier Θ≈1.4–1.8, a broad scenario is approximately **1,737–2,234 CPU-days**; the project's particular identity-labelling grid model suggests **1,800–1,900 CPU-days**. These are engineering estimates, not rigorous upper bounds on wall-clock work. The old bootstrap interval is associated with a superseded cost definition and must not be reused as a confidence interval for the new estimate.

**The mathematical certificates are unaffected.** This correction concerns a computational cost proxy only. It does not change the PAQC theorem, local box certificates, or their completed independent replay evidence; nor does it establish any result for all convex pentagons.
