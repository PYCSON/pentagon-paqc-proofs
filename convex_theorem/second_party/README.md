# second_party — review by the second party (GPT, an AI system), and the producer's review of its checker

**Not peer review.** `replay_status.json` is the precise record of what has been checked independently and what has
not. The independent run on all 1,947,824 targets of Theorems N and F finished on the author's computer, every stored
witness was replayed, and the producer's acceptance program agrees with the replay shard by shard (`final_reports/`).
Theorems N, F and the main theorem are signed by the second party (PENT-26, 2026-11-28).

| path | content |
|---|---|
| `replay_status.json` | status of every part |
| `reports/` | the reviewer's reports (in Chinese): PENT-22 (far region, mathematics and sample structure), PENT-23 (near region and assembly), PENT-24 (independent certificate program, full structure), PENT-24b (long-run package, label criteria), PENT-24c (resumable replay, supplements), PENT-26a (outside subboxes for the 13 hardest near boxes, KS constant), PENT-26b (report comparator, final checklist), **PENT-26 (final review of the full run and signatures on Theorems N, F and the main theorem: `GPT-REPORT-PENT26-FINAL.md`, with `PENT26-FINAL-HANDOFF.md`)**; the producer's acceptance notes of PENT-26a/26b |
| `KS-PROOF.md` | the reviewer's independent rational enclosure of 2π j₀,₁² used by Theorem N |
| `LABEL_ROUTES.md` | the reviewer's ledger: statement and proof of every label and table criterion of its program |
| `independent_checker/` | the reviewer's programs (package r3 with the replay of PENT-24c, and the outside route of PENT-26a: `code/outside26a.py`, `code/supplement26a.py`, `code/replay26a.py`, hashes in `PENT26a_UPDATE_SHA256SUMS`): generation of witnesses, exact replay bound to the frozen covers, structure checks; its own README |
| `independent_checker/PENT24_bundle_b8ab7c1.zip` | the frozen input of the run, byte for byte (SHA-256 `1aba5f0a693313b515d9b3b0952747cd61ced2b4b5fb2d3c9c5783333ff2f958`); the programs of r3 refuse any other input. Two members of an inner archive contain the user name of the author's computer in a path. |
| `final_reports/` | the reports of the full run on the author's computer: the reviewer's replay report `PENT26a_replay_report.tgz`, the producer's acceptance report `claude_accept_report.tgz`, the log of the second supplementary round, the shard-by-shard comparison of the two reports (`compare_reports.log`, verdict CONSISTENT) and the producer's acceptance note (in Chinese) |
| `evidence/` | the reviewer's evidence archives (logs, witnesses of the samples; `PENT26-GPT-EVIDENCE.zip`, SHA-256 `11a8ad9bbd9c37cd1819a8871fb38f785554a52eec4e0b8cc30eee36dbc45380`, is the evidence of the final review); see `evidence/SHA256-OF-ARCHIVES.txt` for the earlier archives. They are kept byte for byte as delivered; paths such as `/mnt/data/…` in them are the reviewer's former working directories and are not meant to be executed as absolute paths |
| `producer_review/` | the producer's review of the reviewer's programs (README.md), its pointwise cross-check and its acceptance tool, with logs |
