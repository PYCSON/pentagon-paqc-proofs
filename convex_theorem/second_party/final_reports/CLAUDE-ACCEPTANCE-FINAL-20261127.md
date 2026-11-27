# Claude 对用户机器上正式全量结果的验收（2026-11-27 晚）

## 输入（用户 11-27 20:57 附上，SHA256 与用户机器上 `sha256sum` 的输出逐个一致）

| 文件 | SHA256 | 字节 |
|---|---|---|
| `PENT26a_replay_report.tgz`（GPT 重放报告） | `a9cd3b826d39b0e098aacebe1f192d7bb7fd71f2911a82f94671883ca1e48736` | 4,250,528 |
| `claude_accept_report.tgz`（Claude 收货工具 v2.1 在用户机器上的报告） | `4f77091685bd1c4ca8753aeddfe6096838f30d4d4867da6f034bceb20ef42b50` | 648,067 |
| `supplement_2.log` | `c3f9a7d38ff5c3a2b64ab633e1b23447053b39b4e6c3f01ad7a1dc7330579d7b` | 1,308 |
| `PENT24_bundle_b8ab7c1.zip`（原冻结输入，放在 `release/v2-inputs/`） | `1aba5f0a693313b515d9b3b0952747cd61ced2b4b5fb2d3c9c5783333ff2f958` | 14,401,842 |

## 运行经过（用户终端记录，`probes/PENT/local/run_logs/`）

- **主跑**：11-22 19:45 – 11-24 09:35。1,947,824 个任务里 1,947,811 个通过，剩下 13 个近区 LG 盒子。
- **补跑第 1 轮**（N + 2，budget 128）：11-25 22:2x，证出 1 个，是 near/132012/26。
- **补跑第 2 轮**（outside，budget 256）：11-26 11:41–11:42，剩下 12 个全部 PROVED，共 86 个 LG 叶子、109 个 outside 叶子。
- **重放**：11-26 11:42 开始。
  - WSL 在 11-26 17:35 前停过一次，17:50 从断点续跑；补跑第 2 轮重跑时选中 12 个，没有新的结果行。
  - 11-27 08:46 `CHAIN DONE`。
- **收货**：08:50 – 09:00。exit code 0，结论 ACCEPT。

## 1. 比较器 v2（`compare_reports_v2_20261127.log`）：**VERDICT: CONSISTENT**

| 阶段 | 目标 | 分片清单（GPT / Claude） | 结果（两边相同） |
|---|---|---|---|
| near | 38,280 | 316 / 316，逐条相同（main 300 + supplement_1 8 + supplement_2 8） | PROVED 38,280 |
| q_lg | 66,138 | 517 / 517，相同 | PROVED 66,138 |
| q_other | 118,962 | 930 / 930，相同 | PROVED 118,962 |
| far_std | 293,846 | 2,296 / 2,296，相同 | PROVED 293,846 |
| far_E2 | 47,804 | 374 / 374，相同 | PROVED 47,804 |
| far_E3 | 35,865 | 281 / 281，相同 | PROVED 35,865 |
| B_std | 934,801 | 7,304 / 7,304，相同 | PROVED 634,807 + CONDITIONAL_TABLE 299,994 |
| B_E2 | 185,823 | 1,452 / 1,452，相同 | PROVED 185,823 |
| B_E3 | 226,305 | 1,769 / 1,769，相同 | PROVED 226,305 |

合计 1,947,824 个目标，没有 UNPROVED、ERROR、MISSING、REJECT。

- **GPT 的闸门**：`structure`、`all_shards_enumerated`、`q_table_closed`、`N_computational_gate`、`F_computational_gate` 全是布尔 true；`attempt_reject_stages` 为 `[]`；`sample` 为 0；`unscanned_shards.json` 九个阶段全是空表。
- **补跑计划的绑定**（`supplement_plans.json`）：
  - 两轮都绑定同一个 main plan（`2d41b5b3…`）和同一个冻结包 SHA；
  - 第 1 轮代码 ID 是 `8ab8cb73…`，第 2 轮是 `beb45023…`（就是 PENT-26a 的 supplement 更新 ID），第 2 轮只允许那 13 个 ID。
- **结构日志**：GEOMETRY exit 0，errors 0。

## 2. 我自己的收货工具（v2.1，在用户机器上读原始分片）

- 九个阶段全部 PROVED（B_std 另有 299,994 个 CONDITIONAL_TABLE，表是完整的：q_lg 和 q_other 都全部 PROVED）。
- **叶子数**：near 主跑 121,786；far_std 348,623；far_E2 88,427；far_E3 48,842；q_lg 69,065。
- **近区补跑**：
  - 第 1 轮 8 个分片、11 个 LG 叶子、0 个 outside 叶子；
  - 第 2 轮 8 个分片、86 个 LG 叶子、109 个 outside 叶子。这和 `supplement_2.log` 逐行相加的结果一致：LG 1+1+4+15+21+1+22+6+1+4+7+3 = 86，outside 6+6+8+12+11+6+10+12+6+10+12+10 = 109。
  - 每个 outside 叶子都用定理 N 的 C 和 (3/5 + 2⁻⁴⁰)² 精确复核过；主跑里出现 outside 叶子会被拒绝，实际没有出现。
- **逐点抽查**：在 676,840 个已证叶子里均匀抽 400 个，每个查 2 个点，失败 0，最小相对余量 3.119 × 10⁻⁴。
- 最后一行唯一：`ACCEPT: coverage and binding complete for the requested stages`。

## 3. 和 PENT-26a 测试见证的区别

GPT 在自己环境里测试时，13 个盒子用了 94 个 LG 叶子加 118 个 outside 叶子。那是测试数据，不进正式结论。正式结论用用户机器上的数据：第 1 轮 1 个盒子、11 个 LG 叶子；第 2 轮 12 个盒子、86 个 LG 叶子加 109 个 outside 叶子。

## 结论

- 🤝 **用户机器上的正式全量结果，GPT 的独立重放和 Claude 的收货逐分片一致**：1,947,824 / 1,947,824，九个阶段的分片清单逐条相同，五个闸门全为 true。
- 按签字总表，计算这一侧（A01–A12）已经齐了；N、F、总定理的签字还要等 GPT 的最终审查（PENT-26c 预审加这些真实报告）。在那之前，对外措辞仍然是 candidate。
- 本节只说明两边读到的是同样的字节、得出的计数相同；它不是数学重放本身。数学重放是 GPT 的 `replay26a.py` 和我的逐点抽查各自做的。原始分片留在用户机器上（约 15 GB）。
