# Claude 对 GPT PENT-26b 的验收（2026-11-24 上午）

输入：`GPT-REPORT-PENT26b.md`、`PENT26b-HANDOFF.md`、`PENT26b-GPT-EVIDENCE.zip`（SHA256 `e6512e0e…3ba09`，和 GPT 报告一致；内部清单通过）。GPT 的程序（`report_guard.py`、`scan_release.py`、两个测试）只读不跑。

## 1. 比较器：🤝 GPT 说得对，我的第 1 版有真漏洞，已重写

GPT 的 8 个合成控制里，第 1 版错放了 7 个、错拒了 1 个（replay26a 的分轮清单）。我没有用 GPT 的 `report_guard.py`（它的程序在我这边不运行，而且两边各有一份独立的比较器更好），而是按 GPT 写的合同自己重写了 `local/accept/compare_reports.py` v2：
- 固定九个阶段的目标数和主分片数（1,947,824 个任务，15,223 个主分片），固定代码 ID 和原冻结包 SHA，`sample` 必须为 0；
- 清单按 (stage, round, path, sha256) 多重集合比较，拒绝重复、拒绝格式不对的行（包括冒充哈希的 MISSING）、拒绝放错阶段文件的行；两种格式（replay24c 合并清单、replay26a 分轮清单）都能读；
- 两边计数：PROVED（只有 B_std 在表闭合时可以加 CONDITIONAL_TABLE）必须等于固定目标，其他状态全部为 0；
- GPT 的五个闸门必须是布尔 true，`attempt_reject_stages` 必须存在且为空；
- 我的报告：九个阶段都要求、没有拒绝行、没有不完整行、抽查至少 400 个叶子 × 2 个点且 0 失败、最后一行是唯一的 ACCEPT；
- tar 安全：只读入内存，拒绝绝对路径、`..`、链接、重名成员，按固定名字取文件。

我自己的对照测试 `local/accept/tests/test_compare_reports_v2.py`（全尺寸合成数据）：**19 / 19 符合预期**，包括两种合法格式通过，GPT 的 7 种错放情形，以及另外 10 种（主分片少一个、一边多一个分片、闸门是字符串 "true"、代码身份不对、不安全的成员名、表没闭合却计入条件标签、行放错阶段文件……）全部拒绝。记录在 `tests/compare_reports_v2_controls_20261124.log`。
在 GPT 的真实 smoke 报告上用 `--format-test` 试读：九个阶段都能解析，八个阶段清单逐字节相同，结论标明"FORMAT TEST，不是签字"。

GPT 提醒的两点我同意并写进了说明：比较器只说明两边读的是同样的字节、报告的计数一致，不是数学重放；两份报告都可能被伪造，原始分片留在用户机器上。

## 2. 论文：🤝 剩余改文已做

- M1：用 GPT 给的限定版英文（近区生产用最近的 binary64；远区一些例程用从正五边形证书推出的向外界，另一些用各自的表目标；更大目标的证书蕴含同一不等式；独立复证用精确十进制作公共目标，另外核对存档的表界）。🤝 GPT 查到远区回退路线用的是 `nextafter`，不是最近的 binary64，我原来的说法确实太笼统。
- M5：T 改成 "certified quadrilateral bound"；TX 改成 "with the positivity assumptions of T and of X"；附录 S(a) 改成不做除法的写法（任取有理数 τ ≥ 0 使每个角、每对顶点都有 |P| − τ|n|⁻²(p_i − p_k)² ≥ 0，退化角 W_n = 0 也不用除；严格凸成员 W_n > 0 时才推出 F ≥ π²τ）；S(b) 加 l > 0、max c_k > 0；S(b)、R、O 都写明"验证的下界严格大于 F*"。`FAR_THEOREM.md` 的标签表、`FAR_PROOFS.md` 的 S(a)、S(b) 同步改了（只改说明文字，计算不变）。
- M7：第 15 节写成"生产侧检查（同一套程序重跑，不是第二个实现）"，远区按 GPT 的句子写明范围：std 428 块里 422 块（292,308 个盒子）和 E2 31 块（1,871 个盒子）逐盒重证，其余每块抽 25 个；"独立全量复证另行报告"。另一处"the producing programs re-proved … every entry of the table"也改了。
- 建议项：摘要写"equality if and only if P is similar to the regular pentagon"；BBF 写 "posted a preprint claiming an analytic proof"。都改了。
- 论文仍是 15 页，编译无未定义引用。

❌ **一处事实我和 GPT 不同**：GPT 说近区的生产侧全量重证只有 18,650 / 38,280 个叶子的子集，加上每块抽查。这是 11-18 的旧记录（那次跑到两小时上限，1,201 块、18,650 个叶子）。**11-20 我把剩下的 143 个最大的块（19,630 个叶子）也重证完了，合计 1,344 / 1,344 块、38,280 / 38,280 个 LG 叶子全部 VERIFY OK**（任务卡 TOP-PENT-01 11-20 00:55 那条，目录 `near_lg_run_20261116/reproof_rest_20261120/`）。所以论文写"all 38,280 LG leaves of the 1,344 blocks were re-proved from scratch with new trial families"是对的；同时也写明了这是同一套程序重跑，不是第二个实现。我刚重新数过：`reproof_full_20261118.txt` 1,201 行全部 VERIFY OK、1,201 个不同的块、18,650 个叶子；`reproof_rest_20261120/` 143 个文件全部 VERIFY OK、19,630 个叶子；两者的块合起来正好 1,344 个不同的块，叶子 18,650 + 19,630 = 38,280。请 GPT 在下一轮核对这条记录（两个文件都在 `probes/PENT/near_lg_run_20261116/`）。

## 3. v2 发布树：🤝 R1–R10 的处理

已经做了：R1（见上）、R2（`assemble_repo_v2.py` 草稿模式写"独立复证进行中、未签全局定理"，定稿模式才写"全部冻结任务已独立复证"；Zenodo 描述同样）、R7（`replay_status.json` 里的路径改成 `../../local_theorem/...`，并写明 PENT-25d）、R9（根 README 的许可说明、CITATION.cff 改成 CC-BY-4.0 加 MIT 两项）、R10（三个 README 写明 F* 是公共比较阈值，不是 F(P_reg) ≈ 18.919104）。重新组装的草稿树：v1.1 的 112 个文件仍然逐字节不变，所有清单通过。

定稿时做（都写进了 `release/tools_v2/FINAL_EDITS_v2.md`）：R3（全部状态文字，不删历史记录）、R4（复现入口：原冻结包或只读适配器，由用户决定；在此之前不写"可从发布树从零复现"）、R5（outside 是否实际使用）、R6（真实报告和签字的索引）、R8（私人路径：用经授权的公开副本，记录每个改动成员的新旧哈希，原件留在私有仓库）。R4 的原冻结包 `PENT24_bundle_b8ab7c1.zip` 我正在做含内部归档的递归隐私扫描，结果补在下面。

## 4. 重放时间和命令：🤝

- GPT 估计全量重放 30–50 小时，稳妥按 48–72 小时安排；
- 生成和重放共用一把排他锁，同一工作目录不能同时跑，不要绕过；
- 收货工具要在重放结束后再跑（16 GB 内存）。
用户步骤改成第 4 版（`local/accept/USER_STEPS_after_run_20261123.md`）：第 4 步先确认锁空闲再二选一启动重放（26a 加 `--task-timeout 7200`，两种都是 all 加 summary）；第 5 步重放结束后再跑收货，按补跑情况给了三条完整命令，不用用户改变量；第 6 步发回两份小报告。所有命令都过了 `bash -n` 语法检查。

时间表因此要更新：生成大约 11-25 结束，加上补跑和 2–3 天的重放，独立复证大约 11-27 到 11-28 结束；PENT-26 签字和发布树定稿在 11-29 到 12-02。离 12-07 的期限还有几天余量，但已经不多。

## 5. 最终签字清单：🤝

GPT 的 `PENT26-FINAL-CHECKLIST.md`（12 项数据和报告核对、32 条数学义务、7 项发布闸门）我会照着准备材料。

## 结论

GPT 这轮的意见除了第 2 节那一处事实以外我都接受，并且都已改完或者已经写进定稿清单。N、F、总定理仍然 NOT SIGNED，等真实全量报告。
