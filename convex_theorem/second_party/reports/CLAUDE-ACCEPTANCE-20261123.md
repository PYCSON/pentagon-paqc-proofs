# Claude 对 GPT PENT-25d + PENT-26a 的验收（2026-11-23 晚）

输入：`GPT-REPORT-PENT25d-26a-20261123.md`、`PENT25d-HANDOFF.md`、`PENT26a-HANDOFF.md`、`PENT26a-GPT-EVIDENCE.zip`（SHA256 `3ddb467a3140537dc01bf0cebbbd44a2b124491e80ab75acfd3b7ab9620a50f1`，和 GPT 报告的一致；内部清单 40 / 40 通过）。GPT 的程序只读不跑；它的数据（见证、日志、报告包）用我自己的工具检查。

## 1. PENT-25d：🤝 接受 SIGN-OFF (final tree)

- 最小差异 134 / 134 行一致，9 个文件变化，0 增删，证书和生成代码不变；Zenodo JSON 可解析；PDF 8 页无 DRAFT。和我的 `CHANGES-vs-PENT25c.diff` 一致。
- 签字的树存档为 `release/v1.1-final-signed/PENT25d_bundle_cc3a52a.zip`（SHA256 `967b2c08…e11d587`），取代 PENT25c；`RELEASE_STEPS_v1.1.md` 已指向它。
- 本地试做公开提交（克隆公开仓库、以 Puyi Sun 身份做一个本地提交，**没有推送**）：在 v1.0 之上 49 个文件变化，v1.0 原有文件里只改了 `.zenodo.json`、`CITATION.cff`、`MANIFEST-SHA256.txt`、`README.md`、`tools/make_manifest.py`，树里 116 个文件。和 PENT-25c 时一样。
- GPT 没有拿到 Wang 的 PDF，如实说明了只读了我的阅读记录，没有当作原文通读。🤝
- GPT 建议的 v2 定位句（"independent computer-assisted proof candidate … does not rely on Wang's claimed proof … no priority claim"，闭合后去掉 candidate，不写 "verification of Wang's proof"）和我已经写进论文和 README 的一致。🤝

## 2. PENT-26a 任务 1：近区 13 个盒子

GPT 在自己的环境里用 outside 叶子把 13 个全部重新证出：94 个 LG 叶子（N = 10）加 118 个 outside 叶子，生成加重放 53 秒（3 进程）。

**我的独立检查**（`accept_pent24.py` v2.1，只读 GPT 的见证数据 `logs/near13_all.jsonl`）：
- 用这 13 行做一个补跑目录 `supplement_2`（主结果里这 13 个记为 UNPROVED），跑 near 阶段：**13 个全部 PROVED，0 拒绝**；叶子计数 94 个 LG、118 个 outside，和 GPT 的一样；每个 outside 叶子都用同一个 C 和 (3/5 + 2⁻⁴⁰)² 精确复核通过；树的二分、叶子盒子、目标值、族的大小全部绑定通过。
- 94 个 LG 叶子全部做了逐点精确检查（`lg_point_check.py`：从顶点公式重建多边形，u 和 χ 写成物理坐标下的多项式，精确积分，再用 Lemma LG），每个叶子 4 个点（一个角加 3 个随机点）：见下面第 6 节的结果。
- outside 判据：GPT 的 `outside26a.py`（读过）用同样的六个十六进制 C、同样的 `OUT2 = (3/5 + 2⁻⁴⁰)²`、同样的 `≥`，只允许那 13 个 ID、只在 near；比我的工具更严，不扩大接受范围。🤝
- 读过 `supplement26a.py`：只写新的 `work/supplement_2`，先核对主结果和第 1 轮补跑每个分片都完整、绑定冻结任务，只对第 1 轮之后仍未证出的这 13 个里的目标做；`replay26a.py` 只写 `work/replay26a`，报告格式（`replay/replay_summary.json`、`replay/<stage>.shards.sha256`）和我的比较器 `compare_reports.py` 对得上。
- 这些见证还不在用户机器的正式结果里；按 GPT 的设计，要在主跑和第 1 轮补跑之后，在用户机器上用 `supplement26a.py --number 2` 重新生成。用户步骤已经写进 `local/accept/USER_STEPS_after_run_20261123.md` 第 3b 步。

**比较器的实测**：GPT 这次附了 9 个阶段的真实小样本重放报告（`logs/smoke26_report.tgz`，每阶段 5 个任务）。我的工具此前对 PENT-24c 证据里的同一批样本跑过一次；两边的分片清单在 near、q_other、far_std、far_E2、far_E3、B_std、B_E2、B_E3 八个阶段**逐字节相同**。q_lg 不同，是因为我那次用的样本带着补跑目录（另一批 q_lg 样本）。这说明比较器的格式和两边的分片哈希约定完全一致。

## 3. PENT-26a 任务 2：签字总表预审 🤝

32 条：大多数同意；G1、G3、N6、N10、F13 有条件；N9、F9–F11 等全量结果。GPT 写的"outside 叶子加进来后 N10 的严格证明"（用 z 本身作简单成员的锚点、闭盒覆盖、outside 子叶不交目标集）我读过，和我在签字总表 N9 下的论证一致，而且更完整。已记进签字总表。

## 4. PENT-26a 任务 3：论文 CHANGES-REQUIRED 🤝（已改）

- **M1（F* 的精确含义）**：GPT 说得对，我把 F* 写成了 "a binary64 number"，而 GPT 的 r3 用的是精确十进制；binary64 比十进制大 715956287/8589934592000000000000000 ≈ 8.33 × 10⁻¹⁷（我用 Fraction 重算过）。我查了生产程序：`near_lg.py`、`cutcell.py`、`qtable.py`、`far_block.py` 在传给检查函数时用的是 binary64（更大），所以生产证书蕴含十进制的结论。改为：F* 是精确十进制；生产侧用更大阈值的证书蕴含同一不等式；binary64 输入按其编码的有理数理解。NEAR_THEOREM.md、FAR_THEOREM.md 加了同样的注。
- **M2**：Lemma S 的证明补上"两个不同顶点"和"相邻边交于另一点必然重叠"两步，用 z 本身作锚点；定理 N 的陈述加上"关于一个内点星形"。
- **M4**：单个多边形的 LG 段落加上"简单多边形、扇形行列式全正、绕数为 1"的前提。
- **M5**：标签 T（LB > F*、τ 的上界为正）、TX（严格分离 θ_T > sup g）、X（分母为正）、S(a)（W_n 和 |n|² 的写法与附录一致，并要求 π²τ > F*）、R（C_m > 0 且 C_{m+2} > 0）、N（分母严格为正、向外舍入、C_p 是对应标号的那一对）。
- **M6**：纤维引理补上横截线零延拓、C¹ 拼接无边界项；计算段落说明 φ 在每个开区间上都被证明为正（双曲或仿射段至多一个零点，三角段还要短于半个周期；我核对了 `fibre.py` 确实这样做），并指向 GPT 的 LABEL_ROUTES。
- **M7**：写明生产侧重放的范围（6 个 std 块的 1,538 个盒子只抽查过）和"第二方全量独立重证是另一回事"。全量结果的数字等最后再填。
- 建议项：主定理写成"简单、正面积的凸五边形，等号当且仅当相似于正五边形"；"three gaps" 改成 "gaps"；非凸那句改成"数值实验表明……这里不对一般非凸五边形作任何断言"。都改了。
- **M3（outside 叶子的论证）**：只有正式补跑用到 outside 叶子时才加；记在 `release/tools_v2/FINAL_EDITS_v2.md`。

论文现在 15 页，编译无未定义引用。

## 5. PENT-26a 任务 4：KS 常数 🤝

GPT 用纯有理数（Machin 公式求 π，J₀ 交错级数加严格余项，J₁ > 0 证明首根唯一）得到 2.40482555769 < j₀,₁ < 2.40482555770，36.3368290709000… < 2πj₀,₁² < 36.3368290712022…，KS = 36.33682907 的余量 ≥ 9.0 × 10⁻¹⁰。和我 11-22 的球算术结果（2πj₀,₁² = 36.336829071074…）一致，两条路线互相独立。N8 闭合。

## 6. 逐点检查的结果

`accept_pent24.py` v2.1，`--sample 94 --points 4 --seed 20261123`（报告 `claude_accept_near13_report.txt`）：
- 94 个 LG 叶子全部抽到（94 / 94），每个 4 个点（一个角加 3 个随机点），共 376 次精确逐点检查：**失败 0**；逐点下界相对目标 F* 的最小余量 **6.0 %**；
- 行的 SHA-256 `c0380f13…16b027b0`，补跑分片清单 SHA-256 `b86812ec…e922a9b6`；
- 结论：`ACCEPT`（near 阶段这 13 个目标的覆盖和绑定完整）。
盒子上多项式的正性（PD/PC）是 GPT 重放器的工作，这里不重复；逐点检查从另一条路线（物理坐标下的多项式、精确积分、Lemma LG）确认了每个叶子的见证在所查的点上成立。

## 结论

- v1.1：双方签字，**可以发布**，等用户决定日期、是否附勘误、由谁推送。
- v2：13 个盒子有了独立证出的路线；论文的必改项除 M3 外都已改；N、F、总定理仍等全量结果和 PENT-26。
