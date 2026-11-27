# PENT-24c — GPT 独立重放、路线证明与补跑交接（2026-11-22）

**结论：实现及样本验收 PASS；全量数学验证仍 PARTIAL；Theorem N / F NOT SIGNED。** 本轮先完成 PENT-25c final tree SIGN-OFF，再处理 PENT-24c。零公开、零外联、没有改仓库；PENT-24c 执行的是 GPT 上轮自己的代码及本轮新增代码。未运行 Claude 的生产者程序。

## 0. 冻结身份与交付

- 原样退回的 PENT24b 用户包 SHA256 `691ea281c6dd13dcd5af2369ef768403302dbbdd1a20c1b8a96d77adbf5cb175`。
- 冻结 PENT24 输入 ZIP SHA256 `1aba5f0a693313b515d9b3b0952747cd61ced2b4b5fb2d3c9c5783333ff2f958`。
- r3 code_id `7cf79beefbc81c796969279b7c4d8bf35c7ecb0525e41f2fe3d88d2c629e2d3b`。新版 ZIP 与原样退回 ZIP 的 driver.py、proofs24b.py、lg_runner.py、lg_exact.py、lg_search.py、lg_stdlib.py、label_exact.py、q_fibre_runner.py **8/8 逐字节相同**；另有八文件 SHA 清单。
- 最终 replay_id `0863771d072ea6e0db750948c0152fb21ad9c2011f3a3e71f98d847cdf5a28af`；supplement_id `8ab8cb73fb880cb4a9af0cc895bd40913ab12722f91658d12e36dc9f342eb619`。
- 用户 ZIP `PENT24c_user_20261122.zip`，38 文件，89,100 bytes，SHA256 `2ae625ea555f01e9654beddd553be92e740f4b7e1693d4e7aa3b487ba1f6b7a6`。ZIP CRC、内部清单、八文件字节比对均通过。不含冻结大数据。只复制 README-24c.md 指定的新增文件，绝不能替换正在用的生成器或重新初始化 work。

## 1. 可续跑重放器

新增 `code/replay24c.py`、`replay_stage.sh`、`replay_progress.sh`。固定接受 r3 plan 和实际磁盘八文件 code_id；冻结输入从 ZIP 独立解包并验证；重建 task 清单及 SHA，绑定 chart/j/盒子/目标。支持 near,q_lg,q_other,far_std,far_E2,far_E3,B_std,B_E2,B_E3,A,B,all,structure,summary，按分片追加 fsync journal，记录源文件 SHA、非 PROVED ID/原因、计数、CPU/worker wall/RSS、原生成耗时。重启时同时比较 plan、重放代码及源分片 SHA；任一源变更即重放相应分片。共享 work/run.lock；缺失源列 MISSING。journal 尾部损坏自动截断。STOP / 磁盘 / 内存保护与生成相同。20 worker×512 MiB；E2/E3（含标签）10 worker×1024 MiB，默认保留 2 GiB 内存、20 GiB 磁盘。沙箱 4-vCPU/4GiB 不测试生产 20 worker；测试时降低保护阈值，用户端保持默认。资源保护不是对 OS 总 RSS 的形式化保证；父进程和系统缓存另有开销。

LG 见证要求 d=6（near/far）或 d=4（q），N∈[1,12]，族大小 1 或 d+1；主生成 N 等于冻结任务 N。检查指数逐变量≤31；对 packed 5-bit exponent 增加乘法前的各轴次数上界（避免溢出进位后解码貌似合法）、结果多项式逐项检查。对现有精确多项式内核进行临时包装，仅重放侧使用，不改变八个文件。结构由独立 GPT `structure.py`、`near_robust.py`、`q_audit.py` 和 `replay_all_impl.geometry` 重证；结构记录与日志绑定 SHA。

summary 只汇总已经记录的 journal/snapshot，与原结果 plan/summary 快照一起归档 `work/PENT24b_replay_report.tgz`，并提供逐阶段源分片 SHA 清单、状态计数、全部已扫描非 PROVED ID/原因以及此前失败尝试。仅跑部分阶段时，未扫描分片以索引及缺失计数给出（未扫描 ID 未枚举）；先运行 all 方可取得完整 ID 列表。summary 不重新哈希实物，若新见证到达，需要先重新运行相应阶段。小报告仅证明本地记录的执行快照，**不能从摘要和哈希自行推断远端实物的数学正确性**，需要可信本地执行及保留原见证供抽查。

N/F 计算闸门只在完整非样本数据上开放，REJECT 尝试需审计；表 conditional 依赖 q 全闭合及结构根界链接。F 闸门仅指 FT §7 的远区 disjunction/near handover；整个凸五边形证明另需 N、L 和解析路线认可。本轮不打开任何定理签字闸门。

### 1.1 真实样本及中断恢复

初版新重放器（数学内核相同，后续只调整编排/summary/补跑身份，初版文件随证据附上）对保存的上轮 smoke 样本先跑 all。near/q_lg/q_other/far_std 四阶段各 5 件落盘后，于 far_E2 写 STOP，退出码 75；恢复退出码 0，前四个分片 SHA 缓存命中，无重复数学重放。九阶段合计 45/45 PROVED，零 REJECT、零 conditional；结构亦通过。日志 `smoke_all.log`、`smoke_resume.log`，完整 journal 随包。

最终候选的结构实际重跑：远区 1,260 blocks、196,798 根列、1,724,444 终态盒，0 error；三图卡候选根格 487,424，必须根无遗漏（std 保留额外 1,592 格）；近区 1,344/1,344 roots，38,280 LG、13,118 far、36 inner，582 touch，0 errors/uncertain。q 表 57,093 rows 与根 LB 链接逐条匹配，0 bad。结构耗时 72.709 s（在受限沙箱，非用户实测）。q_audit 中 `requires_more` 是待数学认证的部分，不当成通过。

### 1.2 回归基准：原 P24b 存有实际见证的 1066 项

|阶段|实际/目标|PROVED|CONDITIONAL|MISSING|重放 CPU s|worker wall s|原生成单件时间之和 s|wall ratio|峰值单 worker RSS KiB|
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|near|17/200|17|0|183|17.219|21.549|40.689|53.0%|114496|
|q_lg|200/200|200|0|0|7.114|9.197|14.899|61.7%|111144|
|q_other|200/200|200|0|0|2.386|2.676|13.167|20.3%|111400|
|far_std|17/200|17|0|183|16.447|16.772|29.899|56.1%|170380|
|far_E2|16/200|16|0|184|51.067|61.006|117.877|51.8%|212204|
|far_E3|16/200|16|0|184|82.601|122.536|149.313|82.1%|186072|
|B_std|200/200|131|69|0|0.089|0.089|0.114|78.8%|30688|
|B_E2|200/200|200|0|0|0.324|0.326|1.336|24.4%|70368|
|B_E3|200/200|200|0|0|0.328|0.329|1.629|20.2%|70368|

实际 997 PROVED、69 CONDITIONAL_TABLE、734 MISSING，零 REJECT，与上轮旧重放记录相符；A/B 退出 0。计时为 worker 单件墙钟累计（不是全流程墙钟），分母来自上一轮保存的生成记录，在不同进程争用条件下仅作相对估计。重放精确稀疏多项式仍很重；E3 尾部风险明显，保留原内存并行上限。没有完成可选四类慢 LG 各200 新生成基准，不外推为经过验证的总工期。此表用初版新重放器（id `af3e795c5d73f047f46a691cb45cf81cdf1cf6265ab697b3257a71c8a3fae119`）；最终版同数学核心，最终控制/补跑/结构另运行。真实 benchmark journals 和初版源码随证据，非虚构数字。

### 1.3 篡改与保护

最终代码真实 pristine PASS；13/13 MUST_REJECT：冻结身份、missing witness、非通过记录携带见证、digest、降低目标、几何替换、子盒错误、N=0、N=13、主 N 不匹配、族大小、维数、重复记录。journal 尾部恢复、磁盘/内存/STOP、共享锁冲突拒绝全部实际 PASS，`final_controls.log`。共享锁冲突安全退出 75。初版中断/续跑测试如上；并未在用户实际设备模拟 20 worker 或电源故障。

## 2. LABEL_ROUTES.md

独立路线文档完整包含 (a)–(j)，逐项注明函数、行号、FP/FT/FL 章节，区分原规范与 GPT 新增的相关三角形、条带阈值、纤维断点和截断等。重点：`q_strip/rect_strip` 的 >= 与 any(a>0) 不直接推出严格性，q 表只需 ≥；B T/TX 使用明确抬高的目标 `(T/den)*(1+10^-10)` 给严格结论。零面积角点无需除宽度，依靠各角点的分别凹不等式。文档说明纤维 OD/PS 下界、外凸包、面积下界、cap 与转移正性、内外切分覆盖。classical PS/Makai 用冻结规范的先验定理，不作外部文献核验。表界、N/边图卡交接的条件性明确。等待 Claude 对 B/q_other 每条解析路线独立验收，不能因为样本全过宣布全域结论。

## 3. UNPROVED/ERROR 补跑与整合

`code/supplement24c.py` / `supplement_stage.sh` 只读主结果 UNPROVED/ERROR，独立 `work/supplement_<n>` 与独立 plan，原生成结果不动；LG N→min(N+2,12)，外层切分预算默认128、最多256；标签及 q_other 提升切分预算；可选 stdlib 搜索。保留原 task_sha256。每次配置改变选新 number；同配置续跑不重复已提交任务。主和补跑分别绑定原 plan/新代码及原冻结 task 身份；仅补跑放开 N 等于原提案次数，保留 1..12、几何/目标/族/严格重放。一个合法补跑见证可补回主 UNPROVED/ERROR；无合法见证不标 PASS，conditional 不提前闭合。损坏的尝试记录保留 REJECT，即使另有合法证书也不自动开门。

**实际控制试验（不是新增生产失败数）**：复制真实 q_lg 样本，故意把 `qstd/13340/0` 主记录改为 UNPROVED 且清空见证。独立补跑生成 N=8（原 N=6）严格 LG 见证一件；原 task_sha256 `c1f9a256d373ba4cc019a64b70cc7d2ef97e0e43e866fbb661345a9b63535acf` 不变。补跑 exit 0，新重放 exit 0，合并阶段 5/5 PROVED、0 REJECT。保存主测试副本、补跑 plan+证书、日志和 journal；没有把这个合成控制当作新的生产未证出样本。小报告实际生成 2423 bytes，SHA256 `c80ad92f37188023acca1e078bd43a8d42b5eb21542ca608d3f07929ec0c6a5e`，未扫描阶段明确缺失40项，N/F false。

## 4. 用户执行与签字边界

用户生成中不要覆盖任何 r3 文件，不要强制停止，仅复制新增文件。生成结束或安全停止后，再运行重放。完整说明见 README-24c.md。命令 `bash replay_stage.sh --stage A`，再 B，structure，或 all；`bash replay_progress.sh`；安全暂停 `touch work/STOP`，确认退出后移除 STOP 再续跑；补跑 `bash supplement_stage.sh --number 1 --stage A --budget 128`；最后 `bash replay_stage.sh --stage all` 与 `bash replay_stage.sh --stage summary`。保留完整结果，只回传 tgz、SHA 和日志。不能并行占用同一 work 锁。当前用户40件试跑与生产全量均未回传，不声称其成功。

**PENT-24c: IMPLEMENTED + SAMPLE TESTED; PARTIAL (full computation and analytic review pending).** Theorem N/F NOT SIGNED. 不宣称全域证明完成，不撤销历史 Theorem L 的局部签字。

### 封装前最终版本附加回归

最终 replay_id `0863771d...a28af` 在真实旧 smoke 样本上从空 journal 开始执行 `--stage all`，退出 0；九阶段分别 5/5 PROVED，共 45/45，零缺失、拒绝及条件性结果。结构同次通过。`--stage summary` 实际退出 0，输出 4,966 bytes 的 `PENT24b_replay_report.tgz`，SHA256 `8372f70466660b2681c208c03351832fa19fad6545efdcb205c45249dad07f42`；大小仅为样本结果，不能保证全量报告同样小。两种小报告（全样本、部分阶段）随证据包。所有测试进程退出，未留下后台计算。
