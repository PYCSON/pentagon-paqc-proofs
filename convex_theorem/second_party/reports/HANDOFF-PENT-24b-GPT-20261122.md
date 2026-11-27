# PENT-24b | 独立长跑程序、补强、基准及可信重放交接

**2026-11-22；总判定：PARTIAL / LONG-RUN CANDIDATE FOR CLAUDE REVIEW。Theorem N/F 均 NOT SIGNED。**

完成了用户要求的主要程序、相同标签样本补强、指定四边形样本、旧见证回归、全九阶段 smoke、保护和续跑测试。未完成四类 LG 的随机 200 件基准，不得称任务书全部 PASS；未运行用户机器全量，不得签字。用户机器包供 Claude 先审，不应未经审查直接启动。

## 0. 身份和信任边界

- PENT25b 先审，报告单独交付，判定 SIGN-OFF（局部 v1.1 发布树）；PENT24b 后做。
- 冻结输入 `PENT24_bundle_b8ab7c1.zip` SHA256 `1aba5f0a693313b515d9b3b0952747cd61ced2b4b5fb2d3c9c5783333ff2f958`；返回 GPT 证据 ZIP SHA256 `8f0cbae220edeb24998df822666d1cdd4a92f430e24b326df6b1dee977ac3e33`，内部 62/62 PASS。读取生产冻结数据与任务书；P24b 数学执行完全基于旧 GPT 独立内核和本轮新代码，未执行/导入生产程序。没有外联、公开或改动仓库。
- 代码版本 `PENT24b-20261122-r3`，核心 code_id `7cf79beefbc81c796969279b7c4d8bf35c7ecb0525e41f2fe3d88d2c629e2d3b`。实际沙箱 Linux x86_64、Python 3.13.5、NumPy 2.3.5、SciPy 1.17.0，cgroup CPU quota 4、内存 4 GiB。用户 Python 3.14 / Windows+WSL2 D 盘未实测，不宣称跨平台已验收。
- 新程序接受条件由精确整数/Fraction 判据决定，SciPy 只用于 Ritz/χ 候选。独立重放用 `python -S`，不需要第三方数值库。不使用用户回传的代码。

### 目标库存（实际冻结输入重枚举）

| 类 | 数量 |
|---|---:|
| near | 38,280 |
| q_lg | 66,138 |
| q_other | 118,962 |
| far_std | 293,846 |
| far_E2 | 47,804 |
| far_E3 | 35,865 |
| B_std | 934,801 |
| B_E2 | 185,823 |
| B_E3 | 226,305 |
| A 小计 | 600,895 |
| B 小计 | 1,346,929 |
| 总计 | 1,947,824 |

q_other 包含全部非 LG parts 与没有 parts 的根格（不仅限 fibre）。far C 合计 377,515。B+far C=1,724,444。没有把结构记录数误作数学认证数。

## 1. 并行和规模：实现完成，用户平台待验

新 `code/driver.py` 是统一长跑入口。请求 `--jobs 20` 默认；普通阶段每 worker RLIMIT_AS 512 MiB，far_E2/E3 默认 1024 MiB，调度按 `min(jobs,10240//worker_mib)` 限流，故后二者最多 10 个。父进程软限额 768 MiB、硬限额 2 GiB（spawn 子进程需能提高软限额）；总正常地址空间预算 10.75 GiB，低于 12 GB。操作系统页缓存等不在该逐进程限额内，额外以 MemAvailable / cgroup 余额作保护。不能声称实测了 20 进程长期峰值；沙箱仅 4 vCPU。`--jobs` 1..20；`--worker-mib` / `--edge-worker-mib` 可调，但保持合计预算，请勿盲降内存限额。旧代码 CLI 保留为回归记录，新统一入口不走其 4 进程限制。

逐阶段独立 run，按冻结输入流式枚举，128 记录分片，每条最多 4 MiB（否则 UNPROVED），故单个结果分片 ≤512 MiB；进度每约 5 秒更新，含 PID、时间、计数。记录以 canonical JSON 绑定 task SHA256、冻结 ZIP、code_id、引擎与 budget；append+fsync，异常尾行仅截断不完整部分，同一 ID 重复或身份冲突拒绝，不混用配置。独占 flock；STOP/Ctrl-C/资源阈值安全止，完整记录可续，未完成任务重试。已记录 UNPROVED/ERROR 不是通过、不会在同配置重复搜索；新的补强需独立目录。单任务超时保护默认 900 秒（包含排队），可放宽同计划的 timeout。归档前额外检查空间，写 SHA256SUMS，再原子改名打印 tgz SHA256。D 盘 Linux 挂载的 fsync/flock 行为需用户 smoke 检查；不能保证实际 I/O 速度。

已实际遇到并修复：fork 继承地址空间导致内存失败，改 spawn；旧父进程硬上限导致 1GiB edge worker 无法 setrlimit，已修正父 soft/hard；2GiB 内存阈值正确主动停机。日志完整保留 `smoke_fork_failed.log`、`smoke_memory_guard.log`、`smoke_hardlimit.log`。最终沙箱为了在 4GiB 内存测试九阶段，显式用 `--mem-min 0.5` 恢复；用户默认仍 2GiB。实际测试日志 `smoke.log`，不是模拟。

## 2. 补强：实现及指定样本 PASS（非全量）

### 2(a) LG 子盒重拟合

`code/proofs24b.py:refit_lg` 在每一子盒重新求一阶仿射 Ritz/χ family，重新 assemble 原始精确 det、PD、PC，检查有理数 μ/严格正性；原盒失败后依据失败多项式变化量选择切分轴（仅启发式），子盒重新拟合，不沿用父盒系数。预算最多 64 终叶、尝试上限 2*budget-1，深度限制，只有两子树完整闭合才接受盒；不把部分成功当 PASS。对于 q 的非单调充分条件，允许严格重算更高 L'≥存档 L。二阶族是可选项，本轮未实现。`lg_bound` 重放绑定原 chart/j/N/L、准确有理 lo/hi，禁止自定义 vertices，重算多项式 digest 与严格正性，看完整 dyadic 切树，不依赖存档 PROVED。

### 2(b) 标签问题、具体修正和同批复验

阈值：原 T 用 `L*(1-cut*tau)>F*`，其中 cut 为盒内 t u 上界，tau 为耳面积占比上界；TX 用最小有效 untouched 阈值，若真实 t u 超阈值，则由剩余耳判据退出，否则用谱界。没有发现原阈值方向写反的证据。旧判据主要损失来自独立区间端点的相关性、四边形谱界较粗、纤维的步长和方向；未发现数学反例。不宣称 E3 X 齐次角面问题全部解决，本轮重点是 T/TX/O。

第一处补强：真实三角形外包关系，令双面积 `A_Q=y+D`，被切耳双面积记 E，`P=A_Q-cut*E`；合格凸成员的剩余面积 ≥P。同一个四边形外包三角形双面积 `A_Q+A B/C`，C>0。直接检查 `TRI*P*C-F*up*(A_Q*C+A*B)>0`，不拆成独立最坏分子分母。图卡四参数的耳不变量是多仿射，乘积为每变量次数≤2，按 tensor Bernstein degree 2 系数分组（除正的 binomial 因子不改变符号），全部严格正即整盒严格正。`TRI` 为严格 4π²/√3 下界，用 Machin π 有理下界与 √3 上界。保留旧严格闭式、依赖表界的 CONDITIONAL_TABLE，增加独立 q-strip、P-strip、P-triangle 等可用谱界替代。

第二处：纤维定向真实凸包。独立从全部顶点角得到整数凸包，投影方向候选为固定方向、边/对角方向及凸包边方向；候选排序不作接受判据。把凸包投影到 [0,1] 后，宽度函数分段仿射，在所有真实凸包折点和均匀 N=16..256 节点切分；每一段用两端真实宽度最大值严格上界 W，势函数下界 π²_lo/W²；用目标稍抬高后的 m，令 v=min(π²_lo/W²,m+4/h²)。在每段计算 `y''=(v-m)y` 的精确向外区间转移（80 位定点，C(z)=Σ z^k/(2k)!，S(z)=Σ z^k/(2k+1)!，|z|≤4，截断 k=20，余项 <2^-80，附加向外 1 ulp）。初值 y(0)=0,y'(0)=1，各后续节点 y>0；振荡段长度<π、非振荡段无二零点，故整区间内 y>0；正解/ground-state transform 给出一维 PS 谱下界，再通过包络凸包与面积下界还原 F。保留 rectangle 和旧严格 fibre 退路。失败才 dyadic split；每一子盒重新构造各路线；所有子盒闭合才通过。实际 route 存入见证；重放以相同精确公式和原盒重算，不能伪造 route 名称。

同批重跑按上一轮 `B_spot_fixed.jsonl` 已完成的 890 块，以旧每块每标签 seed=20261122 的 shuffle 首项重建身份；校对旧失败记录行号。全 1,447 个指定标签都重跑 + 精确重放。旧失败优先只是调度，不改样本集合。

| 标签 | 原样本 | 原 UNPROVED | 新 PROVED | 新 UNPROVED |
|---|---:|---:|---:|---:|
| E2/T | 384 | 118 | 384 | 0 |
| E2/TX | 384 | 13 | 384 | 0 |
| E2/O | 242 | 10 | 242 | 0 |
| E3/T | 78 | 5 | 78 | 0 |
| E3/TX | 78 | 2 | 78 | 0 |
| E3/O | 36 | 1 | 36 | 0 |
| std/O | 245 | 3 | 245 | 0 |
| 总 | 1447 | 152 | 1447 | 0 |

旧失败的实际新路线：E2/T 中 q_fibre 74、P strip(S) 20、q_strip 20、triangle_joint 4，均未再切分；E2/TX 中 q_strip 8、triangle_joint 4、S 1；E3/T 5 和 E3/TX 2 走 q_fibre；std/O 3 走 P_fibre；E2/O 10 中 2 盒用切分（终叶 P_fibre 23、R 3）；E3/O 1 经切分（4 P_fibre 终叶）。故确实主要受到原 fibre/几何谱界松弛，不是阈值方向错误；不能把该同批样本结论推到未指定盒。

### 2(c) q table 非 LG 入口和 1000 随机样本

顶层 `q_fibre_runner.py` 使用同一断点引擎，运行 `q_other` 全部 118,962 项。旧 q_audit 留下 96,817 fibre parts + 2,583 root_fibre 需要更多证明（共 99,400）；从这个集合排除旧已闭式证出的 ID，以 reservoir seed=20261122 均匀选 1000：971 fibre parts、29 root_fibre，**1000/1000 PROVED、0 UNPROVED、0 ERROR**；10 个盒使用切分；终叶路线 q_fibre 924、q_strip 88、empty 1。本轮重放同样 1000/1000。这里“全量闭式入口”是严格的解析/有理纤维/切分的入口，不是说全部 99,400 已运行。所有 1000 记录含任务、见证、秒数在 fastbench/q_fibre1000.jsonl。

### 2(d) 回归

旧 `replay.py` / `lg_exact.py` 以旧见证完整重放，582 inner-touch + 2000 q LG = 2582/2582 PASS，wall 134.522924715s；额外严格绑定 2582 个旧 ID 到冻结任务 chart/j/N/lo/hi/L（各自唯一）。日志 regression_old.log、regression_binding.log。

## 3. 基准：五类 200/200，四类仅部分，外推非完成承诺

执行 `benchmark.py` reservoir seed 20261122，实际使用 SciPy；队列按组轮转，180 秒混合基准达到限时，另快速标签/q 和 qLG 基准补齐。真实 `fastbench` 全部 3247 条完成，92.152s（2 jobs，与其它任务共享 4 vCPU）；qLG 200 条 17.121s（1 job）；慢 LG 组仅下列实际完成数，其余 NOT RUN/未完成（没有计入 UNPROVED 或 PASS）。测试存在并行 CPU 争用，时长为本机实际 per-task wall，非独占的 CPU time。所有具体样本 ID、原始 witness、耗时、RSS 保留在 `bench`, `fastbench`, `qbench`；`bench_stats.json` 去重按 ID，完整组用完整运行结果。P95 nearest-rank。见证字节为 canonical witness，不含 task；磁盘估计用整条实际结果行的平均字节。

| 阶段 | 实际/200 | PROVED/条件/未证出/错误 | 中位 s | P95 s | 见证字节中位/P95 | 理想20进程 h | 配置进程 h | 原始 GiB |
|---|---:|---|---:|---:|---:|---:|---:|---:|
| near | 17/200 | 17/0/0/0 | 1.4667 | 9.2661 | 49702/280199 | 1.273 | 1.273 | 2.761 |
| q_lg | 200/200 | 200/0/0/0 | 0.0567 | 0.1736 | 8250/15301 | 0.068 | 0.068 | 0.547 |
| q_other | 200/200 | 200/0/0/0 | 0.0741 | 0.1359 | 74/79 | 0.109 | 0.109 | 0.030 |
| far_std | 17/200 | 17/0/0/0 | 1.3059 | 4.1881 | 14308/44817 | 7.178 | 7.178 | 4.940 |
| far_E2 | 16/200 | 16/0/0/0 | 4.2000 | 37.9623 | 20604/73792 | 4.891 | 9.783 | 1.176 |
| far_E3 | 16/200 | 16/0/0/0 | 5.3416 | 28.2989 | 14694/60334 | 4.649 | 9.297 | 0.770 |
| B_std | 200/200 | 131/69/0/0 | 0.0002 | 0.0009 | 42/55 | 0.007 | 0.007 | 0.225 |
| B_E2 | 200/200 | 200/0/0/0 | 0.0002 | 0.0491 | 45/77 | 0.017 | 0.017 | 0.044 |
| B_E3 | 200/200 | 200/0/0/0 | 0.0002 | 0.0495 | 45/76 | 0.026 | 0.026 | 0.053 |

按上述已完成样本均值的条件性算术外推：理想 20 worker 18.22 h，实际配置（far E2/E3 10）27.76 h，原始结果约 10.55 GiB。

**重要：上述 27.76 h 不是可信完成时间承诺。** 四类 LG 只有 16–17 个完成观测、长尾被截尾，估计可能严重偏低；完整 smoke 中 E2/3_4/77 单盒耗时 80.95 s、8 子盒（15 次独立拟合），说明长尾明显。用户具体 CPU 性能未知；20 进程的内存与磁盘竞争未实测。B 的纯数学判据耗时外推只有约 0.05h，但 1.35M 条的 IPC、写盘和每16条 fsync 尤其在 WSL `/mnt/d` 上可能主导实际时间；不能将它当成 3 分钟运行承诺。建议先做 A；2–3 倍压力情景 56–84 h，超 3 天完全可能。先对用户机器完整代表块计时，再调整预期；B 设独立第二阶段。磁盘按样本原始约10.55GiB，另外归档空间和 20GiB 安全量，200GB 看似足够，但极端长尾见证、drvfs 和压缩率不确定；保留保护。不要以外推的 0 未证出率签字。q_fibre1000 中位 0.07518s, P95 0.15526s, 无条件/失败见上。

## 4. 用户包和实际 smoke

最终代码实际完成 `bash run_all.sh --smoke --mem-min 0.5`；这是在沙箱先触发真实 2GiB 内存阈值后，用同一 r3 计划重启，之前没有丢失已记录项。9 阶段每类 5，**45/45 PROVED**、0 conditional/error/unproved；实际进度 100% 各类，归档 404,635 bytes，sha256 `f7f53d32e03c5f459f40a28407b94436a868517fde3f72864b3a7d4c84abee2c`。本轮最终 r3 各阶段日志：near 6.545s, q_lg 1.186s, q_other 0.938s, far_std 8.892s, far_E2 95.592s，其他见 raw smoke.log（包含 5/5 全通过）。不将沙箱降低阈值推广到用户机器；README 默认仍 2GiB。

另纯标准库 `python -S` 实测 `qstd/13340/0` 0.3488s 和 `near/023111/127` 13.6735s，两项均生成见证且精确重放 PROVED，未做九阶段标准库全 smoke。SciPy wheels 下载大小 50–100MB 仅流量规划，未联网测量。用户 ZIP 不重复携带冻结数据，用户已有原 ZIP。中文 README 有解压路径、校验、依赖、WSL 命令、启动/停止/重启/分 A/B、进度及结果文件，先经 Claude 审查。

本地启动 `setsid nohup bash run_all.sh > run_all.log 2>&1 < /dev/null &`，进度 `bash progress.sh`；若优先 A `--stage A --partial-pack`，之后 `--stage B`。进度显示 PROVED、CONDITIONAL_TABLE、UNPROVED/ERROR，字面 PASS 不被可信复验器直接信任。已测试原始尾行损坏修复、磁盘阈值模拟（把余盘要求设很大）、实际内存停机与重启。没有实际 20-worker、16GB、长时间 WSL2 环境测试。每个 p24b_user 包核心身份可复现，根清单逐文件 SHA256。

## 5. 可信重放器

`replay_all.py RESULTS_DIR --data-zip ... --jobs 4` 调用保存的 `code/replay_all_impl.py`，用冻结原 ZIP 重新解包，不信任回传的数据、程序、summary 或 PROVED。SHA 清单严格文件集合；每个阶段重新枚举 ID、canonical SHA 以及每片唯一记录，拒绝额外/重复/错位/非法几何，缺项为 MISSING；精确重放 LG 见证与切分树、离散标签 route 和所有子盒；四边形表条件性见证只在 q_LG、q_other 两全量闭合及根表谱界精确绑定后可计入 F，根格结构与近区鲁棒几何也独立重验。样本/部分包允许报告 MISSING，不签字。完整模式拒绝 `--skip-structure`；仅样本可选择跳过。数值候选引擎构造和独立算术接受分离。回传代码绝不执行。重放每批最多 2*jobs 条大见证，避免整片解析导致堆内存过大；重放器默认 4、最多10 worker 每个1024MiB。

旧回归另见 2(d)。新重放焦点实际 2447/2447 PROVED (1447+1000), 17.0965s；来自冻结重新枚举，0 REJECT。篡改对照 8/8 拒绝：错 task SHA、无 witness、UNPROVED 夹带见证、假多项式 digest、降低 L、几何 override、错 child box、重复记录；原例接受。损坏尾行被安全修复、磁盘保护触发。日志 test_controls24b.log。完整结构的本轮远区结果：1,260块、196,798根列、1,724,444终态盒、0错误，三图卡 487,424 候选根格必须保留却漏存为0；std保守多存1592。近区1344块全部匹配、38,280 LG、13,118 far、36 inner、582内球面触片、0错误。进一步的独立 replay smoke、bench 的实际终态见本报告末尾补记。

## 6. 签字口径与剩余清单

- Theorem L：历史 FULL-PASS；P25b SIGN-OFF 与本轮 N/F 无相互替代关系。
- N：只有 near 38,280 全部严格通过，近区结构、鲁棒中心以及局部 L 拼接等数学义务全闭合，方可单独签。
- F：只有 far C 377,515、q LG 66,138、q_other 118,962、B 1,346,929 全量重放完成且 0 未证出/0 error/0 missing/0 reject，CONDITIONAL_TABLE 的依赖表界全闭合，几何/覆盖及局部—近区—远区解析拼接均通过，方可签。严格凸五边形范围；不扩大至任意非凸简单五边形。自动打印 computational gate 不是自动签字。
- **剩余**：四类 LG 200 随机基准未达目标；用户机器实测未做；剩余全域数学证书尚未运行；尤其不能以该批样本推断 E3/X 齐次角面余项等其它标签全部闭合。建议 Claude 先审 proof/程序，用户先 smoke，先 A、报告实际进度与未证出，后 B；最终向 GPT 交回完整结果重放。不存在本轮后台未完成全量任务。

## 7. 最终补记：全部本轮样本重放与封装

- `replay_all.py smoke_work/sample_5 --jobs 2`：从冻结 ZIP 全新提取并重建完整结构，近区鲁棒中心、四边形根覆盖及 57,093/57,093 根 LB 链接全部通过；九类 45/45 数学见证重放 PASS、0 MISSING、0 REJECT、0 CONDITIONAL。实际 wall 123.4912 sec（包含结构检查、冻结提取、重新枚举和证明重放），退出 0。
- `replay_all.py ../bench_replay_results --jobs 2 --skip-structure`：基准 200 计划每类，下列实际 1,066 条（997 PROVED、69 CONDITIONAL_TABLE）全部状态与严格重放一致，0 REJECT；明确 734 MISSING（near 183、std far 183、E2 far 184、E3 far 184）。限时未计算的部分不作真值判定。实际 87.7616 sec、退出 0（sample-only），完整结构已在上一项实际验证。表界依赖、N/F computational gate 均 false；此为正确的拒绝自动签字行为。
- `replay_bench_focus.py`：2,447/2,447（1,447 标签与 1,000 q）实际重放，0 REJECT；17.0965 sec。
- 用户 ZIP 新空目录解包，逐文件清单验证 27/27 OK，启动入口及重放入口 `--help` 均成功。用户 ZIP 61,014 bytes，SHA256 `691ea281c6dd13dcd5af2369ef768403302dbbdd1a20c1b8a96d77adbf5cb175`。无重复冻结数据。证据包包含代码、上述真实 logs、原始 smoke 归档、bench 日志与可重放样本、报告及 SHA256SUMS；不含任何生产程序执行产物。全部本轮启动的计算均已结束。

**最终：长跑候选可交 Claude 静态/技术审查，用户平台与剩余 734 个基准项目仍未完成；全域 N/F not signed。**
