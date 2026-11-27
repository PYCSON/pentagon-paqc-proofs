# PENT-24c：在正在使用的 r3 长跑目录中增加重放与补跑功能

**不要替换正在运行的生成器、不要删除 work 或重新初始化 plan。** 这次只新增下文列出的运行、文档、测试和校验文件。生成的八个文件保持 r3 逐字节不变，`code_id` 必须为 `7cf79beefbc81c796969279b7c4d8bf35c7ecb0525e41f2fe3d88d2c629e2d3b`。新包不含冻结数据，直接使用原来的 PENT24 数据 ZIP。先交 Claude 审查路线文档。零公开、零外联、不要把任何用户结果放入公开仓库。

## 1. 更新（推荐复制新文件，不动原目录）

Windows 解压新 ZIP 到另一个位置，复制以下相对路径到**现有** `/mnt/d/PENT24b_user_20261122`，保持 `code/` 路径不变；不要覆盖原生成文件、`run_all.sh` 或 `work`：

- `code/replay24c.py`
- `code/supplement24c.py`
- `code/test_controls24c.py`（可选测试）
- `replay_stage.sh`、`replay_progress.sh`、`supplement_stage.sh`
- `LABEL_ROUTES.md`、`README-24c.md`、`R3_GENERATOR_SHA256SUMS`、`UPDATE24c.sha256`

WSL 进入现有目录，验证：

```bash
cd /mnt/d/PENT24b_user_20261122
sha256sum -c R3_GENERATOR_SHA256SUMS
sha256sum -c UPDATE24c.sha256
python3 -S -c 'import sys;sys.path.insert(0,"code");import driver;print(driver.code_id())'
```

八项必须全 OK，最后一行必须精确匹配上述 r3 code_id。旧包的 `MANIFEST-SHA256.txt` 仍可用于旧文件；新增文件另有清单，新 ZIP 自带全文件清单 `MANIFEST-SHA256.txt`。**不需要、更不应该为使用重放器而重跑 r3 已完成的生成。**

## 2. 停止生成后才启动重放

同一个 `work/run.lock` 是排他锁：生成、重放、补跑不能同时运行。可以先让生成阶段完成，或者 `touch work/STOP`，等待 `SAFE STOP` 且所有生成 worker 退出，再自己移除 STOP。不要在生成仍活跃时 `rm work/STOP`。默认按 20 进程请求，near/q/std/B_std 每 worker 512 MiB，E2/E3（含标签）每 worker 1024 MiB、最多 10 个。地址空间总限额不超过 10 GiB（另有父进程、系统缓存开销）。磁盘至少 20 GiB、可用内存至少 2 GiB，超限保护退出码 75。默认不要降低保护阈值；沙箱测试才使用小阈值适应 4GiB 限制。新重放不需 SciPy；`-S` 纯标准库。每个 worker 重放已有见证，只做严格整数/有理数接受判定，不求解候选函数。

生成退出后，按需：

```bash
rm -f work/STOP
setsid nohup bash replay_stage.sh --stage A > replay_A.log 2>&1 < /dev/null &
```

只重放 near: `bash replay_stage.sh --stage near`。允许 q_lg、q_other、far_std、far_E2、far_E3、B_std、B_E2、B_E3、A、B、all、structure。`all` 包含结构；`A`、`B` 不自动跑结构，必须单独 `bash replay_stage.sh --stage structure`。推荐先 A，再 B，再 structure，或者最后执行 `--stage all`（已验分片按 SHA256 命中，只核对而不重新计算；补齐未扫描分片，并重算变化分片）。源结果未生成的分片会列为 MISSING，不偷填 PASS。须运行 `--stage all` 才能在报告中列出每一项缺失 ID；只跑部分阶段时，未扫描的分片以范围/编号和 MISSING 总数呈现，报告会明确标注不完整。

若要试跑现有 `smoke_work/sample_5`，不要与全量共用 replay 目录：

```bash
bash replay_stage.sh --work smoke_work --results smoke_work/sample_5 --stage all --jobs 4
```

40件试跑例：将两个路径换为试跑的实际 work 和其 `sample_40` 子目录。`--work` 必须是实际 r3 生成目录，`--results` 必须是该目录内原名 `results` 或 `sample_N`；这样无法绕开原锁。脚本在每次启动从冻结 ZIP 重新提取到临时目录，避免缓存清单被替换；需要约几百 MiB 临时磁盘空间，不运行生产者代码。

## 3. 进度、暂停、续跑

```bash
bash replay_progress.sh
```

输出每阶段当前已有源结果的重放分片数/总分片数、扫描数、六状态、REJECT 尝试数、最新更新时间。暂停 `touch work/STOP`，看到 `SAFE STOP` 再关机。安全退出只丢失尚未提交的当前分片；完整分片日志 `work/replay/<stage>.jsonl` 逐条 fsync。恢复前移除 STOP，再用原命令。若末行断裂自动截断；任一源分片 SHA 变化会重放该分片。如果 replay 代码或 plan 身份变了，不重用旧 journal，移动旧 `work/replay` 为归档目录，重新开始。不要自行编辑日志。重放器记录每次源的 SHA256；若用户生成了更多结果，必须重跑受影响阶段，最终推荐 `--stage all`，才能取得更新的检查快照。

## 4. 补跑（只在主生成已结束，保留原 r3 输出）

```bash
setsid nohup bash supplement_stage.sh --number 1 --stage A --budget 128 > supplement_1.log 2>&1 < /dev/null &
```

仅选主结果的 UNPROVED/ERROR，不选缺失盒子（缺失要继续 r3 生成），不重跑已经有通过见证的主任务。LG 试探次数 N→min(N+2,12)，切分预算默认 128，比 r3 64 更大；B/q_other 加切分预算（未更改数学判据）；`--engine stdlib` 可离线，不需要 SciPy，但搜索效率不同。补跑记录在 `work/supplement_1/`，按原分片序号，仍使用原冻结任务 SHA256，独立 plan 绑定原 plan 和新代码。默认 128，最多 256；提高预算或改变引擎时选 `--number 2` 等新目录。暂停/续跑/锁同上，完成的补跑行不会再次尝试。旧主结果不修改。若某个任务主结果失败而补跑见证有效，就算合并 PROVED；有条件的表界只能 CONDITIONAL_TABLE。重放仍强制 chart、j、父子盒精确位置、L≥冻结目标、N=1..12、族大小 1 或 d+1。仅补跑放开与原 N 的相等，试探函数的次数不是定理前提。出现伪造/损坏记录单独标 REJECT，即使另一个合法见证通过也保留该尝试，阻止自动计算闸门直到问题已审计。

补跑后：`bash replay_stage.sh --stage all`（会发现新补跑目录并核查），全部源仍按 SHA 绑定。

## 5. 只回传小包（不用传约10 GiB 见证）

完成所需重放阶段后，最后执行：

```bash
bash replay_stage.sh --stage summary
ls -lh work/PENT24b_replay_report.tgz
sha256sum work/PENT24b_replay_report.tgz
```

这个报告是已记录 SHA 的**历史快照**：summary 不重新扫真实见证；若源结果有变，先执行 `--stage all` 再 summary。包中只有 replay journals、结构日志、plan/summary 原样快照、每个阶段源分片 SHA256 清单，以及所有已扫描任务的非 PROVED ID/原因与旧尝试记录。未扫描阶段列出所有未扫描分片序号和 MISSING 数量；不能说已完成全量。完整阶段都扫描后，每个非 PROVED ID 都逐条列出。报告可能随失败数增大，但不会带 10 GiB 试探系数。手机热点只需回传这个 tgz、SHA256 和必要日志。请保留本地真实 results/supplements 和冻结 ZIP，供之后抽查。哈希与计数可以远程审计一致性，**小报告不能替代实际见证的远端数学重放**；计算在用户机以已验可信代码执行，需要信任该本地执行环境。样本绝不开放 N/F 闸门。Theorem F 闸门对应 FT §7 的严格凸成员结论/near 移交，最终全域结论仍需 Theorem N 和 L 以及全部解析路线经第二方签字。
