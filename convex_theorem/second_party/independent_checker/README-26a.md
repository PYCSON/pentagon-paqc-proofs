# PENT26a 备用补跑：用户机器使用说明

**现在不停止 r3，不改 r3，不删除已有结果。** 本包仅新增文件，需已有 PENT24c 用户包（完整 r3 + 24c）、同一冻结数据 ZIP、原环境 SciPy（或可用 --engine stdlib）。等待九阶段主跑完全结束，再运行常规补跑第 1 轮 budget 128（N+2）；如仍有指定 13 个未证出的 near 目标，才启动本次备用补跑。不要同时运行生成、补跑和重放。它们共用 work/run.lock。先保留 work 全部原文件。以下在原用户包根目录执行，`/path/to/evidence` 指本证据 ZIP 解压目录。

```bash
cp /path/to/evidence/code/{outside26a,supplement26a,replay26a}.py code/
cp /path/to/evidence/PENT26a_UPDATE_SHA256SUMS .
sha256sum -c R3_GENERATOR_SHA256SUMS && sha256sum -c PENT26a_UPDATE_SHA256SUMS
```

主跑结束以后，如果没有运行常规补跑：

```bash
setsid nohup bash supplement_stage.sh --number 1 --stage all --budget 128 --task-timeout 7200 > supplement_1.log 2>&1 < /dev/null &
```

等第 1 轮结束后，新的备用一行启动（默认 supplement_2，4 进程，不修改旧结果）：

```bash
setsid nohup python3 -u code/supplement26a.py --number 2 --jobs 4 --task-timeout 7200 > supplement_2.log 2>&1 < /dev/null &
```

进度：`tail -n 30 supplement_2.log`。若冻结数据文件不是默认名 `PENT24_bundle_b8ab7c1.zip`，启动时加 `--data-zip /实际路径`。程序检查主跑及第 1 轮每个分片是否完整，核对冻结 ID，选择仍未完成的指定任务，已证出的不会重跑。只在实际解证可信且完整后才接受。暂停 `touch work/STOP`；等程序安全停止后 `rm work/STOP` 再用原命令续跑。预设内存可用量≥2 GiB、磁盘可用≥20 GiB、每 worker 1 GiB，最多10个 worker。13 个盒子的本轮真实计算+重放 3 worker 墙钟 53.36 秒；用户机器计算段为分钟级参考，遍历百万级旧分片/读盘开销未实测，整条命令耗时不能承诺。

如未使用 outside，可保持原 24c 重放器。如使用了 outside，则最终用下列新增重放器（旧 journal 不动，新目录 work/replay26a；注意先确保旧生成/补跑已经停止）：

```bash
setsid nohup bash -c 'python3 -S -u code/replay26a.py --stage all && python3 -S -u code/replay26a.py --stage summary' > replay26a.log 2>&1 < /dev/null &
```

查询：`python3 -S code/replay26a.py --stage progress`。也支持 --stage near、A、B、q_lg、structure 等分阶段重放。传回 `work/PENT26a_replay_report.tgz` 和 `sha256sum work/PENT26a_replay_report.tgz` 输出；用户本机保存主跑、各补充轮全部原始见证和 journal。Claude 收货必须用 PENT26a 冻结包内 **v2.1** 的 accept_pent24.py，不能用单独上传的旧 v2 ZIP；补充目录列表含 work/supplement_1、work/supplement_2。原始大文件不需要为小报告全部上传，但抽查时必须可取；小报告本身不是完整数学证明。

`code/test26a.py`、`code/test_near13_all.py`、`code/check_accept_agreement.py` 是本轮离线回归工具，包含测试输入位置（见源码），不是用户生产入口。实际生成的 13 个见证为 logs/near13_all.jsonl；它们与冻结盒子身份绑定，可在另一侧复验，不代表已写入你正在跑的生产结果。

## 独立测试的复现命令（原 r3+24c 根目录）
先复制新增三个文件及测试工具到 `code/`，并将证据中 `inputs/near13_boxes.jsonl` 指向实际路径（下例 `E` 和 `DATA` 请设绝对路径）。输出放在另外的测试目录，不要写到生产 work：

```bash
python3 code/test26a.py "$DATA" "$E/inputs/near13_boxes.jsonl" /tmp/near132012_7.jsonl
python3 code/test_near13_all.py "$DATA" "$E/inputs/near13_boxes.jsonl" /tmp/near13_all.jsonl
python3 code/check_accept_agreement.py "$DATA" "$E/inputs/near13_boxes.jsonl" /path/to/PENT26a_31e2f15/accept/accept_pent24.py /tmp/near13_all.jsonl
python3 -S "$E/code/ks_exact.py"
```

这些测试会重新拟合和核对，所用候选搜索依赖原环境；纯数学接受仍是有理数。真实证据行见 logs，不将它直接覆盖/伪装成主跑分片。
