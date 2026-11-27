# PENT-24b 独立检查器：Windows / WSL2 使用说明

**先交 Claude 审查，再由用户启动。** 本包只运行 GPT 独立代码；不运行生产程序、不联网、不修改仓库。不包含冻结数据 ZIP，避免重复下载。此版本是长跑候选，不能保证所有盒子都证出；只有严格复验全部闭合后才可签字。请先看随附报告的实际 smoke/基准边界。

## 1. 放置文件

用 Windows 资源管理器把 `PENT24b_user_20261122.zip` 解压到 D 盘，保证路径是 `D:\PENT24b_user_20261122`。把之前收到的 PENT24 冻结数据 ZIP **复制**进同一个文件夹，文件名改成 `PENT24_bundle_b8ab7c1.zip`（不是 GPT-EVIDENCE）。保留原 ZIP。打开 PowerShell，输入：

```powershell
wsl
```

进入 Ubuntu 后，逐条粘贴：

```bash
cd /mnt/d/PENT24b_user_20261122
sha256sum -c MANIFEST-SHA256.txt
sha256sum PENT24_bundle_b8ab7c1.zip
python3 --version
```

全部 OK 且数据哈希为 `1aba5f0a693313b515d9b3b0952747cd61ced2b4b5fb2d3c9c5783333ff2f958` 才继续。实际沙箱测试 Python 3.13；用户的 3.14 未在本机实测。脚本需要 Linux/WSL2，不在 PowerShell 内直接运行。

## 2. 依赖与先试跑

```bash
python3 -c 'import numpy, scipy; print(numpy.__version__, scipy.__version__)'
```

若上一行成功，不必下载。若失败，正常联网时可以由你自己决定是否运行：

```bash
python3 -m pip install --user numpy scipy
```

两包 wheel 常为数十 MB，按约 50–100 MB 手机流量预算（不是本次实测，架构、Python 版本、依赖会改变；若必须源码编译请取消）。Ubuntu 若提示 externally-managed-environment，不要破坏系统 Python；可用已安装的 venv：

```bash
python3 -m venv .venv
.venv/bin/python -m pip install numpy scipy
```

之后在所有命令前设置 `export PYTHON="$PWD/.venv/bin/python"`。若没有网络、wheel 不可用、pip/venv 未安装，直接选择纯标准库引擎：无需安装，但候选搜索更慢，未证出比例可能不同。精确接受判据本身不依赖 SciPy。

先做各阶段 5 件独立 smoke（单独目录，不计入全量）：

```bash
bash run_all.sh --smoke
bash progress.sh --work smoke_work --sample 5
```

标准库试跑用：

```bash
bash run_all.sh --smoke --engine stdlib --work smoke_stdlib
bash progress.sh --work smoke_stdlib --sample 5
```

烟雾测试是程序线路检查，允许 UNPROVED；ERROR、内存保护/超时、缺阶段都需要保留日志，不可把 smoke 当全量认证。若用户机器 smoke 遇故障，先发日志给 GPT/Claude，不盲目启动长跑。

## 3. 一行启动、查看进度

正常 SciPy 引擎，一行启动全量：

```bash
setsid nohup bash run_all.sh > run_all.log 2>&1 < /dev/null &
```

默认 `--jobs 20`，按 near → q LG → q 非 LG（含全部 fibre/closed）→ far std → far E2 → far E3 → B std/E2/E3 顺序。每阶段 worker 内存限额 512 MiB，far E2/E3 默认 1024 MiB 因此实际最多 10 个进程；请求 20 不会强行超过内存预算。进程组的工作进程地址空间总限额 10 GiB，父进程软限额 768 MiB；系统与文件缓存另有消耗。默认可用内存低于 2 GiB 或余盘低于 20 GiB 时安全停止。WSL2 的 16 GB 设置需用户独立确认。不要同时启动两个全量任务，更不要关闭保护。每个结果分片最大 512 MiB。

纯标准库退路（全量新目录，不可与 scipy 同一个 work 混用）：

```bash
setsid nohup bash run_all.sh --engine stdlib --work work_stdlib > run_stdlib.log 2>&1 < /dev/null &
```

进度（随时另开 PowerShell，先 wsl，再 cd 到上述目录）：

```bash
bash progress.sh
```

查看最后日志：

```bash
tail -n 30 run_all.log
```

STDlib 使用 `bash progress.sh --work work_stdlib`。PROVED=严格证出；CONDITIONAL=依赖四边形表；UNPROVED=本箱未证出，不是反例；ERROR=程序/资源失败。进度更新约 5 秒或在分片完成后。关闭终端通常不影响 nohup，关机、睡眠或 `wsl --shutdown` 会停止。

## 4. 暂停与继续

```bash
touch work/STOP
```

等日志出现 `SAFE STOP`，进度中 active=0 后再关机。不要删除 work。继续时再用同样的命令，原始结果按盒子 ID 与输入哈希自动跳过：

```bash
setsid nohup bash run_all.sh >> run_all.log 2>&1 < /dev/null &
```

标准库改用对应的 work_stdlib 路径和原命令。进程数可调整，例如 `--jobs 12`，退出码 75 表示可恢复保护停止。未完成的尾行会被安全截断；已写完整的结果保留。UNPROVED/ERROR 是已经完成的一次尝试，默认续跑不会重复尝试；补强它们需下一版或新 work 目录，不能混入不同版本。不要手改 JSON 或 SHA。若需缩小进程数，保持 engine、budget 与 work 不变，如 `bash run_all.sh --jobs 8`。若某盒超时可用 `--task-timeout 1800` 重试。程序版本、切分 budget 或引擎改变必须新目录。

## 5. 分两阶段运行（推荐先 A，B 留以后）

A 包含近区、四边形表和远区 C，B 是全部远区标签。若总量太大或想先评估，请启动：

```bash
setsid nohup bash run_all.sh --stage A --partial-pack > run_A.log 2>&1 < /dev/null &
```

完成 A 后会生成 `work/PENT24b_partial_results.tgz`（只表示 A 的计算完成，不表示定理成立）。B 可以之后接着运行：

```bash
setsid nohup bash run_all.sh --stage B > run_B.log 2>&1 < /dev/null &
```

独立运行一个阶段可指定 `--stage near`、`--stage q_lg`、`--stage q_other`、`--stage far_std`、`--stage far_E2`、`--stage far_E3`、`--stage B_std`、`--stage B_E2`、`--stage B_E3`。纯标准库独立 fibre 入口 `python3 q_fibre_runner.py --engine stdlib --work work_stdlib`。同一 work 必须串行使用，程序有锁防止并发。第一次会扫描并核对全部数据清单，稍等后才出现具体阶段进度。

## 6. 跑完交回什么

如果所有阶段的尝试都已记录，自动生成 `work/PENT24b_results.tgz`（即使有 UNPROVED，也会明确保留）；日志最后打印 SHA256。打包前需要额外余盘，保守估计约等于原始结果再加 20 GiB 余量。若打包因磁盘锁停止，释放其他无关空间后同样启动即可；不要删 work 的记录。

```bash
ls -lh work/PENT24b_results.tgz
sha256sum work/PENT24b_results.tgz
```

把 tgz、SHA256、`run_all.log`（或 `run_A.log`、`run_B.log`）作为私有附件交回；不要上传公开仓库。标准库把路径换为 `work_stdlib`。若文件太大，不要用手机热点强行上传，先告知文件大小，让 GPT/Claude 设计分阶段交接。独立审查端先解包到空目录，再使用**自己保留的原版可信代码**和冻结数据执行：

```bash
python3 replay_all.py /path/to/unpacked/results --data-zip PENT24_bundle_b8ab7c1.zip --jobs 4
```

不要从用户回传结果里执行任何程序。默认核对冻结输入哈希、全量盒绑定、逐盒覆盖与数学见证、结构和表界依赖。只有 0 REJECT、0 UNPROVED、0 ERROR、0 MISSING，所需条件性表界全部闭合且解析拼接复核通过，才具备 Theorem N / F 分别签字的条件。样本程序退出成功不等于定理成立。当前交付**不签署 N/F**。

`code/` 内部分旧版脚本保留用于回归与追踪（包括旧 CLI 的 4 进程限制），并非本次长跑入口。用户只用顶层 `run_all.sh` / `progress.sh` / `q_fibre_runner.py`；新 `driver.py` 不调用这些旧 CLI。重放端必须使用审查者自己保留的本次代码，不能从回传 ZIP 导入程序。
