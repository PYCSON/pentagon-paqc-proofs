# run3：严格几何 T1 格子扫描（本地长跑）

- 做什么：对 run2 的同一批 480 个中心，测"严格版 T1"（分支定界，球算术）在每个中心最大能放多大的格子 ρ。结果用来算**严格全覆盖**真正要多少格子、多少 CPU。
- 性质：几何部分（凸性、扇形雅可比、方差界、分支定界）全是严格的球算术；**角点特征值用浮点**（再扣掉 2e-5 的相对余量，模拟严格角点证书的损失）。所以这仍然**不是证明**，只用来估成本。
- 每个中心做 3–6 次测试，每次：64 个浮点角点（每个约 0.5 秒）+ 快速浮点预检 + 分支定界（几秒到十分钟，上限 6,000 个盒子；ε 取 Efron–Stein 与精确 e 中较小的上界）。预计每个中心 5–30 分钟。
- 8 个工作进程，nice 19（只用空闲 CPU），每个进程约 150–300 MB 内存。预计 8–20 小时跑完（和 RSN 抢 CPU 时更久）；可以随时停、随时续。

## 依赖
- Python 3.14、numpy、scipy（run2 已装）；
- 新增 **python-flint 0.9.0**（球算术库，PyPI 有 Python 3.14 的现成安装包）。

## 断点续跑
结果逐行追加到 `out_run3/results.jsonl`。停了以后用同一行命令重新启动，做完的中心自动跳过。

## 安全停机
可用内存低于 1.5 GB 或磁盘剩余低于 2 GB 时，写 `out_run3/STOPPED_guard.txt` 并安全退出。

## 跑完
生成 `out_run3/DONE.txt`、`summary.json`、`SHA256SUMS`。把 `results.jsonl`、`summary.json`、`SHA256SUMS` 三个文件发给 Claude。

## 文件
- `run_t1sweep.py`：扫描主程序；
- `rig_t1.py`、`rig_paqc.py`、`rig_box.py`、`mps.py`：与仓库 `probes/PENT/scripts/` 和 `MEAS-01/claude/scripts/C456/` 中同名文件相同（打包时复制）；
- `fstar_exact.json`：正五边形 λ₁·面积 的严格上界；
- `centres_run2.jsonl`：run2 的结果，用作中心列表和起始 ρ；
- `params.json`：参数。
