# HANDOFF PENT-26a-GPT — 2026-11-23

**总判定：四项预审完成；13 个指定备用盒子的独立 LG 见证 13/13 PASS；N8 PASS；签字总表预审有条件通过；全文草稿 CHANGES-REQUIRED（常数及表述前提）。Theorem N/F/总定理 NOT SIGNED。** 用户正在长跑，本轮不停止、不改变原生成器，不提交或修改任何仓库，零公开、零外联。

## 0. 冻结输入及信任边界
- PENT26a_bundle_31e2f15.zip SHA256 `11a1c55768604e606571db2fd6ebe502174e91e8e7c2173966c4209d13721d14`，24/24 根哈希通过。
- PENT25d_bundle_cc3a52a (1).zip SHA256 `967b2c08cff6255a12581c77db441462617a8518848483e3fa7d2ce10e11d587`；PENT25d 另行 SIGN-OFF (final tree)，见同包 `PENT25d-HANDOFF.md`。
- 额外 `claude_accept_v2 (1).zip` SHA256 `c38759bc7a9456760e867cbe9066c527005d17a4a503d51f78eb324bbd7f1943`，4/4 清单通过；**它是旧 v2，不支持 outside。正式收货要用 PENT26a 内 accept/accept_pent24.py v2.1**。用 v2.1 的真实 walk 例程核对几何规则，不把供方程序当本轮独立 LG 证明。
- 上轮 r3+24c 用户包与 PENT24 冻结数据用于本轮测试，数据 SHA256 `1aba5f0a693313b515d9b3b0952747cd61ced2b4b5fb2d3c9c5783333ff2f958`。原 r3 的八文件 SHA256 8/8，代码 ID `7cf79beefbc81c796969279b7c4d8bf35c7ecb0525e41f2fe3d88d2c629e2d3b`。无改动；原 `replay24c.py` 及 `supplement24c.py` 也保留。报告的用户机器主跑截至23日08:56的 near UNPROVED 13、q_lg/q_other 全过、far_std 68% 是供方状态记录，实际主跑分片没有附上，不是本轮独立全量验证结果。

## 1. near13 备用路线
新增 `code/outside26a.py`、`code/supplement26a.py`、`code/replay26a.py`。只处理指定的 13 个 ID；读完整主跑+第1轮普通补跑（budget 128，N+2）并核对每片冻结任务绑定；仅在二者仍未证出时才重新搜索，其他已证出任务不重跑。使用原 task_sha256，同步原运行锁、状态保护、追加 fsync、续跑。默认4 worker，每 worker 1 GiB、总≤10 worker、2 GiB可用内存及20 GiB可用磁盘保护。全量用户工作目录未交付，未实际执行完整 CLI 扫描用户的百万级主跑；仅独立枚举冻结数据、真实跑目标及重放测试。本包为新增代码，不重新分发八个 r3 文件。

数学：C 用任务书指定的六个 hex binary64 的精确有理值；对闭盒 [l,h]，min||z−C||²=Σ max(l_i−C_i,C_i−h_i,0)²；严格使用 ≥(3/5+2⁻⁴⁰)²。只许 `PENT26a-outside-v1` 的 near 指定补充树出现 `{\"outside\":true}`，主跑/原普通补跑禁止，远区及四边形禁止，整棵至少一个 LG 叶。树的每次中点二分等覆盖父盒，最长相对边优先，先测 outside；第6层起对 retained subbox 自行拟合新 N=min(12,N_frozen+2)=10 并严格重放；LG 叶 cap256、最大深度14。LG 叶数不包含 outside。拟合可用 SciPy，但接受是精确整数/有理数内核而非浮点正性判断。

真实结果：冻结目标与独立枚举 13/13 匹配，测试 `logs/near13_all.jsonl` 保存全部逐盒原始见证，生成/重放 13/13，0 失败，94 LG+118 outside。重放同时严格检查原盒子任务哈希、几何及试探次数。

| ID（near/省略） | LG | outside | 生成秒 | 重放秒 |
|---|---:|---:|---:|---:|
| 023112/23 | 1 | 6 | 1.130 | 0.421 |
| 023121/14 | 1 | 6 | 1.143 | 0.333 |
| 132012/151 | 21 | 11 | 18.017 | 7.149 |
| 132012/26 | 8 | 9 | 8.249 | 2.828 |
| 132012/47 | 15 | 12 | 12.976 | 5.260 |
| 132012/7 | 4 | 8 | 3.043 | 1.413 |
| 132022/157 | 22 | 10 | 22.965 | 7.501 |
| 132022/72 | 1 | 6 | 1.032 | 0.338 |
| 132102/10 | 6 | 12 | 6.147 | 3.821 |
| 133022/30 | 1 | 6 | 1.049 | 0.439 |
| 133022/40 | 4 | 10 | 3.179 | 1.525 |
| 232012/19 | 7 | 12 | 5.282 | 2.681 |
| 233021/2 | 3 | 10 | 4.076 | 1.864 |
| 合计 | 94 | 118 | 88.288 | 35.573 |

时间合计为分项墙钟之和（并行3 worker），整批真实墙钟53.357s；最大单行见证772,227字节。机器配置、磁盘读全部旧结果及 guard 开销不同，不能承诺用户全程同速。

### 对照测试（实际 stdout 的结论；详见 logs/outside_test.log）
| 测试 | 输出 |
|---|---|
| a 132012/7 实际独立拟合（第一次） | PROVED；4 LG +8 outside，3.2147s；独立重放 PROVED 1.3371s |
| b 擦边错误 outside | REJECT invalid outside leaf |
| c 主跑 outside | REJECT outer split structure |
| 额外 普通第1轮补跑 outside | REJECT outer split structure |
| d far_std outside | REJECT outside supplement restricted to near13 |
| d q_lg outside | REJECT outside supplement restricted to near13 |
| 额外 全 outside、零 LG | REJECT witness/status mismatch REJECT |

原供方 26+41 是另一固定细分方案；本轮自适应且 exact outside 较早截停，4+8 并不意味着少铺盖。实际过树证明全部子盒覆盖，严格 LG 通过。与 Claude **v2.1** 交叉核对：同一 hex C，同一有理平方半径 `10880332376564647921188889/30223145490365729367654400`、同一 ≥ 方向、同一 near 且补充许可。独立 411 个节点精确判据一致，13/13 几何接收、全部作 main 被拒绝；见 logs/accept_agreement.log。v2.1 对 near 补充一般允许，本 GPT 额外限制固定13 ID及明确版本，较严格、不扩大接受范围。

### 新代码身份与用户命令
`PENT26a_UPDATE_SHA256SUMS`（与 R3_GENERATOR_SHA256SUMS 格式一致）；new replay ID `4bcad88762ad3853b72997aac34f1b139c68c78d3fb4c33c78dd8d3eccb8f141`；outside supplement ID `beb45023db067a2d2b38fa66746345cc0cdd3f5d5dee99086beb2c57da294e9b`。原 24c 用户工作路径中仅复制新增三文件，先核对 hash。详细先后顺序、进度/暂停/恢复及新重放命令见 `README-26a.md`。第1轮常规补跑完成之后仍有13目标未过时，执行：

```bash
setsid nohup python3 -u code/supplement26a.py --number 2 --jobs 4 --task-timeout 7200 > supplement_2.log 2>&1 < /dev/null &
```

进度 `tail -n 30 supplement_2.log`。如果用了 outside，须用新增 `replay26a.py` 重放、生成 `work/PENT26a_replay_report.tgz`。旧 journal 及原结果不改，新增 `work/replay26a`；主分片 SHA 列表格式保持，只另列各轮补跑 SHA。

## 2. 签字总表预审
见 `SIGNOFF-MATRIX-REVIEW.md` 共32条逐项意见及完整 outside 拼接证明。主体逻辑成立有条件；G3 须统一目标精确语义。N6/N10 明确 near-only 子盒点为简单性锚点、全扇形正性与绕数，F13 需表界所有数值、标签、C 盒、图卡覆盖及来源链均闭合。N9/F9–F11 等待实际主跑+补跑+数学重放；不将用户进度或本次13盒通过写成总定理签字。N 的范围可包含近区简单非凸，但总定理只严格凸（平角可由已知三/四边形结论另行处理）。

## 3. 论文预审
见 `PAPER-REVIEW.md`。14页PDF全文和 TeX；列 M1–M7 必改及可选措辞、英文替换句。关键 M1：草稿把 F* 称 binary64，但 r3 使用 exact decimal rational 18.919134195733914，binary64 比它大约8.3348e−17；统一声明为精确十进制即可，不必修改 r3。Lemma S、N 的 star-shaped、outside 子叶绕数、LG admissibility、标签分母/宽度方向/严格性、纤维正性及生产/独立重放范围需明确。Wang PDF 未附，按任务书用供方摘要，不作原文已核对声明。v2 应定位为严格凸五边形独立计算辅助证明候选（闭合后去候选），不主张优先、不称验证 Wang 的证明。

## 4. KS 独立包络
**PASS**，纯 Fraction，Machin 40项及 J0 30项交错误差，J1>0 证明首根唯一定位。`KS-PROOF.md` 的完整推导和 `code/ks_exact.py` 可重放。实际输出：

- 2.40482555769 < j₀,₁ < 2.40482555770；
- 36.3368290709000122596833 < 2πj₀,₁² < 36.3368290712022115518922；
- 对 κ₀=36.33682907 的安全裕量 ≥9.000122596833×10⁻¹⁰。

N8 数学常数独立闭合（以 Hong–Krahn–Szegő 已知定理为前提）。

## 5. 下一步与签字闸门
1. PENT25d final tree SIGN-OFF，可由用户独立决定 v1.1 发布。
2. 继续原 r3；等九阶段结束再跑第1轮普通补跑。仅仍有指定 near 目标 UNPROVED/ERROR 时安装并运行本次备用；不得覆盖既有证据。若其他目标失败，另外逐条诊断。
3. 如采用 outside，运行新增重放器并用 Claude v2.1 收货，传回带 SHA 的小报告；保留完整原始分片。全量认证必须严格核对 38,280 near、66,138 q_lg、118,962 q_other、377,515 C 及1,346,929 B（总1,947,824），0 MISSING、0 最终未证、0 REJECT、表条件闭合、结构 PASS。
4. 完成论文 M1–M7 的必要调整、标明实际路线目录及机器证据等级。PENT26 最终验收再对签字总表与真实报告检查；当下 N、F、total NOT SIGNED。

## 6. 最终实际回归、断点与原主分片 SHA
新重放器对 PENT24c 原真实 smoke（9 阶段各5件）进行全结构+数学重放，第一次由于人为外层 110 s 限时到 E2 被终止（exit124，保留原始日志）；第二次从持久化断点恢复，结构 cache 命中，前4分片共20件跳过，实际完成剩余25件。**最终45/45 PROVED，0 REJECT，0 MISSING；续跑 exit0，summary exit0**。全结构实际运行 wall74.873s（含 q_audit、near_robust），远区1260块/196798根列/1724444盒、近区1344块/51434叶无错误，57093四边形根LB链接一致。结构仅确认证书身份及根界关系，不替代 q_other/LG 的全量数学认证。小报告实际4963字节，SHA256 `4bfb00e70f35ceb7ef90287a572625eacb45e8f723b3840d2ed4c6692a9082c1`，随证据 `logs/smoke26_report.tgz`。原 24c 九个主分片原始 SHA 重算 9/9；新 summary 生成的九份 `*.shards.sha256` 与旧报告 **9/9 逐字节一致**（logs/main_hashes.log）。无 outside 的旧主分支不受影响。本次新 outside 数学重放和供方 v2.1 一致性已另做真实13盒；尚未接收用户全量结果或执行其完整 CLI。
