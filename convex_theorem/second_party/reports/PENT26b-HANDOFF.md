# PENT-26b — GPT 独立预审交接（2026-11-23）

**总判定：预审完成，CHANGES-REQUIRED（论文少量改文、发布草稿、比较器）；N/F/总定理 NOT SIGNED。无需重启或修改 r3。**

输入 PENT26b_bundle_16640e2.zip SHA256 `b8f95940232265f1761a2af1baca44e598e722c3f58d9337491d0969ce34d7c7`，20/20 清单通过。以包内任务书为合同，读取旧签字树及旧GPT证据进行字节对照，所有运算本地，只读生产文件，不执行生产生成程序；仅执行哈希、扫描、PDF本地渲染/抽取、比较器合成控制、候选报告校核器样本及排他锁测试。无外联、公开、仓库修改。实际用户全量结果尚未回传，本报告不验证新全量数学证书。

## 1 论文 M1–M7
完整逐项见 **PAPER-M1-M7.md**（英文替换句、源码位置、M1严格算术）。diff应用到31e2f15旧稿后与新稿逐字节一致，树内TeX/PDF一致，15页通读。M1共同精确十进制目标通过，生产实际阈值口径需限定；near 的 nearest binary64 略大，far fallback可用 nextafter更大阈值，四边形表独立LB，不可把每项都视作同一个F*。M2/M4/M6及三个建议项PASS；M3条件性延后同意（若实际使用outside，置于N生产证明后的独立复证段）；M5主体修复但附录S(a)角点宽度可能为零，TX前提及证明副本须同步；M7仍有未限定“re-proved”表述，与抽样范围冲突。小范围改文，无须重跑r3。

## 2 v2公开树
**V2-RELEASE-REVIEW.md** 给R1–R10必改及建议。249 vs 116文件，112旧文件逐字节不变，4允许根文件变化，133新增、0删；根248/248、局部42/42、convex129/129、near18/18、far36/36、checker37/37哈希通过。嵌套清单覆盖口径见修正日志。八个r3源文件和九个24c更新文件不变，双数据归档原样。主要缺口：候选及完成措辞冲突、FINAL_EDITS覆盖不足、原冻结ZIP入口无法由当前树直接复现、26a若采用outside需附代码、最终证据尚待实际填入、路径隐私及许可。递归扫描9683对象，除已知Windows路径与13个探针外，还有生产复审和归档中的home账户/工作路径及GPT /mnt/data日志等；不直接修改原证书或v1.1历史签字文件，不未经许可公开私有原包。外部文献依据随包材料，不声称独立外部核验。

## 3 比较器
**COMPARATOR-REVIEW.md**：8个合成控制，7个应拒被接受，1个合法分轮格式被误拒。原比较器不能用于签字。提供 `code/report_guard.py` 候选只读校核器；实际GPT smoke与合成Claude格式夹具九阶段元数据测试通过，默认拒样本和重复记录。未有双方全量真实报告，候选亦不签字/不取代数学重放。最终按固定目标、轮次、source SHA、计划、journal、新鲜摘要、完整性及抽查状态联合核对。

## 4 重放
**REPLAY-TIME-AND-COMMANDS.md**：实际24c样本算术理想17.482h；结合用户机std/E2速度及E3情景，工程估计30–50h，稳健安排48–72h（非上限）。16GB/20进程条件下E2/E3自动10进程。现有driver与两版重放器共用排他锁，不能在同一work同时生成B与重放已完成阶段；实测锁竞争拒绝。不要中断主跑、绕过锁、删锁或复制活动数据；待主跑与所有补跑结束后 all+summary，收货串行。26a需补7200 timeout。完整安全命令及暂停/续跑在专文。

## 5 最终签字清单
**PENT26-FINAL-CHECKLIST.md**：12项收货/重放检查、32条原签字义务、7项最终发布闸门，每条有输入、字段/通过条件，以及一次回传材料清单。总计1,947,824任务（A600,895/B1,346,929）、15,223主分片，正式补跑另外计；尚未实跑。L既有签字保留，N/F及总定理仍待完成真实认证、严谨文案、最终签字。PENT25d仍为历史SIGN-OFF；V2尚不可作正式证明发布。

## 复核索引
logs/input_sha.log、paper_diff.log、paper_extracted.txt、threshold_verify.log、tree_manifests.log、manifest_coverage.log、frozen_link.log、r3_tree.log、privacy_scan.log、comparator_controls.log、guard_smoke.log、replay_time.log、lock_control.log 均为本轮检查原始输出；tree_diff.log保留早期路径前缀误报，后由manifest_coverage.log纠正。code含本轮扫描及比较器测试代码，不修改r3。Python3标准库执行。报告ZIP保留隐私扫描日志仅供私有审查，不自动作为公开树材料。
