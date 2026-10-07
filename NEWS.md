# monopolyAtlas NEWS

## v3.1.0 (2026-10-07) — 口径修正: 技术基因分母政策对齐管线

**行为变更 (breaking for absolute TMI values)**

- `top_share()` / `compute_tmi()`: 技术基因 `^(MT-|MTRNR|RPL|RPS)` 现在只剔出排名、
  **保留在分母** (分母=样本全部基因总和), 与 mono_v2 手稿管线 step1 完全一致。
  v3.0 的分母为剔除后的非技术基因总和, 导致 TMI 系统性偏高
  (真实 TCGA 10,156 样本上比值中位 1.36x, q90 1.81x; 组内 p95 垄断判定一致率 97.7%)。
- 遗留行为可通过 `top_share(..., renormalize = TRUE)` 重现 (仅供对照, 勿用于主分析)。
- `score_monopoly()` 的基因级 top-k recurrence 现在同样剔技术基因出排名 (此前包含)。
- `carrier_dominance()` 新增 `exclude_technical = TRUE` 参数: top1 排名剔技术基因,
  分母保持全基因总和。scRNA 场景建议用户自行再剔 MALAT1/NEAT1 (项目政策)。
- 绝对 TMI 值与 v3.0 不可直接比较;  shipped atlas 数据 (atlas_*) 不受影响
  (均为管线离线计算结果)。

**验证**: 53+4 testthat 断言全绿; 与管线口径的数值锚点 (手算玩具矩阵 + 真实队列锚点
IGHG1 x SKCM = 0.625) 钉死在 test-compute.R / test-atlas.R。
