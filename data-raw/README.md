# data-raw 版本化说明

atlas 数据与手稿版本绑定（投稿可复现性声明用）：

- **当前版本**: v2026.10（随 atlas_meta.rda 声明，`atlas_version()` 可读）
- **源数据**: E:\cancer_mg\mono_v2\data\upgrade\（MANIFEST 登记的数据包）+ step1/step2 管线表
- **重建**: `Rscript data-raw/build_atlas_data.R`（内置硬锚点断言：
  IGHG1×SKCM mono_freq=0.625、LUAD Myeloid/MAST FTL top1=0.573、74=36+38、74×32、74×52；
  任何口径漂移会在构建期直接 stop）
- **版本升级纪律**: 手稿口径冻结后更新 atlas_version 字段（vYYYY.MM），
  升数据必升版本号，build_test_log 同步入库。
