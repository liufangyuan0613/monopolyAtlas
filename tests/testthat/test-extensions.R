# 拓展模块测试

test_that("export_pdf writes PDF and PNG", {
  skip_on_cran()
  df <- data.frame(cancer = paste0("C", 1:10), baseline = runif(10, 0, 0.5), mono_rate = runif(10))
  p <- plot_headroom(df)
  d <- file.path(tempdir(), "mono_export_test")
  out <- export_pdf(p, "tmp_panel", dir = d, width = 5, height = 4)
  expect_true(file.exists(out["pdf"])); expect_true(file.exists(out["png"]))
  expect_gt(file.size(out["pdf"]), 1000)
  unlink(d, recursive = TRUE)
})

test_that("simulate_share_null gives sane null TMI distribution", {
  sim <- simulate_share_null(n_genes = 5000, n_samples = 50, seed = 1)
  expect_equal(dim(sim$matrix), c(5000, 50))
  expect_true(all(sim$tmi5_null > 0 & sim$tmi5_null < 0.5))  # 零模型 TMI5 应很低
  expect_lt(mean(sim$tmi5_null), 0.05)
})

test_that("carrier_dominance reproduces group-level carrier rates", {
  # 按列构造: col1=(80,5,5) col2=(10,70,5) col3=(10,5,5) col4=(10,5,80)
  m <- matrix(c(80, 5, 5, 10, 70, 5, 10, 5, 5, 10, 5, 80), nrow = 3)
  rownames(m) <- c("WFDC2", "FTL", "OTHER")
  colnames(m) <- paste0("cell", 1:4)
  groups <- c("Malignant", "Malignant", "Myeloid", "Myeloid")
  dom <- carrier_dominance(m, carriers = c("WFDC2", "FTL"), groups = groups)
  expect_equal(dom$carrier_top1_rate[dom$group == "Malignant"], 1)   # WFDC2/FTL 主导
  expect_equal(dom$carrier_top1_rate[dom$group == "Myeloid"], 0.5)   # 一半 WFDC2 一半 OTHER
})

test_that("morans_i_perm detects clustered vs random patterns", {
  set.seed(1)
  coords <- expand.grid(x = 1:12, y = 1:12)
  # 聚集信号: 左侧高
  x_clu <- as.numeric(coords$x <= 6) + rnorm(nrow(coords), 0, 0.05)
  r1 <- morans_i_perm(x_clu, as.matrix(coords), n_perm = 99)
  expect_gt(r1$I, 0.3); expect_lt(r1$p_emp, 0.05)
  # 随机信号
  x_rnd <- rnorm(nrow(coords))
  r2 <- morans_i_perm(x_rnd, as.matrix(coords), n_perm = 99)
  expect_lt(abs(r2$I), 0.2)
})
