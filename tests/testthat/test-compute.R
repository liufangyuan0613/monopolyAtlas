# 核心计算锚点测试 (口径钉死, 勿放松)

test_that("top_share computes TMI5 correctly on toy matrix", {
  m <- matrix(c(50, 20, 10, 10, 5, 5, 1, 1, 1, 1, 1, 1), nrow = 6)
  rownames(m) <- paste0("G", 1:6); colnames(m) <- c("S1", "S2")
  # S1 = (50,20,10,10,5,5) sum=100, top5=95 -> 0.95; S2 = 6×1 sum=6, top5=5 -> 5/6
  expect_equal(unname(top_share(m, k = 5, exclude_technical = FALSE)[1]), 0.95, tolerance = 1e-9)
  expect_equal(unname(top_share(m, k = 5, exclude_technical = FALSE)[2]), 5 / 6, tolerance = 1e-9)
})

test_that("technical genes excluded per pipeline policy", {
  m <- matrix(c(90, 10, 5, 5), nrow = 2)
  rownames(m) <- c("MT-CO1", "IGHG1"); colnames(m) <- c("S1", "S2")
  # 剔除 MT- 后只剩 IGHG1: share=1
  expect_equal(unname(top_share(m, k = 5, exclude_technical = TRUE)[1]), 1, tolerance = 1e-9)
})

test_that("monopoly_status uses 95th percentile within cohort", {
  set.seed(1); x <- runif(200)
  st <- monopoly_status(x, pctile_cut = 0.95)
  expect_equal(sum(st), 10)
  expect_true(all(x[st] >= quantile(x, 0.95, type = 8)))
})

test_that("compute_tmi handles sparse matrices", {
  skip_if_not_installed("Matrix")
  m <- Matrix::Matrix(matrix(c(50, 20, 10, 10, 5, 5), nrow = 3), sparse = TRUE)
  rownames(m) <- c("A", "B", "C"); colnames(m) <- c("S1", "S2")
  # S1=(50,20,10) sum=80, top2=70 -> 0.875
  expect_equal(unname(compute_tmi(m, k = 2, exclude_technical = FALSE)[1]),
               0.875, tolerance = 1e-9)
})

test_that("gini/hhi sane bounds", {
  m <- matrix(runif(500 * 20), nrow = 500)
  rownames(m) <- paste0("G", 1:500)
  g <- compute_gini(m); h <- compute_hhi(m)
  expect_true(all(g >= 0 & g <= 1)); expect_true(all(h > 0 & h <= 1))
})
