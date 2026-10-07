# 核心计算锚点测试 (口径钉死, 勿放松)

test_that("top_share computes TMI5 correctly on toy matrix", {
  m <- matrix(c(50, 20, 10, 10, 5, 5, 1, 1, 1, 1, 1, 1), nrow = 6)
  rownames(m) <- paste0("G", 1:6); colnames(m) <- c("S1", "S2")
  # S1 = (50,20,10,10,5,5) sum=100, top5=95 -> 0.95; S2 = 6×1 sum=6, top5=5 -> 5/6
  expect_equal(unname(top_share(m, k = 5, exclude_technical = FALSE)[1]), 0.95, tolerance = 1e-9)
  expect_equal(unname(top_share(m, k = 5, exclude_technical = FALSE)[2]), 5 / 6, tolerance = 1e-9)
})

test_that("technical genes excluded from ranking, kept in denominator (pipeline policy)", {
  m <- matrix(c(90, 10, 5, 5), nrow = 2)
  rownames(m) <- c("MT-CO1", "IGHG1"); colnames(m) <- c("S1", "S2")
  # 管线口径: 分母=全基因总和; MT- 只剔出排名
  # S1: cs=100, IGHG1 share=0.1; S2: cs=10, IGHG1 share=0.5
  expect_equal(unname(top_share(m, k = 5, exclude_technical = TRUE)[1]), 0.1, tolerance = 1e-9)
  expect_equal(unname(top_share(m, k = 5, exclude_technical = TRUE)[2]), 0.5, tolerance = 1e-9)
  # 遗留口径 (renormalize=TRUE): 分母=非技术基因总和 -> share=1
  expect_equal(unname(top_share(m, k = 5, exclude_technical = TRUE, renormalize = TRUE)[1]), 1,
               tolerance = 1e-9)
})

test_that("compute_tmi matches manual pipeline computation (denominator = total sum)", {
  # 模拟真实构成: 30% 技术基因 + 5 个非技术基因
  m <- matrix(c(300, 100, 80, 60, 40, 20,
                150, 50, 40, 30, 20, 10), nrow = 6)
  rownames(m) <- c("RPL1", "A", "B", "C", "D", "E"); colnames(m) <- c("S1", "S2")
  # 管线手算: 分母 S1=600, top5 非技术=(100+80+60+40+20)/600=0.5; S2: 分母=300, top5=150/300=0.5
  expect_equal(unname(compute_tmi(m, k = 5)[1]), 0.5, tolerance = 1e-9)
  expect_equal(unname(compute_tmi(m, k = 5)[2]), 0.5, tolerance = 1e-9)
  # k=3: S1 (100+80+60)/600=0.4; S2 (50+40+30)/300=0.4
  expect_equal(unname(compute_tmi(m, k = 3)[1]), 0.4, tolerance = 1e-9)
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
