test_that("TMI works", {set.seed(42); m <- matrix(runif(1000*20), 1000); rownames(m) <- paste0("G", 1:1000); t <- compute_TMI(m, 5); expect_length(t, 20); expect_true(all(t >= 0 & t <= 1))})
test_that("Gini works", {m <- matrix(abs(rnorm(500*10)), 500); g <- compute_Gini(m); expect_length(g, 10)})
test_that("HHI works", {m <- matrix(runif(500*10), 500); h <- compute_HHI(m); expect_length(h, 10); expect_true(all(h >= 0 & h <= 1))})
test_that("score_monopoly works", {set.seed(42); m <- matrix(runif(500*30), 500); rownames(m) <- paste0("G", 1:500); colnames(m) <- paste0("S", 1:30); r <- score_monopoly(m, 5); expect_s3_class(r, "monopolyResult")})
