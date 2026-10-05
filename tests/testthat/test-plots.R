# 可视化函数 smoke tests (对象可创建, 图层非空)

test_that("plot_headroom builds with permutation annotation", {
  df <- data.frame(cancer = paste0("C", 1:20), baseline = runif(20, 0, 0.6),
                   mono_rate = runif(20, 0, 0.8))
  p <- plot_headroom(df, perm_null = rnorm(1000, 0, 0.2))
  expect_s3_class(p, "ggplot")
  expect_gt(length(p$layers), 1)
})

test_that("plot_tmi_ridge builds", {
  df <- data.frame(sample = paste0("S", 1:200), tmi5 = runif(200),
                   cohort = rep(c("LUAD", "OV", "KIRC"), length.out = 200))
  expect_s3_class(plot_tmi_ridge(df), "ggplot")
})

test_that("plot_state_curve builds", {
  df <- data.frame(cohort = paste0("C", 1:25), baseline = runif(25, 0, 0.5),
                   mono_rate = runif(25))
  expect_s3_class(plot_state_curve(df), "ggplot")
})

test_that("plot_carrier_matrix falls back to ggplot without ComplexHeatmap", {
  df <- data.frame(gene = rep(paste0("G", 1:8), 4), cancer = rep(c("A", "B", "C", "D"), each = 8),
                   mono_freq = runif(32))
  p <- plot_carrier_matrix(df, use_complexheatmap = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_pathway_atlas and dependency boundary build", {
  df <- data.frame(cancer = rep(c("OV", "LUAD"), 5), pathway = rep(paste0("P", 1:5), each = 2),
                   delta = rnorm(10, -0.2), fdr = runif(10, 0, 0.1),
                   delta_excl = rnorm(10, -0.1), fdr_excl = runif(10, 0, 0.2))
  expect_s3_class(plot_pathway_atlas(df), "ggplot")
  d2 <- data.frame(gene = paste0("G", 1:15), group = "carrier", frac_le_neg1 = runif(15, 0, 0.8))
  expect_s3_class(plot_dependency_boundary(d2), "ggplot")
})

test_that("plot_competence with paired diffs", {
  perf <- data.frame(model = c("Demand", "Supply", "Joint"), ROC_AUC = c(0.93, 0.76, 0.93),
                     ROC_CI_low = c(0.85, 0.68, 0.87), ROC_CI_high = c(0.99, 0.84, 0.97))
  diffs <- data.frame(comparison = "Joint - Demand", delta_auc = 0.001, ci_low = -0.04, ci_high = 0.05)
  res <- plot_competence(perf, diffs)
  expect_type(res, "list"); expect_s3_class(res$performance, "ggplot")
})

test_that("plot_gene_card works on shipped atlas", {
  card <- plot_gene_card("WFDC2")
  expect_s3_class(card$demand_plot, "ggplot")
  expect_match(card$info, "cancer_specific")
})
