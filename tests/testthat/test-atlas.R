# Atlas 数据一致性与正文锚点 (与冻结数字对齐)

test_that("atlas_genes is the 74-gene canonical set", {
  data("atlas_genes", package = "monopolyAtlas")
  expect_equal(nrow(atlas_genes), 74)
  expect_equal(sum(atlas_genes$group == "cancer_specific"), 36)
  expect_equal(sum(atlas_genes$group == "tissue_shared"), 38)
  expect_equal(sum(duplicated(atlas_genes$gene)), 0)
})

test_that("atlas_gene_cancer reproduces IGHG1 x SKCM anchor (freq = 0.625)", {
  data("atlas_gene_cancer", package = "monopolyAtlas")
  row <- atlas_gene_cancer[atlas_gene_cancer$gene == "IGHG1" & atlas_gene_cancer$cancer == "SKCM", ]
  expect_equal(nrow(row), 1)
  expect_equal(row$mono_freq, 0.625, tolerance = 1e-6)
})

test_that("atlas_cancer covers 32 cancers with sane rates", {
  data("atlas_cancer", package = "monopolyAtlas")
  expect_gte(nrow(atlas_cancer), 30)
  expect_true(all(atlas_cancer$mono_rate >= 0 & atlas_cancer$mono_rate <= 1))
})

test_that("atlas_demand covers 74 genes x >=50 tissues", {
  data("atlas_demand", package = "monopolyAtlas")
  expect_equal(length(unique(atlas_demand$gene)), 74)
  expect_gte(length(unique(atlas_demand$tissue)), 50)
  expect_true(all(atlas_demand$pctile >= 0 & atlas_demand$pctile <= 1))
})

test_that("query_gene returns annotated entry; unknown gene errors", {
  q <- query_gene("WFDC2")
  expect_equal(q$group, "cancer_specific")
  expect_error(query_gene("NOT_A_GENE_XYZ"), "not found")
})
