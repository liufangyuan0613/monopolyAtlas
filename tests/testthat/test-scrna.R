# scRNA 层与版本声明测试

test_that("atlas_version declares manuscript binding", {
  v <- atlas_version()
  expect_s3_class(v, "atlasMeta")
  expect_equal(v$n_genes, 74)
  expect_equal(unname(v$gene_split), c(36, 38))
  expect_match(v$atlas_version, "^v20")
})

test_that("atlas_scrna covers 12 datasets with LUAD myeloid FTL anchor", {
  data("atlas_scrna", package = "monopolyAtlas")
  expect_gte(length(unique(atlas_scrna$dataset)), 12)
  row <- atlas_scrna[atlas_scrna$dataset == "LUAD" & atlas_scrna$cell_type == "Myeloid/MAST", ]
  expect_equal(nrow(row), 1)
  expect_equal(row$ftl_top1_rate, 0.573, tolerance = 0.005)
})

test_that("query_scrna filters and plot_scrna_dominance builds", {
  d <- query_scrna("LUAD", min_cells = 100)
  expect_gt(nrow(d), 3)
  expect_true(all(d$dataset == "LUAD"))
  expect_s3_class(plot_scrna_dominance("LUAD"), "ggplot")
})
