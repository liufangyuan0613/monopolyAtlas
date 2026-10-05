# 数据版本与 scRNA 层

#' Atlas data version declaration (投稿可复现性一句话)
#' @return list: atlas_version, build_date, manuscript, n_genes, anchors, regenerate
#' @export
atlas_version <- function() {
  utils::data("atlas_meta", package = "monopolyAtlas", envir = environment())
  class(atlas_meta) <- c("atlasMeta", "list")
  atlas_meta
}

#' @export
print.atlasMeta <- function(x, ...) {
  cat("monopolyAtlas data", x$atlas_version, "(built", x$build_date, ")\n")
  cat(" manuscript:", x$manuscript, "\n")
  cat(" genes:", x$n_genes, sprintf("(%d specific + %d shared)", x$gene_split[1], x$gene_split[2]), "\n")
  cat(" anchors:", paste(x$anchors, collapse = "; "), "\n")
  cat(" regenerate:", x$regenerate, "\n")
  invisible(x)
}

#' 单细胞垄断摘要 (Fig3 层, 12 数据集 x celltype)
#' @param dataset 可选过滤 (如 "LUAD")
#' @param min_cells 每组最小细胞数 (默认 50, 防小样本噪声)
#' @return data.frame: dataset, cell_type, n_cells, tmi5_median, carrier_top1_rate, ftl_top1_rate, malignant_frac
#' @export
query_scrna <- function(dataset = NULL, min_cells = 50) {
  utils::data("atlas_scrna", package = "monopolyAtlas", envir = environment())
  d <- atlas_scrna[atlas_scrna$n_cells >= min_cells, ]
  if (!is.null(dataset)) d <- d[d$dataset == toupper(dataset), ]
  rownames(d) <- NULL
  d
}

#' scRNA carrier dominance 条形图 (Fig3 视觉语言)
#' @param dataset 目标数据集 (默认 "LUAD")
#' @param top_n 展示组数
#' @return ggplot
#' @export
plot_scrna_dominance <- function(dataset = "LUAD", top_n = 10) {
  d <- query_scrna(dataset, min_cells = 50)
  if (!nrow(d)) stop("no groups for dataset ", dataset)
  d <- d[order(-d$carrier_top1_rate), ][seq_len(min(top_n, nrow(d))), ]
  d$cell_type <- factor(d$cell_type, levels = rev(d$cell_type))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$carrier_top1_rate, y = .data$cell_type)) +
    ggplot2::geom_col(ggplot2::aes(fill = .data$malignant_frac), width = 0.65, alpha = 0.92) +
    ggplot2::scale_fill_gradient(low = "#56B4E9", high = "#D55E00", name = "malignant frac") +
    ggplot2::labs(x = "Carrier top-1 rate (74-gene set)", y = NULL,
                  title = paste0(toupper(dataset), " — carrier dominance by cell type"),
                  subtitle = paste0("n_groups=", nrow(d), " | min_cells=50")) +
    theme_monopoly(base_size = 11)
}
