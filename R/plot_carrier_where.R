# CARRIER / WHERE 层 (Fig2/Fig4 视觉语言)

#' Carrier monopoly matrix: genes x cancers (ComplexHeatmap 优先, 否则 ggplot 热图)
#' @param mat data.frame long: gene, cancer, mono_freq; 或 wide matrix (genes x cancers)
#' @param cluster 行列聚类 (默认 TRUE)
#' @param use_complexheatmap 强制使用 ComplexHeatmap (缺省=可用即用)
#' @return ComplexHeatmap::Heatmap 或 ggplot
#' @export
plot_carrier_matrix <- function(mat, cluster = TRUE, use_complexheatmap = NULL) {
  if (is.data.frame(mat)) {
    stopifnot(all(c("gene", "cancer", "mono_freq") %in% names(mat)))
    w <- reshape(mat[, c("gene", "cancer", "mono_freq")],
                 idvar = "gene", timevar = "cancer", direction = "wide")
    rownames(w) <- w$gene; w$gene <- NULL
    mat <- as.matrix(w)
    colnames(mat) <- sub("^mono_freq\\.", "", colnames(mat))
  }
  mat[is.na(mat)] <- 0
  use_ch <- if (is.null(use_complexheatmap)) requireNamespace("ComplexHeatmap", quietly = TRUE) else use_complexheatmap
  if (use_ch) {
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE))
      stop("ComplexHeatmap required: BiocManager::install('ComplexHeatmap')")
    col_fun <- grDevices::colorRamp2(c(0, 0.25, 0.5, 0.75, 1),
                                     c("#F7F7F7", "#56B4E9", "#0072B2", "#E69F00", "#D55E00"))
    return(ComplexHeatmap::Heatmap(
      mat, name = "mono_freq", col = col_fun,
      cluster_rows = cluster, cluster_columns = cluster,
      row_names_gp = grid::gpar(fontsize = 7), column_names_gp = grid::gpar(fontsize = 9),
      column_title = "Carrier monopoly frequency (gene x cancer)"))
  }
  df <- as.data.frame(as.table(mat))
  names(df) <- c("gene", "cancer", "mono_freq")
  ggplot2::ggplot(df, ggplot2::aes(x = .data$cancer, y = .data$gene, fill = .data$mono_freq)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradientn(colors = c("#F7F7F7", "#56B4E9", "#0072B2", "#E69F00", "#D55E00"),
                                  limits = c(0, 1)) +
    ggplot2::labs(x = NULL, y = NULL, title = "Carrier monopoly frequency") +
    theme_monopoly(base_size = 9) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   axis.text.y = ggplot2::element_text(size = 5))
}

#' Gene card: demand tissue profile + cancer occurrence + structure (单基因一页)
#' @param gene_symbol gene symbol
#' @param demand data.frame: gene, tissue, pctile (默认=内置 atlas_demand)
#' @param gc data.frame: gene, cancer, mono_freq (默认=内置 atlas_gene_cancer)
#' @param genes_meta data.frame: gene 注释 (默认=内置 atlas_genes)
#' @param top_n 展示的组织/癌种数
#' @return patchwork-style list of ggplots (或合并 grob 若 patchwork 可用)
#' @export
plot_gene_card <- function(gene_symbol, demand = NULL, gc = NULL, genes_meta = NULL, top_n = 10) {
  if (is.null(demand)) { utils::data("atlas_demand", package = "monopolyAtlas", envir = environment()); demand <- atlas_demand }
  if (is.null(gc)) { utils::data("atlas_gene_cancer", package = "monopolyAtlas", envir = environment()); gc <- atlas_gene_cancer }
  if (is.null(genes_meta)) { utils::data("atlas_genes", package = "monopolyAtlas", envir = environment()); genes_meta <- atlas_genes }
  d <- demand[demand$gene == gene_symbol, ]
  if (!nrow(d)) stop("gene not in atlas_demand")
  d <- d[order(-d$pctile), ][seq_len(min(top_n, nrow(d))), ]
  d$tissue <- factor(d$tissue, levels = rev(d$tissue))
  g1 <- ggplot2::ggplot(d, ggplot2::aes(x = .data$pctile, y = .data$tissue)) +
    ggplot2::geom_col(fill = "#0072B2", width = 0.65, alpha = 0.9) +
    ggplot2::geom_vline(xintercept = 0.9, linetype = 2, color = "#D55E00", linewidth = 0.5) +
    ggplot2::labs(x = "GTEx expression percentile", y = NULL,
                  title = paste0(gene_symbol, " — host demand (WHERE)")) +
    theme_monopoly(base_size = 10)
  c2 <- gc[gc$gene == gene_symbol & gc$mono_freq > 0, ]
  c2 <- c2[order(-c2$mono_freq), ][seq_len(min(top_n, nrow(c2))), ]
  c2$cancer <- factor(c2$cancer, levels = rev(c2$cancer))
  g2 <- ggplot2::ggplot(c2, ggplot2::aes(x = .data$mono_freq, y = .data$cancer)) +
    ggplot2::geom_col(fill = "#D55E00", width = 0.65, alpha = 0.9) +
    ggplot2::labs(x = "Monopoly frequency", y = NULL, title = "Cancer occurrence (STATE)") +
    theme_monopoly(base_size = 10)
  meta <- genes_meta[genes_meta$gene == gene_symbol, ]
  info <- if (nrow(meta)) paste0("group: ", meta$group, " | module: ", meta$module,
                                 " | family: ", meta$family) else ""
  list(demand_plot = g1, cancer_plot = g2, info = info,
       meta = if (nrow(meta)) as.list(meta[1, ]) else list())
}

#' Headroom panel: baseline x mono_rate + permutation null band (Fig4 语言)
#' @param df data.frame: cancer, baseline (normal median), mono_rate, archetype (可选)
#' @param perm_null numeric vector of permutation rhos (可选, 画零分布注释)
#' @return ggplot
#' @export
plot_headroom <- function(df, perm_null = NULL) {
  stopifnot(all(c("cancer", "baseline", "mono_rate") %in% names(df)))
  rho_obs <- stats::cor(df$baseline, df$mono_rate, method = "spearman")
  sub <- paste0("Spearman rho = ", round(rho_obs, 3))
  if (!is.null(perm_null)) {
    p_emp <- max(mean(abs(perm_null) >= abs(rho_obs)), 1 / (length(perm_null) + 1))
    sub <- paste0(sub, " | permutation P = ", signif(p_emp, 2), " (n=", length(perm_null), ")")
  }
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$baseline, y = .data$mono_rate))
  if ("archetype" %in% names(df)) {
    p <- p + ggplot2::geom_point(ggplot2::aes(color = .data$archetype), size = 2.6, alpha = 0.9) +
      ggplot2::scale_color_manual(values = c("new_occupancy (low baseline)" = "#D55E00",
                                             "program_continuation (mid baseline)" = "#0072B2",
                                             "locked (high baseline)" = "grey55"), name = NULL)
  } else {
    p <- p + ggplot2::geom_point(color = "#0072B2", size = 2.6, alpha = 0.9)
  }
  p + ggplot2::geom_smooth(method = "loess", se = FALSE, color = "grey40", linewidth = 0.7,
                           linetype = 2, formula = y ~ x) +
    ggplot2::labs(x = "Normal baseline TMI5 (headroom = 1 - baseline)",
                  y = "Monopolized tumor fraction",
                  title = "Tissue permissiveness atlas", subtitle = sub) +
    theme_monopoly() + ggplot2::theme(legend.position = "top")
}
