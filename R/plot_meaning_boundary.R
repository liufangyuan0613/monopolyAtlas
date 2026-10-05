# WHICH / MEANING / BOUNDARY 层 (Fig5-7 视觉语言)

#' Competence decomposition: model AUC bars + paired delta forest (Fig5 语言)
#' @param perf data.frame: model, ROC_AUC, ROC_CI_low, ROC_CI_high (PR 可选)
#' @param diffs data.frame (可选): comparison, delta_auc, ci_low, ci_high
#' @return ggplot (或 list 若给 diffs)
#' @export
plot_competence <- function(perf, diffs = NULL) {
  stopifnot(all(c("model", "ROC_AUC") %in% names(perf)))
  perf$model <- factor(perf$model, levels = perf$model[order(perf$ROC_AUC)])
  g1 <- ggplot2::ggplot(perf, ggplot2::aes(x = .data$ROC_AUC, y = .data$model)) +
    ggplot2::geom_col(fill = "#0072B2", width = 0.6, alpha = 0.9) +
    {if (all(c("ROC_CI_low", "ROC_CI_high") %in% names(perf)))
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$ROC_CI_low, xmax = .data$ROC_CI_high),
                              height = 0.25, linewidth = 0.5)} +
    ggplot2::geom_vline(xintercept = 0.5, linetype = 3, color = "grey60") +
    ggplot2::xlim(0.4, 1) +
    ggplot2::labs(x = "ROC AUC", y = NULL, title = "Carrier competence: Demand dominates") +
    theme_monopoly()
  if (is.null(diffs)) return(g1)
  diffs$comparison <- factor(diffs$comparison, levels = diffs$comparison)
  g2 <- ggplot2::ggplot(diffs, ggplot2::aes(x = .data$delta_auc, y = .data$comparison)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$ci_low, xmax = .data$ci_high),
                            height = 0.25, color = "#0072B2", linewidth = 0.6) +
    ggplot2::geom_point(color = "#D55E00", size = 2.6) +
    ggplot2::labs(x = "paired delta AUC", y = NULL, title = "Increment audit") +
    theme_monopoly()
  list(performance = g1, paired_differences = g2)
}

#' Pathway effect atlas: original vs top5-excluded 双层 (Fig6 语言)
#' @param df data.frame: cancer, pathway, delta, fdr, delta_excl, fdr_excl
#' @return ggplot (两层 facet)
#' @export
plot_pathway_atlas <- function(df) {
  stopifnot(all(c("cancer", "pathway", "delta", "delta_excl") %in% names(df)))
  d1 <- data.frame(cancer = df$cancer, pathway = df$pathway, delta = df$delta,
                   sig = df$fdr < 0.05, layer = "original")
  d2 <- data.frame(cancer = df$cancer, pathway = df$pathway, delta = df$delta_excl,
                   sig = df$fdr_excl < 0.05, layer = "top5-excluded")
  dd <- rbind(d1, d2)
  dd$layer <- factor(dd$layer, levels = c("original", "top5-excluded"))
  lim <- max(abs(dd$delta), na.rm = TRUE)
  ggplot2::ggplot(dd, ggplot2::aes(x = .data$cancer, y = .data$pathway)) +
    ggplot2::geom_point(ggplot2::aes(size = abs(.data$delta), color = .data$delta,
                                     alpha = .data$sig)) +
    ggplot2::scale_color_gradient2(low = "#0072B2", mid = "grey90", high = "#D55E00",
                                   midpoint = 0, limits = c(-lim, lim)) +
    ggplot2::scale_alpha_manual(values = c(`TRUE` = 0.95, `FALSE` = 0.25), guide = "none") +
    ggplot2::facet_wrap(~layer) +
    ggplot2::labs(size = "|rank-biserial|", color = "effect",
                  title = "Program meaning: pathway effects and top-5 attenuation") +
    theme_monopoly(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}

#' Dependency boundary: Chronos distribution by group + exceptions (Fig7 语言)
#' @param df data.frame: gene, group, chronos (long) 或 summary (gene, group, frac_le_neg1)
#' @param exceptions 标注的基因 (默认 c("IFITM3","WFDC2"))
#' @return ggplot
#' @export
plot_dependency_boundary <- function(df, exceptions = c("IFITM3", "WFDC2")) {
  if ("chronos" %in% names(df)) {
    summ <- do.call(rbind, lapply(split(df, df$gene), function(s) {
      data.frame(gene = s$gene[1], group = s$group[1],
                 frac_le_neg1 = mean(s$chronos <= -1, na.rm = TRUE),
                 median_chronos = stats::median(s$chronos, na.rm = TRUE))
    }))
  } else {
    stopifnot(all(c("gene", "group", "frac_le_neg1") %in% names(df)))
    summ <- df
  }
  summ <- summ[order(summ$frac_le_neg1), ]
  summ$gene <- factor(summ$gene, levels = summ$gene)
  summ$lab <- summ$gene %in% exceptions
  ggplot2::ggplot(summ, ggplot2::aes(x = .data$frac_le_neg1, y = .data$gene)) +
    ggplot2::geom_vline(xintercept = 0.9, linetype = 2, color = "#D55E00", linewidth = 0.6) +
    ggplot2::geom_col(ggplot2::aes(fill = .data$lab), width = 0.7, alpha = 0.9) +
    ggplot2::scale_fill_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "#0072B2"), guide = "none") +
    ggplot2::annotate("text", x = 0.91, y = 1, label = "common-essential threshold (90%)",
                      hjust = 0, vjust = -0.5, size = 3, color = "#D55E00") +
    ggplot2::labs(x = "Fraction of cell lines with Chronos <= -1", y = NULL,
                  title = "Marker is not dependency",
                  subtitle = paste0(sum(summ$frac_le_neg1 >= 0.9, na.rm = TRUE),
                                    "/", nrow(summ), " genes pass strict common-essential")) +
    theme_monopoly(base_size = 10) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 6))
}
