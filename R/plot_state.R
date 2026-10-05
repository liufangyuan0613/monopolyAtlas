# STATE 层 (Fig1 视觉语言)

#' TMI5 distribution by cohort (ridge/violin)
#' @param df data.frame: sample, tmi5, cohort (cancer/tissue), state (mono/nonmono/normal 可选)
#' @param order_by "median" (默认) 或 "none"
#' @return ggplot
#' @export
plot_tmi_ridge <- function(df, order_by = c("median", "none")) {
  order_by <- match.arg(order_by)
  stopifnot(all(c("sample", "tmi5", "cohort") %in% names(df)))
  if (order_by == "median") {
    ord <- stats::aggregate(tmi5 ~ cohort, df, stats::median, na.rm = TRUE)
    df$cohort <- factor(df$cohort, levels = ord$cohort[order(ord$tmi5)])
  }
  fill_col <- if ("state" %in% names(df)) "state" else NULL
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$tmi5, y = .data$cohort)) +
    ggplot2::geom_violin(ggplot2::aes(fill = .data[[fill_col %||% "cohort"]]),
                         scale = "width", trim = TRUE, linewidth = 0.2, alpha = 0.85) +
    ggplot2::geom_boxplot(width = 0.12, outlier.size = 0.2, alpha = 0.6) +
    ggplot2::labs(x = "TMI5 (top-5 expression share)", y = NULL,
                  title = "Tissue-level TMI5 landscape") +
    theme_monopoly() +
    ggplot2::theme(legend.position = "none")
  if (!is.null(fill_col)) {
    p <- p + ggplot2::scale_fill_manual(values = .mono_group_colors, name = NULL) +
      ggplot2::theme(legend.position = "top")
  }
  p
}

#' State curve: excess share vs baseline (U-shape / state relation)
#' @param df data.frame: cohort, baseline (normal median TMI5), tumor_median, mono_rate
#' @param label_top 标注前 N 个极端 cohort (默认 6)
#' @return ggplot
#' @export
plot_state_curve <- function(df, label_top = 6) {
  stopifnot(all(c("cohort", "baseline", "mono_rate") %in% names(df)))
  df <- df[order(df$baseline), ]
  df$ext <- rank(-abs(stats::fitted(stats::loess(mono_rate ~ baseline, df)) - df$mono_rate))
  lab <- df[df$ext <= label_top, ]
  ggplot2::ggplot(df, ggplot2::aes(x = .data$baseline, y = .data$mono_rate)) +
    ggplot2::geom_smooth(method = "loess", se = TRUE, color = "#0072B2",
                         fill = "#56B4E9", alpha = 0.5, linewidth = 0.8) +
    ggplot2::geom_point(ggplot2::aes(color = .data$mono_rate), size = 2.4, alpha = 0.9) +
    ggplot2::geom_text(data = lab, ggplot2::aes(label = .data$cohort),
                       size = 3, vjust = -0.8, check_overlap = TRUE) +
    ggplot2::scale_color_gradient(low = "#56B4E9", high = "#D55E00", guide = "none") +
    ggplot2::labs(x = "Normal baseline TMI5 (GTEx median)",
                  y = "Monopolized tumor fraction",
                  title = "Baseline determinism: headroom decides where monopoly happens",
                  subtitle = paste0("Spearman rho = ",
                                    round(stats::cor(df$baseline, df$mono_rate, method = "spearman"), 3))) +
    theme_monopoly()
}

`%||%` <- function(a, b) if (is.null(a)) b else a
