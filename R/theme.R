# 主题与调色板 — 与 mono_v2 正文视觉语言一致
# (ComplexHeatmap 成对模块/气泡矩阵配色; 主图禁火山图纪律不在此层约束)

#' monopoly palette (Okabe-Ito 扩展, 正文主色)
#' @param n 需要的颜色数 (1-8)
#' @return character vector of hex colors
#' @export
monopoly_palette <- function(n = 8) {
  pal <- c("#D55E00", "#0072B2", "#009E73", "#E69F00", "#56B4E9",
           "#CC79A7", "#F0E442", "#999999")
  pal[seq_len(min(n, length(pal)))]
}

#' ggplot2 theme for monopoly figures
#' @param base_size base font size
#' @export
theme_monopoly <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      axis.line = ggplot2::element_line(linewidth = 0.4),
      axis.ticks = ggplot2::element_line(linewidth = 0.3),
      axis.title = ggplot2::element_text(face = "bold"),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0),
      plot.subtitle = ggplot2::element_text(color = "grey30"),
      legend.title = ggplot2::element_text(size = ggplot2::rel(0.9)),
      legend.background = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", color = NA),
      strip.text = ggplot2::element_text(face = "bold")
    )
}

# 分组配色 (74 基因三组 / 三态)
.mono_group_colors <- c(cancer_specific = "#D55E00", tissue_shared = "#0072B2",
                        background = "grey70", mono = "#D55E00", nonmono = "grey60",
                        normal = "#009E73")
