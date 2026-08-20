# Monopoly Atlas query + visualization

#' Query a gene in the Monopoly Atlas
#' @param gene_symbol Gene symbol (e.g. "EEF1A1")
#' @param atlas_data Optional data frame
#' @return Named list with monopoly dimensions
#' @export
query_gene <- function(gene_symbol, atlas_data = NULL) {
  if (is.null(atlas_data)) {
    data("monopoly_atlas", package = "monopolyAtlas", envir = environment())
    atlas_data <- monopoly_atlas
  }
  row <- atlas_data[atlas_data$gene == gene_symbol, ]
  if (nrow(row) == 0) stop("Gene not found in atlas")
  as.list(row[1, ])
}

#' Get top monopoly genes
#' @param n Number of genes
#' @return data.frame
#' @export
get_top_monopoly <- function(n = 10) {
  data("monopoly_atlas", package = "monopolyAtlas", envir = environment())
  head(monopoly_atlas[order(-monopoly_atlas$MonopolyScore), ], n)
}

#' Plot monopoly gene profile
#' @param gene_symbol Gene symbol
#' @param atlas_data Optional data frame
#' @return A ggplot object
#' @export
plot_gene <- function(gene_symbol, atlas_data = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("ggplot2 is required for plotting. Install with install.packages('ggplot2')")

  if (is.null(atlas_data)) {
    data("monopoly_atlas", package = "monopolyAtlas", envir = environment())
    atlas_data <- monopoly_atlas
  }
  row <- atlas_data[atlas_data$gene == gene_symbol, ]
  if (nrow(row) == 0) stop("Gene not found in atlas")

  df <- data.frame(
    Feature = c("MonopolyScore", "CCLE top-1", "Breadth/33", "Compact", "Low TE"),
    Value = c(
      max(row$MonopolyScore, 0, na.rm = TRUE),
      min(max(row$CCLE_Cf, 0, na.rm = TRUE) * 5, 1),
      max(row$n_cancers_mono, 0, na.rm = TRUE) / 33,
      1 - min(max(row$gene_len_kb, 0, na.rm = TRUE) / 100, 1),
      1 - min(max(row$gene_TE_All, 0, na.rm = TRUE) * 5, 1)
    )
  )

  cls <- if (is.null(row$conservation_class) || is.na(row$conservation_class))
    "Unknown" else row$conservation_class

  ggplot2::ggplot(df, ggplot2::aes(x = Feature, y = Value, fill = Feature)) +
    ggplot2::geom_col(width = 0.5, alpha = 0.85) +
    ggplot2::scale_fill_manual(values = c(
      "MonopolyScore" = "#D55E00", "CCLE top-1" = "#E69F3B",
      "Breadth/33" = "#56B4E9", "Compact" = "#009E73", "Low TE" = "#CC79A7"
    )) +
    ggplot2::ylim(0, 1.05) +
    ggplot2::labs(
      title = paste("Monopoly Profile:", gene_symbol),
      subtitle = paste0(
        "CCLE top-1: ", round(max(row$CCLE_Cf,0,na.rm=TRUE)*100, 1), "% | ",
        max(row$n_cancers_mono,0,na.rm=TRUE), "/33 cancers | ", cls
      ),
      x = "", y = "Score (0-1)"
    ) +
    ggplot2::theme_minimal(14) +
    ggplot2::theme(
      legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", color = "#D55E00")
    )
}

#' Launch interactive Shiny app
#' @export
launch_shiny <- function() {
  app_dir <- system.file("shiny", package = "monopolyAtlas")
  if (app_dir == "") stop("Shiny app not found in package")
  shiny::runApp(app_dir)
}
