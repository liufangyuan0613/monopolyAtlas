# Monopoly Atlas query + visualization

query_gene <- function(gene_symbol, atlas_data = NULL) {
  if (is.null(atlas_data)) {
    data("monopoly_atlas", package = "monopolyAtlas", envir = environment())
    atlas_data <- monopoly_atlas
  }
  row <- atlas_data[atlas_data$gene == gene_symbol, ]
  if (nrow(row) == 0) stop("Gene not found")
  as.list(row[1, ])
}

get_top_monopoly <- function(n = 10) {
  data("monopoly_atlas", package = "monopolyAtlas", envir = environment())
  head(monopoly_atlas[order(-monopoly_atlas$MonopolyScore), ], n)
}

plot_gene <- function(gene_symbol, atlas_data = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  if (is.null(atlas_data)) {
    data("monopoly_atlas", package = "monopolyAtlas", envir = environment())
    atlas_data <- monopoly_atlas
  }
  row <- atlas_data[atlas_data$gene == gene_symbol, ]
  if (nrow(row) == 0) stop("Gene not found")
  
  df <- data.frame(
    Feature = c("Monopoly\nScore", "CCLE\ntop-1%", "Breadth\n/33", "Compact\n(1/kb)", "Low\nTE"),
    Value = c(
      row$MonopolyScore,
      min(row$CCLE_Cf * 5, 1),
      row$n_cancers_mono / 33,
      1 - min(row$gene_len_kb / 100, 1),
      1 - min(row$gene_TE_All * 5, 1)
    )
  )
  ggplot2::ggplot(df, ggplot2::aes(x = Feature, y = Value, fill = Feature)) +
    ggplot2::geom_col(width = 0.5, alpha = 0.85) +
    ggplot2::scale_fill_manual(values = c("#D55E00","#E69F3B","#56B4E9","#009E73","#CC79A7")) +
    ggplot2::ylim(0, 1.05) +
    ggplot2::labs(
      title = paste("Monopoly Profile:", gene_symbol),
      subtitle = paste0("CCLE top-1: ", round(row$CCLE_Cf*100,1), "% | ",
                        row$n_cancers_mono, "/33 cancers | ", row$conservation_class),
      x = "", y = "Score"
    ) +
    ggplot2::theme_minimal(14) +
    ggplot2::theme(legend.position = "none",
      plot.title = ggplot2::element_text(face = "bold", color = "#D55E00"))
}

launch_shiny <- function() {
  shiny::runApp(system.file("shiny", package = "monopolyAtlas"))
}
