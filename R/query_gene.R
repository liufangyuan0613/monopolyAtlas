# Monopoly Atlas query functions
query_gene <- function(gene_symbol, atlas_data=NULL) {
  if(is.null(atlas_data)) { data("monopoly_atlas", package="monopolyAtlas", envir=environment()); atlas_data <- monopoly_atlas }
  row <- atlas_data[atlas_data$gene == gene_symbol, ]
  if(nrow(row) == 0) stop("Gene not found")
  as.list(row[1,])
}
get_top_monopoly <- function(n=10) {
  data("monopoly_atlas", package="monopolyAtlas", envir=environment())
  head(monopoly_atlas[order(-monopoly_atlas$MonopolyScore), ], n)
}
launch_shiny <- function() {
  shiny::runApp(system.file("shiny", package="monopolyAtlas"))
}
