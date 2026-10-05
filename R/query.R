# Atlas 数据查询

#' Query a gene in the monopoly atlas (v3: 74-gene 全注释)
#' @param gene_symbol gene symbol
#' @param genes_meta optional data.frame (默认=内置 atlas_genes)
#' @return named list
#' @export
query_gene <- function(gene_symbol, genes_meta = NULL) {
  if (is.null(genes_meta)) {
    utils::data("atlas_genes", package = "monopolyAtlas", envir = environment())
    genes_meta <- atlas_genes
  }
  row <- genes_meta[genes_meta$gene == gene_symbol, ]
  if (nrow(row) == 0) stop("Gene not found in atlas (74 monopoly genes)")
  as.list(row[1, ])
}

#' Top monopoly genes by cancer occurrence
#' @param n number
#' @param gc optional data.frame (默认=内置 atlas_gene_cancer)
#' @return data.frame
#' @export
get_top_monopoly <- function(n = 10, gc = NULL) {
  if (is.null(gc)) {
    utils::data("atlas_gene_cancer", package = "monopolyAtlas", envir = environment())
    gc <- atlas_gene_cancer
  }
  agg <- stats::aggregate(mono_freq ~ gene, gc, function(x) c(n_hit = sum(x > 0), max_freq = max(x)))
  out <- data.frame(gene = agg$gene, n_cancers = agg$mono_freq[, "n_hit"],
                    max_freq = agg$mono_freq[, "max_freq"])
  out <- out[order(-out$n_cancers, -out$max_freq), ]
  utils::head(out, n)
}
