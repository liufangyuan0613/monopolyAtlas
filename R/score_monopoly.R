compute_TMI <- function(expr, k = 5, normalize = FALSE) {
  if (is.data.frame(expr)) expr <- as.matrix(expr)
  if (!is.matrix(expr) || !is.numeric(expr)) stop("expr must be a numeric matrix")
  if (normalize) expr <- sweep(expr, 2, colSums(expr)/1e6, "/")
  n <- ncol(expr); tmi <- numeric(n)
  for (i in seq_len(n)) {
    s <- sort(expr[, i], decreasing = TRUE)
    tmi[i] <- sum(s[seq_len(min(k, length(s)))]) / sum(s)
  }
  names(tmi) <- colnames(expr); tmi
}

compute_Gini <- function(expr) {
  if (is.data.frame(expr)) expr <- as.matrix(expr)
  if (!is.matrix(expr) || !is.numeric(expr)) stop("expr must be a numeric matrix")
  n <- ncol(expr); gini <- numeric(n)
  for (i in seq_len(n)) {
    x <- sort(expr[, i]); x <- x[x > 0]
    if (length(x) < 2) { gini[i] <- 0; next }
    nx <- length(x)
    gini[i] <- (2 * sum(seq_len(nx) * x) - (nx + 1) * sum(x)) / (nx * sum(x))
  }
  names(gini) <- colnames(expr); gini
}

compute_HHI <- function(expr) {
  if (is.data.frame(expr)) expr <- as.matrix(expr)
  if (!is.matrix(expr) || !is.numeric(expr)) stop("expr must be a numeric matrix")
  n <- ncol(expr); hhi <- numeric(n)
  for (i in seq_len(n)) { shares <- expr[,i]/sum(expr[,i]); hhi[i] <- sum(shares^2) }
  names(hhi) <- colnames(expr); hhi
}

score_monopoly <- function(expr, k_top = 5, min_samples = 1) {
  if (is.data.frame(expr)) expr <- as.matrix(expr)
  if (!is.matrix(expr) || !is.numeric(expr)) stop("expr must be numeric matrix (genes x samples)")
  if (is.null(rownames(expr))) stop("expr must have rownames (gene symbols)")
  message("Scoring ", ncol(expr), " samples x ", nrow(expr), " genes")
  tmi <- compute_TMI(expr, k = k_top)
  gini <- compute_Gini(expr)
  hhi <- compute_HHI(expr)
  scores <- data.frame(sample = colnames(expr), TMI = tmi, Gini = gini, HHI = hhi)
  ns <- ncol(expr); ng <- nrow(expr)
  top_genes <- matrix("", ns, k_top)
  colnames(top_genes) <- paste0("top", seq_len(k_top)); rownames(top_genes) <- colnames(expr)
  gene_top_count <- setNames(rep(0, ng), rownames(expr))
  gene_frac_sum <- setNames(rep(0, ng), rownames(expr))
  gene_frac_n <- setNames(rep(0, ng), rownames(expr))
  for (i in seq_len(ns)) {
    frac <- expr[, i] / sum(expr[, i])
    idx <- order(frac, decreasing = TRUE)[seq_len(min(k_top, ng))]
    top_genes[i, ] <- rownames(expr)[idx]
    for (j in seq_along(idx)) {
      g <- rownames(expr)[idx[j]]
      gene_top_count[g] <- gene_top_count[g] + 1
      gene_frac_sum[g] <- gene_frac_sum[g] + frac[idx[j]]
      gene_frac_n[g] <- gene_frac_n[g] + 1
    }
  }
  gs <- data.frame(gene = rownames(expr), top_count = gene_top_count,
    top_freq = gene_top_count / ns, mean_frac = gene_frac_sum / pmax(gene_frac_n, 1))
  gs <- gs[order(-gs$top_freq, -gs$mean_frac), ]
  mg <- gs$gene[gs$top_count >= min_samples]
  res <- list(scores = scores, gene_stats = gs, monopoly_genes = mg,
    top_genes = top_genes, params = list(k_top = k_top, min_samples = min_samples))
  class(res) <- "monopolyResult"; res
}

print.monopolyResult <- function(x, ...) {
  cat("Monopoly Analysis: ", ncol(x$top_genes), " samples\n")
  cat("Top genes (k=", x$params$k_top, "): ", length(x$monopoly_genes), " monopoly genes\n", sep="")
  print(head(x$scores, 5))
  cat("\nTop monopoly genes:\n")
  print(head(x$gene_stats, 8))
}

plot.monopolyResult <- function(x, ...) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("ggplot2 required")
  p <- ggplot2::ggplot(x$scores, ggplot2::aes(x = TMI)) +
    ggplot2::geom_histogram(bins = 30, fill = "#D55E00", alpha = 0.8) +
    ggplot2::labs(title = "TMI Distribution",
      subtitle = paste("Mean TMI =", round(mean(x$scores$TMI), 3),
        "| Monopoly genes =", length(x$monopoly_genes)),
      x = "TMI", y = "Count") +
    ggplot2::theme_minimal(12)
  print(p); invisible(p)
}
