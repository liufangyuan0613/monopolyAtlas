# Core monopoly metrics (v3) — vectorized, sparse-aware.
# 口径与 mono_v2 管线 step1/step2 一致 (v3.1.0 起钉死):
#   share_i = expr_i / sum(expr)  (列=样本; 分母=样本全部基因总和)
#   TMI(k) = sum(top-k shares)
#   技术基因 ^(MT-|MTRNR|RPL|RPS) 只剔出排名, 保留在分母 (管线 step1 政策)
#   renormalize=TRUE = v3.0 遗留口径 (分母=非技术基因总和), 仅供对照, 勿用于主分析

.as_matrix <- function(expr) {
  if (methods::is(expr, "sparseMatrix") || methods::is(expr, "Matrix")) return(expr)
  if (is.data.frame(expr)) expr <- as.matrix(expr)
  if (!is.matrix(expr) || !is.numeric(expr)) stop("expr must be a numeric matrix (genes x samples)")
  expr
}

.col_sums <- function(m) {
  if (methods::is(m, "sparseMatrix")) return(Matrix::colSums(m))
  base::colSums(m)
}

#' Top-k expression share per sample
#' @param expr numeric matrix (genes x samples), sparse supported
#' @param k number of top genes (default 5, 与正文 TMI5 一致)
#' @param exclude_technical drop ^(MT-|MTRNR|RPL|RPS) from RANKING only
#'   (管线 step1 政策: 技术基因不参与排名, 但保留在分母)
#' @param renormalize legacy v3.0 behavior: 剔除后按剩余(非技术)基因总和重归一.
#'   默认 FALSE = 管线口径 (分母=全部基因总和); TRUE 仅供与 v3.0 对照
#' @return numeric vector of top-k shares per sample (0-1)
#' @export
top_share <- function(expr, k = 5, exclude_technical = TRUE, renormalize = FALSE) {
  expr <- .as_matrix(expr)
  if (is.null(rownames(expr))) stop("expr must have rownames (gene symbols)")
  cs <- .col_sums(expr)
  cs[cs <= 0] <- NA_real_
  keep <- rep(TRUE, nrow(expr))
  if (exclude_technical) {
    keep <- !grepl("^(MT-|MTRNR|RPL|RPS)", rownames(expr))
    if (renormalize) {
      # legacy: 分母改为非技术基因总和
      expr <- expr[keep, , drop = FALSE]
      cs <- .col_sums(expr)
      cs[cs <= 0] <- NA_real_
      keep <- rep(TRUE, nrow(expr))
    }
  }
  share <- expr
  if (methods::is(share, "sparseMatrix")) {
    share <- Matrix::t(Matrix::t(share) / cs)
  } else {
    share <- sweep(expr, 2, cs, "/")
  }
  share <- share[keep, , drop = FALSE]
  k <- min(k, nrow(share))
  out <- vapply(seq_len(ncol(share)), function(j) {
    x <- share[, j]
    if (methods::is(x, "sparseVector") || methods::is(x, "sparseMatrix")) x <- as.numeric(x)
    x <- x[is.finite(x)]
    if (!length(x)) return(NA_real_)
    sum(sort(x, decreasing = TRUE)[seq_len(k)])
  }, numeric(1))
  names(out) <- colnames(expr)
  out
}

#' Transcriptomic Monopoly Index (TMI)
#' @inheritParams top_share
#' @param normalize library-size normalize to pseudo-CPM first (对原始 counts 使用)
#' @return numeric vector, TMI per sample
#' @details v3.1.0 起分母恒为样本全部基因总和 (管线口径);
#'   技术基因仅剔出排名. 与 v3.0 的绝对值不可直接比较 (v3.0 分母偏小, TMI 偏高).
#' @export
compute_tmi <- function(expr, k = 5, normalize = FALSE, exclude_technical = TRUE) {
  expr <- .as_matrix(expr)
  if (normalize) {
    cs <- .col_sums(expr); cs[cs <= 0] <- NA_real_
    expr <- if (methods::is(expr, "sparseMatrix")) Matrix::t(Matrix::t(expr) / (cs / 1e6)) else sweep(expr, 2, cs / 1e6, "/")
  }
  top_share(expr, k = k, exclude_technical = exclude_technical)
}

#' Monopoly status per sample (管线 step2 判定)
#' @param tmi5 TMI5 vector (from compute_tmi k=5)
#' @param pctile_cut within-cohort percentile cutoff (默认 0.95)
#' @param cohort optional factor; 在其内部按分位判定 (缺省=全体)
#' @return logical vector, TRUE = monopolized
#' @export
monopoly_status <- function(tmi5, pctile_cut = 0.95, cohort = NULL) {
  if (is.null(cohort)) {
    cut <- stats::quantile(tmi5, pctile_cut, na.rm = TRUE, names = FALSE, type = 8)
    return(tmi5 >= cut)
  }
  cohort <- as.factor(cohort)
  out <- rep(NA, length(tmi5))
  for (lv in levels(cohort)) {
    idx <- which(cohort == lv)
    cut <- stats::quantile(tmi5[idx], pctile_cut, na.rm = TRUE, names = FALSE, type = 8)
    out[idx] <- tmi5[idx] >= cut
  }
  out
}

#' Gini coefficient per sample
#' @export
compute_gini <- function(expr) {
  expr <- .as_matrix(expr)
  vapply(seq_len(ncol(expr)), function(i) {
    x <- sort(as.numeric(expr[, i])); x <- x[x > 0]
    nx <- length(x)
    if (nx < 2) return(0)
    (2 * sum(seq_len(nx) * x) - (nx + 1) * sum(x)) / (nx * sum(x))
  }, numeric(1)) |> stats::setNames(colnames(expr))
}

#' Herfindahl-Hirschman index per sample
#' @export
compute_hhi <- function(expr) {
  expr <- .as_matrix(expr)
  cs <- .col_sums(expr); cs[cs <= 0] <- NA_real_
  vapply(seq_len(ncol(expr)), function(i) {
    s <- as.numeric(expr[, i]) / cs[i]
    sum(s^2, na.rm = TRUE)
  }, numeric(1)) |> stats::setNames(colnames(expr))
}

#' Score a matrix: per-sample TMI/Gini/HHI + per-gene top-k recurrence
#' @param expr numeric matrix (genes x samples)
#' @param k_top top-k for TMI/recurrence (默认 5)
#' @return list(samples=data.frame, genes=data.frame)
#' @export
score_monopoly <- function(expr, k_top = 5) {
  expr <- .as_matrix(expr)
  if (is.null(rownames(expr))) stop("expr must have rownames (gene symbols)")
  message("Scoring ", ncol(expr), " samples x ", nrow(expr), " genes")
  tmi <- compute_tmi(expr, k = k_top, exclude_technical = TRUE)
  samples <- data.frame(sample = colnames(expr), TMI = tmi,
                        Gini = compute_gini(expr), HHI = compute_hhi(expr),
                        mono = monopoly_status(tmi), row.names = NULL)
  cs <- .col_sums(expr); cs[cs <= 0] <- NA_real_
  share <- if (methods::is(expr, "sparseMatrix")) Matrix::t(Matrix::t(expr) / cs) else sweep(expr, 2, cs, "/")
  # 技术基因不参与 top-k 排名 (管线口径, v3.1.0 起); 分母保持全基因总和
  keep <- !grepl("^(MT-|MTRNR|RPL|RPS)", rownames(share))
  share <- share[keep, , drop = FALSE]
  cnt <- setNames(integer(nrow(share)), rownames(share))
  frac_sum <- setNames(numeric(nrow(share)), rownames(share))
  for (j in seq_len(ncol(share))) {
    x <- as.numeric(share[, j])
    idx <- order(x, decreasing = TRUE)[seq_len(min(k_top, length(x)))]
    cnt[idx] <- cnt[idx] + 1L
    frac_sum[idx] <- frac_sum[idx] + x[idx]
  }
  genes <- data.frame(gene = rownames(share), top_count = as.integer(cnt),
                      mean_share_when_top = ifelse(cnt > 0, frac_sum / pmax(cnt, 1), 0),
                      row.names = NULL)
  list(samples = samples, genes = genes)
}
