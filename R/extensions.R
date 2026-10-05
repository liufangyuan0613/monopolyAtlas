# 拓展模块: 出图导出 / 零模型模拟 / 单细胞垄断 / 空间 Moran's I

#' 每 panel 单独导出 PDF + PNG (项目图件纪律)
#' @param plot ggplot (或 ComplexHeatmap 对象, 走 pdf/png 设备)
#' @param name 文件基名 (不含扩展名)
#' @param dir 输出目录
#' @param width,height 尺寸 (英寸)
#' @return  invisible 输出路径向量
#' @export
export_pdf <- function(plot, name, dir = ".", width = 7, height = 5) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  f_pdf <- file.path(dir, paste0(name, ".pdf"))
  f_png <- file.path(dir, paste0(name, ".png"))
  grDevices::pdf(f_pdf, width = width, height = height)
  print(plot); grDevices::dev.off()
  grDevices::png(f_png, width = width, height = height, units = "in", res = 300)
  print(plot); grDevices::dev.off()
  invisible(c(pdf = f_pdf, png = f_png))
}

#' 零模型表达矩阵 + TMI 零分布 (检验统计是否由定义几何自动产生)
#' @param n_genes,n_samples 维度
#' @param shape gamma 形状参数 (越小越集中; 0.5 接近真实重尾)
#' @param seed 随机种子
#' @return list(matrix, tmi5_null)
#' @export
simulate_share_null <- function(n_genes = 20000, n_samples = 200, shape = 0.5, seed = 42) {
  set.seed(seed)
  m <- matrix(stats::rgamma(n_genes * n_samples, shape = shape, rate = 1),
              nrow = n_genes)
  rownames(m) <- paste0("G", seq_len(n_genes))
  colnames(m) <- paste0("S", seq_len(n_samples))
  list(matrix = m, tmi5_null = compute_tmi(m, k = 5, exclude_technical = FALSE))
}

#' 单细胞/单点垄断摘要: 每组 carrier 主导率 (Fig3 FTL-myeloid 型分析)
#' @param expr genes x cells 矩阵 (稀疏支持)
#' @param carriers carrier 基因集 (如 atlas_genes$gene)
#' @param groups 每细胞分组 (celltype 等), 长度=ncol(expr)
#' @return data.frame: group, n_cells, carrier_top1_rate, carrier_share_median
#' @export
carrier_dominance <- function(expr, carriers, groups) {
  expr <- if (is.data.frame(expr)) as.matrix(expr) else expr
  if (is.null(rownames(expr))) stop("expr must have rownames")
  stopifnot(length(groups) == ncol(expr))
  keep <- intersect(carriers, rownames(expr))
  if (!length(keep)) stop("no carrier genes found in expr")
  cs <- if (methods::is(expr, "sparseMatrix")) Matrix::colSums(expr) else colSums(expr)
  cs[cs <= 0] <- NA_real_
  share <- if (methods::is(expr, "sparseMatrix")) Matrix::t(Matrix::t(expr) / cs) else sweep(expr, 2, cs, "/")
  gidx <- match(keep, rownames(share))
  top1 <- apply(share, 2, function(x) rownames(share)[which.max(x)])
  cshare <- if (methods::is(share, "sparseMatrix")) Matrix::colSums(share[gidx, , drop = FALSE]) else colSums(share[gidx, , drop = FALSE])
  df <- data.frame(group = groups, top1 = top1, carrier_share = as.numeric(cshare))
  do.call(rbind, lapply(split(df, df$group), function(s) {
    data.frame(group = s$group[1], n_cells = nrow(s),
               carrier_top1_rate = mean(s$top1 %in% keep),
               carrier_share_median = stats::median(s$carrier_share, na.rm = TRUE))
  }))
}

#' Moran's I + 置换检验 (空间 spot, 小鼠 HCC 口径: kNN6, 499 置换)
#' @param x 数值向量 (如 spot 的 TMI5)
#' @param coords n x 2 坐标矩阵
#' @param k kNN 邻居数 (默认 6)
#' @param n_perm 置换次数 (默认 499)
#' @param seed 种子
#' @return list(I, z, p_emp, expected)
#' @export
morans_i_perm <- function(x, coords, k = 6, n_perm = 499, seed = 42) {
  stopifnot(length(x) == nrow(coords))
  n <- length(x)
  dd <- as.matrix(stats::dist(coords))
  diag(dd) <- Inf
  W <- matrix(0, n, n)
  for (i in seq_len(n)) {
    nb <- order(dd[i, ])[seq_len(k)]
    W[i, nb] <- 1; W[nb, i] <- 1
  }
  s0 <- sum(W)
  moran <- function(v) {
    z <- v - mean(v)
    (n / s0) * sum(W * (z %o% z)) / sum(z^2)
  }
  i_obs <- moran(x)
  set.seed(seed)
  null <- replicate(n_perm, moran(sample(x)))
  mu <- mean(null); sd_ <- stats::sd(null)
  z_sc <- (i_obs - mu) / sd_
  p_emp <- max(mean(null >= i_obs), 1 / (n_perm + 1))
  list(I = i_obs, z = z_sc, p_emp = p_emp, expected = mu,
       note = paste0("kNN", k, ", ", n_perm, " permutations"))
}
