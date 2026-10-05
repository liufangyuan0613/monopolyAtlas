# data-raw/build_atlas_data.R
# 从 mono_v2 upgrade 数据包生成随包 .rda (可复现; 需要本机 upgrade 目录)
# 用法: Rscript data-raw/build_atlas_data.R
library(data.table)

UP <- "E:/cancer_mg/mono_v2/data/upgrade"
AN <- "E:/cancer_mg/mono_v2/data/analysis"
OUT <- "data"
dir.create(OUT, showWarnings = FALSE)

# ---- atlas_genes: 74 垄断基因注释 ----
g5 <- fread(file.path(UP, "Fig5_WHICH/fig5_gene_master_v2.tsv"))
g5 <- g5[label_monopoly == 1]
f74 <- fread(file.path(UP, "Fig5_WHICH/fig5F_full74.tsv"))
host <- fread(file.path(UP, "Fig4_WHERE/fig4C_host_demand.tsv"))
atlas_genes <- merge(
  g5[, .(gene, group, module, family, secreted_flag)],
  f74[, .(gene, best_pctile, mean_pctile, n_tissue_ge90, pred_demand, pred_joint)],
  by = "gene", all.x = TRUE)
atlas_genes <- merge(atlas_genes, host[, .(gene, host_tissue, gtex_percentile)],
                     by = "gene", all.x = TRUE)
setnames(atlas_genes,
         c("best_pctile", "mean_pctile", "n_tissue_ge90", "pred_demand", "pred_joint"),
         c("demand_best_pctile", "demand_mean_pctile", "demand_n_tissue_ge90",
           "competence_demand_score", "competence_joint_score"))
atlas_genes <- as.data.frame(atlas_genes)
stopifnot(nrow(atlas_genes) == 74)

# ---- atlas_gene_cancer: 74 x 32 垄断频率 (严格复刻提取器口径) ----
ex <- fread("E:/cancer_mg/mono_v2/data/step2_excess_per_sample.tsv")
pt <- ex[sample_type == "Primary Tumor"]
mono_ids <- pt[pctile5 >= 0.95, .(file_id, project)]
t50 <- fread("E:/cancer_mg/mono_v2/data/step1_top50_tcga.tsv.gz")
setorder(t50, file_id, -share)
t50[, rank := seq_len(.N), by = file_id]
top5 <- t50[rank <= 5]
G74 <- atlas_genes$gene
sub <- top5[file_id %in% mono_ids$file_id & gene %in% G74]
sub <- merge(sub, mono_ids, by = "file_id")
nmono <- mono_ids[, .(n_mono_samples = .N), by = project]
freq <- sub[, .(n_mono = .N), by = .(gene, project)]
freq <- merge(freq, nmono, by = "project")
freq[, mono_freq := n_mono / n_mono_samples]
atlas_gene_cancer <- merge(
  CJ(gene = G74, project = unique(pt$project)),
  freq[, .(gene, project, n_mono, n_mono_samples, mono_freq)],
  by = c("gene", "project"), all.x = TRUE)
atlas_gene_cancer[is.na(mono_freq), `:=`(n_mono = 0L, mono_freq = 0)]
atlas_gene_cancer[, cancer := sub("TCGA-", "", project)]
atlas_gene_cancer <- as.data.frame(atlas_gene_cancer[, .(gene, cancer, mono_freq, n_mono, n_mono_samples)])

# 锚点硬断言: IGHG1 x SKCM = 0.625
anchor <- atlas_gene_cancer[atlas_gene_cancer$gene == "IGHG1" & atlas_gene_cancer$cancer == "SKCM", "mono_freq"]
stopifnot(abs(anchor - 0.625) < 1e-6)

# ---- atlas_cancer: 32 癌种汇总 ----
byc <- pt[, .(n_tumor = .N, n_mono = sum(pctile5 >= 0.95), tmi5_median = median(tmi5)),
          by = project]
byc[, mono_rate := n_mono / n_tumor]
byc[, cancer := sub("TCGA-", "", project)]
topc <- atlas_gene_cancer[atlas_gene_cancer$mono_freq > 0, ]
topc <- topc[order(topc$cancer, -topc$mono_freq), ]
topc <- aggregate(gene ~ cancer, topc, function(x) paste(utils::head(x, 3), collapse = "|"))
setnames(topc, "gene", "top_carriers")
atlas_cancer <- merge(as.data.frame(byc[, .(cancer, n_tumor, n_mono, mono_rate, tmi5_median)]),
                      topc, by = "cancer", all.x = TRUE)

# ---- atlas_demand: 74 x 52 GTEx ----
dm <- fread(file.path(UP, "Fig5_WHICH/fig5D_demand_matrix.tsv.gz"))
atlas_demand <- as.data.frame(dm[gene %in% G74, .(gene, tissue, pctile, med_tpm)])
stopifnot(length(unique(atlas_demand$gene)) == 74)

# ---- atlas_scrna: 12 单细胞数据集 dataset x celltype 垄断摘要 (Fig3 层) ----
cells <- fread(file.path(UP, "Fig3_SOURCE/fig3_umap_all.tsv.gz"))
G74v <- atlas_genes$gene
cells[, is_carrier_top1 := top1 %in% G74v]
cells[, is_ftl_top1 := top1 == "FTL"]
sc <- cells[, .(n_cells = .N,
                tmi5_median = median(TMI5, na.rm = TRUE),
                carrier_top1_rate = mean(is_carrier_top1, na.rm = TRUE),
                ftl_top1_rate = mean(is_ftl_top1, na.rm = TRUE),
                malignant_frac = mean(malignant == TRUE, na.rm = TRUE)),
            by = .(dataset, cell_type)]
setorder(sc, dataset, -n_cells)
atlas_scrna <- as.data.frame(sc)
# 锚点硬断言: LUAD Myeloid/MAST ftl_top1_rate = 0.573 (fig3C 冻结值)
# 注意: stopifnot 对空向量会静默通过——必须先断言行数
anch <- atlas_scrna[atlas_scrna$dataset == "LUAD" & atlas_scrna$cell_type == "Myeloid/MAST", "ftl_top1_rate"]
stopifnot(length(anch) == 1, abs(anch - 0.573) < 0.005)

# ---- atlas_meta: 数据版本声明 ----
atlas_meta <- list(
  atlas_version = "v2026.10",
  build_date = as.character(Sys.Date()),
  manuscript = "mono_v2 (proof_jtm_v8 line, Cancer Research/Cell Reports target)",
  n_genes = nrow(atlas_genes),
  gene_split = c(cancer_specific = 36, tissue_shared = 38),
  anchors = c("IGHG1xSKCM mono_freq=0.625", "LUAD Myeloid/MAST FTL top1=0.573"),
  source_package = "E:/cancer_mg/mono_v2/data/upgrade (MANIFEST-registered)",
  regenerate = "data-raw/build_atlas_data.R")

save(atlas_genes, file = file.path(OUT, "atlas_genes.rda"), compress = "xz")
save(atlas_cancer, file = file.path(OUT, "atlas_cancer.rda"), compress = "xz")
save(atlas_gene_cancer, file = file.path(OUT, "atlas_gene_cancer.rda"), compress = "xz")
save(atlas_demand, file = file.path(OUT, "atlas_demand.rda"), compress = "xz")
save(atlas_scrna, file = file.path(OUT, "atlas_scrna.rda"), compress = "xz")
save(atlas_meta, file = file.path(OUT, "atlas_meta.rda"), compress = "xz")
cat("atlas data built:", nrow(atlas_genes), "genes,", nrow(atlas_gene_cancer), "gene-cancer rows,",
    nrow(atlas_demand), "demand rows,", nrow(atlas_scrna), "scrna summary rows\n")
