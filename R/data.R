# 数据集文档
#' 74 monopoly genes annotation (36 cancer-specific + 38 tissue-shared)
#' @format data.frame: gene, group, module, family, secreted_flag, competence scores
"atlas_genes"

#' 32 cancer types monopoly summary
#' @format data.frame: cancer, n_tumor, n_mono, mono_rate, tmi5_median, top_carriers
"atlas_cancer"

#' 74 genes x 32 cancers monopoly frequency (long)
#' @format data.frame: gene, cancer, mono_freq, n_mono, n_mono_samples
"atlas_gene_cancer"

#' 74 genes x 52 GTEx tissues demand (expression percentile, long)
#' @format data.frame: gene, tissue, pctile, med_tpm
"atlas_demand"
