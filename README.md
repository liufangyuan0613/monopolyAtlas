# Transcriptomic Monopoly Atlas (monopolyAtlas v3)

Pan-cancer transcriptomic monopoly: compute, query, visualize, and explore.
与 mono_v2 手稿同一口径 (TMI5 / top-5 share / 95th-percentile monopoly / 74-gene atlas)。

## Install

```r
remotes::install_github("liufangyuan0613/monopolyAtlas")
```

## Quick start

```r
library(monopolyAtlas)

# 1) 从表达矩阵算垄断指标 (bulk / scRNA 稀疏矩阵均可)
tmi  <- compute_tmi(expr, k = 5)
mono <- monopoly_status(tmi)

# 2) 查 atlas (74 垄断基因 x 32 癌种 x 52 GTEx 组织)
query_gene("WFDC2")
get_top_monopoly(10)

# 3) 正文视觉语言一键出图
plot_carrier_matrix(atlas_gene_cancer[atlas_gene_cancer$mono_freq > 0, ])
plot_headroom(headroom_df, perm_null = my_null)
card <- plot_gene_card("WFDC2"); card$demand_plot

# 4) 交互式 atlas
launch_shiny()
```

## Function map

| 层 | 函数 | 对应正文 |
|---|---|---|
| Compute | `compute_tmi` `top_share` `monopoly_status` `compute_gini` `compute_hhi` `score_monopoly` | pipeline step1/2 |
| Query | `query_gene` `get_top_monopoly` | atlas 数据 |
| STATE | `plot_tmi_ridge` `plot_state_curve` | Fig1 |
| CARRIER | `plot_carrier_matrix` `plot_gene_card` | Fig2 |
| WHERE | `plot_headroom` | Fig4 |
| WHICH | `plot_competence` | Fig5 |
| MEANING | `plot_pathway_atlas` | Fig6 |
| BOUNDARY | `plot_dependency_boundary` | Fig7 |
| Theme | `theme_monopoly` `monopoly_palette` | 全部 |

## Data shipped

- `atlas_genes` — 74 monopoly genes (36 cancer-specific + 38 tissue-shared) 全注释
- `atlas_gene_cancer` — 74 x 32 垄断频率 (严格管线口径; IGHG1xSKCM=0.625 锚点测试)
- `atlas_cancer` — 32 癌种汇总 (n_mono / mono_rate / top carriers)
- `atlas_demand` — 74 x 52 GTEx tissue expression percentile (WHERE 层)

## Citation

Liu F et al. "Transcriptomic monopoly reveals an evolutionarily conserved architecture
of gene-expression resource allocation in cancer." (2026, in preparation)
