# Transcriptomic Monopoly Atlas

Pan-cancer transcriptomic monopoly gene query tool and interactive Shiny dashboard.

## Install

```r
remotes::install_github("liufangyuan0613/monopolyAtlas")
```

## Quick Start

```r
library(monopolyAtlas)
query_gene("EEF1A1")  # Query a gene
launch_shiny()         # Launch interactive app
```

## Dimensions

| Dimension | Description |
|-----------|-------------|
| MonopolyScore | Composite (0-1) |
| CCLE top-1% | Cell line top expression frequency |
| n_cancers | 33 cancer types |
| Gene structure | Length, exons, GC, TE |
| Cross-species | 5 species conservation |

## Citation

Liu F et al. "Transcriptomic monopoly reveals an evolutionarily conserved architecture of gene-expression resource allocation in cancer." (2026)
