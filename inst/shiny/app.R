library(shiny)
library(ggplot2)

genes_meta <- get0("atlas_genes", envir = asNamespace("monopolyAtlas"))
if (is.null(genes_meta)) { data("atlas_genes", package = "monopolyAtlas"); genes_meta <- atlas_genes }
data("atlas_gene_cancer", package = "monopolyAtlas")
data("atlas_demand", package = "monopolyAtlas")
data("atlas_cancer", package = "monopolyAtlas")

gchoices <- sort(genes_meta$gene)
cchoices <- sort(unique(atlas_gene_cancer$cancer))

ui <- fluidPage(
  titlePanel(h2("Transcriptomic Monopoly Atlas v3", style = "text-align:center")),
  sidebarLayout(
    sidebarPanel(
      selectizeInput("g", "Gene:", choices = gchoices, selected = "WFDC2"),
      selectInput("c", "Cancer:", choices = c("ALL", cchoices), selected = "ALL"),
      hr(), tableOutput("info"), width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Demand (WHERE)", plotOutput("demand", height = 480)),
        tabPanel("Occurrence (STATE)", plotOutput("occur", height = 480)),
        tabPanel("Cancer atlas", plotOutput("catlas", height = 480)),
        tabPanel("Table", DT::DTOutput("tbl"))
      ), width = 9
    )
  )
)

server <- function(input, output, session) {
  gd <- reactive({ genes_meta[genes_meta$gene == input$g, ] })

  output$info <- renderTable({
    g <- gd()
    if (!nrow(g)) return(data.frame())
    data.frame(
      Dim = c("group", "module", "family", "secreted"),
      Val = c(g$group, g$module, g$family, as.character(g$secreted_flag))
    )
  }, rownames = FALSE, colnames = FALSE)

  output$demand <- renderPlot({
    d <- atlas_demand[atlas_demand$gene == input$g, ]
    d <- d[order(-d$pctile), ][seq_len(min(12, nrow(d))), ]
    d$tissue <- factor(d$tissue, levels = rev(d$tissue))
    ggplot(d, aes(pctile, tissue)) +
      geom_col(fill = "#0072B2", width = 0.65) +
      geom_vline(xintercept = 0.9, linetype = 2, color = "#D55E00") +
      labs(x = "GTEx expression percentile", y = NULL,
           title = paste0(input$g, " — host demand")) +
      theme_classic(base_size = 13)
  })

  output$occur <- renderPlot({
    d <- atlas_gene_cancer[atlas_gene_cancer$gene == input$g & atlas_gene_cancer$mono_freq > 0, ]
    d <- d[order(-d$mono_freq), ][seq_len(min(12, nrow(d))), ]
    d$cancer <- factor(d$cancer, levels = rev(d$cancer))
    ggplot(d, aes(mono_freq, cancer)) +
      geom_col(fill = "#D55E00", width = 0.65) +
      labs(x = "Monopoly frequency", y = NULL,
           title = paste0(input$g, " — cancer occurrence")) +
      theme_classic(base_size = 13)
  })

  output$catlas <- renderPlot({
    d <- atlas_gene_cancer
    if (input$c != "ALL") d <- d[d$cancer == input$c, ]
    d <- d[d$mono_freq > 0, ]
    top <- stats::aggregate(mono_freq ~ gene, d, max)
    top <- top[order(-top$mono_freq), ][seq_len(min(15, nrow(top))), ]
    dd <- d[d$gene %in% top$gene, ]
    dd$gene <- factor(dd$gene, levels = rev(top$gene))
    ggplot(dd, aes(mono_freq, gene, fill = cancer)) +
      geom_col(position = "dodge", alpha = 0.9) +
      scale_fill_manual(values = monopolyAtlas::monopoly_palette(8)) +
      labs(x = "Monopoly frequency", y = NULL,
           title = ifelse(input$c == "ALL", "Top carriers (pan-cancer)", paste0("Top carriers — ", input$c))) +
      theme_classic(base_size = 13) + theme(legend.position = "none")
  })

  output$tbl <- DT::renderDT({
    atlas_gene_cancer[atlas_gene_cancer$gene == input$g, ]
  })
}

shinyApp(ui, server)
