library(shiny); library(ggplot2); library(plotly); library(DT)

atlas <- readRDS("data/demo_atlas.rds")
atlas <- atlas[!is.na(atlas$gene), ]
rownames(atlas) <- atlas$gene
gchoices <- sort(atlas$gene)

ui <- fluidPage(
  titlePanel(h2("Transcriptomic Monopoly Atlas", style = "text-align:center")),
  sidebarLayout(
    sidebarPanel(
      selectizeInput("g", "Gene Symbol:", choices = gchoices, selected = "EEF1A1",
        options = list(placeholder = 'Search...')),
      hr(), tableOutput("info"), width = 3
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Radar", plotlyOutput("radar", height = 450)),
        tabPanel("Structure", plotOutput("struct", height = 450)),
        tabPanel("Species", plotOutput("species", height = 450)),
        tabPanel("Table", DTOutput("tbl"))
      ), width = 9
    )
  )
)

server <- function(input, output, session) {
  gd <- reactive({ req(input$g); atlas[input$g, ] })
  
  output$info <- renderTable({
    g <- gd(); cls <- if(is.null(g$conservation_class)) "?" else g$conservation_class
    data.frame(
      Dim = c("MonopolyScore","CCLE%","Cancers","Length","GC%","TE","Class"),
      Val = c(sprintf("%.3f",g$MonopolyScore), sprintf("%.1f%%",g$CCLE_Cf*100),
              sprintf("%d/33",g$n_cancers_mono), sprintf("%.1fkb",g$gene_len_kb),
              sprintf("%.1f%%",g$GC_content*100), sprintf("%.3f",g$gene_TE_All), cls)
    )
  }, rownames = FALSE, colnames = FALSE)

  output$radar <- renderPlotly({
    g <- gd()
    v <- c(MS=g$MonopolyScore, CCLE=min(g$CCLE_Cf*5,1), Breadth=g$n_cancers_mono/33,
           Compact=1-min(g$gene_len_kb/100,1), LowTE=1-min(g$gene_TE_All*5,1),
           Cons=(g$human_monopoly+g$pig_monopoly+g$cattle_monopoly+g$chicken_ortholog+g$mouse_ortholog)/5)
    plot_ly(type="scatterpolar",r=v,theta=names(v),fill="toself",marker=list(color="#D55E00"))%>%
      layout(polar=list(radialaxis=list(range=c(0,1))), title=input$g)
  })

  output$struct <- renderPlot({
    g <- gd()
    df <- data.frame(F=c("Gene(kb)/10","Exons","GC%","TE*100"),
                     V=c(g$gene_len_kb/10,g$n_exon,g$GC_content*100,g$gene_TE_All*100))
    ggplot(df,aes(F,V,fill=F))+geom_col(width=.5,alpha=.85)+
      scale_fill_manual(values=c("#D55E00","#E69F3B","#56B4E9","#009E73"))+
      labs(title=paste(input$g,"Structure"),x="",y="")+theme_minimal(13)+theme(legend.position="none")
  })

  output$species <- renderPlot({
    g <- gd()
    df <- data.frame(S=c("Human","Pig","Cattle","Chicken","Mouse"),
                     P=c(g$human_monopoly,g$pig_monopoly,g$cattle_monopoly,g$chicken_ortholog,g$mouse_ortholog))
    cls <- if(is.null(g$conservation_class)) "?" else g$conservation_class
    ggplot(df,aes(S,P,fill=P))+geom_col(width=.5,alpha=.85)+
      scale_fill_gradient(low="grey90",high="#D55E00")+
      labs(title=paste(input$g,"-",cls),x="",y="")+theme_minimal(13)+theme(legend.position="none")
  })

  output$tbl <- renderDT({
    datatable(atlas[,c("gene","MonopolyScore","CCLE_Cf","gene_len_kb","gene_TE_All","conservation_class")],
      rownames=FALSE,options=list(pageLength=15))
  })
}
shinyApp(ui,server)
