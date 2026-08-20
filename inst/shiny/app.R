library(shiny);library(ggplot2);library(plotly);library(DT)
atlas=as.data.frame(readRDS("data/demo_atlas.rds"));rownames(atlas)=atlas$gene
ui=fluidPage(titlePanel(h2("Transcriptomic Monopoly Atlas")),sidebarLayout(
sidebarPanel(selectizeInput("g","Gene:",sort(atlas$gene),selected="EEF1A1"),hr(),tableOutput("info")),
mainPanel(tabsetPanel(tabPanel("Radar",plotlyOutput("radar",height=480)),
tabPanel("Structure",plotOutput("struct",height=480)),tabPanel("Species",plotOutput("species",height=480)),
tabPanel("Table",DTOutput("tbl"))))))
server=function(input,output){
gd=reactive({req(input$g);atlas[input$g,]})
output$info=renderTable({g=gd();data.frame(Dim=c("MonopolyScore","CCLE%","Cancers","Length","GC%","TE","Class"),Val=c(sprintf("%.3f",g$MonopolyScore),sprintf("%.1f%%",g$CCLE_Cf*100),sprintf("%d/33",g$n_cancers_mono),sprintf("%.1fkb",g$gene_len_kb),sprintf("%.1f%%",g$GC_content*100),sprintf("%.3f",g$gene_TE_All),g$conservation_class))},rn=F,colnames=F)
output$radar=renderPlotly({g=gd();v=c(MS=g$MonopolyScore,CCLE=min(g$CCLE_Cf*5,1),Breadth=g$n_cancers_mono/33,Compact=1-pmin(g$gene_len_kb/100,1),LowTE=1-pmin(g$gene_TE_All*5,1),Cons=(g$human_monopoly+g$pig_monopoly+g$cattle_monopoly+g$chicken_ortholog+g$mouse_ortholog)/5);plot_ly(type="scatterpolar",r=v,theta=names(v),fill="toself",marker=list(color="#D55E00"))%>%layout(polar=list(radialaxis=list(range=c(0,1))),title=input$g)})
output$struct=renderPlot({g=gd();df=data.frame(F=c("Gene","Exon","Intron","nExon"),V=c(g$gene_len_kb,g$exon_total/1000,g$intron_total/1000,g$n_exon));ggplot(df,aes(F,V,fill=F))+geom_col(alpha=.85,width=.5)+scale_fill_manual(values=c("#D55E00","#E69F3B","#56B4E9","#009E73"))+labs(title=paste(input$g,"Structure"),x="")+theme_minimal(13)+theme(legend.position="none")})
output$species=renderPlot({g=gd();df=data.frame(S=c("Human","Pig","Cattle","Chicken","Mouse"),P=c(g$human_monopoly,g$pig_monopoly,g$cattle_monopoly,g$chicken_ortholog,g$mouse_ortholog));ggplot(df,aes(S,P,fill=P))+geom_col(alpha=.85,width=.5)+scale_fill_manual(values=c("TRUE"="#D55E00","FALSE"="grey80"))+labs(title=paste(input$g,"-",g$conservation_class),x="")+theme_minimal(13)+theme(legend.position="none")})
output$tbl=renderDT({datatable(atlas[,c("gene","MonopolyScore","CCLE_Cf","n_cancers_mono","gene_len_kb","gene_TE_All","conservation_class")],options=list(pageLength=15))})}
shinyApp(ui,server)
