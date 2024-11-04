library(shiny)
library(plotly)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ComplexHeatmap)
library(ggpubr)
library(poolr)
library(reshape2)
library(rstatix)
library(gginnards)
library(GSVA)
library(edgeR)
library(limma)
library(Biobase)
library(GSVAdata)
library(EnhancedVolcano)
library(limma)



ui <- fluidPage(
  tabsetPanel(
    tabPanel(
      h4("SSGSEA Plot", style = "color: #800080;"),
      
      sidebarLayout(
        sidebarPanel(width = 3,
                     fileInput('file1', 'Upload Data File for SSGSEA Plots:',
                               accept=c('text/csv','text/comma-separated-values,text/plain')),
                     
                     checkboxInput('header', 'Data File has Variable Names as Column Headers.', TRUE),
                     
                     selectInput("type", "Cancer Subtype:",
                                 c("AML" = "AML",
                                   "Synovial" = "Synovial")),
                     
                     #Data file seperator 
                     
                     radioButtons('sep', 'Data File Separator Value:',
                                  c(Comma=',',
                                    Semicolon=';',
                                    Tab='\t')
                     )
        ),
        mainPanel(
          h4("SSGSEA Plot", align = "center"),
          plotlyOutput("volcanoPlot", height = "500px"),
          dataTableOutput("Table"),
          verbatimTextOutput("summary")
        )
        
      )
      
      
      
    )
    
    
  )
  
  
)

server <- function(input, output, session) {
  
  
  read1 <- reactive({
    inFile1 <- input$file1
    if (is.null(inFile1))
      return(NULL)
    
    corum_complexes <- read.csv(inFile1$datapath,
                                header=input$header,
                                na.strings = input$na.strings,
                                sep=input$sep)
    
    corum_complexes <- corum_complexes %>% as.data.frame()
    
    corum_complexes <- corum_complexes %>%
      dplyr::select("ComplexID", "ComplexName", "Gene_name") %>%
      as.data.frame()
    
    return(corum_complexes)
    
  })
  
  read2 <- reactive({
  
    list <- strsplit(as.character(read1()$Gene_name), ";")
    names(list) = read1()$ComplexID
    
    l2 = lapply(names(list), function(name){
      l2 = list[[name]] %>%
        as.vector()
      l2[!(l2 %in% Genes)]
      
    })
    l1 <- list[sapply(l2, function(x) length(x) >= 5)]
    return(l1)
    
    
  })
  
  
  read3 <- reactive({
    
    common_ess = read.csv("./CommonEssentials.csv")
    
    Genes = common_ess$gene
    
    depmap <- readRDS("./depmap_Q3_screen_data.rds")
    depmap_crispr_aml_nonaml = depmap$gene_fitness_crispr_depmap_q3
    
    ## Function to subset data 
    subset_object <- function(
    data_object, 
    meta_selection, 
    meta_id="DepMap_ID"){
      
      new_data_object <- list()
      new_data_object$genes <- data_object$genes
      new_data_object$meta_data <- data_object$meta_data[data_object$meta_data[[meta_id]] %in% meta_selection, ]
      new_data_object$mrna <- data_object$mrna[, meta_selection]
      new_data_object$scna <- data_object$snca[, meta_selection]
      new_data_object$gene_fitness_crispr_depmap_q3 <- data_object$gene_fitness_crispr_depmap_q3[, meta_selection]
      new_data_object$gene_fitness_rnai_depmap <- data_object$gene_fitness_rnai_depmap[, meta_selection]
      
      return(new_data_object)
      
    }
    
    ## Process AML data
    aml_cell_lines <- rownames(depmap$meta_data[grep(input$type, depmap$meta_data$subtype),])
    aml_depmap <- subset_object(depmap, meta_selection = aml_cell_lines)
    aml_crispr_fitness = aml_depmap$gene_fitness_crispr_depmap_21q3 %>% as.data.frame()
    
    aml_crispr_fitness = aml_crispr_fitness %>%
      dplyr::mutate(Gene = row.names(aml_crispr_fitness))%>%
      dplyr::filter(!(Gene %in% Genes))
    
    
    aml_crispr_fitness$Gene = NULL
    
    ## Process Non AML data
    non_aml_cell_lines <- rownames(depmap$meta_data[-grep(input$type, depmap$meta_data$subtype),])
    non_aml_depmap <- subset_object(depmap, meta_selection = non_aml_cell_lines)
    
    non_aml_crispr_fitness = non_aml_depmap$gene_fitness_crispr_depmap_21q3 %>% as.data.frame()
    
    
    non_aml_crispr_fitness = non_aml_crispr_fitness %>%
      dplyr::mutate(Gene = row.names(non_aml_crispr_fitness))%>%
      dplyr::filter(!(Gene %in% Genes))
    
    non_aml_crispr_fitness$Gene = NULL
    
    
    aml_gene_format = cbind(aml_crispr_fitness, non_aml_crispr_fitness)
    rownames(aml_gene_format) = rownames(aml_crispr_fitness)
    
    ##return(aml_gene_format)
    
    ##GSVA create expression set
    
    aml_gene_format_exprs = as.matrix(aml_gene_format)
    
    
    aml_gene_format_exprs[is.na(aml_gene_format_exprs)] <- 0
    
    
    pData = data.frame("cancer_type" = c(rep("AML",53), rep("nonAML", 1743)))
    rownames(pData) = colnames(aml_gene_format_exprs)
    
    phenoData <- new("AnnotatedDataFrame",
                     data=pData)
    
    aml_noaml_exprs_eset = ExpressionSet(assayData=aml_gene_format_exprs,
                                         phenoData=phenoData)
    
    aml_noaml_exprs_es <- gsva(aml_noaml_exprs_eset, read2(), method = "ssgsea")
    
    ##corum_complexes_names_matched = (read2()[match(rownames(aml_noaml_exprs_es), read2()$ComplexID),])$ComplexName
    
    
    ##mod <- model.matrix(~ factor(aml_noaml_exprs_es$cancer_type))
    mod1 <- model.matrix(~ factor(aml_noaml_exprs_es$cancer_type) -1)
    colnames(mod1) <- c("AML", "NonAML")
    fit2 <- lmFit(aml_noaml_exprs_es, mod1)
    contrast.matrix <- makeContrasts("AML-NonAML", levels = mod1)
    contrast.matrix
    fit2C <- contrasts.fit(fit2, contrast.matrix)
    fit2C <- eBayes(fit2C)
    ##topTable(fit2C)
    return(fit2C)
    
  })
  
  read4 <- reactive({
    
    res <- decideTests(read3(), p.value=0.10)
    
    return(res)
  })
  
  
  
  
  read5 <- reactive({
    
    tt <- topTable(fit2C, n=Inf)
    
    corum_complexes_names_matched = (read2()[match(rownames(tt), read2()$ComplexID),])$ComplexName
    
    
    tt$ComplexName = corum_complexes_names_matched
    tt$ComplexID = rownames(tt)
    return(tt)
    
  })
  
  
  
  
  read6 <- reactive({
    
    ## --------------------------
    ## Fix inputs for volcano plot
    ## --------------------------
    ## fold-change
    fc1 = read5()$logFC
    fc1[fc1 < -0.2] <- -0.2
    fc1[fc1 > 0.2] <- 0.2
    
    # adjusted p-value
    qval1 = read5()$adj.P.Val
    
    # gene names
    genes1 <- make.unique(rownames(read5()))
    
    # Create a data frame using gene names, fold changes, and q-values
    
    df1 <- data.frame(gene_name = as.character(genes1), lfc = fc1, q = qval1, complexName = read5()$ComplexName, stringsAsFactors = F)
    rownames(df1) <- genes1
    return(df1)
    
  })
  
  output$summary<- renderPrint({
    
    read6()
    
  })
  
  read7 <- reactive({  
    
    # adjusted p-value
    qval1 = read5()$adj.P.Val
    
    ## Plot volcano
    # y-axis upper limit
    max_pval1 <- -log10(min(qval1)) + 0.1
    
    
    
    
  })
  
  output$volcanoPlot <- renderPlotly({
    
    
    labels1 = (read5()$ComplexName)[1:15]
    plot <- EnhancedVolcano(read6(), x = 'lfc', y = 'q', lab = read6()$complexName,
                            pCutoff = 0.10, FCcutoff = 1.0,
                            gridlines.major = FALSE, gridlines.minor = FALSE,
                            drawConnectors = T, legendLabSize = 12, colAlpha = 0.75,
                            cutoffLineCol = "red", cutoffLineType = "dashed",
                            border = "full",
                            ##colCustom = cols1,
                            legendPosition = "right",
                            pointSize = 2, cutoffLineWidth = 0.4,
                            labFace = "plain", subtitle = "",
                            ylim = c(0, read7()), xlim = c(-0.4, 0.4),
                            axisLabSize = 12, captionLabSize = 12,
                            xlab = "SSGSEA enrichment score difference", ylab = "-Log10 adjusted p-value", title = "",
                            caption = paste0('Total = ', nrow(df1), ' genes'),typeConnectors = "closed", legendIconSize = 2,
                            selectLab = labels1,
                            borderWidth = 1.5, boxedLabels = TRUE)
  })
  
 
  
}

shinyApp(ui, server)










