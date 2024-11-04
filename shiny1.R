# --------------------------------------------------
# Shiny app for mining        
# --------------------------------------------------

require(shiny)
require(shinythemes)
require(ggplot2)
require(plotly)

load("data/DepMap_Q2.RData")

####################################
# User interface                   #
####################################
ui <- fluidPage(theme = shinytheme("readable"),
                navbarPage(
                  "DepMap",
                  
                # -------------------------------
                # Tab 1. Single gene DepMap CRISPR score
                # -------------------------------
                tabPanel("Single Gene Score",
                         sidebarPanel(
                           textInput("singlegene", "Search for a gene:", "")
                         ), 
                         
                         mainPanel(
                           h2("DepMap 22Q2, CRISPR screen, Achilles & SCORE"),
                           h3("Gene Effect score (Chronos)"),
                           plotOutput(outputId = "plot1.1", height = 1200)
                        )
                ), 
                
                # -------------------------------------
                # Tab 2. Protein complexes DepMap CRISPR score
                # -------------------------------------               
                tabPanel("Protein Complexes Score", 
                         selectInput("epifactors", "Choose an epigenetic protein complex (Epifactors):",
                                     list(`Epifactors epigenetics complex` = list("a", "b", "c")
                                      )),
                         selectInput("corum", "Choose a CORUM protein complex:",
                                     list(`Complex group name` = list("A", "B", "C")
                                     )),                         
                         ),
                
                # -------------------------------
                # Tab 3. Multi-gene DepMap CRISPR score
                # -------------------------------              
                tabPanel("Multi-gene Score", 
                         sidebarLayout(
                           sidebarPanel(
                             textAreaInput("multigene", "Enter a list of genes", "", 
                                           width = "200px"),
                             textOutput("multigeneout")
                           ),
                           
                         mainPanel(
                           tabsetPanel(
                             tabPanel("Plot", plotOutput("plot3.1")), 
                             tabPanel("Summary", verbatimTextOutput("tab3_empty")), 
                             tabPanel("Table", tableOutput(NULL))
                           )
                         )
                )),
                
                # -------------------------------
                # Tab 4. Cancer-specific score
                # -------------------------------              
                tabPanel("Cancer-specific score", 
                         sidebarLayout(
                           sidebarPanel(
                             textAreaInput("tab4_gene", "Enter gene(s)", "", 
                                           width = "200px"),
                             textOutput("tab4_gene_out"),
                             selectInput("Lineage_sub_subtype", "Choose a cancer sub-type:",
                                         list(`Cancer lineage sub-type` = list("Synovial sarcoma", "AML", "CML"))
                                         )
                           ),
                           
                           mainPanel(
                             tabsetPanel(
                               tabPanel("Plot", plotOutput("plot4.1")), 
                               tabPanel("Summary", verbatimTextOutput("tab4_empty")), 
                               tabPanel("Table", tableOutput(NULL))
                             )
                           )
                         ))
  )
)

####################################
# Server                           #
####################################
server <- function(input, output) {
  
  # -------------------------------
  # Tab 1. Single gene DepMap CRISPR score
  # ------------------------------- 
  
  # plot1.1
  output$plot1.1 <- renderPlot({
    crispr_eff <- unlist(geneeffect[input$singlegene, ])
    df1.1 <- data.frame(Sample = names(crispr_eff),
                     Lineage = sampleinfo$lineage[match(names(crispr_eff), sampleinfo$DepMap_ID)],
                     CRISPR_eff = crispr_eff)
    # sort by median
    medians <- aggregate(df1.1$CRISPR_eff, list(df1.1$Lineage), median)
    medians <- medians[order(medians$x), ]
    df1.1$Lineage <- factor(df1.1$Lineage, levels = rev(medians$Group.1))
    if (input$singlegene != ""){
    g <- ggplot(df1.1, aes(x=Lineage, y=CRISPR_eff, color=Lineage)) +
      geom_hline(yintercept=0, color = "gray70", size=1, linetype="dashed")+
      geom_boxplot(outlier.shape = NA, size=1)+
      geom_jitter(size=1.5, stroke=NA, alpha=0.60, aes(color=Lineage),position = position_jitter(width = .1)) +
      coord_flip() +  theme_bw() + 
        theme(legend.position = "None", axis.text=element_text(size=14),
              panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
      ylab("CRISPR Gene Effect (Chronos)") + xlab("")
    
    print(g)
    }
  })
  
  # -------------------------------
  # Tab 3. Multi-gene DepMap CRISPR score
  # -------------------------------

  # Side bar panel
    output$multigeneout <- renderText({ 
      multigene_fixed <- strsplit(input$multigene, "\n|,| ")[[1]]
      multigene_fixed <- multigene_fixed[multigene_fixed != ""]
    output <- paste(multigene_fixed, collapse = ", ") 
    paste0("Selected genes: ", output)
    })  

  # Plot tab
  output$plot3.1 <- renderPlot({
    plot(1:10)
  }) 
  
  # Summary tab
  output$tab3_empty <- renderText({ 
  paste("This panel is empty")
  }) 
  
  # Table tab
  
  # -------------------------------
  # Tab 4. Cancer-specific score
  # -------------------------------
  # Side bar panel
  output$tab4_gene_out <- renderText({ 
    tab4_gene_fixed <- strsplit(input$tab4_gene, "\n|,| ")[[1]]
    tab4_gene_fixed <- tab4_gene_fixed[tab4_gene_fixed != ""]
    output <- paste(tab4_gene_fixed, collapse = ",") 
    paste0("Selected genes: ", output)
  })  
  
  # Plot tab
  output$plot4.1 <- renderPlot({
    plot(1:10)
  }) 
  
  # Summary tab
  output$tab4_empty <- renderText({ 
    paste("This panel is empty")
  }) 
  
  # Table tab
}

####################################
# Create the shiny app             #
####################################
shinyApp(ui = ui, server = server)
