# hTF analysis user interface

# initialize environment
rm(list=ls())
options(scipen=999)
options(digits=22)

# load analysis libraries
library(ggplot2)
library(ggseqlogo)
library(gridExtra)

# load shiny app libraries
library(shiny)
library(shinythemes)
library(DT)

# scale transformation function for rounding axes
scaleFUN <- function(x) sprintf("%.1f", x)

# function to obtain SNV value matrix
zscore_to_SNV_matrix <- function(PBM, PBM_col, target_seq_col, target_seq) {
  
  # initialize a matrix to hold the SNV probe z-score values
  nucs <- c("A", "C", "G", "T")
  SNV <- matrix(nrow=4, ncol=nchar(as.character(target_seq)))
  
  # collect z-score values for each position in the matrix
  for (i in 1:nrow(SNV)) {
    for (j in 1:ncol(SNV)) {
      # construct the current sequence of interest using the current SNV
      SNV_seq <- paste(substr(target_seq, 1, j-1), nucs[i], substr(target_seq, j+1, ncol(SNV)), sep="")
      stopifnot(SNV_seq %in% as.character(PBM[, target_seq_col]))
      
      # collect the z-score for that given SNV probe sequence and condition
      SNV[i, j] <- PBM[which(as.character(PBM[, target_seq_col])==SNV_seq), PBM_col]
    }
  }
  
  # return the SNV energy matrix to the parent function
  rownames(SNV) <- nucs
  return(SNV)
}

# function to transform a SNV energy matrix using the column-wise median to plot as energy logo
SNV_to_SNV_trans <- function(SNV) {
  
  # transform z-score values to reflect deviation from column-wise median
  SNV_trans <- matrix(nrow=4, ncol=ncol(SNV))
  rownames(SNV_trans) <- c("A", "C", "G", "T")
  for (i in 1:nrow(SNV_trans)) {
    for (j in 1:ncol(SNV_trans)) {
      SNV_trans[i, j] <- SNV[i, j] - median(SNV[,j])
    }
  }
  
  # return the PWM matrix to be plotted as a sequence logo
  return(SNV_trans)
}

# declare version-independent annotation columns
annot_cols <- c("TF_name", "ID", "family", "TF_class", "cons_seq")

# shiny ui section
ui <- navbarPage("hTF array analysis (v1.0.0)",
                 
  # wrap a layout theme around the app
  theme = shinytheme("flatly"),
  
  # "Data explorer" tab panel
  tabPanel("Data explorer",
           sidebarLayout(
             
             # sidebar panel to help user select their data
             sidebarPanel(
               
               # user-selected file input
               fileInput("hTF_file", "Upload an hTF data file",
                         multiple = FALSE,
                         accept = c("text/plain",
                                    ".dat")),
               width = 2),
             
             # main panel for data exploration
             mainPanel(
               tabsetPanel(type = "tabs",
                           tabPanel("consensus probes only", DT::dataTableOutput("datatable_hTF")),
                           tabPanel("consensus + single variants")),
               width = 10)
             )
  ),
  
  # "Scatterplots" tab panel
  tabPanel("Scatterplots",
           sidebarLayout(
             sidebarPanel(
               
               # user-selected experiment comparison
               uiOutput("selectScatterY"),
               uiOutput("selectScatterX"),
               
               # user-selected TF group to compare
               uiOutput("selectTFGroup")
               
             ),
             mainPanel(
               tabsetPanel(type = "tabs",
                           tabPanel("TF-level comparison", plotOutput("scatter1")),
                           tabPanel("Family-level comparison", plotOutput("scatter2")),
                           tabPanel("Class-level comparison", plotOutput("scatter3")))
             )
           )
  ),
  
  # "Motifs" tab panel
  tabPanel("Motif grid",
           sidebarLayout(
             sidebarPanel(
               
               # user-selected PBM experiments to compare
               uiOutput("selectMotifExp"),

               # user-selected TF group to compare
               uiOutput("selectMotifTF")
               
             ),
             mainPanel(
               tabsetPanel(type = "tabs",
                           tabPanel("Energy logos", plotOutput("motif1", height="auto")),
                           tabPanel("Position-weight logos"))
             )
           )
  ),
  
  # "Heatmaps" tab panel
  tabPanel("Recruitment heatmap",
           sidebarLayout(
             sidebarPanel(
               # TODO
             ),
             mainPanel(
               tabsetPanel(type = "tabs",
                           tabPanel("Tab 1"),
                           tabPanel("Tab 2"),
                           tabPanel("Tab 3"))
             )
           )
  )
  
)

# shiny server section
server <- function(input, output) {
  
  ### READ INPUT DATA FILE AND PREPROCESS DESIGN-LEVEL VARIABLES
  
  # enforce file size limit of 30mb
  options(shiny.maxRequestSize=30*1024^2)
  
  # function to build a reactive data table from the user-specified file input
  df_hTF <- eventReactive(input$hTF_file, {
    read.table(input$hTF_file$datapath,
               header = T,
               sep = "\t",
               stringsAsFactors = F)
  })
  
  # collect names of PBM experiments performed in current data file
  PBM_exp <- reactive({
    names(df_hTF())[grepl("br", names(df_hTF()))]
  })
  
  # collect names of TFs in the current design
  TF_choices <- reactive({
    TF_temp <- unique(df_hTF()[, "TF_name"])
    TF_temp[!is.na(TF_temp)]
  })
  
  ### DATA TABLE CONSTRUCTION FOR EXPLORER
  
  # construct a consensus probe datatable to display (seed probes only, rounded z-scores)
  cons_dtab_hTF <- reactive({
    cons_dtab_temp <- df_hTF()[which(df_hTF()$SNV_pos_offset==0), c(annot_cols, PBM_exp())]
    cons_dtab_temp[,-1*(1:length(annot_cols))] <- round(cons_dtab_temp[,-1*(1:length(annot_cols))], digits = 3)
    names(cons_dtab_temp) <- gsub("v[0-9]_a[0-9]_run[0-9]_br_", "", names(cons_dtab_temp))
    cons_dtab_temp
  })
  
  # output the consensus probe data table
  output$datatable_hTF <- DT::renderDataTable({
    cons_dtab_hTF()
  })
  
  ### CONDITIONAL SCATTERPLOT UI OPTIONS
  
  # user-specified y variable
  output$selectScatterY <- renderUI({
    selectInput("y_hTF", "Select Y variable", PBM_exp())
  })
  
  # user-specified x variable
  output$selectScatterX <- renderUI({
    selectInput("x_hTF", "Select X variable", PBM_exp())
  })
  
  # user-specified TF group
  output$selectTFGroup <- renderUI({
    selectInput("TFs_hTF", label="Choose TFs to compare", TF_choices(), multiple=T)
  })
  
  ### BUILD SCATTERPLOTS
  
  # build dataframe for scatterplot
  df_scatter <- reactive({
    df_hTF()[which(df_hTF()$TF_name %in% input$TFs_hTF), c(input$y_hTF, input$x_hTF, "TF_name", "family", "TF_class")]
  })
  
  # construct scatterplot to output (TF-level)
  output$scatter1 <- renderPlot({
    ggplot(df_scatter(), aes_string(y=input$y_hTF, x=input$x_hTF, color="TF_name")) +
      geom_point(size=3) +
      scale_x_continuous(labels=scaleFUN) + scale_y_continuous(labels=scaleFUN) +
      ylab(input$y_hTF) + xlab(input$x_hTF) +
      labs(color = "TF name") +
      theme_classic(base_size = 16) +
      theme(legend.title = element_text(size=20)) +
      theme(legend.text = element_text(size=18))
  })
  
  # construct scatterplot to output (Family-level)
  output$scatter2 <- renderPlot({
    ggplot(df_scatter(), aes_string(y=input$y_hTF, x=input$x_hTF, color="family")) +
      geom_point(size=3) +
      scale_x_continuous(labels=scaleFUN) + scale_y_continuous(labels=scaleFUN) +
      ylab(input$y_hTF) + xlab(input$x_hTF) +
      labs(color = "TF family") +
      theme_classic(base_size = 16) +
      theme(legend.title = element_text(size=20)) +
      theme(legend.text = element_text(size=14))
  })
  
  # construct scatterplot to output (Class-level)
  output$scatter3 <- renderPlot({
    ggplot(df_scatter(), aes_string(y=input$y_hTF, x=input$x_hTF, color="TF_class")) +
      geom_point(size=3) +
      scale_x_continuous(labels=scaleFUN) + scale_y_continuous(labels=scaleFUN) +
      ylab(input$y_hTF) + xlab(input$x_hTF) +
      labs(color = "TF class") +
      theme_classic(base_size = 16) +
      theme(legend.title = element_text(size=20)) +
      theme(legend.text = element_text(size=14))
  })
  
  ### CONDITIONAL MOTIF GRID UI OPTIONS
  
  # user-specified PBM experiment group
  output$selectMotifExp <- renderUI({
    selectInput("PBM_exp_motif", "Choose experiments to compare", PBM_exp(), multiple = T)
  })
  
  # user-specified TF group
  output$selectMotifTF <- renderUI({
    selectInput("TFs_motif", label="Choose TFs to compare", TF_choices(), multiple=T)
  })
  
  ### PROCESS AND PLOT MOTIF GRID
  
  # expand out the combinations of TFs and PBM experiments
  TF_PBM_comb <- eventReactive({input$PBM_exp_motif
    input$TFs_motif}, {
    TF_PBM_comb_temp <- expand.grid(input$PBM_exp_motif, input$TFs_motif, stringsAsFactors = F)
    names(TF_PBM_comb_temp) <- c("PBM_exp", "TFs")
    TF_PBM_comb_temp
  })
  
  # collect a motif plot list
  motif_list <- reactive({
    lapply(seq(nrow(TF_PBM_comb())), function(m) {
      
      # determine seed sequence for the current TF
      seed_seq <- df_hTF()[which(df_hTF()$TF_name==TF_PBM_comb()[m, "TFs"] & df_hTF()$SNV_pos_offset==0), "target_seq"]
      
      # determine the seed score for the current TF, PBM_exp combo
      seed_score <- df_hTF()[which(df_hTF()$TF_name==TF_PBM_comb()[m, "TFs"] & df_hTF()$SNV_pos_offset==0), TF_PBM_comb()[m, "PBM_exp"]]
      
      # collect SNV matrix for the current combination
      SNV_mat <- zscore_to_SNV_matrix(df_hTF()[which(df_hTF()$TF_name==TF_PBM_comb()[m, "TFs"]), ], TF_PBM_comb()[m, "PBM_exp"], "target_seq", seed_seq)
      
      # transform the SNV matrix
      SNV_mat_trans <- SNV_to_SNV_trans(SNV_mat)

      # plot the transformed SNV matrix
      ggp <- ggseqlogo(SNV_mat_trans, method="custom", seq_type="dna")
      my.ggp.yrange <- ggplot_build(ggp)$layout$panel_params[[1]]$y.range
      ggp <- ggp +
        ggtitle(paste(TF_PBM_comb()[m, "TFs"], TF_PBM_comb()[m, "PBM_exp"], round(seed_score, digits = 3), sep="\n")) +
        scale_x_continuous(expand = c(0.01,0.01)) +
        scale_y_continuous(expand = c(0.01,0.01), labels=scaleFUN) +
        annotate('rect', xmin = 0.5, xmax = ncol(SNV_mat_trans)+0.5, ymin = my.ggp.yrange[1], ymax = 0.0, alpha = 0.75, col='white', fill='white')
        
      # place a transparent window over the logo if the seed score is below threshold
      ggp <- ggp +
        if (seed_score < 1.5) {
          annotate('rect', xmin = 0.5, xmax = ncol(SNV_mat_trans)+0.5, ymin = my.ggp.yrange[1], ymax= my.ggp.yrange[2], alpha = 0.75, col='white', fill='white')
        }  
      
      # finixh plotting axis labels and borders
      ggp <- ggp +
        annotate('rect', xmin = 0.5, xmax = ncol(SNV_mat_trans)+0.5, ymin = 0, ymax = my.ggp.yrange[2], col="lightgray", fill=NA, size=1.5) +
        annotate('rect', xmin = 0.5, xmax = ncol(SNV_mat_trans)+0.5, ymin = my.ggp.yrange[1], ymax = 0, col="lightgray", fill=NA, size=1.5) +
        ylab(expression(paste(Delta, "z-score", sep=""))) +
        theme(axis.text.x=element_blank()) +
        theme(axis.title.x=element_blank()) +
        theme(axis.text.y=element_text(size=10, face="plain"), axis.title.y=element_text(size=10, face="plain")) +
        theme(plot.title = element_text(size=10))
      
      # return the full plot
      ggp
      
    })
  })
  
  # declare grid height function (150px x number of TFs)
  grid_height_fn <- eventReactive({input$PBM_exp_motif
    input$TFs_motif}, {
    150 * length(input$TFs_motif)
  })
    
  # output a motif grid
  output$motif1 <- renderPlot({
    do.call(gridExtra::grid.arrange, c(motif_list(), ncol=length(input$PBM_exp_motif)))
  }, height=grid_height_fn)
  
}

# call shiny app using ui and server
shinyApp(ui, server)
