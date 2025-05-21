#
# This is the UI for R shiny web application.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
source('functions.R')

#-----------------------------Libraries-----------------------------------------
#import necessary libraries
library(shiny)
library(bslib)
library(readxl)
library(plotly)
library(ggplot2)
library(colourpicker)
library(RColorBrewer)
library(DT)
library(igraph)
library(shinycssloaders)
library(shinyWidgets)
library(shinyjs)

required_ui_packages <- c("shiny", "bslib", "readxl", "plotly", "ggplot2", 
                          "colourpicker", "RColorBrewer", "DT", "igraph", 
                          "shinycssloaders", "shinyWidgets")
install_and_load(required_ui_packages)

#----------
bs_theme_custom <- bs_theme(
  bg = "#FFFFFF",         #overall background (light white)
  fg = "#333333",         #foreground/text color (dark gray)
  primary = "#1a73e8",    #primary color (used for buttons, etc.)
  base_font = font_google("Roboto"),  #customize font if needed
  "navbar-bg" = "#e3f2fd",    #navbar background color (light blue)
  "navbar-color" = "#333333", #navbar text color (dark gray)
  "navbar-link-color" = "#333333",  #navbar link color (dark gray)
  "navbar-link-hover-color" = "#1a73e8"  #link hover color (blue)
)
#--------------------------------UI--------------------------------------------
ui <- fluidPage(
  #enable Javascript functions
  useShinyjs(),
  #theme = bslib::bs_theme(bootswatch = "flatly"),
  theme = bs_theme_custom,
  #----------------------------- Title------------------------------------------
  #Application title and introductions
  tags$head(
    tags$title("RNA-Seq Analysis"), #window title
    tags$style(HTML("
      .title-container {
        display: flex;
        align-items: center;
        justify-content: space-between;
        width: 100%;
        padding: 20px;
      }
      .title-text {
        font-size: 36px; 
        font-weight: bold;
        text-align: center;
      }
      .title-img {
        height: 120px; 
        width: auto;
      }
    "))),
  
    div(class = "title-container",
        img(src = "rna_seq_logo.png", class = "title-img"),  
        h4("RNA-Seq Analysis", class = "title-text"),     
        img(src = "rna_logo.png", class = "title-img") 
    ),

  p("Developed by: Krupa Sampat"),

  #------------------------------Tabs-------------------------------------------
  #Nested Tabs
  navbarPage("Hi!",
             #---------------------------Home page---------------------------------------
             tabPanel(
               #Tab title 
               "Home Page", icon = icon("home"),
               fluidRow(
                 column(
                   width = 8,
                   #title = "Welcome to the RNA-Seq Analysis App", width = 12, status = "primary", solidHeader = TRUE,
                   h2("Welcome to the RNA-Seq Analysis App"),
                   br(),
                   h4("Introduction"),
                   p("This application allows users to perform RNA-Seq data analysis. 
                   You can perform differential expression and functional analysis using the app."),
                   br(),
                   h4("How to use this application:"),
                   tags$div(tags$ol(
                     tags$li("Use the 'Samples' tab to load and view sample metadata."),
                     tags$li("Use the 'Counts' tab to upload count data and perform exploratory data analysis."),
                     tags$li("The 'DE Analysis' tab allows for differential expression analysis using `DESeq2` package."),
                     tags$li("In the 'Functional Analysis' tab, perform gene ontology and pathway analysis."),
                     tags$li("The 'GSEA' tab is to perform the gene set enrichment analysis using `fgsea` package.")
                  )),
                   br(), br(),
                   h5("Citations:"),
                   tags$div(tags$ul(
                     tags$li("DESeq2 (v1.42.0): Love, M.I., Huber, W., Anders, S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2 Genome Biology 15(12):550 (2014)"),
                     tags$li("ClusterProfiler (v4.10.0): T Wu, E Hu, S Xu, M Chen, P Guo, Z Dai, T Feng, L Zhou, W Tang, L Zhan, X Fu, S Liu, X Bo, and G Yu. clusterProfiler 4.0: A
                     universal enrichment tool for interpreting omics data. The Innovation. 2021, 2(3):100141"),
                     tags$li("fGsea (v1.28.0), G. Korotkevich, V. Sukhov, A. Sergushichev. Fast gene set enrichment analysis. bioRxiv (2019), doi:10.1101/060012")
                  ))
                ),
                column(
                  width = 4,
                  align ='right',
                  tags$img(src="rnaseq_image.png", height = '500px', width = '500px')
                ))
             ),
             #------------------------------------------------------------------
             tabPanel("Analyze", icon = icon("wand-magic-sparkles"),
                      tabsetPanel(
                        #---------------------------Sample------------------------------------------
                        tabPanel(
                          #tab title
                          "Samples", icon = icon("table"),
                          sidebarLayout(
                            #side bar panel
                            sidebarPanel(
                              #load sample data (metadata)
                              fileInput("data", "Load sample metadata", accept = c(".xlsx", ".csv"),
                                        buttonLabel = "Browse...", placeholder = "No file selected"), 
                              tags$div(tags$p("Note: Ensure that the metadata contains the following columns:"),
                                       tags$ul(
                                         tags$li("`sample_id`: eg: RNA0000 (ideally one used for downloading the FASTQ files)"),
                                         tags$li("`replicate`/`date`/`batch`"),
                                         tags$li("`condition`: Ensure this contains ONLY the condition/treatment group; if details (such as concentration, etc) are required, add another column named `condition_description`"),
                                         tags$li("`cell_type`: The cell line used for study"),
                                         tags$li("`ercc_spike_in`: The ERCC spike-in mix (1 or 2) used; if details (such as concentration, etc) are required, add another column named `ercc_spikein_description`")
                                       )),
                              
                              #select columns to subset from metadata
                              uiOutput("select_columns"),
                              br(),
                              #select the column to be used as "condition" column
                              uiOutput("select_conditionCol"),
                              #select columns to plot and group by for histogram 
                              # uiOutput("histogram_columns"),
                            ),
                            #       ),
                            #main panel
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Metadata", 
                                         br(), 
                                         uiOutput("summary_table_title"),
                                         br(), 
                                         tableOutput("summary_table"), 
                                         br(), 
                                         uiOutput("metadata_title"), br(),
                                         dataTableOutput("data_table")),
                                #tabPanel("Table", br(), uiOutput('metadata_title), br(), 
                                #dataTableOutput("data_table")),
                                #tabPanel("Plot", "Please wait a few seconds for the plot to load...", plotly::plotlyOutput("data_plot"))
                              )
                            )
                          )),
                        
                        #---------------------------Counts------------------------------------------
                        tabPanel(
                          #tab title
                          "Counts", icon = icon("chart-bar"),
                          sidebarLayout(
                            #side bar panel
                            sidebarPanel(
                              #load the counts file
                              fileInput("counts_file", "Load the count data file", accept = c(".txt", ".csv"),
                                        buttonLabel = "Browse...", placeholder = "No file selected"),
                              br(),
                              
                              #select samples for analysis
                              uiOutput("sample_selector"),
                            ),
                            
                            #main panel
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Counts Summary",
                                         br(),
                                         uiOutput("counts_summary_title"), br(),
                                         tableOutput("counts_summary_table"), br(),
                                         dataTableOutput("metadata_subset"), br(),
                                         dataTableOutput("counts_subset")),
                                tabPanel("Diagnostic Plots",
                                         tabsetPanel(tabPanel("Exploratory Data Analysis (EDA)", "Please wait for a few seconds for the plot to load...", 
                                                              withSpinner(plotly::plotlyOutput("library_size", height="700px"), type = 4), 
                                                              uiOutput('downloadButton_libSize'), 
                                                              br(),
                                                              withSpinner(plotly::plotlyOutput("sample_distribution_boxplot", height="700px"), type = 4),
                                                              uiOutput('downloadButton_sampleDist'),
                                                              br()
                                                              ),
                                                     tabPanel("Heatmap", "Please wait for a few seconds for the plot to load...", 
                                                              withSpinner(plotOutput("dist_heatmap", height="700px"), type = 5), 
                                                              br(), uiOutput('downloadButton_dist_heatmap'), br()
                                                              ), 
                                                     #withSpinner(plotOutput("gene_clust_heatmap"))),
                                                     tabPanel("Multi-dimensional Scaling (MDS) Plot", "Please wait for a few seconds for the plot to load...", 
                                                              withSpinner(plotOutput("mds_plot")),
                                                              uiOutput('downloadButton_mds_plot')
                                                              ))),
                                tabPanel("Principal Component Analysis (PCA)", 
                                         tabsetPanel(tabPanel("PCA", "Select the PCs to plot and click the 'Plot PCA' button. 
                                          Please wait for a few seconds for the plot to load...",
                                                              sidebarLayout(
                                                                sidebarPanel(
                                                                  #Select which principal components to plot - by default PC1 and PC2 are selected
                                                                  selectInput("pc_x", "X axis principal component to plot:",
                                                                              choices = paste0("PC", 1:10), selected="PC1"),
                                                                  selectInput("pc_y", "Y axis principal component to plot:",
                                                                              choices = paste0("PC", 1:10), selected = "PC2"),
                                                                  br(),
                                                                  actionButton("plot_pca", "Plot PCA"),
                                                                  br(), br(), #),
                                                                p("To generate PCA using `plotPCA()` function from DESeq2: "),
                                                                actionButton("deseq2_pca","Plot PCA (DESeq2)")),
                                                                
                                                                mainPanel(plotOutput("pca_plot", height="900px", width="900px"),
                                                                          uiOutput('downloadButton_pca_plot'), #)
                                                                          br(), br(),
                                                                          plotOutput("deseq2_pcaPlot", height="900px", width="900px"),
                                                                          uiOutput('downloadButton_deseq2_pcaPlot'))
                                                              )),
                                                     tabPanel("Normalized Counts PCA", "This PCA is generated using DESeq2 normalized counts. 
                                                              Please wait for a few seconds for the plot to load...", 
                                                              withSpinner(plotOutput("norm_pca", height="900px", width="1100px")),
                                                              uiOutput('downloadButton_norm_pca'))
                                         )),
                                tabPanel("Normalized Counts",
                                         sidebarLayout(
                                           sidebarPanel(
                                             #select the normalization method
                                             actionButton("help_norm", label = icon("question-circle"), class="btn btn-info"),
                                             radioButtons("norm_method", "Choose Normalization Method:",
                                                          choices = c("CPM" = "CPM", "FPKM/RPKM" = "FPKM", "TPM" = "TPM", "DESeq2"="DESeq2"#), 
                                                          ), selected=NULL),
                                             #choice to round the counts
                                             checkboxInput("round_normCounts", "Round Normalized Counts", value = FALSE),
                                             conditionalPanel(condition= "input.round_normCounts == true",
                                                              numericInput("round_digits", "Decimal Places:", value=2, min=0, max=10)),
                                             actionButton("normalize", "Normalize"),
                                             uiOutput("downloadButtonUI")
                                           ),
                                           mainPanel("Normalized Counts", withSpinner(dataTableOutput("normCounts")))
                                         )),
                                tabPanel("ERCC Spike in Analysis", "Click the 'Analyze ERCC spike ins to generate the plot.
                     Please wait for a few seconds for the table and plots to load...", 
                                         sidebarLayout(
                                           sidebarPanel(
                                             textAreaInput('ercc_dilution', "Please enter the ERCC dilution factor used: (e.g. 0.02 for 2uL in 100uL (i.e. 1:100 dilution))", 
                                                           value='0.02'),
                                             br(),
                                             #select type of ERCC plot
                                             radioButtons('ercc_plot_type', "Select ERCC plot type", 
                                                          choices = c("Sample-wise" = 'sample_wise_ercc',
                                                                      "Entire ERCC dataset" = 'entire_ercc_dataset')),
                                             radioButtons('ercc_mix_type', 'Select ERCC Mix Type: ', 
                                                          choices = c("Mix 1" = 'mix1', "Mix 2" = 'mix2'), selected = 'mix1'),
                                             br(),
                                             actionButton("generate_ercc_plot", "Analyze ERCC spike ins"),
                                             ),
                                           mainPanel(
                                             plotly::plotlyOutput("ercc_plot", width="1000px", height="900px"), 
                                             br(),
                                             uiOutput('downloadButton_ercc_plot'), br(),
                                             verbatimTextOutput("ercc_model_stats"), br(), 
                                             dataTableOutput("ercc_table"),
                                             #plotOutput('ercc_hist_l2fc')
                                           )
                                         )
                                )
                              )
                            )
                          )
                        ),
                        
                        #-------------------------------DE------------------------------------------
                        tabPanel(
                          #tab title
                          "DE Analysis", icon = icon("flask"),
                          sidebarLayout(
                            #side bar Panel
                            sidebarPanel(
                              #help for design for deseq2
                              actionButton("help_design", label = 'Help?', class="btn btn-info", icon = icon("question-circle")),
                              #input for design formula 
                              textInput("design_formula", "Enter the columns to use in the Design formula for DESeq2", value="~ batch + condition"),
                              #slider input - filter no. of genes - by default 10 is selected
                              sliderInput("filter_slider", "Select the minimum number of counts required for a gene to be included in the analysis. Genes having total counts below this threshold will be filtered out.",
                                          min = 0, max = 100, value = 10, step = 10),
                              #checkboxInput - filter genes in no. of samples (minimum no. of samples that need the counts to be >= genes threshold)
                              checkboxInput('filter_sample', "Include minimum no. of samples that contain the above threshold?", value = TRUE),
                              conditionalPanel(
                                condition = 'input.filter_sample = true',
                                uiOutput('sample_filter')
                              ),
                              #select main factor for comparisons:
                              uiOutput('factor_selector'),
                              #set reference for comparisons
                              uiOutput("select_ref"),
                              #check if ercc spike in was used
                              checkboxInput("ercc_used", "Include ERCC spike in analysis?", value = TRUE),
                              #
                              conditionalPanel(condition = "input.ercc_used == false",
                                               radioButtons("use_genes_dds", "Which genes to use in DESeq2?",
                                                            choices = c("Use all genes (including ERCC spike-in genes)?" = "use_all",
                                                                        "Use genes (excluding ERCC spike-in genes)?" = "use_wo_ercc"),
                                                            selected = "use_wo_ercc")),
                              #allow user to set contrasts for DEG results
                              uiOutput("contrasts_selector"),
                              actionButton("run_deseq2", "Normalize Data with DESeq2")
                            ),
                            
                            #main panel
                            mainPanel(
                              tabsetPanel(
                                #tab for plots
                                tabPanel("Diagnostic Plots", "Please wait for a few seconds for the plots to load...",
                                         withSpinner(plotOutput("disp_plot", height = "700px", width = "900px")), 
                                         br(), 
                                         withSpinner(plotOutput("meansd_rld", height = "700px", width = "900px")), 
                                         br(), 
                                         withSpinner(plotOutput("meansd_vsd", height = "700px", width = "900px"))),
                                tabPanel("Generate DEG Results & Summary", 
                                         sidebarLayout(
                                           sidebarPanel(
                                             #select p adj value
                                             radioButtons("alpha_cutoff", "Select the adjusted p value (FDR) cutoff: ", 
                                                          choices = c(0.05, 0.01), selected = 0.01),
                                             #select l2fc value
                                             #sliderInput("fc_slider", "Select the log2 fold change threshold (for getting significant DEGs): ", 
                                            #             min=0, max=5, step = 0.5, value = 1),
                                            numericInput("fc_cutoff", "Input the log2 fold change threshold (for getting significant DEGs):",
                                                         value=1, step = 0.001),
                                             #select conditions/contrasts to compare
                                             p("Select the 'contrasts'/comparisons to be used for generating results. 
                            The order of the contrasts changes the interpretability of results."),
                                             uiOutput("contrast1_ui"),
                                             uiOutput("contrast2_ui"),
                                             actionButton("generate_deg_res", "Generate DEG Results")
                                           ),
                                           #main panel
                                           mainPanel(withSpinner(dataTableOutput("deg_summary_table")))
                                         )
                                ),
                                tabPanel("DEG Results Table", actionBttn("help_deg_res", label = "", style = "fill", icon = icon("question-circle")),
                                         dataTableOutput("deg_table"), uiOutput('downloadButton_degres')),
                                tabPanel("Volcano Plot & MA plot", 
                                         sidebarPanel(
                                           #choice to annotate the volcano plot
                                           p('Annotate volcano plot with top DEGs?'),
                                           shinyWidgets::prettyToggle('annotate_volcano', label_on = "Yes!", icon_on = icon("thumbs-up"), 
                                                                      label_off = "No..", icon_off = icon("thumbs-down"), outline = TRUE, plain=TRUE),
                                           conditionalPanel(
                                             condition = 'input.annotate_volcano == true',
                                             sliderInput("top_genes_cutoff", "Select the number of DEGs (top) to display: ", min = 0, max = 100, step = 5, value = 10)
                                           ),
                                           actionButton("vol_plot", "Generate Volcano plot"),
                                           p("^--- Click to generate the volcano plot!")),
                                         #main panel
                                         mainPanel(
                                           conditionalPanel(
                                             condition = 'input.vol_plot > 0',
                                             withSpinner(plotlyOutput("volcano_plot", height = "900px", width = "1200px"), type=5)), 
                                           br(), 
                                                   uiOutput('downloadButton_vp'), br(),
                                                   #p("MA plot"),
                                                   plotOutput('ma_plot', height = "700px"),
                                                   uiOutput('downloadButton_ma_plot'))
                                ),
                                tabPanel("Significant Gene Lists", p("Signficant Up-Regulated Genes"), dataTableOutput("up_regTable"),
                                         br(), br(),
                                         p("Significant Down-Regulated Genes"), dataTableOutput("down_regTable"), br(), uiOutput('download_geneLists'))
                              )
                            )
                          )),
                        
                        #---------------------------Functional Analysis-----------------------------
                        tabPanel(
                          #tab title
                          "Functional Analysis", icon = icon("dna"),
                          sidebarLayout(
                            sidebarPanel(
                              #upload gene lists
                              checkboxInput('upload_geneLists', 'Upload Gene lists?', value = FALSE),
                              conditionalPanel(
                                condition = 'input.upload_geneLists == true',
                                fileInput('genelist_file', "Please load a file containing genes for functional analysis.", accept = c(".csv", ".txt", ".xlsx")),
                                p("The gene list should contain a column with Ensembl ids and another with Gene symbols.")
                              ),
                              #or use pre-generated genelists
                              checkboxInput('generated_geneLists', 'Use the Gene lists generated through `DE analysis` tab? (only the significant genes will be used for this)', value = FALSE),
                              conditionalPanel(
                                condition = 'input.generated_geneLists == true',
                                radioButtons('geneList_to_use', 'Which gene list to use?', 
                                             choices = c('Up-regulated Gene list' = 'use_upregList',
                                                         'Down-regulated Gene List' = 'use_downregList'))
                              ),
                              #select the functional analysis methods
                              selectInput('func_analysis_types', 'Select which type(s) of functional analysis would you like to perform:',
                                          choices = c('Gene ontology - Biological Processes' = 'go_bp', 'Gene ontology - Molecular Functions' = 'go_mf', 
                                                      'Gene ontology - Cellular Components' = 'go_cc', 'KEGG' = 'kegg'), 
                                          multiple = TRUE, selectize = TRUE),
                              #no. of categories to display in the plots
                              numericInput('show_categories', 'Number of categories to display?', value = 10),
                              #plot options
                              radioButtons('func_plot_type', 'Choose the Plot type (Visualization):',
                                           choices = c('Bar plot' = 'bar', 'Dot plot' = 'dot', 'Enrichment Plot' = 'emap')),
                              sliderInput("fontSize", "Font size for plots:", min = 8, max = 24, value = 12),
                              actionButton('run_func_analysis', 'Perform Functional Analysis')
                            ),
                            #main panel
                            mainPanel( uiOutput('func_plots_ui')
                              #tabsetPanel(
                                #tabPanel("Plot", "Please wait for the analysis to run and the plots to load..", withSpinner(plotOutput("func_plots", height = '2500px'))),
                                #tabPanel("Plots", "Please wait for the analysis to complete and the plots to load..",
                                        # uiOutput('func_plots_ui')),
                                #tabPanel("Table", withSpinner(uiOutput('dynamic_func_tabs'))) #dataTableOutput("func_table"))
                              #)
                            )
                          )),
                        
                        #----------------------------------GSEA-------------------------------------
                        tabPanel(
                          #tab title
                          "GSEA", icon = icon("database"),
                          sidebarLayout(
                            sidebarPanel(
                              #upload gene sets
                              radioButtons("data_source", "Choose data source:",
                                           choices = c("Upload Gene set?" = "upload_geneSet",
                                                       "Use the DESeq2 results generated from `DE Analysis` tab?" = "use_res"),
                                           selected = "use_res"),
                              conditionalPanel(
                                condition = "input.data_source == 'upload_geneSet'",
                                fileInput('gsea_data', 'Upload Gene Enrichment Data', accept = c(".csv", ".tsv", ".xlsx")),
                                p("Please ensure that the uploaded data contains the following columns: `Geneid` and `log2FoldChange`")
                              ),
                              selectInput("species", "Select Species:",
                                          choices = c("Homo sapiens", "Mus musculus")),
                              radioButtons("geneId_type", "Gene Id type used in your data:", 
                                           choices = c("Gene symbol", "Ensembl"), selected = "Ensembl"),
                              actionBttn("help_gsea", label = 'Help?', style = "fill", icon = icon("question-circle")),
                              br(), br(),
                              selectInput("gsea_category", "Select Category:", 
                                          choices = c("Hallmark Gene sets (H)" = "H",
                                                      "Positional Gene sets (C1)" = "C1",
                                                      "Curated Gene sets (C2)" = "C2",
                                                      "Regulatory Target Gene sets (C3)" = "C3",
                                                      "Computational Gene sets (C4)" = "C4",
                                                      "GO Gene sets (C5)" = "C5",
                                                      "Oncogenic signatures (C6)" = "C6",
                                                      "Immunogenic signatures (C7)" = "C7",
                                                      "Cell type signatures (C8)" = "C8")),
                              sliderInput("top_cat", "Top number of categories to plot", min = 5, max = 50, value = 10),
                              actionButton('run_gsea', 'Run GSEA')
                            ),
                            #main panel
                            mainPanel(
                              tabsetPanel(
                                tabPanel("Plot", "Please wait for the analysis to run and the plots to load..", 
                                         br(),
                                         uiOutput("downloadButton_gseaPlot"),
                                         br(),
                                         withSpinner(plotOutput("gsea_plot"))),
                                tabPanel("Table", 
                                         actionBttn("help_gsea_res", label = "", style = "fill", icon = icon("question-circle")), 
                                         withSpinner(dataTableOutput('gsea_table')), uiOutput("downloadButton_gseaRes")) 
                              )
                            )
                          )
                        )
                      ))
  ))