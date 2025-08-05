# This is the server function for R shiny web application.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

source('functions.R')

#list of required packages
required_packages <- c("tidyverse", "DESeq2", "AnnotationDbi", "org.Hs.eg.db", 
                       "pheatmap", "magrittr", "GenomicFeatures", "vsn", 
                       "viridis", "genefilter", "EnsDb.Hsapiens.v86", 
                       "clusterProfiler", "enrichplot", "patchwork", "fgsea",
                       "RUVSeq")
install_and_load(required_packages)

set.seed(627)
#-----------------------------Libraries-----------------------------------------
#import necessary libraries
library(tidyverse)
library(readxl)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(RUVSeq)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)
library(magrittr)
library(GenomicFeatures)
library(vsn)
library(viridis)
library(genefilter)
library(EnsDb.Hsapiens.v86)
library(clusterProfiler)
library(enrichplot)
library(patchwork)
library(fgsea)
library(grid)
library(gridGraphics)

#-------------------------------SERVER------------------------------------------
server <- function(input, output, session) {
  #exceed the max size for file inputs
  options(shiny.maxRequestSize = 70*1024^2)
  #------------------------------Tabs-------------------------------------------
  showDownloadButton <- reactiveValues()
  #---------------------------Sample--------------------------------------------
  #TAB 1: SAMPLES
  #load sample data
  
  load_metadata <- reactive({
    req(input$data)
    ext <- tools::file_ext(input$data$name)
    
    metadata <- switch(ext,
                       "csv" = read.csv(input$data$datapath, header=TRUE),
                       "xlsx" = read_excel(input$data$datapath),
                       "xls" = read_excel(input$data$datapath),
                       "tsv" = readr::read_tsv(input$data$datapath, col_names = TRUE),
                       validate('Unsupported file type. Please upload a `.csv` or `.xlsx`, `.xls`, or `.tsv` file.'))
    
    #format colnames
    colnames(metadata) <- tolower(gsub("_+$", "", gsub(" +", "_", gsub("[.#]", "", colnames(metadata)))))
    return(metadata)
  })
  
  selected_sampleId_column <- reactiveVal(NULL)
  selected_condition_column <- reactiveVal(NULL)
  update_metadata <- reactive({
    metadata <- load_metadata()
    req(metadata)
    
    if (!is.null(selected_sampleId_column())) {
      names(metadata)[names(metadata) == selected_sampleId_column()] <- "sample_id"
    }
    
    rename_cols <- c("ercc_exfold_spike_in_mix" = "ercc_exfold_spike_in",
                     "ercc_spike-in_mix" = "ercc_spike_in",
                     "conditions" = "condition_description",
                     "1:100_diluted_ercc_spike-in_(ul)" = "ercc_spike_in")
    names(metadata) <- dplyr::recode(names(metadata), !!!rename_cols)
    
    #factor the conditions column
    if ("condition_description" %in% colnames(metadata)) {
      metadata <- metadata %>%
        dplyr::mutate(condition = case_when(
          grepl("0.01% DMSO", condition_description) & grepl("not split.", condition_description) ~ "DMSOnotsplit",
          grepl("0.01% DMSO", condition_description) & grepl("split at", condition_description) ~ "DMSOsplit",
          grepl("500 nM VTP-50649", condition_description) & grepl("not split.", condition_description) ~ "VTPnotsplit",
          grepl("500 nM VTP-50649", metadata$condition_description) & grepl("split at", metadata$condition_description) ~ "VTPsplit",
          grepl("0.01% DMSO", condition_description) ~ "DMSO", #for TC-797
          TRUE ~ "VTP" #for TC-797
        ))
    } else if ("condition" %in% colnames(metadata)) {
      metadata <- metadata %>%
        mutate(condition = gsub("VTP50469", "VTP", condition))
    }
    
    if ("ercc_spike_in_mix" %in% colnames(metadata)) {
      metadata <- metadata %>%
        mutate(ercc_spike_in_mix = gsub("ERCC Spike-in mix 1", "1", ercc_spike_in_mix))
    }
    
    if (!"batch" %in% colnames(metadata)) {
      #metadata$date <- as.Date(metadata$date, format = '%m%d%Y')
      if ("date" %in% colnames(metadata)) {
        uniqueDates <- unique(metadata$date)
        metadata$batch <- as.factor(as.integer(factor(metadata$date, levels=uniqueDates)))
      }}
    
    print(colnames(metadata))
    
    return(metadata)
    
  })
  
  observe({
    metadata <- load_metadata()  # Load the metadata
    req(metadata)  # Ensure metadata exists
    
    if (!"sample_id" %in% colnames(metadata)) {
      # Show a modal if 'sample_id' is missing
      showModal(modalDialog(
        title = "Select Sample-ID Column",
        "It looks like 'sample_id' is not present in the metadata. Please select the appropriate column that represents the sample IDs.",
        selectInput('sampleId_column', "Select column:", choices = c("", colnames(metadata)), selected = NULL, multiple = FALSE),
        footer = modalButton("Confirm")
      ))
    }
  })
  
  observeEvent(input$sampleId_column, {
    
    selected_sampleId_column(input$sampleId_column)
    removeModal()
  })
  
  observeEvent(input$selected_conditionCol, {
    selected_condition_column(input$selected_conditionCol)
  })
  
  #subset the metadata
  subset_metadata <- reactive({
    #metadata <- load_metadata()
    metadata <- update_metadata()
    req(metadata)
    if (is.null(metadata)) {
      showNotification("No data uploaded yet.", type='error')
      return(NULL)
    }
    
    #subset the required columns from metadata
    selected_columns <- input$selected_columns
    
    #add batches to metadata subset according to the dates
    metadata_subset <- metadata %>%
      dplyr::select(!!selected_columns)
    
    print(head(metadata_subset))
    
    return(metadata_subset)
  })
  
  #----------------------------Output-------------------------------------------
  #get columns to subset 
  output$select_columns <- renderUI({
    req(input$data, update_metadata())
    #load metadata and get colnames
    columns <- names(update_metadata())
    
    #update the CheckboxGroupInput choices accordingly
    checkboxGroupInput("selected_columns", "Select columns to Include in the analysis", choices = columns, selected = NULL)
  })
  
  #select the "condition" column
  output$select_conditionCol <- renderUI({
    req(update_metadata())
    metadata <- update_metadata()
    selectInput("selected_conditionCol", "Please select the column to be used as 'condition'/'treatment group'/'to make comparisons'.",
                choices = colnames(metadata), selected = "condition")
  })
  
  # #get columns to plot and group by for histogram
  # output$histogram_columns <- renderUI({
  #   req(subset_metadata())
  #   metadata_subset <- subset_metadata()
  #   tagList(
  #   selectInput("columns_to_plot", "Select column for histogram", choices = names(metadata_subset)),
  #   selectInput("column_to_groupby", "Select column to Group by", choices = names(metadata_subset)))
  # })
  
  #output for `summary_table`
  output$summary_table_title <- renderUI({
    req(subset_metadata())
    h5('Summary of metadata')
  })
  
  output$summary_table <- renderTable({
    req(subset_metadata())
    metadata_subset <- subset_metadata()
    summary_df <- data.frame(Column_name = names(metadata_subset),
                             Type = sapply(metadata_subset, class))
    
    return(summary_df)
  })
  
  #output for `data_table`
  output$data_table <- renderDataTable({
    metadata_subset <- subset_metadata()
    metadata_subset <- datatable(metadata_subset, options = list(paging=TRUE, scrollX=TRUE), 
                                 caption = "Metadata Table (entire) with selected columns")
    if (is.null(metadata_subset)) {
      return("No data uploaded yet")
    }
    isolate({
      return(metadata_subset)
    })
  })
  
  # #output for histogram
  # output$data_plot <- renderPlotly({
  #   req(input$columns_to_plot, input$column_to_groupby)
  #   metadata_subset <- subset_metadata()
  #   
  #   plot_histogram(metadata_subset, input$columns_to_plot, input$column_to_groupby)
  # })
  
  #---------------------------Counts--------------------------------------------
  #TAB 2: COUNTS
  #load counts data
  load_counts_data <- reactive({
    req(input$counts_file)
    
    #check if meatadata is uploaded
    if (is.null(input$data)) {
      showNotification("Metadata file missing!!! Please upload the metadata file first.", type='error', duration = NULL)
      return(NULL)
    }
    
    metadata_subset <- subset_metadata()
    
    counts_data <- read_counts(file = input$counts_file, ext = tools::file_ext(input$counts_file$name))
    
    print("Initial counts data structure: ")
    print(str(counts_data))
    
    #remove unnecessary columns from count data
    #counts_data <- counts_data[, -c(2:6)]
    unnecessary_cols <- c("Chr", "Start", "End", "Strand", "Length")
    counts_data <- counts_data[, !colnames(counts_data) %in% unnecessary_cols]
    
    #rename colnames according to sample IDs - if necessary
    pattern <- "RNA[0-9]+"
    
    if (!all(grepl(pattern, colnames(counts_data)))) {
      extracted_ids <- sapply(colnames(counts_data)[-1], function(name){regmatches(name, regexpr(pattern, name))})
      print(extracted_ids)
      
      extracted_ids <- unlist(extracted_ids)
      if(length(extracted_ids) != length(colnames(counts_data)[-1])) {
        stop("Error in extracting or matching sample IDs")
      }
      
      colnames(counts_data)[-1] <- extracted_ids
    }
    
    head(colnames(counts_data))
    
    counts_data <- na.omit(counts_data)
    
    #error message if sample_id column is not selected from metadata
    if (!("sample_id" %in% colnames(metadata_subset))) {
      msg <- "The `sample_id` column is missing from metadata. Please select all the required columns."
      showNotification(msg, type='error', duration = NULL)
      return()
    } else if (!any(colnames(counts_data)[-1] %in% metadata_subset$sample_id)) {
      #check if all samples are present in the metadata
      missing_ids <- colnames(counts_data)[-1][!colnames(counts_data)[-1] %in% metadata_subset$sample_id]
      showNotification(paste0("The following sample IDs from counts data are missing in metadata: ", paste(missing_ids, collapse = ", ")), type = 'error', duration = NULL)
      return(NULL)
    }
    
    print(head(counts_data))
    
    return(counts_data)
  })
  
  #get gene lengths from featureCounts data - for normalization
  get_gene_lengths <- reactive({
    counts_data <- read_counts(file = input$counts_file, ext = tools::file_ext(input$counts_file$name))
    
    if (!"Length" %in% colnames(counts_data)) {
      showNotification("Gene Length column not found in the uploaded file. Please upload a counts file which contains gene length (named as `Length`) column.", type='error')
      return(NULL)
    }
    
    #extract gene lengths from count data as a vector
    gene_lengths <- as.numeric(counts_data[["Length"]])
    
    print(head(gene_lengths))
    print(str(gene_lengths))
    
    return(gene_lengths)
  })
  
  #summary table
  counts_summary_table <- reactive({
    req(load_counts_data(), load_metadata(), filtered_metadata(), filtered_counts_sub())
    
    counts <- load_counts_data()
    metadata <- load_metadata()
    metadata_sub <- filtered_metadata()
    counts_sub <- filtered_counts_sub()
    
    if(is.null(counts)){
      return("No data uploaded yet")
    }
    
    num_samples_counts <- ncol(counts) - 1 #Assuming 1st column is GeneId
    total_genes <- nrow(counts)
    num_samples <- nrow(metadata)
    num_selected_samples <- nrow(metadata_sub)
    
    counts_summary_table <- data.frame(Measure = c("Number of Samples in Metadata", "Number of Samples in Counts data", "Number of genes", "Number of Selected samples"),
                                       Values = c(num_samples, num_samples_counts, total_genes, num_selected_samples))
    
    return(counts_summary_table)
  })
  
  #convert count data to matrix and remove geneid column for DESeq2
  get_counts_matrix <- reactive({
    req(filtered_counts_sub())
    counts_sub <- filtered_counts_sub()
    
    counts_matrix <- as.matrix(counts_sub[, -which(colnames(counts_sub) == "Geneid")])
    rownames(counts_matrix) <- counts_sub$Geneid
    
    return(counts_matrix)
  })
  
  # ----- RUVSeq ---
  #filter counts data for ruvseq
  filter_ruv <- reactive({
    req(input$ruv_norm)
    req(filtered_counts_sub(), filtered_metadata(), selected_condition_column(), input$ruv_filter_slider, input$ruv_sample_slider)
    
    counts_sub <- subset(filtered_counts_sub())
    metadata_sub <- filtered_metadata()
    conditionCol <- selected_condition_column()
    counts_to_filter <- as.numeric(input$ruv_filter_slider)
    samples_counts_filter <- as.numeric(input$ruv_sample_slider)
    
    metadata_sub[[conditionCol]] <- factor(metadata_sub[[conditionCol]])
    #filter out non-expressing genes - 
    filter_ruv <- apply(counts_sub[-1], 1, function(x) length(x [ x > counts_to_filter]) >= samples_counts_filter)
    filtered_counts_ruv <- counts_sub[filter_ruv, ]
    
    print(head(filtered_counts_ruv))
    
    #set rownames as Geneid column
    rownames(filtered_counts_ruv) <- filtered_counts_ruv$Geneid
    
    return(filtered_counts_ruv)
  })
  
  names_ruv <- reactive({ 
    req(filter_ruv())
    
    ruv_filtered_counts <- filter_ruv()
    column <- if ('Geneid' %in% colnames(ruv_filtered_counts)) {
      print("Using Geneid column-")
      ruv_filtered_counts$Geneid
    } else {
      print("Using rownames-")
      rownames(ruv_filtered_counts)
    }
    
    #extract gene ids 
    genes_ruv <- column[grepl("ENSG", column)]
    print(head(genes_ruv))
    
    #extract spike in gene id
    spikes_ruv <- column[grepl("ERCC", column)]
    print(head(spikes_ruv))
    
    #return list
    list(genes = genes_ruv, spikes = spikes_ruv)
  })
  
  prepare_ruv_obj <- reactive({
    req(filter_ruv(), names_ruv(), filtered_metadata(), selected_condition_column())
    
    ruv_filtered_counts <- filter_ruv()
    ruv_matrix <- as.matrix(ruv_filtered_counts[-1])
    rownames(ruv_matrix) <- ruv_filtered_counts$Geneid
    metadata_sub <- filtered_metadata()
    sample_col <- metadata_sub$sample_id
    condition_col <- metadata_sub[[selected_condition_column()]]
    set_ruv <- newSeqExpressionSet(as.matrix(ruv_matrix),
                                   phenoData = data.frame(condition = condition_col, 
                                                          row.names = sample_col))
    print(set_ruv)
    
    return(set_ruv)
  })
  
  run_ruvseq <- reactive({
    req(names_ruv(), prepare_ruv_obj(), input$ruv_k)
    
    set_ruv <- prepare_ruv_obj()
    spikes <- names_ruv()$spikes
    
    set_ruv <- betweenLaneNormalization(set_ruv, which = "upper")
    
    set1_ruv <- RUVg(set_ruv, spikes, k = input$ruv_k)
    print(pData(set1_ruv))
    
    return(set1_ruv)
  })
  
  observeEvent(input$ruv_norm, {
    req(filter_ruv(), names_ruv(), prepare_ruv_obj(), run_ruvseq())
    print("RUV normalization triggered")
  })
  

  #library size plot
  get_library_sizes <- reactive({
    req(get_counts_matrix())
    counts_matrix <- get_counts_matrix()
    
    filtered_counts <- counts_matrix[rowSums(counts_matrix) >= 10, ]
    
    library_size <- colSums(filtered_counts)
    names(library_size) <- colnames(filtered_counts)
    
    print(library_size)
  })
  
  #PCA
  generate_pca_plot <- eventReactive(input$plot_pca, {
    req(filtered_counts_sub(), filtered_metadata(), input$pc_x, input$pc_y, selected_condition_column())
    
    counts_sub <- filtered_counts_sub()
    metadata_sub <- filtered_metadata()
    conditionCol <- selected_condition_column()
    
    pca_results <- prcomp(t(counts_sub[-1]), center=TRUE)
    
    if (!conditionCol %in% colnames(metadata_sub)) {
      showNotification("The selected `condition` column is missing from the metadata. Please select all required columns before proceeding.")
      return()
    }
    
    conditions <- unique(metadata_sub[[conditionCol]])
    
    color_mapping <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(conditions)), conditions)
    
    pca_plot <- make_pca(pca_results, metadata_sub, conditionCol, color_mapping, 
                         "PCA Plot", input$pc_x, input$pc_y)
    
    return(pca_plot)
  })

  #normalized counts PCA 
  generate_norm_pca <- reactive({
    req(filtered_metadata(), get_deseq2_norm_counts(), input$pc_x, input$pc_y, selected_condition_column())
    norm_counts <- get_deseq2_norm_counts()
    metadata_sub <- filtered_metadata()
    conditionCol <- selected_condition_column()
    
    #identify genes with 0 variance
    zero_var <- apply(norm_counts, 1, function(x) var(x, na.rm = TRUE) == 0)
    #filter norm counts matrix and remove 0 variance genes
    norm_counts <- norm_counts[!zero_var, ]
    
    if (nrow(norm_counts) < 2) {
      showNotification("Not enough variable genes to perform PCA. Try tweaking parameters to get normalized counts again.", type = "error", duration = NULL)
      return(NULL)
    }
    
    pca_results <- prcomp(t(norm_counts), scale. = TRUE)
    
    if (!conditionCol %in% colnames(metadata_sub)) {
      showNotification("The selected `condition` column is missing from the metadata. Please select all required columns before proceeding.")
      return()
    }
    
    conditions <- unique(metadata_sub[[conditionCol]])
    color_mapping <- setNames(colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(length(conditions)), conditions)
    
    pca_plot <- make_pca(pca_results, metadata_sub, conditionCol, color_mapping, 
                         "PCA Plot", input$pc_x, input$pc_y)
    
    return(pca_plot)
  })
  
  #initialize reactive values for metadata and counts subset
  filtered_metadata <- reactiveVal()
  filtered_counts_sub <- reactiveVal()
  #initialize reactive value for normalized counts
  normalized_counts <- reactiveVal()
  #showDownloadButton <- reactiveVal(FALSE)
  showDownloadButton$norm_counts <- FALSE
  
  #get normalized counts according to the selected method
  observeEvent(input$normalize, {
    req(load_counts_data(), get_gene_lengths(), filtered_counts_sub(), input$norm_method)
    counts_data <- load_counts_data()
    counts_sub <- filtered_counts_sub()
    gene_lengths <- get_gene_lengths()
    
    norm_counts_data <- switch(input$norm_method,
                               "CPM" = calculate_cpm(counts_sub[, -which(colnames(counts_sub) == "Geneid")]),
                               "FPKM" = calculate_fpkm(counts_sub[, -which(colnames(counts_sub) == "Geneid")], gene_lengths),
                               "TPM" = calculate_tpm(counts_sub[, -which(colnames(counts_sub) == "Geneid")], gene_lengths),
                                "DESeq2" = {if (!deseq2_run()) {
                                  showModal(modalDialog(title='Error', "Please run DESeq2 from the DE Analysis Tab first to display the DESeq2 normalized counts.", easyClose = TRUE, size='m'))
                                  return(NULL)
                                } else {
                                  req(get_deseq2_norm_counts())
                                  norm_counts_data <- get_deseq2_norm_counts()
                                }
                                }
    )
    
   
    #round the normalized counts if user chooses to 
    if (input$round_normCounts) {
      norm_counts_data <- round(norm_counts_data, digits = input$round_digits)
    }
    
    #covert to data frame
    norm_counts_data <- as.data.frame(norm_counts_data)
    
    #add Gene id column
    if (input$norm_method == "DESeq2") {
      norm_counts_data$Geneid <- rownames(norm_counts_data)
      norm_counts_data <- norm_counts_data[, c("Geneid", setdiff(colnames(norm_counts_data), "Geneid"))]
    } else {
      norm_counts_data <- cbind(Geneid = counts_data$Geneid, norm_counts_data)
    }
    
    #add gene symbol column
    norm_counts_data <- get_gene_symbols(norm_counts_data)
    
    norm_counts_data <- norm_counts_data[, c("Geneid", "symbol", setdiff(colnames(norm_counts_data), c("Geneid", "symbol")))]
    
    print(head(norm_counts_data))
    
    #update reactive value
    normalized_counts(norm_counts_data)
    
    #update and display download button
   isolate({showDownloadButton$norm_counts <- TRUE}) 
  })
  
  #ERCC Spike in Analysis 
  #filter ercc counts from counts data
  get_ercc_counts <- reactive({
    req(filtered_counts_sub())
    counts_sub <- filtered_counts_sub()
    
    rownames(counts_sub) <- counts_sub$Geneid
    ercc_counts <- counts_sub[grep("ERCC", rownames(counts_sub)), ]
    if (nrow(ercc_counts) == 0) stop("No ERCC spike-ins found in the count data.")
    
    return(ercc_counts)
  })
  
  #use the functions `prepare_ercc_info()` and `prepare_combined_ercc_data()` to generate data for the dose response curve plot 
  run_ercc_analysis <- reactive({
    req(get_gene_lengths(), get_ercc_counts(), filtered_metadata(), filtered_counts_sub())
    counts_sub <- filtered_counts_sub()
    row.names(counts_sub) <- counts_sub$Geneid
    gene_lengths <- get_gene_lengths()
    ercc_data <- get_ercc_counts()
    metadata_sub <- filtered_metadata()
    conditionCol <- selected_condition_column()
    
    #calculate fpkm for ERCC
    #fpkm_spike <- calculate_fpkm(ercc_data[-1], gene_lengths)
    fpkm <- calculate_fpkm(counts_sub[-1], gene_lengths)
    row.names(fpkm) <- counts_sub$Geneid
    fpkm_spike <- fpkm[grep('ERCC', rownames(fpkm)), ]
    head(fpkm_spike)
    
    ercc_file_name <- "ercc_analysis.txt" 
    erccPath <- "/mount/eagen/krupa/rnaseq/AnalysisR/"
    ercc_file_path <- file.path(erccPath, ercc_file_name)
    #ercc_file_path <- "/mount/eagen/krupa/rnaseq/AnalysisR/ercc_analysis.txt"
    dilution_factor <- as.numeric(input$ercc_dilution)
    ercc_info <- prepare_ercc_info(ercc_file_path, dilution_factor)
    
    combined_data <- prepare_combined_ercc_data(fpkm_spike, ercc_data, ercc_info, metadata_sub, conditionCol)
    return(combined_data)
  })
  
  #get ercc linear model statistics using the data generated by `run_ercc_analysis()`
  ercc_model_stats <- reactive({
    req(run_ercc_analysis())
    
    plot_data <- run_ercc_analysis()
    stats <- list()
    
    mix_type <- input$ercc_mix_type
    
    if (mix_type == 'mix1') {
      stats <- calculate_stats(plot_data, "transcript_molecules_mix1")
    } else if (mix_type == 'mix2') {
      stats <- calculate_stats(plot_data, "transcript_molecules_mix2")
    } else if (mix_type == 'both_mixes') {
      stats$mix1 <- calculate_stats(plot_data, "transcript_molecules_mix1")
      stats$mix2 <- calculate_stats(plot_data, "transcript_molecules_mix2")
    }
    
    return(stats)
  })
  
  #functions to download the plots
  createDownloadUI <- function(plot_id, output, input, session, showDownloadButton) {
    dropdown_id <- paste0("dropdown_", plot_id)
    filename_id <- paste0("filename_", plot_id)
    filetype_id <- paste0("filetype_", plot_id)
    width_id <- paste0("plot_width_", plot_id)
    height_id <- paste0("plot_height_", plot_id)
    download_id <- paste0("download_", plot_id)
    close_id <- paste0("close_", plot_id)
    
    # Create download dropdown button
    output[[paste0("downloadButton_", plot_id)]] <- renderUI({
      req(showDownloadButton[[plot_id]])
 
      dropdownButton(
        inputId = dropdown_id,  # Assign an ID for JS control
        circle = FALSE,
        icon = icon("download"),
        status = "primary",
        size = "sm",
        tooltip = tooltipOptions(title = "Download Options"),
        inline = TRUE,
        
        textInput(filename_id, "File name:", plot_id),
        radioButtons(filetype_id, "File type:", choices = c('PNG' = 'png', 'PDF' = 'pdf'), selected = 'pdf'),
        numericInput(width_id, "Width (inches):", value = 8, min = 1, max = 30),
        numericInput(height_id, "Height (inches):", value = 6, min = 1, max = 30),
        downloadButton(download_id, "Download"),
        
        div(style = "display: flex; justify-content: space-between;",
            actionButton(close_id, "Close", class = "btn btn-danger")
        )
      )
    })
    
    # Close dropdown when close is clicked
    observeEvent(input[[close_id]], {
      runjs("$('.dropdown-menu').removeClass('show');")  # Hides the dropdown
    })
    
    return(list(download_id = download_id, filename_id = filename_id, filetype_id = filetype_id,
                width_id = width_id, height_id = height_id))
  }

  # #function to create a download handler for the plots
  createDownloadHandler_plots <- function(plot_id, plot_function, data_function, output, input, session, ids) {
    output[[ids$download_id]] <- downloadHandler(
      filename = function() {
        paste0(input[[ids$filename_id]], ".", input[[ids$filetype_id]])
      },
      content = function(file) {
        args <- data_function()
        
        #check if the function returns a static (ggplot) object 
        maybe_gg <- if (is.list(args)) do.call(plot_function, args) else plot_function(args, plotly = FALSE)
  
        if (inherits(maybe_gg, "ggplot")) {
          ggsave(file, maybe_gg, width = input[[ids$width_id]], height = input[[ids$height_id]], dpi = 300, units = "in", device   = input[[ids$filetype_id]])

          #otherwise assume it’s base-graphics: open device, draw, close
        } else {
          w    <- input[[ids$width_id]]
          h    <- input[[ids$height_id]]
          type <- input[[ids$filetype_id]]
          
          if (type == "png") {
            png(file, width = w, height = h, units = "in", res = 300)
          } else {
            pdf(file, width = w, height = h)
          }
          
          if (is.list(args)) do.call(plot_function, args) else plot_function(args, plotly = FALSE)
          
          dev.off()
        }
      }
    )
  }

  
  #----------------------------Output-------------------------------------------
  #output for `sample_selector`
  output$sample_selector <- renderUI({
    req(subset_metadata(), load_counts_data())
    metadata_subset <- subset_metadata()
    counts <- load_counts_data()
    
    all_samples <- str_sort(union(metadata_subset$sample_id, colnames(counts[-1])))
    
    #update and allow user to select samples
    tagList(checkboxGroupInput("selected_samples", "Select samples to include for analysis", choices = all_samples, selected = all_samples),
            actionButton("submit_samples", "Submit Selected Samples"))
  })
  
  output$ruvseq_ui <- renderUI({
    req(input$submit_samples)
    
    tagList(
      #RUVg normalization for ERCC spike ins
      checkboxInput("ruvseq", "Normalize data using RUVSeq?", value = TRUE),
      conditionalPanel(condition="input.ruvseq == true", 
                       div(
                         sliderInput("ruv_filter_slider", "Select the minimum number of counts required for a gene to be included in the analysis. Genes having total counts below this threshold will be filtered out.",
                                     min = 0, max = 100, value = 5, step = 5),
                         checkboxInput('ruv_filter_sample', "Include minimum no. of samples that contain the above threshold?", value = TRUE),
                         conditionalPanel(
                           condition = 'input.ruv_filter_sample = true',
                           uiOutput('ruv_sample_filter')),
                         p("Using RUVg method to normalize ERCC spike in genes."),
                         numericInput("ruv_k", "Number of unwanted factors (k) to correct for:", value = 1),
                         br(),
                         actionButton("ruv_norm", "Normalize - RUVseq")
                       )))
    
  })
  
  observeEvent(input$submit_samples, {
    req(input$submit_samples)
    metadata_subset <- subset_metadata()
    counts <- load_counts_data()
    
    #error message if samples are not present in count data or metadata
    #sample IDs missing in counts
    missing_samples_in_counts <- setdiff(input$selected_samples, colnames(counts[-1]))
    if (length(missing_samples_in_counts) > 0) {
      msg <- paste("The following sample IDs are missing in counts file: ", paste(missing_samples_in_counts, collapse = ", "), ". Please upload a file containing all the desired samples.")
      showNotification(msg, type='error', duration = NULL)
      return() 
    }
    
    #sample IDs missing in metadata
    missing_samples_in_metadata <- setdiff(input$selected_samples, metadata_subset$sample_id)
    if (length(missing_samples_in_metadata) > 0) {
      msg <- paste("The following sample IDs are missing in metadata:", paste(missing_samples_in_metadata, collapse = ", "), ". Please upload a file containing all the desired samples.")
      showNotification(msg, type = "error", duration = NULL)
      return()
    }
    
    #filter counts and metadata according to the selected samples
    metadata_sub <- metadata_subset[metadata_subset$sample_id %in% input$selected_samples, ]
    counts_sub <- counts[, c("Geneid", input$selected_samples)]
    
    #update reactive values
    filtered_metadata(metadata_sub)
    filtered_counts_sub(counts_sub)
    
    print(nrow(counts_sub))
  })
  
  #output for filtered metadata subset 
  
  output$metadata_subset <- renderDataTable({
    req(filtered_metadata())
    metadata_sub <- datatable(filtered_metadata(), options = list(paging=TRUE, scrollX=TRUE), 
                              caption = "Filtered metadata for selected samples")
    return(metadata_sub)
  })
  
  #output for filtered counts subset
  output$counts_subset <- renderDataTable({
    req(filtered_counts_sub())
    counts_sub <- datatable(filtered_counts_sub(), options = list(paging=TRUE, scrollX=TRUE), 
                            caption = "Raw Counts Table for selected samples") 
    return(counts_sub)
  })
  
  #output for sample filter
  output$ruv_sample_filter <- renderUI({
    req(input$filter_sample, filtered_metadata())
    sliderInput('ruv_sample_slider', 'Select the minimum number of samples that should pass the above threshold.', 
                min = 0, max = length(filtered_metadata()$sample_id), value = 2, step = 1)
  })
  
  #output for `summary_table`
  output$counts_summary_title <- renderUI({
    req(counts_summary_table())
    h5("Counts Summary Table")
  })
  
  output$counts_summary_table <- renderTable({
    req(counts_summary_table())
    counts_table <- counts_summary_table()
    
    if (!is.null(counts_table)){
      return(counts_table)
    } else {
      return("No data uploaded yet.")
    }
  })
  
  #output for `library_size` plot
  showDownloadButton$libSize <- FALSE
  output$library_size <- renderPlotly({
    req(get_library_sizes(), filtered_metadata())
    library_sizes <- get_library_sizes()
    metadata_subset <- filtered_metadata()
    library_size_plot <- make_libSize_plot(library_sizes, metadata = metadata_subset, plotly=TRUE)
    showDownloadButton$libSize <- TRUE
    return(library_size_plot)
  })
  
  libSize_ids <- createDownloadUI("libSize", output, input, session, showDownloadButton)
  createDownloadHandler_plots("libSize",  make_libSize_plot,
                              data_function = function() {
                                list(library_sizes = get_library_sizes(), metadata = filtered_metadata())
                              }, output, input, session, libSize_ids)
  
  #output for `sample_distribution` boxplot
  showDownloadButton$sampleDist <- FALSE
  output$sample_distribution_boxplot <- renderPlotly({
    req(get_counts_matrix())
    counts_matrix <- get_counts_matrix()
    sample_dist_plot <- plot_sample_distributions(counts_matrix, plotly=TRUE)
    showDownloadButton$sampleDist <- TRUE
    return(sample_dist_plot)
  })
  
  sampleDist_ids <- createDownloadUI("sampleDist", output, input, session, showDownloadButton)
  createDownloadHandler_plots("sampleDist", plot_sample_distributions, 
                              data_function = get_counts_matrix, output, input, session, sampleDist_ids)
  
  #output for `rle_before_ruv` boxplot
  showDownloadButton$rle_before_ruv <- FALSE
  output$rle_before_ruv <- renderPlot({
    req(prepare_ruv_obj(), selected_condition_column(), filtered_metadata())
    set_ruv <- prepare_ruv_obj()
    metadata_sub <- filtered_metadata()
    conditionCol <- selected_condition_column()
    
    RLE_before <- make_rlePlot(data = set_ruv, metadata = metadata_sub, condition = conditionCol)
    showDownloadButton$rle_before_ruv <- TRUE
    return(RLE_before)
  })
  
  rle_before_ids <- createDownloadUI("rle_before_ruv", output, input, session, showDownloadButton)
  createDownloadHandler_plots("rle_before_ruv", make_rlePlot, 
                              data_function = function() {
                                list(data = prepare_ruv_obj(), metadata = filtered_metadata(), condition = selected_condition_column()) 
                                }, output, input, session, rle_before_ids)
  
  #output for `rle_after_ruv` boxplot
  showDownloadButton$rle_after_ruv <- FALSE
  output$rle_after_ruv <- renderPlot({
    req(run_ruvseq(), selected_condition_column(), filtered_metadata())
    set1_ruv <- run_ruvseq()
    metadata_sub <- filtered_metadata()
    conditionCol <- selected_condition_column()
    
    RLE_after <- make_rlePlot(data = set1_ruv, metadata = metadata_sub, condition = conditionCol)
    showDownloadButton$rle_after_ruv <- TRUE
    return(RLE_after)
  })
  
  rle_after_ids <- createDownloadUI("rle_after_ruv", output, input, session, showDownloadButton)
  createDownloadHandler_plots("rle_after_ruv", make_rlePlot, 
                              data_function = function() {
                                list(data = run_ruvseq(), metadata = filtered_metadata(), condition = selected_condition_column()) 
                              }, output, input, session, rle_after_ids)
  
  
  #output for distance heatmap
  showDownloadButton$dist_heatmap <- FALSE
  output$dist_heatmap <- renderPlot({
    if (!deseq2_run()) {
      showModal(modalDialog(title='Error', "Please run DESeq2 from the DE Analysis Tab first to display the distance heatmap.", easyClose = TRUE, size='m'))
    } else {
      req(prep_dds(), get_vsd_counts(), selected_condition_column())
      vsd <- get_vsd_counts()
      conditionCol <- selected_condition_column()
      condition <- filtered_metadata()[[conditionCol]]
      dist_heatmap <- make_distHeatmap(vsd, condition)
      showDownloadButton$dist_heatmap <- TRUE
      return(dist_heatmap)
    }
  })
  
  dist_heatmap_ids <- createDownloadUI("dist_heatmap", output, input, session, showDownloadButton)
  createDownloadHandler_plots("dist_heatmap", make_distHeatmap, 
                              data_function = function() {
                                list(vsd = get_vsd_counts(), condition = filtered_metadata()[[selected_condition_column()]])
                              }, output, input, session, dist_heatmap_ids)
  
  #output for MDS plot
  showDownloadButton$mds_plot <- FALSE
  output$mds_plot <- renderPlot({
    if (!deseq2_run()) {
      showModal(modalDialog(title='Error', "Please run DESeq2 from the DE Analysis Tab first to display the MDS plot.", easyClose = TRUE, size='m'))
    } else {
      req(prep_dds(), get_vsd_counts(), selected_condition_column())
      vsd <- get_vsd_counts()
      conditionCol <- selected_condition_column()
      conditions <- filtered_metadata()[[conditionCol]]
      #condition <- filtered_metadata()$condition
      mds_plot <- make_mdsplot(vsd, conditions)
      showDownloadButton$mds_plot <- TRUE
      return(mds_plot)
    }
  })
  
  mds_plot_ids <- createDownloadUI('mds_plot', output, input, session, showDownloadButton)
  createDownloadHandler_plots('mds_plot', make_mdsplot,
                              data_function = function() {
                                list(vsd = get_vsd_counts(), condition = filtered_metadata()[[selected_condition_column()]])
                              }, output, input, session, mds_plot_ids)

  #output for PCA plot
  showDownloadButton$pca_plot <- FALSE
  observeEvent(input$plot_pca, {
    output$pca_plot <- renderPlot({
      req(generate_pca_plot())
      pca_plot <- generate_pca_plot()
      showDownloadButton$pca_plot <- TRUE
      return(pca_plot)
    })
  })
  
  createDownloadUI('pca_plot', output, input, session, showDownloadButton)
  
  output$download_pca_plot <- downloadHandler(
    filename = function() {
      paste0(input$filename_pca_plot, '.', input$filetype_pca_plot)
    }, content = function(file) {
      pca_plot <- generate_pca_plot()
      ggsave(filename = file, plot = pca_plot, width = input$plot_width_pca_plot,
      height = input$plot_height_pca_plot, dpi = 300, device = input$filetype_pca_plot)
    }
  )
  
  #output for PCA plot before ruvseq
  showDownloadButton$pca_before_ruv <- FALSE
  output$pca_before_ruv <- renderPlot({
    req(prepare_ruv_obj(), selected_condition_column(), filtered_metadata())
    set_ruv <- prepare_ruv_obj()
    metadata_sub <- filtered_metadata()
    conditionCol <- selected_condition_column()
    
    pca_before_ruv <- make_pca_ruv(data = set_ruv, metadata = metadata_sub, condition = conditionCol)
    showDownloadButton$pca_before_ruv <- TRUE
    return(pca_before_ruv)
  })
  
  pca_before_ruv_ids <- createDownloadUI('pca_before_ruv', output, input, session, showDownloadButton)
  createDownloadHandler_plots('pca_before_ruv', make_pca_ruv,
                              data_function = function() {
                                list(prepare_ruv_obj(), filtered_metadata(), selected_condition_column())
                              }, output, input, session, pca_before_ruv_ids)
  
  #output for PCA plot using DESeq2
  showDownloadButton$deseq2_pcaPlot <- FALSE
  observeEvent(input$deseq2_pca, {
    if (!deseq2_run()) {
      showNotification("Please run DESeq2 from the DE Analysis Tab first to display the PCA plot.", type="error", duration=NULL)
    } else {
      output$deseq2_pcaPlot <- renderPlot({
        req(get_vsd_counts(), filtered_metadata(), selected_condition_column())
        metadata_sub <- filtered_metadata()
        vsd <- get_vsd_counts()
        conditionCol <- selected_condition_column()
        conditions <- metadata_sub[[conditionCol]]
        pca_data <- plotPCA(vsd, intgroup=c(conditionCol), returnData=TRUE)
        pca_plot <- edit_plotPCA(pca_data, metadata_sub, condition = conditionCol, input$pc_x, input$pc_y)
        showDownloadButton$deseq2_pcaPlot <- TRUE
        return(pca_plot)
      })
    }
  })
  
  deseq2_pcaPlot_ids <- createDownloadUI('deseq2_pcaPlot', output, input, session, showDownloadButton)
  createDownloadHandler_plots('deseq2_pcaPlot', edit_plotPCA,
                              data_function = function() {
                                req(get_vsd_counts(), filtered_metadata(), selected_condition_column())
                                metadata_sub <- filtered_metadata()
                                vsd <- get_vsd_counts()
                                conditionCol <- selected_condition_column()
                                conditions <- metadata_sub[[conditionCol]]
                                pca_data <- plotPCA(vsd, intgroup=c(conditionCol), returnData=TRUE)
                                
                                list(pca_data = pca_data, metadata = metadata_sub, condition = conditionCol, input$pc_x, input$pc_y)
                              }, output, input, session, deseq2_pcaPlot_ids)
  
  #output for PCA plot after ruvseq
  showDownloadButton$norm_pca_after_ruv <- FALSE
  output$norm_pca_after_ruv <- renderPlot({
    req(run_ruvseq(), selected_condition_column(), filtered_metadata())
    set1_ruv <- run_ruvseq()
    metadata_sub <- filtered_metadata()
    conditionCol <- selected_condition_column()
    
    norm_pca_after_ruv <- make_pca_ruv(data = set1_ruv, metadata = metadata_sub, condition = conditionCol)
    showDownloadButton$norm_pca_after_ruv <- TRUE
    return(norm_pca_after_ruv)
  })
  
  norm_pca_after_ruv_ids <- createDownloadUI('norm_pca_after_ruv', output, input, session, showDownloadButton)
  createDownloadHandler_plots('norm_pca_after_ruv', make_pca_ruv,
                              data_function = function(){
                                list(run_ruvseq(), filtered_metadata(), selected_condition_column())
                              }, output, input, session, norm_pca_after_ruv_ids)
  
  
  #output for normalized PCA plot
  showDownloadButton$norm_pca <- FALSE
  output$norm_pca <- renderPlot({
    if (!deseq2_run()) {
      showNotification("Please run DESeq2 from the DE Analysis Tab first to display the normalized PCA plot.", type="error", duration=NULL)
    } else {
      req(generate_norm_pca())
      norm_pca <- generate_norm_pca()
      showDownloadButton$norm_pca <- TRUE
      return(norm_pca)
    }
  })
  
  createDownloadUI('norm_pca', output, input, session, showDownloadButton)
  
  output$download_norm_pca <- downloadHandler(
    filename = function() {
      paste0(input$filename_norm_pca, '.', input$filetype_norm_pca)
    }, content = function(file) {
      norm_pca <- generate_norm_pca()
      ggsave(filename = file, plot = norm_pca, width = input$plot_width_norm_pca,
             height = input$plot_height_norm_pca, dpi = 300, device = input$filetype_norm_pca)
    }
  )
  
  #output for description of normalization methods
  observeEvent(input$help_norm, {
    showModal(modalDialog(
      title = "Normalization methods description:",
      tags$ul(
        tags$li(strong("CPM (Counts Per Million): "), "This method normalizes the gene expression measurements by scaling counts to the total number of reads per sample and then adjusting up to per million reads."),
        tags$li(strong("FPKM/RPKM (Fragments/Reads per Kilobase of transcripts per Million mapped reads): "), "This method normalizes gene expression by both the sequencing depth (total number of reads) and the gene length and scales them to per million reads. 
                FPKM is useful for comparing expression levels of different genes within same sample or across samples with similar sequencing depths."),
        tags$li(strong("TPM (Transcripts per Kilobase per Million mapped reads): "), "Similar to FPKM, but TPM first adjusts the counts for gene length, then normalized them by the total number of reads, and then scaled by 1 million reads. 
                TPM is useful when comparing gene expression levels across multiple samples, as it accounts for both sequencing depth and gene length, providing a consistent metric across different conditions."),
        tags$li(strong("DESeq2: "), "Here, the counts are divided by sample-specific size factors determined by median ratio of gene counts relative to geometric mean per gene. 
               This accounts for gene count comparisons between samples and for DE analysis. NOT for within sample comparisons.")
      ),
      easyClose = TRUE,
      size='m'
    ))
  })
  
  #output for the normalized counts table
  output$normCounts <- renderDataTable({
    req(normalized_counts(), input$norm_method)
    
    normCounts <- datatable(normalized_counts(), rownames=FALSE, options = list(paging=TRUE, scrollX=TRUE),
                            caption = paste("Normalized Counts after applying", input$norm_method, "method"))
    return(normCounts)
  })
  
  output$downloadButtonUI <- renderUI({
    req(normalized_counts(), showDownloadButton$norm_counts)
    downloadButton("download_norm_counts", "Download Normalized Counts")
  })
  
  #download the normalized counts
  output$download_norm_counts <- create_download_handler(filename_prefix = reactive({paste0("normalized_counts_", tolower(input$norm_method), "_")}),
                                                         data_function = normalized_counts,
                                                         data_prep = function(data) data)
  
  #output for ERCC counts table
  output$ercc_table <- renderDataTable({
    req(get_ercc_counts())
    ercc_data <- datatable(get_ercc_counts(), options = list(autoWidth=TRUE, paging=TRUE, scrollX=TRUE),
                           caption = "Raw counts for ERCC spike-in Genes")
    return(ercc_data)
  })
  
  #output for ERCC dashboard plot
  showDownloadButton$ercc_plot <- FALSE
  observeEvent(input$generate_ercc_plot, {
    output$ercc_plot <- renderPlotly({
      req(run_ercc_analysis(), ercc_model_stats())
      
      plot_data <- run_ercc_analysis()
      ercc_stats <- ercc_model_stats()
      
      m1 <- make_ercc_plot(plot_data, input$ercc_mix_type, input$ercc_plot_type, ercc_stats, plotly=TRUE)
      showDownloadButton$ercc_plot <- TRUE
      return(m1)
    })
    
    ercc_plot_ids <- createDownloadUI('ercc_plot', output, input, session, showDownloadButton)
    createDownloadHandler_plots("ercc_plot", make_ercc_plot,
                                data_function = function() {
                                  req(run_ercc_analysis(), ercc_model_stats())
                                  
                                  plot_data <- run_ercc_analysis()
                                  ercc_stats <- ercc_model_stats()
                                  list(plot_data, input$ercc_mix_type, input$ercc_plot_type, ercc_stats)
                                }, output, input, session, ercc_plot_ids)
    
    output$ercc_model_stats <- renderText({
      stats <- ercc_model_stats()
      text_out <- ""
      
      if (input$ercc_mix_type == 'mix1'){
        text_out <- paste0(text_out, "Mix 1:\n",
                           "Dynamic Range: ", round(stats$dynamic_range, 2), "\n",
                           "Intercept: ", round(stats$intercept, 2), "\n",
                           "Slope: ", round(stats$slope, 2), "\n",
                           "LLD Concentration: ", round(stats$lld_conc, 2), "\n",
                           "R-squared: ", round(stats$r2, 2), "\n\n")
      } else if (input$ercc_mix_type == 'mix2') {
        text_out <- paste0(text_out, "Mix 2:\n",
                           "Dynamic Range: ", round(stats$dynamic_range, 2), "\n",
                           "Intercept: ", round(stats$intercept, 2), "\n",
                           "Slope: ", round(stats$slope, 2), "\n",
                           "LLD Concentration: ", round(stats$lld_conc, 2), "\n",
                           "R-squared: ", round(stats$r2, 2), "\n\n")
      } 
      
      return(text_out)
    })
  })
  
  # #check if the log2foldchange distribution of ERCC spike ins is centered at 0
  # observeEvent(input$use_genes_dds == "use_all", {
  #   output$ercc_hist_l2fc <- renderPlot({
  #     req(deseq2_run(), run_deseq2(), deg_results())
  #     res <- deg_results()
  #     ercc_res <- res[grep("ERCC", rownames(res)), ]
  #     summary(ercc_res$log2FoldChange)
  #     hist(ercc_res$log2FoldChange, breaks=20, main="ERCC log2FC distribution", xlab="log2 Fold Change")
  #   })
  # })
    
  
  #--------------------------------DE-------------------------------------------
  #update design formula is ruvseq normalization is done
  default_formula <- "~ batch + condition"
  observeEvent(input$ruvseq, {
    req(input$ruv_k)
    
    if (isTRUE(input$ruvseq) && input$ruv_k > 0) {
      w_terms <- paste0("W_", seq_len(input$ruv_k))
      new_formula <- paste("~", paste(c(w_terms, "condition"), collapse = " + "))
      
      updateTextInput(
        session,
        "design_formula",
        value = new_formula
      )
    } else {
      # Reset to default when RUV is turned off
      updateTextInput(
        session,
        "design_formula",
        value = default_formula
      )
    }
  })
  #set reference for comparisons
  debounced_formula <- debounce(reactive(input$design_formula), 3000)
  
  observe({
    req(filtered_metadata(), debounced_formula())
    
    design = as.formula(debounced_formula())
    factor_names <- all.vars(design)
    
    #check if all factor names are present in colData
    #missing_factors <- setdiff(factor_names, colnames(filtered_metadata()))
    
    if (isTRUE(input$ruvseq) && !is.null(input$ruv_k) && input$ruv_k > 0) {
      w_terms <- paste0("W_", seq_len(input$ruv_k))
      non_w_factors <- setdiff(factor_names, w_terms)
      missing_factors <- setdiff(non_w_factors, colnames(filtered_metadata()))
    } else {
      missing_factors <- setdiff(factor_names, colnames(filtered_metadata()))
    }
    
    if (length(missing_factors) > 0){
      msg = paste("The following factors are missing in the metadata: ", paste(missing_factors, collapse = ','), ". Please upload the file containing all the desired columns / select all required columns from the `Samples Tab`.")
      showNotification(msg, type = 'error', duration = NULL)
      output$factor_selector <- renderUI({NULL})
      output$select_ref <- renderUI({NULL})
      output$contrast1_ui <- renderUI({NULL})
      output$contrast2_ui <- renderUI({NULL})
    } else {
      output$factor_selector <- renderUI({
        selectInput("selected_factor", "Select column from design formula to be used for comparisons (ideally the last column - eg: condition)", 
                    choices = factor_names, selected = factor_names[2])
      })
    }
  })
  
  observeEvent(input$selected_factor, {
    req(input$selected_factor)
    selected_factor <- input$selected_factor
    #extract unique factor levels for the selected factor from metadata
    levels <- unique(filtered_metadata()[[selected_factor]])
    
    #update UI for reference selection
    output$select_ref <- renderUI({
      selectInput("reference_level", label = paste("Set reference for", selected_factor), choices = levels, selected = levels[1])
    })
    
    output$contrast1_ui <- renderUI({
      selectInput("contrast1", "Contrast 1 (ideally the treatment group)", choices = levels, selected = levels[2])
    })
    
    output$contrast2_ui <- renderUI({
      selectInput("contrast2", "Contrast 2 (ideally the control/reference group)", choices = levels, selected = levels[1])
    })
    
  })
  
  #reactive Value to keep track of whether DESeq2 has been run (dds created)
  deseq2_run <- reactiveVal(FALSE)
  
  #prepare dds for DESeq2
  prep_dds <- reactive({
    withProgress(message = "Preparing DESeq2 dataset...", value = 0, {
      set.seed(627)
      req(input$design_formula, input$filter_slider, input$selected_factor)
      
      incProgress(0.2, detail = "Loading counts matrix and metadata...")
      
      #use set1_ruv if ruvseq was performed, else filtered counts
      if (isTRUE(input$ruvseq)) {
        req(run_ruvseq())
        set1_ruv <- run_ruvseq()
        counts_matrix <- counts(set1_ruv)
        metadata_subset <- pData(set1_ruv)
      } else {
        req(filtered_counts_sub(), get_counts_matrix(), filtered_metadata())
        counts_matrix <- get_counts_matrix()
        metadata_subset <- filtered_metadata()
      }
      
      design_formula = as.formula(input$design_formula)
      
      dds <- DESeqDataSetFromMatrix(countData = counts_matrix,
                                    colData = metadata_subset,
                                    design = design_formula)
      
      incProgress(0.4, detail = "Checking sample consistency...")
      #check if all samples are present in the counts matrix and are in the same order
      sample_check <- all(rownames(metadata_subset$sample_id %in% colnames(counts_matrix))) && all(rownames(metadata_subset$sample_id == colnames(counts_matrix)))
      
      if (!sample_check) {
        showModal(modalDialog(
          title = "Error", "Samples in the metadata do not match those in the counts matrix or are not in the same order.", easyClose = TRUE, size = 'm'
        ))
        return(NULL)
      }
      
      incProgress(0.6, detail = "Filtering genes based on count threshold...")
      #filter - remove low count genes according to the selected threshold
      if (input$filter_sample == FALSE) {
        print(paste('Using filter threshold as', input$filter_slider))
        filter_threshold <- input$filter_slider
        keep_genes <- rowSums(counts(dds)) >= filter_threshold
      } else {
        print(paste("Using filter threshold as", input$filter_slider, "in sample(s)", input$sample_slider))
        filter_threshold <- input$filter_slider
        sample_threshold <- input$sample_slider
        keep_genes <- rowSums(counts(dds) >= filter_threshold) >= sample_threshold
      }
      
      dds <- dds[keep_genes, ]
      
      #rownames of dds - storing for futher use.
      dds_geneids <- rownames(dds)
      
      incProgress(0.8, detail = "Releveling reference...")
      #relevel reference according to the selected reference
      dds[[input$selected_factor]] <- relevel(dds[[input$selected_factor]], ref = input$reference_level)
      
      
      #estimate size factors 
      if (!isTRUE(input$ruvseq)) {
        if (input$ercc_used == FALSE) {
          print("Not using ERCC spike in genes as control genes..")
          #check which genes to use to estimate size factors
          if (input$use_genes_dds == "use_all") {
            print("Using all genes (including ERCC genes) for DESeq2 analysis.")
            dds <- estimateSizeFactors(dds)
          } else if (input$use_genes_dds == "use_wo_ercc") {
            print("Excluding ERCC spike-in genes from DESeq2 analysis.")
            dds <- dds[!grepl("^ERCC", rownames(dds)), ]
            #check if filtering is done correctly
            print(any(grepl("^ERCC", rownames(dds))))
            print(sum(grepl("^ERCC", rownames(dds))))
            dds <- estimateSizeFactors(dds)
          }
        } else {
          print("Using ERCC spike in genes as control genes to calculate the size factors..")
          ercc_genes <- rownames(counts_matrix)[grepl("ERCC", rownames(counts_matrix))]
          if (length(ercc_genes) == 0) {
            showNotification("No ercc genes found in the counts matrix.", type='error')
            return(NULL)
          }
          #create a logical vector for controlGenes - estimateSizeFactors 
          is_ercc <- rownames(dds) %in% ercc_genes
          print(sum(rownames(dds) %in% ercc_genes))
          dds <- estimateSizeFactors(dds, controlGenes = is_ercc)
        }
      } else {
        dds <- estimateSizeFactors(dds)
      }
 
      incProgress(1, detail = "Estimating dispersions...")
      #estimate dispersions
      print('Estimating dispersions...')
      dds <- estimateDispersions(dds)
      
      #set the DESeq2 flag to TRUE 
      deseq2_run(TRUE)
      
      list(dds = dds, dds_geneids = dds_geneids)
    })
  })
  
  #get variance stabilized counts - varianceStabilizingTransformation
  get_vsd_counts <- reactive({
    withProgress(message = "Generating Variance Stabilized Counts...", value = 0, {
      req(prep_dds())
      dds <- prep_dds()$dds
      vsd <- vst(dds, blind=FALSE)
      print("Variance stabilized counts generated.")
      incProgress(1)
      return(vsd)
    })
  })
  
  #get rlog trasformed values
  get_rld_counts <- reactive({
    withProgress(message = "Generating Rlog Transformed Counts...", value = 0, {
      req(prep_dds())
      dds <- prep_dds()$dds
      rld <- rlog(dds, blind = FALSE)
      print("Rlog Transformed counts generated.")
      incProgress(1)
      return(rld)
    })
  })
  
  #get deseq2 normalized counts
  get_deseq2_norm_counts <- reactive({
    withProgress(message = "Generating DESeq2 Normalized Counts...", value = 0, {
      req(prep_dds())
      dds <- prep_dds()$dds
      norm_counts <- counts(dds, normalized=TRUE)
      print("DESeq2 normalized counts generated.")
      print(head(norm_counts))
      incProgress(1)
      return(norm_counts)
    })
  })
  
  observeEvent(input$run_deseq2, {
    req(prep_dds())
    
    print("Started Analysis...")
    dds <- prep_dds()$dds
    print(dds)
    
    return(dds)
  })
  
  run_deseq2 <- reactive({
    withProgress(message = "Running DEG Analysis with DESeq2...", value = 0, {
      set.seed(627)
      req(prep_dds())
      dds <- prep_dds()$dds
      print("Started DEG analysis..")
      #run DESeq2
      incProgress(0.5, detail = "Running DESeq2...")
      dds <- DESeq(dds)
      print(dds)
      #print result names for deseq
      print('Result Names for dds: ')
      print(DESeq2::resultsNames(dds))
      incProgress(1, detail = "DESeq2 analysis complete.")
      return(dds)
    })
  })
  
  #DEG results for deseq
  deg_results <- reactiveVal()
  deg_res_df <- reactiveVal()
  
  observeEvent(input$generate_deg_res, {
    withProgress(message = "Generating DEG Results...", value = 0, {
      req(deseq2_run(), input$selected_factor, input$contrast1, input$contrast2)
      dds <- run_deseq2()
      
      design <- as.formula(input$design_formula)
      #factor_name <- all.vars(design)[-1][1] #assuming the first factor is used for generating results
      
      if (is.null(input$contrast1) || is.null(input$contrast2)) {
        showNotification("Please select the contrasts for comparisons.", type = 'error')
        return(NULL)
      }
      
      incProgress(0.4, detail = "Calculating DEG results...")
      print("Generating DEG Results..")
      print(paste0("Using alpha cutoff as ", input$alpha_cutoff))
      res <- results(dds, contrast = c(input$selected_factor, input$contrast1, input$contrast2), alpha = as.numeric(input$alpha_cutoff))
      print(head(res))
      print(summary(res))
      
      print(table(res$padj < 0.01))
      
      incProgress(0.7, detail = "Annotating genes with symbols...")
      res <- get_gene_symbols(res)
      
      incProgress(0.9, detail = "Converting results to dataframe...")
      print("Converting DEG results to a dataframe for further analysis..")
      res_df <- as.data.frame(res)
      #update reactive values
      deg_results(res)
      deg_res_df(res_df)
      incProgress(1, detail = "DEG results generation complete.")
    })
  })
  
  #get gene lists 
  #upregulated genes 
  geneLIST_Up <- reactive({
    req(deg_results(), input$alpha_cutoff, input$fc_cutoff)
    sig_up_degs <- get_sig_upreg_genes(res = deg_results(), alpha = as.numeric(input$alpha_cutoff), l2fc = input$fc_cutoff)
    #return list
    list(sig_UP_degs_table = sig_up_degs$sig_UP_degs_table, genelist_UP = sig_up_degs$genelist_UP, sig_up_table = sig_up_degs$sig_UP_table)
  })
  
  #downregulated genes
  geneLIST_Down <- reactive({
    req(deg_results(), input$alpha_cutoff, input$fc_cutoff)
    sig_down_degs <- get_sig_downreg_genes(res = deg_results(), alpha = as.numeric(input$alpha_cutoff), l2fc = input$fc_cutoff)
    #return list
    list(sig_DOWN_degs_table = sig_down_degs$sig_DOWN_degs_table, genelist_DOWN = sig_down_degs$genelist_DOWN, sig_down_table = sig_down_degs$sig_DOWN_table)
  })
  
  #MA plot
  ma_plot <- reactive({
    req(deg_results())
    res <- deg_results()
    ma_plot <- generate_MAPlot(res, paste('MA plot -', input$selected_factor, '-', input$contrast1, 'vs', input$contrast2))
    return(ma_plot)
  })
  
  #show error if trying to run DESeq2 before RUVseq, when RUVSeq is selected
  observeEvent(input$run_deseq2, {
    if (isTRUE(input$ruvseq) && input$ruv_norm == 0) {
      showModal(modalDialog(
        title = "Normalization required",
        "You have selected RUVSeq to normalize the data, but haven't clicked the `Normalize - RUVSeq` button yet. 
        Please do so before running DESeq2.",
        easyClose = TRUE, footer = modalButton("OK")
      ))
    }
  })
  
  #----------------------------Output-------------------------------------------
  #output for description of design formula
  url <- a("DESeq2 Manual", href="https://www.bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html", target="_blank", rel="noopener noreferrer")
  observeEvent(input$help_design, {
    showModal(modalDialog(
      title = "DESeq2 Design formula:",
      tagList(
        p("`Design formula` indicates how to model the samples."),
        p("If you want to use multiple columns, the order in which they are written MATTERS."),
        p("For example: design = ~ batch + condition"),
        p("Here, we measure the effect of the `condition` (treatment) while controlling for `batch` differences"),
        p("Refer: ", url)
      ),
      easyClose = TRUE,
      size='m'
    ))
  })
  
  #output for sample filter
  output$sample_filter <- renderUI({
    req(input$filter_sample, filtered_metadata())
    sliderInput('sample_slider', 'Select the minimum number of samples that should pass the above threshold.', 
                min = 0, max = length(filtered_metadata()$sample_id), value = 0, step = 1)
  })
  
  #output for dispersion estimate plot
  output$disp_plot <- renderPlot({
    req(deseq2_run())
    dds <- prep_dds()$dds
    disp_plot <- plotDispEsts(dds, cex.lab = 1.4, cex.axis = 1.2)
    title("Dispersion Estimates", cex.main = 1.6, font.main = 2)
    return(disp_plot)
  })
  
  #output for meansd rld values
  output$meansd_rld <- renderPlot({
    req(deseq2_run())
    rld <- get_rld_counts()
    meansd_rld_plot <- meanSdPlot(assay(rld), ranks=FALSE)$gg
    meansd_rld_plot <- meansd_rld_plot +
      ggtitle("Mean-SD Plot of Rlog-transformed Counts") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text = element_text(size = 12))
    return(meansd_rld_plot)
  })
  
  #output for meansd vsd values
  output$meansd_vsd <- renderPlot({
    req(deseq2_run())
    vsd <- get_vsd_counts()
    meansd_vsd_plot <- meanSdPlot(assay(vsd), ranks = FALSE)$gg
    meansd_vsd_plot <- meansd_vsd_plot +
      ggtitle("Mean-SD Plot for Variance-Stabilized Tranformed Counts") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
            axis.title.x = element_text(size = 14),
            axis.title.y = element_text(size = 14),
            axis.text = element_text(size = 12))
    return(meansd_vsd_plot)
  })
  
  #output for help button - deg result columns explanation
  observeEvent(input$help_deg_res, {
    showModal(modalDialog(
      title = "DEG Results Table Explained:",
      tags$ul(
        tags$li(tags$b("Gene:"), " The identifier for the gene, such as Ensembl gene ID."),
        tags$li(tags$b("BaseMean:"), " The average of the normalized counts for all samples, providing a meansure of the overall expression level of the gene."),
        tags$li(tags$b("Log2FoldChange:"), " The log2 fold change in expression between the two conditions being compared. 
                A positive values indicates higher expression in the treatment group (Upregulated), while a negative value indicates a higher expression in the control group (Downregulation)."),
        tags$li(tags$b("lfcSE:"), " The standard error of the log2 fold change estimate, providing a measure of variability or uncertainty in the fold change estimate."),
        tags$li(tags$b("Stat:"), " The test statistic for the differential expression test, typically a Wald test statistic."),
        tags$li(tags$b("P value:"), " The p-value for the differential expression test, indicating the statistical significance for the observed difference in expression."),
        tags$li(tags$b("P adj:"), " The Benjamini-Hochberg (BH) adjusted p-value to account for multiple testing."),
        tags$li(tags$b("Symbol:"), " The identifier for the gene.")
      ),
      easyClose = TRUE,
      size = 'm'
    ))
  })
  
  #output for deg summary table
  output$deg_summary_table <- renderDataTable({
    req(deseq2_run(), run_deseq2(), deg_results())
    dds <- run_deseq2()
    res <- deg_results()
    alpha_cutoff <- as.numeric(input$alpha_cutoff)
    summary_df <- generate_deg_summary(input$design_formula, dds, res, input$fc_cutoff, input$contrast1, input$contrast2, alpha_cutoff)
    summary_table <- datatable(summary_df, options = list(paging=TRUE, scrollX=TRUE),
                               caption = "Summary Table for DEG results")
    return(summary_table)
  })
  
  showDownloadButton$degres <- FALSE
  #output for deg results table
  output$deg_table <- renderDataTable({
    req(deg_res_df())
    res_df <- deg_res_df()
    res_df <- res_df[, c(ncol(res_df), 1:(ncol(res_df) - 1))]
    res_df <- datatable(res_df, colnames = c('Geneid' = 1), options = list(autoWidth=TRUE, paging=TRUE, scrollX=TRUE),
                        caption = "DEG results table")
    showDownloadButton$degres <- TRUE
    return(res_df)
  })
  
  #show download button for deg results
  output$downloadButton_degres <- renderUI({
    req(showDownloadButton$degres)
      tagList(downloadButton('download_degres', "Download DEG results"),
              downloadButton('download_sig_degs', "Download Significant DEGs table"))
  })
  
  #download deg results
  output$download_degres <- create_download_handler(filename_prefix = reactive({paste0('deg_results_', tolower(input$contrast1), "_vs_", tolower(input$contrast2), "_")}),
                                                    data_function = deg_results,
                                                    data_prep = function(data) data.frame(GeneId = rownames(data), data))
  
  output$download_sig_degs <- create_download_handler(filename_prefix = reactive({paste0('sig_deg_results_', tolower(input$contrast1), "_vs_", tolower(input$contrast2), "_")}),
                                                      data_function = deg_results,
                                                      data_prep = function(data) {
                                                        sig_data <- get_sig_degs(data, as.numeric(input$alpha_cutoff), input$fc_cutoff)
                                                        sig_data <- data.frame(Geneid = rownames(sig_data), sig_data)
                                                        return(sig_data)
                                                      })
  
  #output for volcano plot
  showDownloadButton$vp <- FALSE
  observeEvent(input$vol_plot, {
    if (is.null(deg_results()) || nrow(deg_results()) == 0) {
      showNotification("Please generate DEG results first from the `Generate DEG Results & Summary` tab.", type='error', duration=NULL)
    } else {
    output$volcano_plot <- renderPlotly({
      req(input$vol_plot > 0, deg_results())
      res <- deg_results()
      top_genes <- ifelse(input$annotate_volcano, input$top_genes_cutoff, 0)
      vp <- generate_volcano_plot(res, input$fc_cutoff, as.numeric(input$alpha_cutoff), 
                                  title = paste('Volcano plot - ', input$selected_factor, "-", input$contrast1, "vs", input$contrast2), 
                                  top_genes, plotly=TRUE)
      showDownloadButton$vp <- TRUE
      return(vp)
    })
  }
  })
  
  output$downloadButton_vp <- renderUI({
    req(showDownloadButton$vp)
    createDownloadUI('vp', output, input, session, showDownloadButton)
  })
  
  output$download_vp <- downloadHandler(
    filename = function() {
      paste0(input$filename_vp, '.', input$filetype_vp)
    }, content = function(file) {
      req(deg_results())
      res <- deg_results()
      top_genes = ifelse(input$annotate_volcano, input$top_genes_cutoff, 0)
      vp <- generate_volcano_plot(res, input$fc_cutoff, input$alpha_cutoff,
                                  title = paste('Volcano plot - ', input$selected_factor, "-", input$contrast1, "vs", input$contrast2), 
                                  top_genes, plotly=FALSE)
      ggsave(filename = file, plot = vp, width = input$plot_width_vp,
             height = input$plot_height_vp, dpi = 300, device = input$filetype_vp)
    }
  )
  
  #output for MA plot
  showDownloadButton$ma_plot <- FALSE
  output$ma_plot <- renderPlot({
    req(deg_results())
    #ma <- generate_MAPlot(res = deg_results(), 
    #                      title = paste('MA plot -', input$selected_factor, '-', input$contrast1, 'vs', input$contrast2))
    showDownloadButton$ma_plot <- TRUE
    replayPlot(ma_plot())
    #return(ma)
  })
  
  createDownloadUI('ma_plot', output, input, session, showDownloadButton)
  output$download_ma_plot <- downloadHandler(
    filename = function() {
      paste0(input$filename_ma_plot, '.', input$filetype_ma_plot)
    }, content = function(file) {
      req(deg_results())
      res <- deg_results()
      if (input$filetype_ma_plot == "png") {
        png(file, width = input$plot_width_ma_plot, height = input$plot_height_ma_plot, res = 300, units="in")
      } else if (input$filetype_ma_plot == "pdf") {
        pdf(file, width = input$plot_width_ma_plot, height = input$plot_width_ma_plot)
      }
      #DESeq2::plotMA(res, ylim=c(-8, 8), alpha=0.01)
      #title(main = "MA plot", col.main='black', font.main=4, cex.main=1.2)
      replayPlot(ma_plot())
      dev.off()
    }
  )
  
  #output for gene lists
  output$up_regTable <- renderDataTable({
    req(geneLIST_Up())
    genesUP <- geneLIST_Up()$sig_up_table
    up_table <- datatable(genesUP, colnames = c('Geneid' = 1), options = list(searching = TRUE, paging=TRUE, scrollX=TRUE, ordering = TRUE),
                          caption = "Significant Up-regulated DEGs")
    return(up_table)
  })
  
  output$down_regTable <- renderDataTable({
    req(geneLIST_Down())
    genesDOWN <- geneLIST_Down()$sig_down_table
    down_table <- datatable(genesDOWN, colnames = c('Geneid' = 1), options = list(searching = TRUE, paging=TRUE, scrollX=TRUE, ordering = TRUE),
                            caption = "Significant Down-regulated DEGs")
    # if (!is.null(genesListDown) && 'symbol' %in% colnames(genesListDown)) {
    #   paste(genesListDown$symbol, collapse = ", ")
    # } else {
    #   "No significant down-regulated genes found."
    # }
  })
  
  #show download button for downloading gene lists
  output$download_geneLists <- renderUI({
    req(geneLIST_Up(), geneLIST_Down())
    if (!is.null(geneLIST_Up()) && !is.null(geneLIST_Down())) {
      tagList(br(), downloadButton('download_upreg', HTML("Download Significant <strong>UP</strong> regulated DEG results table")), br(), br(),
              downloadButton('download_downreg', HTML("Download Significant <strong>DOWN</strong> regulated DEGs results table")), br(), br(),
              downloadButton('download_upreg_list', HTML("Download only the Gene List (Geneids + Symbols) for <strong>UP</strong> regulated DEGs")), br(), br(),
              downloadButton('download_downreg_list', HTML("Download only the Gene List (Geneids + Symbols) for <strong>DOWN</strong> regulated DEGs"))
      )
    }
  })
  
  #download handler for gene lists
  output$download_upreg <- create_download_handler(filename_prefix = "significant_upregulated_results_",
                                                   data_function = function() geneLIST_Up()$sig_UP_degs_table,
                                                   data_prep = function(data) {
                                                     data_df <- as.data.frame(data)
                                                     data_df <- data.frame(Geneid = rownames(data_df), data_df, row.names = NULL)
                                                     return(data_df)
                                                   })
  
  output$download_upreg_list <- create_download_handler(filename_prefix = "significant_upregulated_genelist_",
                                                        data_function = function() geneLIST_Up()$genelist_UP)
  
  output$download_downreg_list <- create_download_handler(filename_prefix = "significant_down_regulated_genelist_",
                                                          data_function = function() geneLIST_Down()$genelist_DOWN)
  
  output$download_downreg <- create_download_handler(filename_prefix = "significant_down_regulated_results_",
                                                     data_function = function() geneLIST_Down()$sig_DOWN_degs_table,
                                                     data_prep = function(data) {
                                                       data_df <- as.data.frame(data)
                                                       data_df <- data.frame(Geneid = rownames(data_df), data_df, row.names = NULL)
                                                       return(data_df)
                                                     })
  
  #-----------------------FUNCTIONAL ANALYSIS-----------------------------------
  load_geneLists <- reactive({
    print("Loading the genelists...")
    if (input$upload_geneLists) {
      req(input$genelist_file)
      #genelist_file <- input$genelist_file
      file_ext <- tools::file_ext(input$genelist_file$name)
      #using gene IDs 
      print("Using and extracting the 1st column (assuming gene IDs) from the file for further analysis..")
      geneList <- switch(file_ext,
                         "csv" = read.csv(input$genelist_file$datapath, stringsAsFactors = FALSE)[,1],
                         "xlsx" = read_excel(input$genelist_file$datapath)[,1],
                         #"txt" = read.table(input$genelist_file$datapath, stringsAsFactors = FALSE)[[1]],
                         read.table(input$genelist_file$datapath, header=FALSE, stringsAsFactors = FALSE)[,1],
                         validate('Unsupported file type. Please upload a `.csv` or `.xlsx` file.'))
    } else if (input$generated_geneLists) {
      req(geneLIST_Up(), geneLIST_Down())
      geneList <- switch(input$geneList_to_use, 
                         "use_upregList" = geneLIST_Up()$genelist_UP$geneid,
                         "use_downregList" = geneLIST_Down()$genelist_DOWN$geneid)
    } else {
      showNotification("No gene list option selected. Please select one..", type='error')
      return(NULL)
    }
    
    print(head(geneList))
    return(geneList)
  })
  
  #initiate reactive value for storing functional analysis results
  func_analysis_results <- reactiveVal()
  func_analysis_plots <- reactiveVal()
  
  observeEvent(input$run_func_analysis, {
    withProgress(message = "Running functional analysis...", detail = "This may take a while to complete.....", value = 0, {
      req(load_geneLists(), input$func_analysis_types)
      gene_ids <- load_geneLists()
      
      if (length(gene_ids) == 0) {
        showNotification("No genes found for the analysis. Please check your input list.")
        return(NULL)
      }
      
      #generate functional analysis data
      incProgress(0.4, detail = "Running enrichment analysis...")
      enrich_results <- tryCatch({
        generate_func_data(
          gene_ids = gene_ids,
          analysis_types = input$func_analysis_types
        )
      }, error = function(e) {
        showNotification(paste("Error during enrichment analysis:", e$message), type = "error")
        return(NULL)
      })
      
      if (is.null(enrich_results) || length(enrich_results) == 0) {
        showNotification("No results from enrichment analysis. Check input genes and analysis parameters.", type = "error")
        return(NULL)
      }
      
      #save enrichment results
      func_analysis_results(enrich_results)
      
      #generate plots
      incProgress(0.8, detail = "Generating plots...")
      plots <- tryCatch({
        generate_func_plot(
          enrich_results = enrich_results,
          plot_type = input$func_plot_type,
          show_categories = input$show_categories,
          font_size = as.numeric(input$font_size)
        )
      }, error = function(e) {
        showNotification(paste("Error during plot generation:", e$message), type = "error")
        return(NULL)
      })
      
      if (is.null(plots) || length(plots) == 0) {
        showNotification("Plot generation failed. Check input parameters and enrichment results. 
                         Or the Enrichment results are empty.", type = "error")
        return(NULL)
      }
      
      # Save plots to reactive
      func_analysis_plots(plots)
      
      incProgress(1, detail = "Functional analysis complete.")
    })
  })
  
  #----------------------------Output-------------------------------------------
  #output for functional plots ui
  #define tab mapping for generating tabs
  tab_mapping <- list(
    'go_bp' = 'GO_BP', 'go_mf' = 'GO_MF', 'go_cc' = 'GO_CC', 'kegg'='KEGG'
  )
  output$func_plots_ui <- renderUI({
    req(input$func_analysis_types)
    selected_func_types <- input$func_analysis_types
    
    tab_panels <- lapply(selected_func_types, function(type){
      tab_name <- tab_mapping[[type]]
      tabPanel(title = tab_name,
               fluidPage(
                 p("Enrichment Plot"),
                 plotOutput(paste0("func_plot_", tab_name), height = "700px"),
                 br(),
                 uiOutput(paste0('downloadButton_', tab_name)),
                 br(),
                 p("Enrichment Table"),
                 dataTableOutput(paste0("table_", tab_name)))
               )
    })
    
    do.call(tabsetPanel, tab_panels)
  })
  
  observe({
    req(func_analysis_results()) 
    
    enrichResults <- func_analysis_results()  
    plot_list <- generate_func_plot(enrichResults, input$func_plot_type, input$show_categories, font_size = as.numeric(input$fontSize))
    
    lapply(names(plot_list), function(type) {
      output[[paste0("func_plot_", type)]] <- renderPlot({
        req(plot_list[[type]])  
        showDownloadButton[[type]] <- TRUE
        plot_list[[type]]
      })
      
      download_info <- createDownloadUI(type, output, input, session, showDownloadButton)
      
      output[[download_info$download_id]] <- downloadHandler(
        filename = function() { 
          paste0(input[[download_info$filename_id]], ".", input[[download_info$filetype_id]]) 
        }, content = function(file) {
          width <- as.numeric(input[[download_info$width_id]])
          height <- as.numeric(input[[download_info$height_id]])
          filetype <- input[[download_info$filetype_id]]
          
          if (filetype == "png") {
            png(file, width = width, height = height, units = "in", res = 300)
          } else if (filetype == "pdf") {
            pdf(file, width = width, height = height)
          }
          
          print(plot_list[[type]])  
          dev.off()
          }
      )
    })
  })
  
  observe({
    req(func_analysis_results())  
    enrichResults <- func_analysis_results()
    
    lapply(input$func_analysis_types, function(type) {
      tab_name <- tab_mapping[[type]]
      output[[paste0("table_", tab_name)]] <- renderDataTable({
        if (!is.null(enrichResults[[tab_name]])) {
          table_data <- as.data.frame(enrichResults[[tab_name]])
          datatable(table_data, options = list(pageLength = 10, autoWidth = TRUE, scrollX=TRUE),
                    caption = paste("Enrichment results for", tab_name))
        }
      })
    })
  })
  
  
  
  #-----------------------------GSEA--------------------------------------------
  #load data for gsea
  load_gsea_data <- reactive({
    if (input$data_source == "upload_geneSet") {
      req(input$gsea_data)
      ext <- tools::file_ext(input$gsea_data$name)
      print("Using the 1st column (assuming gene IDs) from the file as rownames for data..")
      gseaData <- switch(ext,
                         "csv" = read.csv(input$gsea_data$datapath, header=TRUE, row.names = 1),
                         "tsv" = read.delim(input$gsea_data$datapath, header = TRUE, row.names = 1),
                         "xlsx" = read_excel(input$gsea_data$datapath)[[1]],
                         validate("Unsupported file type. Please upload either a `.csv` or `.tsv` file."))
    } else if (input$data_source == 'use_res') {
      req(deg_results(), deg_res_df())
      res_df <- deg_res_df()
      gseaData <- res_df[order(res_df$log2FoldChange, decreasing = TRUE), ]
    }
    
    return(gseaData)
  })
  
  #prepare gene set list
  get_geneSet <- reactive({
    req(load_gsea_data())
    gseaData <- load_gsea_data()
    
    if (!"log2FoldChange" %in% colnames(gseaData)) {
      showNotification("Error: 'log2FoldChange' column is missing from the data. 
                       Please upload a data containing the necessary columns: 'Geneid' and 'log2FoldChange'",
                       type = 'error')
      stop("Error")
    }
    
    #check if gene IDs column exists
    if (any(duplicated(rownames(gseaData))) || all(is.na(rownames(gseaData))) || all(rownames(gseaData) == "")) {
      # If row names are duplicated, NA, or empty, use the 'gene' column if available
      if ("Geneid" %in% colnames(gseaData)) {
        genes <- get_cleaned_geneIds(gseaData$Geneid)
      } else if ("symbol" %in% colnames(gseaData)) {
        genes <- get_cleaned_geneIds(gseaData$symbol)
      } else {
        showNotification("Error: Row names are invalid or missing and no `gene` or `symbol` column found in the data.", 
                         type = "error")
        return(NULL)  
      }
    } else {
      # Use rownames as gene identifiers if they are valid
      genes <- get_cleaned_geneIds(rownames(gseaData))
    }
    
    #extract log2fc values and gene names
    geneSet <- gseaData$log2FoldChange
    #genes <- get_cleaned_geneIds(rownames(gseaData))
    
    # if (is.null(names(geneSet)) || length(names(geneSet)) == 0 || all(is.na(names(geneSet)))) {
    #   showNotification("Error: Gene set log2FoldChange values are not properly named. Please check the gene identifiers.", type = "error")
    #   return(NULL)
    # }
    
    names(geneSet) <- genes
    geneSet <- sort(geneSet, decreasing = TRUE)
    
    print("Gene set for GSEA:")
    print(head(geneSet))
    return(geneSet)
  }) 
  
  #run GSEA
  fgsea_results <- reactiveVal(NULL)
  observeEvent(input$run_gsea, {
    withProgress(message = "Running GSEA...", value = 0, {
      req(get_geneSet(), input$species, input$gsea_category)
      geneSet <- get_geneSet()
      
      incProgress(0.2, detail = "Loading gene sets from msigdbr...")
      
      #extract the msigdbr gene sets for analysis
      print("Loading the gene sets from msigdbr...")
      m_df <- msigdbr::msigdbr(species = input$species, category = input$gsea_category)
      if (input$geneId_type == "Ensembl") {
        pathways <- split(x = m_df$ensembl_gene, f = m_df$gs_name)
      } else {
        pathways <- split(x = m_df$gene_symbol, f = m_df$gs_name)
      }
      
      incProgress(0.4, detail = "Preparing pathways for GSEA...")
      print("Pathways:")
      print(head(pathways))
      
      incProgress(0.6, detail = "Running GSEA analysis...")
      #run fgsea
      print("Running fgsea...")
      fgsea_res <- fgsea(pathways = pathways,
                         stats = geneSet, 
                         minSize = 15, maxSize = 500)
      print("Fgsea results..:")
      print(head(fgsea_res))
      incProgress(1, detail = "GSEA analysis complete.")
      fgsea_results(fgsea_res)
    })
  })
  
  sig_pathways <- reactive({
    withProgress(message = "Extracting significant pathways...", value = 0, {
      req(fgsea_results())
      fgsea_res <- fgsea_results()
      
      print("Extracting significant pathways from fgsea results...")
      fgsea_res_sorted <- fgsea_res[order(fgsea_res$NES, decreasing = TRUE), ]
      
      incProgress(0.6, detail = "Filtering significant pathways...")
      significant_pathways <- fgsea_res_sorted[fgsea_res_sorted$padj < 0.05, ]
      print(head(significant_pathways))
      
      incProgress(1, detail = "Pathway extraction complete.")
      return(significant_pathways)
    })
  })
  
  #----------------------------Output-------------------------------------------
  #output for help button
  url_gse <- a("GSEA", href="https://www.gsea-msigdb.org/gsea/index.jsp", target="_blank", rel="noopener noreferrer")
  url_msd <- a("Molecular Signatures Database (MSigDB)", href="https://www.gsea-msigdb.org/gsea/msigdb/index.jsp", 
               target="_blank", rel="noopener noreferrer")
  observeEvent(input$help_gsea, {
    showModal(modalDialog(
      title = "GSEA Category Descriptions",
      tags$div(
        tags$ul(
          tags$li(tags$b("Hallmark Gene sets (H):"), " 50 gene sets that summarize and represent specific well-defined biological states/processes and display coherent expression."),
          tags$li(tags$b("Positional Gene sets (C1):"), " Gene sets corresponding to each human chromosome and cytogenic band. These reflect the gene architeture as represented on the primary assembly."),
          tags$li(tags$b("Curated Gene sets (C2):"), " Gene sets collected from various sources such as online pathway databases, publications in PubMed, and knownledge of domain experts"),
          tags$li(tags$b("Regulatory Target Gene sets (C3):"), " Gene sets representing potential targets of regulation by transcription factors or microRNAs."),
          tags$li(tags$b("Computational Gene sets (C4):"), " Gene sets defined by mining large collections of cancer-oriented microarray data."),
          tags$li(tags$b("GO Gene sets (C5):"), " Gene sets derived from the Gene Ontology (GO) database, which describes gene functions."),
          tags$li(tags$b("Oncogenic signatures (C6):"), " Gene sets representing signatures of cellular pathways that are dysregulated in cancer."),
          tags$li(tags$b("Immunologic signatures (C7):"), " Gene sets representing cell states and perturbations within the immune system."),
          tags$li(tags$b("Cell type signatures (C8):"), " Gene sets representing cell type-specific gene expression signatures.")
        ),
        tags$p("For more details, refer the following resources:"),
        tags$p("1. ", url_gse, 
               "(Subramanian A, et al. (2005). 
             Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles"),
        tags$p("2. ", url_msd),
        easyClose = TRUE,
        size='m'
      )))
  })
  
  #output for plot
  showDownloadButton$gseaPlot <- FALSE
  output$gsea_plot <- renderPlot({
    req(input$run_gsea, fgsea_results(), sig_pathways())
    significant_pathways <- sig_pathways()
    plot <- generate_gsea_plot(significant_pathways, input$top_cat)
    showDownloadButton$gseaPlot <- TRUE
    return(plot)
  }, height = function() {
    significant_pathways <- sig_pathways()
    if (!is.null(significant_pathways) && nrow(significant_pathways) > 0) {
      #calculate height dynamically but limit it to a maximum value (e.g., 3000px)
      dynamic_height <- 500 + (nrow(significant_pathways) * 10)
      return(min(dynamic_height, 3000))
    } else {
      return(400)  #default height if no pathways are found
    }
  }, width = 1200)
  
  gsea_ids <- createDownloadUI('gseaPlot', output, input, session, showDownloadButton)
  createDownloadHandler_plots('gseaPlot', 
                              plot_function = generate_gsea_plot,
                              data_function = function() {
                                list(sig_pathways(), input$top_cat)
                              }, output, input, session, gsea_ids)
  
  #output for help button - explains gsea results columns
  observeEvent(input$help_gsea_res, {
    showModal(modalDialog(
      title = "GSEA Results Description",
      tags$div(
        tags$ul(
          tags$li(tags$b("Pathway:"), " Name of the Pathway being analyzed."),
          tags$li(tags$b("Pval:"), " The p-value (statistical significance) of the enrichment score of the pathway."),
          tags$li(tags$b("Padj:"), " The adjusted p-value to account for multiple testing. This helps control for the false dicovery rate (FDR)"),
          tags$li(tags$b("log2err:"), " The log2-transformed error estimate for the enrichment score. It provides measure of the variability or uncertainty in the enrichment score."),
          tags$li(tags$b("ES (Enrichment Score):"), " This reflects the degree to which the pathway is overrepresented at the top or bottom of the rank of genes. 
                  A positive ES indicates enrichment at the top of the list, while negative ES indicates enrichment at the bottom."),
          tags$li(tags$b("NES (Normalized Enrichment Score):"), " This accounts for differences in gene set size and other factors. Higher NES = stronger enrichment"),
          tags$li(tags$b("Size:"), " The number of genes in the gene set that are found in the dataset and used in the enrichment analysis."),
          tags$li(tags$b("LeadingEdge:"), " The subset of genes in the gene set that contribute most to the enrichment score.")
        )
      ),
      easyClose = TRUE,
      size = 'm'
    ))
  })
  
  showDownloadButton$gseaRes <- FALSE
  #output for gsea results table
  output$gsea_table <- renderDataTable({
    req(sig_pathways())
    fgsea_res_df <- as.data.frame(sig_pathways())
    fgsea_res_df <- datatable(fgsea_res_df, options = list(pageLength = 10, autoWidth = TRUE, scrollX=TRUE),
                              caption = "GSEA results")
    showDownloadButton$gseaRes <- TRUE
    return(fgsea_res_df)
  })
  
  output$downloadButton_gseaRes <- renderUI({
    req(showDownloadButton$gseaRes)
      downloadButton("download_gseaRes", "Download GSEA Results")
  })
  
  output$download_gseaRes <- create_download_handler(filename_prefix = reactive({paste0("GSEA_Results_", input$gsea_category, "_")}),
                                                     data_function = sig_pathways, data_prep = function(data){
                                                       df <- as.data.frame(data)
                                                       df[] <- lapply(df, function(col) {
                                                         if(is.list(col)) {
                                                           return(sapply(col, toString))
                                                         } else {
                                                           return (col)
                                                         }
                                                       })
                                                       return(df)
                                                     })
  
}