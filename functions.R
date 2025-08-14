#These functions are used in the App server

#----Libraries----
#function to install and load necessary packages - if not already installed
install_and_load <- function(packages) {
  #ensure BiocManager is installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  #install missing packages
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      BiocManager::install(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

#---- Samples ---
#function to plot the histogram
# plot_histogram <- function(metadata_subset, columns_to_plot, column_to_groupby){
#   hist <- ggplot(metadata_subset, aes(x = !!sym(columns_to_plot), fill = !!sym(column_to_groupby))) +
#     geom_histogram(bins=30, alpha=0.5, position = "identity") +
#     labs(title = paste0("Histogram of ", input$columns_to_plot, " by ", input$column_to_groupby),
#          x = input$columns_to_plot, y = "Count") +
#     theme_bw() +
#     theme(legend.position = "bottom")
#   hist <- ggplotly(hist)
#   
#   return(hist)
# }

#----Counts----- 
#function to read the counts data
read_counts <- function(file, ext){
  counts_data <- switch(ext,
                        "csv" = read.csv(file$datapath, header=TRUE), #row.names=1),
                        "txt" = read_table(file$datapath, comment = '#'),
                        validate("Unsupported file type. Please upload a `.csv` or `.txt` file."))
  return(counts_data)
}

#function to get EMSEMBL IDs (cleaned from Gencode ids) (i.e. remove version numbers)
get_cleaned_geneIds <- function(geneIds){
  if (any(grepl('\\..+', geneIds))) {
    cleaned_ids <- gsub('\\..+', '', geneIds)
  } else {
    cleaned_ids <- geneIds
  }
  return(cleaned_ids)
}

#function to get gene symbols using ENSEMBL gene ids
get_gene_symbols <- function(counts) {
  #get gene-ids
  if ('Geneid' %in% colnames(counts)) {
    gene_ids <- counts$Geneid
  } else {
    gene_ids <- rownames(counts)
  }
  
  #split the ENSEMBL and ERCC ids for annotation
  ens_ids <- gene_ids[grepl("ENSG", gene_ids)]
  ercc_ids <- gene_ids[grepl("ERCC", gene_ids)]
  
  #convert GENOCODE to ENSEMBL IDs 
  keys <- get_cleaned_geneIds(ens_ids)
  
  #add gene symbols
  ens_symbols <- AnnotationDbi::mapIds(EnsDb.Hsapiens.v86, keys = keys,
                                       column = 'SYMBOL', keytype = 'GENEID', multiVals = 'first')
  
  gene_ids <- c(setNames(keys, ens_ids), setNames(ercc_ids, ercc_ids))
  symbols <- c(ens_symbols, setNames(ercc_ids, ercc_ids))
  
  counts$symbol <- symbols[gene_ids]
  
  #return the results
  return(counts)
}

#function to get ENTREZ ids from ENSEMBL gene ids
get_entrez_ids <- function(gene_ensembl_ids, key_type) {
  cleaned_ids <- sub("\\..*", "", gene_ensembl_ids)
  entrez_ids <- bitr(cleaned_ids, fromType = key_type, toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  return(entrez_ids$ENTREZID)
}

#----Visualizations----
#function to make library sizes plot
make_libSize_plot <- function(library_sizes, metadata, plotly=FALSE) {
  valid_metadata <- metadata[metadata$sample_id %in% names(library_sizes), ]
  
  possible_cell_cols <- c("cell_type", "cell_line", "cell_lines")
  fill_col <- intersect(possible_cell_cols, colnames(valid_metadata))[1]
  
  df <- data.frame(sample = valid_metadata$sample_id, 
                   millions_reads = library_sizes[valid_metadata$sample_id]/1e6) 
  #cell_type = factor(valid_metadata$cell_type))
  
  if (!is.na(fill_col)) {
    df$cell <- factor(valid_metadata[[fill_col]])
  }
  
  print(df)
  
  plot_aes <- if (!is.na(fill_col)) {
    aes(x = millions_reads, y = sample, fill = cell)
  } else {
    aes(x = millions_reads, y = sample)
  }
  
  library_size_plot <- ggplot(df, plot_aes) +
    geom_bar(stat = "identity") +
    scale_x_continuous() +
    labs(x = 'Samples', y = 'Number of Reads', title = 'Library Size (Million reads per sample)') +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 15))
  
  if (!is.na(fill_col)) {
    library_size_plot <- library_size_plot +
      scale_fill_brewer(palette = "Paired") +
      labs(fill=fill_col)
  }
  
  if (plotly) {
    return(ggplotly(library_size_plot))
  } else {
    return(library_size_plot)
  }
}

#function to make sample distributions boxplot
plot_sample_distributions <- function(counts_matrix, plotly=FALSE) {
  filtered_counts <- counts_matrix[rowSums(counts_matrix) >= 10, ]
  
  #convert to tibble - if not
  if(!is_tibble(filtered_counts)) {
    filtered_counts <- as_tibble(filtered_counts)
  }
  
  data_long <- filtered_counts %>% pivot_longer(cols = starts_with("R"), names_to = 'sample', values_to = 'counts')
  sample_dist_plot <- ggplot(data_long, aes(x=sample, y=counts, color = sample)) +
    geom_boxplot() + scale_y_log10() +
    labs(x = "Samples", y = "Counts") +
    ggtitle("Sample Distributions (log10 (Counts))") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          text = element_text(size=15))
  
  if (plotly) {
    return(ggplotly(sample_dist_plot))
  } else {
    return(sample_dist_plot)
  }
}

#function to make RUVSeq - RLE boxplots
make_rlePlot <- function(data, metadata, condition) {
  #color mapping 
  sample_conditions <- metadata$condition
  unique_conditions <- unique(sample_conditions)
  colors <- brewer.pal(max(3, length(unique_conditions)), "Set2")
  color_map <- setNames(colors[1:length(unique_conditions)], unique_conditions)
  sample_colors <- color_map[sample_conditions]
  
  #calulate ylim range
  log_expr <- log2(exprs(data) + 1)
  medians <- apply(log_expr, 1, median, na.rm = TRUE)
  rle_values <- sweep(log_expr, 1, medians, FUN = "-")
  
  #trim extremes: use 1st and 99th percentiles for ylim
  lower <- quantile(rle_values, 0.01, na.rm = TRUE)
  upper <- quantile(rle_values, 0.99, na.rm = TRUE)
  ylim_range <- c(max(-3, lower), min(3, upper))
  
  rle_plot <- EDASeq::plotRLE(data, outline=FALSE, ylim = ylim_range, col= sample_colors,
                              main = "Relative Log Expression (RLE) Plot",
                              ylab = "RLE")
  legend("topright", legend = unique_conditions, fill = color_map[unique_conditions],
         border = NA, bty = "n", cex = 1.5)
  
  return(rle_plot)
}

#function to make a sample-to-sample heatmap
make_distHeatmap <- function(vsd, condition) {
  #calculate sample distances
  sample_dist <- dist(t(assay(vsd)))
  #convert to distance matrix
  sample_distmat <- as.matrix(sample_dist)
  #set the rownames and colnames
  rownames(sample_distmat) <- condition
  colnames(sample_distmat) <- colnames(vsd)
  #set the colors for the heatmap
  colors <- colorRampPalette(rev(brewer.pal(20, "Blues")))(255)
  
  #generate the heatmap
  heatmap_grab <- grid.grabExpr(
    pheatmap(sample_distmat,
             clustering_distance_rows=sample_dist,
             clustering_distance_cols=sample_dist,
             col=colors,
             fontsize=15,
             main="Sample-to-sample Distance Heatmap")
  )
  
  heatmap_ggplot <- ggplot() +
    annotation_custom(heatmap_grab) +
    theme_void()
  
  return(heatmap_ggplot)
}

#function to make MDS plot
make_mdsplot <- function(vsd, condition) {
  #calculate sample distances
  sample_dist <- dist(t(assay(vsd)))
  #convert to distance matrix
  sample_distmat <- as.matrix(sample_dist)
  #generate the mds plot
  msd <- as.data.frame(colData(vsd)) %>%
    cbind(cmdscale(sample_distmat))
  msd$condition <- factor(condition)
  mds <- ggplot(msd, aes(x=`1`, y=`2`)) +
    geom_point(aes(color=condition), size=3) +
    geom_text(aes(label=colnames(vsd)), vjust = 1.5, color = "black", size = 4) +
    theme_classic() +
    ggtitle('MDS plot') + xlab('MDS 1') + ylab('MDS 2') +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 15)) 
  #return the plot
  print(mds)
}

#function to generate PCA plot
make_pca <- function(pca_results, metadata, condition_col, color_mapping, title, pc_x, pc_y) {
  #extract pca data from `prcomp()` results
  pca_data <- as.data.frame(pca_results$x)
  
  #get which and how many PCs are available 
  available_pcs <- colnames(pca_data)
  num_pcs <- length(available_pcs)
  
  pca_data$sample_id <- rownames(pca_data)
  pca_data <- merge(pca_data, metadata, by='sample_id')
  
  if (!(pc_x %in% available_pcs) || !(pc_y %in% available_pcs)) {
    # Show an error notification with the number of available PCs
    showNotification(paste("The selected principal component (PC) does not exist.",
                           "The available PCs for this data:", paste(available_pcs, collapse = ", "),
                           "Please select valid PCs between PC1 and PC", num_pcs, "."),
                     type = 'error')
    return(NULL)
  }
  
  
  var <- pca_results$sdev^2
  var_explained <- var / sum(var) * 100
  
  pca_data[[condition_col]] <- factor(pca_data[[condition_col]], levels=names(color_mapping))
  
  aes_mapping <- aes(x = !!sym(pc_x), y = !!sym(pc_y), color = !!sym(condition_col))
  if ("batch" %in% colnames(pca_data)) {
    pca_data$batch <- factor(pca_data$batch)
    aes_mapping <- modifyList(aes_mapping, aes(shape = batch))
  } else {
    showNotification("Batch column is missing. Proceeding without batch information.", type = 'warning', duration = 30)
  }
  
  pca_plot <- ggplot(pca_data, aes_mapping) +
    geom_point(size = 3) +
    geom_text(aes(label = sample_id), vjust = 1.5, color = 'black', size = 4) +
    scale_color_manual(name = 'Condition', values = color_mapping) +
    theme_classic() +
    ggtitle(title) +
    xlab(paste(pc_x, "- ", round(var_explained[as.numeric(gsub("PC", "", pc_x))], 2), "% variance explained")) +
    ylab(paste(pc_y, "- ", round(var_explained[as.numeric(gsub("PC", "", pc_y))], 2), "% variance explained")) +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 15))
  
  return(pca_plot)
}

#function to edit the PCA plot obtained from DESeq2
edit_plotPCA <- function(pca_data, metadata, condition, pc_x, pc_y) {
  #get percent variance explained 
  percentVar <- round(100*attr(pca_data, "percentVar"))
  #merge the metadata and pca data by sample id
  metadata <- metadata %>% select(-all_of(condition)) 
  pca_data <- merge(pca_data, metadata, by.x="name", by.y="sample_id")
  
  print("Column names in the merged PCA data:")
  print(colnames(pca_data))
  
  #condition_x <- paste0(condition, ".x")
  
  # if (condition_x %in% colnames(pca_data)) {
  #   pca_data$condition_col <- factor(pca_data[[condition]])
  # } else {
  #   stop(paste("The expected condition column", condition, "was not found in the merged PCA data."))
  # }
  pca_data[[condition]] <- factor(pca_data[[condition]])
  
  aes_mapping <- aes(x = !!sym(pc_x), y = !!sym(pc_y), color = condition)
  if ("batch" %in% colnames(pca_data)) {
    pca_data$batch <- factor(pca_data$batch)
    aes_mapping <- modifyList(aes_mapping, aes(shape = batch))
  } else {
    showNotification("Batch column is missing. Proceeding without batch information.", type = 'warning', duration = 30)
  }
  #generate the PCA plot
  pca_plot <- ggplot(pca_data, aes_mapping) +
    geom_point(size=3) +
    geom_text(aes(label=name), vjust = 1.5, color = "black", size = 4) +
    ggtitle("PCA plot using `plotPCA()` from DESeq2") +
    xlab(paste(pc_x, " -", percentVar[1], "% variance explained")) + 
    ylab(paste(pc_y, " -", percentVar[2], "% variance explained")) +
    theme_classic() + 
    theme(plot.title = element_text(hjust=0.5),
          text = element_text(size = 15))
  #return the plot
  print(pca_plot)
}

make_pca_ruv <- function(data, metadata, condition) {
  #color mapping
  sample_conditions <- metadata$condition
  unique_conditions <- unique(sample_conditions)
  colors <- brewer.pal(max(3, length(unique_conditions)), "Set2")
  color_map <- setNames(colors[1:length(unique_conditions)], unique_conditions)
  sample_colors <- color_map[sample_conditions]
  
  pca_plot <- EDASeq::plotPCA(data, col=sample_colors, cex=1.1)
  legend("topright", legend = unique_conditions, fill = color_map[unique_conditions],
         border = NA, bty = "n", cex = 1.5)
  
  return(pca_plot)
}

#function to generate ERCC dose response curve plot
make_ercc_plot <- function(data, mix_type, ercc_plot_type, ercc_stats, plotly=FALSE){
  conc_column <- ifelse(mix_type == 'mix1', "transcript_molecules_mix1", "transcript_molecules_mix2")
  title <- paste("ERCC", mix_type, "Dose Response", ifelse(ercc_plot_type == "sample_wise_ercc", "(Sample-wise)",""))
  
  m1 <- ggplot(data, aes(x = log2(!!sym(conc_column)), y = log_fpkm)) +
    geom_smooth(method = 'lm', se=TRUE, color='black', level=0.99) +
    geom_hline(yintercept = log2(1), linetype = "dotted") +
    annotate("text", x = 80, y = log2(1), label = "1 RPKM", hjust = -0.5, vjust = 1.5) +
    labs(title = title,
         x = "Log2 Transcript Molecules",
         y = "Log2 FPKM") +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5),
          text = element_text(size = 15))
  
  #plot type specifics 
  if (ercc_plot_type == 'sample_wise_ercc') {
    m1 <- m1 + facet_wrap(~ sample_id)
  } else {
    m1 <- m1 + annotate("text", x = 70, y = 20, label = ercc_stats$annotate_text, hjust = 1, vjust = 1, size = 4, color = "black")
  }
  
  if (!plotly) {
    m1 <- m1 + 
      geom_point(aes(color=sample_id, shape=subgroup),size = 3, alpha = 0.6)
  } else {
    m1 <- m1 +
      geom_point(aes(color=sample_id, shape=subgroup,
                     text = paste0(
                       "ERCC ID: ", ercc_id
                     )),
                 size = 3, alpha = 0.6)
  }
}

#function to generate volcano plot
generate_volcano_plot <- function(res, l2fc, alpha, title, xlab, ylab, top_genes, plotly=FALSE) {
  #ensure l2fc and alpha are numeric
  alpha <- as.numeric(alpha)
  l2fc <- as.numeric(l2fc)
  
  #convert the results obj to a dataframe
  labelled_res <- as.data.frame(res) %>%
    rownames_to_column(var = 'gene') %>%
    dplyr::filter(!is.na(log2FoldChange) & !is.na(padj)) %>%
    mutate(
      #set the padj == 0 values to a very small number to avoid Inf -log10(padj)
      padj = ifelse(padj == 0, .Machine$double.xmin, padj),
      #labels
      sig = case_when(
        padj < alpha & log2FoldChange > l2fc  ~ "Upregulated",
        padj < alpha & log2FoldChange < -l2fc ~ "Downregulated",
        TRUE                                  ~ "Not Significant"
      )
    )
  
  #get the top up-regulated DEGs for plotting
  top_upreg <- labelled_res %>%
    filter(sig == "Upregulated") %>%
    arrange(padj, desc(log2FoldChange)) %>%
    slice_head(n = top_genes)
  
  #get the top down-regulated DEGs for plotting
  top_downreg <- labelled_res %>%
    filter(sig == "Downregulated") %>%
    arrange(padj, log2FoldChange) %>%
    slice_head(n = top_genes)
  
  label_all <- paste("padj <", alpha, "& abs(l2fc) >", l2fc)
  label_padj <- paste("padj <", alpha)
  label_l2fc <- paste("abs(l2fc) >", l2fc)
  
  labelled_res <- labelled_res %>%
    mutate('DEG_label' = factor(case_when(
      padj < alpha & abs(log2FoldChange) > l2fc ~ label_all,
      padj < alpha ~ label_padj,
      abs(log2FoldChange) > l2fc ~ label_l2fc,
      TRUE ~ "NS"
    ), levels = c(label_all, label_padj, label_l2fc, "NS")))
  
  print(str(labelled_res))
  
  #generate volcano plot
  vp <- ggplot(labelled_res, aes(x = log2FoldChange, y = -log10(padj),
                                 text = paste("Gene:", gene, "<br>Symbol:", symbol, "<br>log2FC:", round(log2FoldChange, 2),
                                              "<br>-log10(Padj):", round(-log10(padj), 2)))) +
    geom_point(aes(color=DEG_label)) +
    labs(x = xlab, y = ylab, title = title) +
    theme_classic() +
    geom_hline(yintercept = -log10(alpha), linetype = 'dashed', color='grey46', linewidth = 0.6) +
    geom_vline(xintercept = c(-l2fc, l2fc), linetype = "dashed", color = "grey46", linewidth = 0.6) +
    theme(plot.title = element_text(hjust=0.5),
          text = element_text(size = 15))
  
  if (plotly==FALSE) {
    if (top_genes > 0) {
      vp <- vp +
        geom_text_repel(data = top_upreg, aes(x = log2FoldChange, y = -log10(padj), label = symbol), color="black", size = 3, max.overlaps = Inf) +
        geom_text_repel(data = top_downreg, aes(x = log2FoldChange, y = -log10(padj), label = symbol), color='black', size = 3, max.overlaps = Inf)
    }
    
    return(vp)
    
  } else {
    vp <- ggplotly(vp, tooltip = "text")
    
    vp <- vp %>%
      style(hoverlabel = list(bgcolor = "lightblue", font = list(color = "black", size = 12)))
    
    if (top_genes > 0) {
      annotations <- list()
      for (i in 1:nrow(top_upreg)) {
        annotations <- c(annotations, list(
          list(
            x = top_upreg$log2FoldChange[i], y = -log10(top_upreg$padj[i]),
            text = top_upreg$symbol[i], xref = "x", yref = "y", showarrow = TRUE, arrowhead = 2, ax = ifelse(i %% 2 == 0, 20, -20),
            ay = ifelse(i %% 2 == 0, -40, 40), yanchor="bottom", font = list(size = 10, color = "grey70")
          )
        ))
      }
      for (i in 1:nrow(top_downreg)) {
        annotations <- c(annotations, list(
          list(
            x = top_downreg$log2FoldChange[i], y = -log10(top_downreg$padj[i]),
            text = top_downreg$symbol[i], xref = "x", yref = "y", showarrow = TRUE, arrowhead = 2, ax = ifelse(i %% 2 == 0, 20, -20),
            ay = ifelse(i %% 2 == 0, -40, 40), yanchor="bottom", font = list(size = 12, color = "grey70")
          )
        ))
      }
      vp <- vp %>% layout(annotations = annotations)
      # if (top_genes > 0) {
      #   annotations <- list()
      #   
      #   for (i in 1:nrow(top_upreg)) {
      #     annotations <- c(annotations, list(
      #       list(
      #         x = top_upreg$log2FoldChange[i],
      #         y = -log10(top_upreg$padj[i]),
      #         text = top_upreg$symbol[i],
      #         xref = "x", yref = "y",
      #         showarrow = TRUE,
      #         arrowhead = 2,
      #         standoff = 5,
      #         font = list(size = 10, color = "grey20"),
      #         arrowcolor = "grey60",
      #         bgcolor = "white",
      #         bordercolor = "grey80"
      #       )
      #     ))
      #   }
      #   
      #   for (i in 1:nrow(top_downreg)) {
      #     annotations <- c(annotations, list(
      #       list(
      #         x = top_downreg$log2FoldChange[i],
      #         y = -log10(top_downreg$padj[i]),
      #         text = top_downreg$symbol[i],
      #         xref = "x", yref = "y",
      #         showarrow = TRUE,
      #         arrowhead = 2,
      #         standoff = 5,
      #         font = list(size = 10, color = "grey20"),
      #         arrowcolor = "grey60",
      #         bgcolor = "white",
      #         bordercolor = "grey80"
      #       )
      #     ))
      #   }
      #   
      #   vp <- vp %>% layout(annotations = annotations)
      # 
    }
    
    #return the plot
    return(vp)
  }
}

#function to generate the MA plot
generate_MAPlot <- function(res, title) {
  ma <- DESeq2::plotMA(res, ylim = c(-8,8), alpha=0.01)
  title(main = title, col.main='black', font.main=4, cex.main=1.2)
  return(recordPlot())
}

#----Normalized Counts----
#functions to compute normalized counts
#CPM
calculate_cpm <- function(counts){
  total_counts <- colSums(counts)
  cpm <- t(t(counts)/total_counts) * 1e6
  return(cpm)
}

#FPKM/RPKM
calculate_fpkm <- function(counts, gene_lengths){
  #convert gene lengths from bases to kilobases
  gene_lengths_kilobases <- gene_lengths/1000
  
  #calculate reads per kilobase
  rpk <- counts/gene_lengths_kilobases
  
  #sum total reads to normalize each library size to one million reads
  total_counts <- colSums(counts)
  fpkm <- t(t(rpk)/total_counts ) * 1e9
  return(fpkm)
}

#TPM 
calculate_tpm <- function(counts, gene_lengths){
  #convert gene lengths from bases to kilobases
  gene_lengths_kilobases <- gene_lengths/1000
  
  #calculate reads per kilobase
  rpk <- counts/gene_lengths_kilobases
  
  #calculate the sum of RPKs per sample (column-wise)
  sum_rpk <- colSums(rpk)
  
  #calculate TPM - scale the sum to one million
  tpm <- t(t(rpk)/sum_rpk) * 1e6
  return(tpm)
}

#-----ERCC Analysis ----
#function to prepare the ercc info, i.e. adjust the concentration according to the dilution factor used.
prepare_ercc_info <- function(ercc_file_path, dilution_factor) {
  ercc_info <- read_delim(ercc_file_path, delim='\t')
  colnames(ercc_info) <- tolower(gsub(" ", "_", colnames(ercc_info)))
  #get new concentrations according to the dilution factor
  ercc_info <- ercc_info %>%
    mutate(adj_conc_mix1 = `concentration_in_mix_1_(attomoles/ul)`*dilution_factor,
           adj_conc_mix2 = `concentration_in_mix_2_(attomoles/ul)`*dilution_factor,
           transcript_molecules_mix1 = `adj_conc_mix1`*6.02214179*10^23,
           transcript_molecules_mix2 = `adj_conc_mix2`*6.02214179*10^23)
  #transcript_molecules_mix1 = `adj_conc_mix1`*6.02214179*10^5,
  #transcript_molecules_mix2 = `adj_conc_mix2`*6.02214179*10^5)
  
  head(ercc_info)
  return(ercc_info)
}

#function to get data to generate the ERCC analysis - dose response plot
prepare_combined_ercc_data <- function(fpkm_spike, ercc_data, ercc_info, metadata, conditionCol) {
  fpkm_spike_df <- as.data.frame(fpkm_spike)
  fpkm_spike_df$ercc_id <- ercc_data$Geneid
  head(fpkm_spike_df)
  
  #convert to long format
  fpkm_long <- pivot_longer(fpkm_spike_df, cols=-ercc_id, names_to = "sample_id", values_to = "fpkm", names_pattern = "(RNA\\d+)")
  
  #merge fpkm data with ercc info data 
  combined_data <- left_join(fpkm_long, ercc_info, by = "ercc_id")
  
  #prepare experiment conditions data frame
  if (conditionCol %in% colnames(metadata)){
    exp_cond <- metadata %>%
      dplyr::select(sample_id, !!sym(conditionCol))
  } else {
    stop(paste("Column", conditionCol, "does not exist in the metadata."))
  }
  #exp_cond <- metadata %>%
  # dplyr::select(sample_id, conditionCol) 
  
  #adjust FPKM values for log transformation and merge with conditions
  combined_data <- combined_data %>%
    mutate(fpkm_adj = fpkm + 1,  #avoid log(0) by adding 1
           log_fpkm = log2(fpkm_adj))
  
  combined_data <- merge(combined_data, exp_cond, by='sample_id')
  
  return(combined_data)
}

#function to compute ERCC dose-response curve statistics
calculate_stats <- function(data, molecule_col) {
  dynamic_range <- max(log2(data[[molecule_col]])) - min(log2(data[[molecule_col]]))
  model <- lm(log_fpkm ~ log2(data[[molecule_col]]), data=data)
  intercept <- coef(model)[1]
  slope <- coef(model)[2]
  lld_conc <- -intercept/slope
  r2 <- summary(model)$r.squared
  
  annotate_text <- paste0(#"n: ", round(dynamic_range, 2), "\n",
    #    "Slope: ", round(slope, 2), "\n",
    "R²: ", round(r2, 2))
  
  list(dynamic_range = dynamic_range, intercept = intercept, slope = slope, lld_conc = lld_conc, r2 = r2, model_summary = model, annotate_text = annotate_text)
}

#---- DEG----
#function to generate data for summary table
generate_deg_summary <- function(design, dds, res, l2fc, contrast1, contrast2, alpha) {
  total_dds_genes <- nrow(dds)
  total_upreg_genes <- nrow(res[res$padj < alpha & res$log2FoldChange > 0 & !is.na(res$padj), ])
  total_downreg_genes <- nrow(res[res$padj < alpha & res$log2FoldChange < 0 & !is.na(res$padj), ])
  total_degs <- total_upreg_genes + total_downreg_genes
  total_sig_genes <- nrow(res[res$padj < alpha & abs(res$log2FoldChange) > l2fc & !is.na(res$padj), ])
  sig_upreg <- nrow(res[res$padj < alpha & res$log2FoldChange > l2fc & !is.na(res$padj), ])
  sig_downreg <- nrow(res[res$padj < alpha & res$log2FoldChange < -l2fc & !is.na(res$padj), ])
  
  summary_df <- data.frame(
    Metric = c("Design formula used", "Contrast", "Total of genes", "Total no. of DEGs (without log2FC threshold)", "Total no. of Up-regulated Genes", "Total no. of Down-regulated Genes", 
               paste("Total no. of Significant DEGs (threshold: adjusted p value: ", alpha, "& abs(log2FC) > ", l2fc),
               "No. of Significant Up-regulated DEGs", "No. of Significant Down-regulated DEGs"),
    Value = c(design, paste0(contrast2, ' vs ', contrast1), total_dds_genes, total_degs, total_upreg_genes, total_downreg_genes, total_sig_genes, sig_upreg, sig_downreg))
  
  return(summary_df)
}

#function to get significant DEGs table 
get_sig_degs <- function(res, alpha, l2fc){
  print("Generating Significant DEGs Table...")
  print(paste("Using the threshold padj <", alpha, "& |log2 fold change| > ", l2fc))
  sig_degs_table <- res[res$padj < alpha & abs(res$log2FoldChange) > l2fc & !is.na(res$padj), ]
  print(summary(sig_degs_table))
  return(sig_degs_table)
}

#function to extract significant up-regulated genes
get_sig_upreg_genes <- function(res, alpha, l2fc){
  sig_UP_degs_table <- res[res$padj < alpha & res$log2FoldChange > l2fc & !is.na(res$padj), ]
  print(summary(sig_UP_degs_table))
  sig_UP_degs_df <- as.data.frame(sig_UP_degs_table)
  sig_UP_table <- sig_UP_degs_df[, c('symbol', 'padj', 'log2FoldChange')]
  #genelist along with gene symbols
  genelist_UP <- data.frame(geneid = get_cleaned_geneIds(rownames(sig_UP_degs_df)), 
                            symbol = sig_UP_degs_df$symbol, row.names = NULL)
  print(head(genelist_UP))
  #geneLIST_Up(TRUE)
  
  #return both the table and the gene list
  list(sig_UP_degs_table = sig_UP_degs_table, genelist_UP = genelist_UP, sig_UP_table = sig_UP_table)
}

#function to extract significant down-regulated genes
get_sig_downreg_genes <- function(res, alpha, l2fc){
  sig_DOWN_degs_table <- res[res$padj < alpha & res$log2FoldChange < -l2fc & !is.na(res$padj), ]
  print(summary(sig_DOWN_degs_table))
  sig_DOWN_degs_df <- as.data.frame(sig_DOWN_degs_table)
  sig_DOWN_table <- sig_DOWN_degs_df[, c('symbol', 'padj', 'log2FoldChange')]
  #genelist along with gene symbols
  genelist_DOWN <- data.frame(geneid = get_cleaned_geneIds(rownames(sig_DOWN_degs_df)), 
                              symbol = sig_DOWN_degs_df$symbol, row.names = NULL)
  print(head(genelist_DOWN))
  #geneLIST_Down(TRUE)
  
  #return both the table and the gene list
  list(sig_DOWN_degs_table = sig_DOWN_degs_table, genelist_DOWN = genelist_DOWN, sig_DOWN_table = sig_DOWN_table)
}

#--- Functional Analysis ----
#function to generate plots according to the selected visualization type
#function to get data for selected functional analysis type
generate_func_data <- function(gene_ids, analysis_types, show_notification = TRUE) {
  enrich_results <- list()
  
  if ('go_bp' %in% analysis_types) {
    if (show_notification) message("Running GO BP enrichment analysis...")
    enrich_results[['GO_BP']] <- enrichGO(gene = gene_ids, OrgDb = org.Hs.eg.db,
                                          keyType = "ENSEMBL", ont = "BP")
  }
  if ('go_mf' %in% analysis_types) {
    if (show_notification) message("Running GO MF enrichment analysis...")
    enrich_results[['GO_MF']] <- enrichGO(gene = gene_ids, OrgDb = org.Hs.eg.db,
                                          keyType = "ENSEMBL", ont = "MF")
  }
  if ('go_cc' %in% analysis_types) {
    if (show_notification) message("Running GO CC enrichment analysis...")
    enrich_results[['GO_CC']] <- enrichGO(gene = gene_ids, OrgDb = org.Hs.eg.db,
                                          keyType = "ENSEMBL", ont = "CC")
  }
  if ('kegg' %in% analysis_types) {
    if (show_notification) message("Running KEGG pathway enrichment analysis...")
    gene_ids_entrez <- bitr(gene_ids, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    if (is.null(gene_ids_entrez) || nrow(gene_ids_entrez) == 0) {
      stop("No gene IDs could be mapped to ENTREZ IDs.")
    }
    enrich_results[['KEGG']] <- enrichKEGG(gene = gene_ids_entrez$ENTREZID,
                                           organism = 'hsa', keyType = "kegg")
  }
  
  #print(str(enrich_results))
  return(enrich_results)
}

#function to generate plots according to the selected type - using the func analysis data
generate_func_plot <- function(enrich_results, plot_type, show_categories, font_size) {
  plot_list <- list()
  
  generate_plot <- function(result, plot_type, show_categories) {
    if (!is.null(result) && nrow(as.data.frame(result)) > 0) {
      #if (!is.null(result) || length(result@result) == 0) {
      #font_size <- ifelse(is.null(font_size) || !is.numeric(font_size) || font_size <= 0, 12, font_size)
      
      if (plot_type == 'bar') {
        return(barplot(result, showCategory = show_categories)) #+
        #theme(text = element_text(size=font_size)))
      } else if (plot_type == 'dot') {
        return(dotplot(result, showCategory = show_categories))# +
        #theme(text = element_text(size=font_size)))
      } else if (plot_type == 'emap') {
        return(emapplot(pairwise_termsim(result), showCategory = show_categories)) # +
        # theme(text = element_text(size=font_size)))
      }
    } else {
      showNotification("Result is NULL or empty, skipping plot generation..", type='warning')
      return(NULL)
    }
  }
  
  for (type in names(enrich_results)) {
    plot_obj <- generate_plot(enrich_results[[type]], plot_type, show_categories)
    if (!is.null(plot_obj)) {
      title_text <- paste("Top", show_categories, "Categories for", type)
      plot_list[[type]] <- func_add_title(plot_obj, title_text, font_size)
    }
  }
  
  return(plot_list)
}


#--- GSEA ----
generate_gsea_plot <- function(significant_pathways, top_categories) {
  top_categories <- as.numeric(top_categories)
  significant_pathways$pathway <- stringr::str_wrap(significant_pathways$pathway, width = 30)
  
  upreg <- significant_pathways %>% 
    filter(NES > 0) %>%
    arrange(desc(NES)) %>%
    head(top_categories)
  
  downreg <- significant_pathways %>%
    filter(NES < 0) %>%
    arrange(NES) %>%
    head(top_categories)
  
  top_pathways <- rbind(upreg, downreg)
  
  top_pathways$direction <- ifelse(top_pathways$NES > 0, "Upregulated", "Downregulated")
  
  p <- ggplot(top_pathways, aes(reorder(pathway, NES), NES)) +
    geom_bar(stat = "identity", aes(fill = direction)) +
    coord_flip() +
    theme_minimal() +
    labs(title = "Top Significantly Enriched Pathways",
         x = "Pathways", y = "Normalized Enrichment Score (NES)") +
    scale_fill_manual(values = c("Downregulated" = "red", "Upregulated" = "blue")) + 
    theme(legend.title = element_blank(),
          axis.text.x = element_text(size = 16),  
          axis.text.y = element_text(size = 12),  
          axis.title.x = element_text(size = 20), 
          axis.title.y = element_text(size = 20), 
          plot.title = element_text(size = 25, face = "bold"),  
          legend.text = element_text(size = 18))
  
  return(p)
}


#---- Additional----
#function to add title to the enrichment plots
func_add_title <- function (plot, title, font){
  if (is.null(plot)) {
    showNotification("Plot object is NULL, skipping title addition.")
    return(NULL)
  }
  
  plot <- plot + labs(title = title) +
    theme(plot.title = element_text(hjust = 0.5, size = font),
          axis.title.x = element_text(size = font),
          axis.title.y = element_text(size = font),
          axis.text.x = element_text(size = font),
          axis.text.y = element_text(size = font),
          legend.title = element_text(size = font),
          legend.text = element_text(size = font),
          strip.text = element_text(size = font),
          plot.caption = element_text(size = font),
          plot.subtitle = element_text(size = font),
          text = element_text(size = font))
  
  return(plot)
}

#function to create download handler for data frames
create_download_handler <- function(filename_prefix, data_function, data_prep = identity){
  downloadHandler(
    filename = function() {
      if (is.function(filename_prefix)) {
        prefix <- filename_prefix()
      } else {
        prefix <- filename_prefix
      }
      paste0(prefix, Sys.Date(), ".csv")
    },
    content = function(file) {
      data <- data_function()
      prepared_data <- data_prep(data)
      write.csv(prepared_data, file, row.names = FALSE)
    }
  )
}