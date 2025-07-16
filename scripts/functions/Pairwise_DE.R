# Functions to perform a pairwise differential expression analysis between 2 given cell types in my expression matrix
# Pairwise_DE_table returns a table overview of differential gene expression
# Pairwise_DE_VP returns a volcano plot visualising results


#########################################################################################################################
############################################ Pairwise_DE_table() ########################################################

Pairwise_DE_table <- function(expr_matrix, sample_info, cell_type_A, cell_type_B){

# Subset to samples of interest (r1 and r2)
selected_samples <- sample_info$cell_type %in% c(cell_type_A, cell_type_B)
expr_subset <- expr_matrix[selected_samples, ]
sample_info_subset <- sample_info[selected_samples, ]

# Create design matrix for differential expression
# cell_type should be a factor
sample_info_subset$cell_type <- factor(sample_info_subset$cell_type, levels = c(cell_type_B, cell_type_A))
design <- model.matrix(~ cell_type, data = sample_info_subset)

# Fit the linear model
fit <- lmFit(t(expr_subset), design)  # transpose so genes are rows
fit <- eBayes(fit)

# Extract results for the contrast: r2 vs r1
results <- topTable(fit, coef = paste0("cell_type", cell_type_A), adjust.method = "BH", number = Inf)

return(results)
}


#########################################################################################################################
############################################## Pairwise_DE_VP() #########################################################

Pairwise_DE_VP <- function(expr_matrix, sample_info, cell_type_A, cell_type_B, logFC_threshold, pval_threshold, n_top_genes){
  
  # Subset to samples of interest (r1 and r2)
  selected_samples <- sample_info$cell_type %in% c(cell_type_A, cell_type_B)
  expr_subset <- expr_matrix[selected_samples, ]
  sample_info_subset <- sample_info[selected_samples, ]
  
  # Create design matrix for differential expression
  # cell_type should be a factor
  sample_info_subset$cell_type <- factor(sample_info_subset$cell_type, levels = c(cell_type_B, cell_type_A))
  design <- model.matrix(~ cell_type, data = sample_info_subset)
  
  # Fit the linear model
  expr_subset_t <- t(expr_subset)
  rownames(expr_subset_t) <- colnames(expr_subset)  # ensure gene names
  
  # Fit the linear model
  fit <- lmFit(expr_subset_t, design)
  fit <- eBayes(fit)
  
  # Volcano plot
  # Extract topTable from limma fit
  tt <- topTable(fit, coef = paste0("cell_type", cell_type_A), number = Inf, sort.by = "P")
  tt$gene <- rownames(tt)
  tt$log10P <- -log10(tt$P.Value)
  
  # Define log10P_threshold 
  log10P_threshold <- -log10(pval_threshold)
  
  # Define cell type labels
  label_A <- paste("Up in ", cell_type_A)
  label_B <- paste("Up in ", cell_type_B)
  
  # Classify significance
  tt <- tt %>%
    mutate(sig = case_when(
      logFC > logFC_threshold & log10P > log10P_threshold ~ label_A,
      logFC < -logFC_threshold & log10P > log10P_threshold ~ label_B,
      TRUE ~ "NS"
    ))
  
  # Get top 15 up and down genes for labeling
  top_up <- tt %>% filter(sig == label_A) %>% arrange(P.Value) %>% head(15)
  top_down <- tt %>% filter(sig == label_B) %>% arrange(P.Value) %>% head(15)
  top_genes <- bind_rows(top_up, top_down)
  
  # Define color mapping dynamically
  color_map <- setNames(c("red", "blue", "grey"), c(label_A, label_B, "NS"))
  
  # Plot
  ggplot(tt, aes(x = logFC, y = log10P)) +
    geom_point(aes(color = sig), alpha = 0.6, size = 1.5) +
    scale_color_manual(values = color_map) +
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
    geom_hline(yintercept = log10P_threshold, linetype = "dashed", color = "black") +
    geom_text(aes(x = max(logFC) * 0.8, y = log10P_threshold, 
                  label = paste0("p = ", format(pval_threshold, scientific = TRUE))),
              vjust = -0.5, hjust = 0, size = 5) +
    geom_text_repel(data = top_genes, aes(label = ID), size = 3, max.overlaps = 100) +
    labs(title = paste("Volcano plot:", cell_type_A, " vs ", cell_type_B),
         x = "log2 Fold Change",
         y = "-log10(P-value)",
         color = "Expression") +
    theme_minimal(base_size = 14)
}
