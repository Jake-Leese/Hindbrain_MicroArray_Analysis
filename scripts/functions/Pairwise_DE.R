# Functions to perform a pairwise differential expression analysis between 2 given cell types in my expression matrix
# Pairwise_DE_table returns a table overview of differential gene expression
# Pairwise_DE_VP returns a volcano plot visualising results

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


Pairwise_DE_VP <- function(expr_matrix, sample_info, cell_type_A, cell_type_B){
  
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
  
  # Volcano plot
  volcanoplot(fit, coef = paste0("cell_type", cell_type_A), highlight = 10, names = colnames(expr_matrix))
}
