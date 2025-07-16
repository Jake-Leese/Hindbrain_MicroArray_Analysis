setwd("/data/Hindbrain_microarray")
.libPaths("/R/libs/R_scvi_integration_v3.4/")

Sys.setenv(OMP_NUM_THREADS = "1")

library(GEOquery)
library(oligo)
library(limma)
library(Biobase)
library(pd.mogene.1.0.st.v1)
library(mogene10sttranscriptcluster.db)
library(annotate)
library(dplyr)
library(ggrepel)
library(ggplot2)
library(pheatmap)
library(scatterplot3d)
library(affyPLM)  # Had to install with configure.args = "CFLAGS='-fopenmp -DOPENMP=1'" to overcome multithreading error

source("./scripts/functions/Pairwise_DE.R")

raw_data_path <- "./raw_files/"
QC_path <- "./QC/"

# Download GEO series and extract raw CEL files
getGEOSuppFiles("GSE48359", baseDir = raw_data_path)
#untar(paste0(raw_data_path, "GSE48359/GSE48359_RAW.tar"), exdir = paste0(raw_data_path, "GSE48359/CEL"))

# Read in CEL files
cel_files <- list.files(paste0(raw_data_path, "GSE48359/CEL"), pattern = "\\.CEL\\.gz$", full.names = TRUE)
sapply(cel_files, gunzip, overwrite = TRUE)

raw_data <- ReadAffy(celfile.path = paste0(raw_data_path, "GSE48359/CEL"))

# Create sample_info dataframe to accompany microarray
sample_names <- sampleNames(raw_data)

sample_info <- data.frame(
  sample = sample_names,
  gsm_id = gsub("(_.*)|\\.CEL$", "", sample_names),
  cell_type = rep(c("m", "r1", "r2", "r3", "r4", "r5"), times = 2)
)

saveRDS(sample_info, "./processed_data/sample_info.rds")

############################################################################################################
########################################## QUALITY CONTROL #################################################

# Raw data plots
# Boxplot of raw log2 intensities. Not overly drastic differences between samples
boxplot(log2(exprs(raw_data)),
        main = "Raw Data (Log2 Intensity)", 
        las = 2, col = "lightblue")

plm <- fitPLM(raw_data)

# RLE Plot
RLE(plm, main = "RLE - Raw Data")

# NUSE Plot
NUSE(plm, main = "NUSE - Raw Data")


# NORMALIZE
norm_data <- rma(raw_data)

# Post-norm Boxplot
boxplot(exprs(norm_data),
        main = "RMA Normalized Data",
        las = 2, col = "lightgreen")

# Convert to expression matrix
expr_matrix <- exprs(norm_data)

# Global normalization to 50th percentile (quantile normalization)
expr_matrix <- normalizeBetweenArrays(expr_matrix, method = "quantile")

# Scaling: center each gene to median across samples
expr_matrix <- sweep(expr_matrix, 1, apply(expr_matrix, 1, median), "-")

# REMOVE NON-EXPRESSED GENES

# Use MAS5 calls to get detection flags
calls <- mas5calls(raw_data)
present_calls <- exprs(calls)  # "P", "M", "A"

# Convert to logical matrix: TRUE = Present
present_matrix <- present_calls == "P"

# Keep only genes that are present in at least one replicate pair
# (You can customize this further by grouping sample_info$cell_type)

# For simplicity: keep genes present in at least 2 samples overall
present_count <- rowSums(present_matrix)
keep_genes <- present_count >= 2
expr_matrix <- expr_matrix[keep_genes, ]

# Remove probes starting with "AFFX"
keep_non_affx <- !grepl("^AFFX", rownames(expr_matrix))
expr_matrix <- expr_matrix[keep_non_affx, ]


############################################################################################################
############################################### PCA ########################################################

# Caclulate PCA on normalised and filtered expr_matrix
pca <- prcomp(t(expr_matrix), scale. = TRUE)

# Calculate proportion of variance explained
variance_explained <- pca$sdev^2 / sum(pca$sdev^2)

# Convert to percent
percent_variance <- round(variance_explained * 100, 2)

# Scree plot
barplot(percent_variance[1:10],
        names.arg = paste0("PC", 1:10),
        las = 2,
        ylab = "Percent Variance Explained",
        main = "Scree Plot")

# Create a color palette
cell_types <- unique(sample_info$cell_type)
cell_colors <- rainbow(length(cell_types))  # Or use RColorBrewer for nicer palettes

# Convert cell types to a factor for coloring
col_vector <- cell_colors[as.factor(sample_info$cell_type)]

# Plot PCA

plot(pca$x[,1:2], col = col_vector,
     main = "PCA of Samples", pch = 19)

legend("right",
       inset = c(-0.15, 0),
       legend = cell_types,
       col = cell_colors,
       pch = 19,
       title = "Cell Type")

# 3D plot of PCA with 3 principle components
scatterplot3d(
  x = pca$x[, 1],
  y = pca$x[, 2],
  z = pca$x[, 3],
  color = col_vector,
  pch = 19,
  xlab = paste0("PC1 (", percent_variance[1], "%)"),
  ylab = paste0("PC2 (", percent_variance[2], "%)"),
  zlab = paste0("PC3 (", percent_variance[3], "%)"),
  main = "3D PCA Plot"
)
legend("topright", legend = unique(sample_info$cell_type),
       col = unique(col_vector), pch = 19, title = "Cell Type")

# Set margin to allow space on the right
par(mar = c(5, 4, 4, 8), xpd = TRUE)


############################################################################################################
################################# Measure distance between replicates ######################################

expr_t <- t(expr_matrix)  # samples x genes

# Remove genes with zero variance across samples
expr_t_clean <- expr_t[, apply(expr_t, 2, function(x) sd(x, na.rm = TRUE) > 0)]

# Compute Pearson correlation matrix
cor_matrix <- cor(t(expr_t_clean), method = "pearson", use = "pairwise.complete.obs")

# Distance = 1 - correlation
dist_matrix <- as.dist(1 - cor_matrix)

# Confirm validity, Check if there are no NA values in dist_matrix
all(!is.na(dist_matrix))

# hierarchical clusterng using average linkage
hc <- hclust(dist_matrix, method = "average")

plot(hc,
     main = "Condition Tree (Average Linkage, Pearson Distance)",
     xlab = "",
     sub = "",
     sample_info$cell_type)


############################################################################################################
################################### Add gene symbols to expr_matrix ########################################

# Load array annotation file
annot <- read.delim("GPL3213-9900.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# Extract current probe IDs
probe_ids <- colnames(expr_t_clean)

# Subset annotation to relevant probes
annot_filtered <- annot[match(probe_ids, annot$ID), ]

# Start with the original probe IDs
new_colnames <- colnames(expr_t_clean)

# Match probe IDs to annotation table
matched <- match(colnames(expr_t_clean), annot_filtered$ID)

# Replace with gene symbol *only if it's non-empty and non-NA
gene_symbols <- annot_filtered$Gene.Symbol[matched]

# Identify where gene symbols are available
has_symbol <- !is.na(gene_symbols) & gene_symbols != ""

# Replace those rownames
new_colnames[has_symbol] <- gene_symbols[has_symbol]

# Assign the new rownames
colnames(expr_t_clean) <- new_colnames

saveRDS(expr_t_clean, file = "./processed_data/expr_matrix.rds")


#############################################################################################################
################################### Differential expression analysis ########################################

expr_t_clean <- readRDS("./processed_data/expr_matrix.rds")
sample_info <- readRDS("./processed_data/sample_info.rds")

# Retransform data back to default setting, genes x samples
expr_clean <- t(expr_t_clean)

# Build the design matrix for a linear model
design <- model.matrix(~ 0 + factor(sample_info$cell_type))
colnames(design) <- levels(factor(sample_info$cell_type))

# Fitting a linear model
fit <- lmFit(expr_clean, design)
contrast.matrix <- makeContrasts(
  m - (r1 + r2 + r3 + r4 + r5)/5,
  r1 - (m + r2 + r3 + r4 + r5)/5,
  r2 - (m + r1 + r3 + r4 + r5)/5,
  r3 - (m + r1 + r2 + r4 + r5)/5,
  r4 - (m + r1 + r2 + r3 + r5)/5,
  r5 - (m + r1 + r2 + r3 + r4)/5,
  levels = design
)

# Apply contrast matrix to original model fit
fit2 <- contrasts.fit(fit, contrast.matrix)

# Perform empirical Bayes moderation of standard errors
fit2 <- eBayes(fit2)

# Repeat for each contrast, or loop:
marker_lists <- lapply(colnames(contrast.matrix), function(coef) {
  topTable(fit2, coef = coef, number = 100, adjust = "fdr", sort.by = "P")
})
names(marker_lists) <- c("m", "r1", "r2", "r3", "r4", "r5")

# Write each element of the list to a separate CSV
for (cell_type in names(marker_lists)) {
  write.csv(marker_lists[[cell_type]],
            file = paste0("./processed_data/markers_", cell_type, ".csv"),
            row.names = TRUE)
}

# Make a heatmap
top_genes <- unique(unlist(lapply(marker_lists, function(x) head(x$ID, 20))))
pheatmap(expr_clean[top_genes, ], 
         annotation_col = data.frame(row.names = sample_info$sample, cell_type = sample_info$cell_type), 
         show_rownames = TRUE)


#############################################################################################################
############################################# Pairwise DEGs #################################################

m_vs_r1 <- Pairwise_DE_table(expr_t_clean, sample_info, "m", "r1")
r1_vs_r2 <- Pairwise_DE_table(expr_t_clean, sample_info, "r1", "r2")
r2_vs_r3 <- Pairwise_DE_table(expr_t_clean, sample_info, "r2", "r3")
r3_vs_r4 <- Pairwise_DE_table(expr_t_clean, sample_info, "r3", "r4")
r4_vs_r5 <- Pairwise_DE_table(expr_t_clean, sample_info, "r4", "r5")

write.csv(m_vs_r1, file = "./processed_data/markers_m_vs_r1.csv", row.names = TRUE)
write.csv(r1_vs_r2, file = "./processed_data/markers_r1_vs_r2.csv", row.names = TRUE)
write.csv(r2_vs_r3, file = "./processed_data/markers_r2_vs_r3.csv", row.names = TRUE)
write.csv(r3_vs_r4, file = "./processed_data/markers_r3_vs_r4.csv", row.names = TRUE)
write.csv(r4_vs_r5, file = "./processed_data/markers_r4_vs_r5.csv", row.names = TRUE)

# Generating volcano plots for pairwise comparisons
Pairwise_DE_VP(expr_t_clean, sample_info, "m", "r1", logFC_threshold = 2, pval_threshold = 0.005, n_top_genes = 15)
Pairwise_DE_VP(expr_t_clean, sample_info, "r1", "r2", logFC_threshold = 2, pval_threshold = 0.005, n_top_genes = 15)
Pairwise_DE_VP(expr_t_clean, sample_info, "r2", "r3", logFC_threshold = 2, pval_threshold = 0.005, n_top_genes = 15)
Pairwise_DE_VP(expr_t_clean, sample_info, "r3", "r4", logFC_threshold = 2, pval_threshold = 0.005, n_top_genes = 15)
Pairwise_DE_VP(expr_t_clean, sample_info, "r4", "r5", logFC_threshold = 2, pval_threshold = 0.005, n_top_genes = 15)

