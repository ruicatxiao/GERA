# ==============================================================================
# ROBUST METHYLATION ANALYSIS SCRIPT
# ==============================================================================

# Load libraries
library(tidyverse)
library(ggpubr)
library(pheatmap)
library(UpSetR)

# 1. Load and Prepare Data
# ==============================================================================
input_file <- "combined_methylation_matrix.tsv"
df <- read.table(input_file, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Reshape to long format
long_df <- df %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(cols = -Gene, names_to = "Full_Name", values_to = "Methylation") %>%
  # Fix: Use separate (works reliably) or manually parse if needed, but separate is usually fine
  separate(Full_Name, into = c("Sample", "Region", "Type"), sep = "_", remove = TRUE)

# Detect Samples
samples <- unique(long_df$Sample)
s1 <- sort(samples)[1]
s2 <- sort(samples)[2]
cat("Analysis:", s2, "vs", s1, "\n")

# ==============================================================================
# 2. Global DNA Methylation Profile (Violin Plot)
# ==============================================================================
ggplot(long_df, aes(x = Sample, y = Methylation, fill = Sample)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  facet_grid(Type ~ Region, scales = "free_y") +
  theme_bw() +
  labs(title = "Global Methylation Profile", y = "Methylation (%)")


# ==============================================================================
# 3. Gene Body vs Promoter Correlation (Grid)
# ==============================================================================
cor_data <- long_df %>%
  pivot_wider(names_from = Region, values_from = Methylation)

ggplot(cor_data, aes(x = gb, y = promoter)) +
  geom_point(alpha = 0.4, color = "darkblue") +
  geom_smooth(method = "lm", se = TRUE, color = "red") +
  stat_cor(method = "pearson", size = 3) +
  facet_grid(Sample ~ Type) +
  labs(title = "Promoter vs Gene Body Correlation", 
       x = "Gene Body (%)", y = "Promoter (%)") +
  theme_bw()


# ==============================================================================
# 4. Differential Analysis (Identify Significant Genes Only)
# ==============================================================================
diff_thresh <- 30
min_meth <- 1

diff_df <- long_df %>%
  pivot_wider(names_from = Sample, values_from = Methylation) %>%
  mutate(
    Diff = !!sym(s2) - !!sym(s1),
    Mean = ( !!sym(s1) + !!sym(s2) ) / 2,
    Status = case_when(
      Diff > diff_thresh & Mean > min_meth ~ paste0("Hyper_", s2),
      Diff < -diff_thresh & Mean > min_meth ~ paste0("Hypo_", s2),
      TRUE ~ "NS"
    )
  )

# Extract ONLY significant genes
sig_diff_df <- diff_df %>% filter(Status != "NS")
sig_genes_list <- unique(sig_diff_df$Gene)
sig_diff_df
write.table(sig_diff_df, 
            file = "Significant_Differential_Genes.tsv", 
            sep = "\t", 
            quote = FALSE, 
            row.names = FALSE)
cat("Significant genes written to: Significant_Differential_Genes.tsv\n")


cat("Total Significant Genes Found:", length(sig_genes_list), "\n")

if(length(sig_genes_list) == 0) {
  stop("No significant genes found with current thresholds. Adjust diff_thresh or min_meth.")
}

# Create a subset of the LONG dataframe containing ONLY significant genes
# This is used for Heatmaps and Violin plots
long_df_sig <- long_df %>% filter(Gene %in% sig_genes_list)

# ==============================================================================
# 5. Heatmaps (Only for Significant Genes)
# ==============================================================================
for (meth_type in unique(long_df_sig$Type)) {
  
  cat(paste("Processing Heatmap:", meth_type, "\n"))
  
  # Create Matrix
  mat_data <- long_df_sig %>%
    filter(Type == meth_type) %>%
    mutate(Condition = paste(Sample, Region, sep = "_")) %>%
    select(Gene, Condition, Methylation) %>%
    pivot_wider(names_from = Condition, values_from = Methylation) %>%
    column_to_rownames("Gene") %>%
    as.matrix()
  
  # Handle NAs
  mat_data[is.na(mat_data)] <- 0
  
  # Filter Zero Variance rows to prevent pheatmap crash
  row_vars <- apply(mat_data, 1, var, na.rm = TRUE)
  mat_data <- mat_data[!is.na(row_vars) & row_vars > 0, ]
  
  if(nrow(mat_data) > 0) {
    pheatmap(mat_data, 
             scale = "none", 
             clustering_distance_rows = "euclidean",
             clustering_distance_cols = "euclidean",
             main = paste("Significant Genes:", meth_type),
             show_rownames = TRUE,  # --- REQUESTED: Show Gene Names ---
             fontsize_row = 16,       # Small font to fit names
             breaks = seq(0, 100, length.out = 100), 
             color = colorRampPalette(c("blue", "white", "red"))(100))
  } else {
    cat("  Skipped: No variable genes in this group.\n")
  }
}

# ==============================================================================
# 6. UpSet Plot (Differential Genes)
# ==============================================================================
# Create list from the significant diff dataframe
upset_list <- sig_diff_df %>%
  mutate(Category = paste(Type, Region, Status, sep = "_")) %>%
  select(Gene, Category) %>%
  split(.$Category) %>%
  map(~ .x$Gene)

# Basic upset call to avoid argument errors
upset(fromList(upset_list), order.by = "freq")

# ==============================================================================
# 7. Violin Plot for Significant Genes Only
# ==============================================================================
ggplot(long_df_sig, aes(x = Sample, y = Methylation, fill = Sample)) +
  geom_violin(trim = TRUE, alpha = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  facet_grid(Type ~ Region, scales = "free_y") +
  theme_bw() +
  labs(title = "Methylation Levels of Significant Genes Only", 
       subtitle = paste("n =", length(sig_genes_list)),
       y = "Methylation (%)")

# ==============================================================================
# COMPREHENSIVE HEATMAP: Significant Genes (0-100% Scale)
# ==============================================================================

# 1. Prepare the Data Matrix
# We use the 'long_df' (full data) filtered by the 'sig_diff_df' (significant genes)
# We include ALL regions, samples, and types.

heatmap_matrix <- long_df %>%
  # Filter for ONLY the significant genes we identified earlier
  filter(Gene %in% sig_diff_df$Gene) %>%
  # Create a combined column name for the x-axis: Sample_Region_Type
  mutate(Condition = paste(Sample, Region, Type, sep = "_")) %>%
  select(Gene, Condition, Methylation) %>%
  # Pivot to wide format (Rows = Genes, Columns = Conditions)
  pivot_wider(names_from = Condition, values_from = Methylation) %>%
  column_to_rownames("Gene") %>%
  as.matrix()

# 2. Data Cleaning
# Replace any remaining NAs with 0 (though unlikely after filtering)
heatmap_matrix[is.na(heatmap_matrix)] <- 0

# 3. Plotting with pheatmap
pheatmap(heatmap_matrix,
         
         # --- SCALE ---
         scale = "none",  # CRITICAL: Use raw 0-100 values, do not normalize/Z-score
         
         # --- CLUSTERING ---
         cluster_rows = TRUE,  # Cluster genes (Y-axis)
         cluster_cols = TRUE,  # Cluster conditions/samples (X-axis)
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         
         # --- COLORS ---
         # Define strict breaks from 0 to 100 so color is absolute
         breaks = seq(0, 100, length.out = 101),
         color = colorRampPalette(c("blue", "white", "red"))(100),
         
         # --- LABELS ---
         show_rownames = TRUE,
         fontsize_row = 12,        # Small font to fit many names
         show_colnames = TRUE,
         fontsize_col = 12,
         angle_col = "315",        # Rotate column names for readability
         
         # --- TITLES ---
         main = "Significant Genes Methylation")
