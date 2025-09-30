################################################################################

###Figure3A###
library(tidyverse)       # loads dplyr, readr, ggplot2, etc.
library(pheatmap)
library(viridis)

# Load and filter function
process_run <- function(data, sample_name) {
  df <- data %>%
    filter(nt_percent_identity >= 85, nt_alignment_length > 500) %>%
    select(tax_id, nt_bpm, tax_level, name, category) %>%
    filter(!is.na(nt_bpm), nt_bpm > 50)
  
  names(df)[names(df) == "nt_bpm"] <- sample_name
  return(df)
}

# Load CSVs
run1.raw <- read.csv("22_pier_calls_755619_taxon_report.csv")
run2.raw <- read.csv("22_espanola_calls_755613_taxon_report.csv")
run3.raw <- read.csv("22_baquerizo_calls_755612_taxon_report.csv")
run4.raw <- read.csv("22_mann_calls_755618_taxon_report.csv")
run5.raw <- read.csv("23_carola_pipe_calls_755614_taxon_report.csv")
run6.raw <- read.csv("23_mann_calls_755615_taxon_report.csv")
run7.raw <- read.csv("23_marinos_calls_755616_taxon_report.csv")
run8.raw <- read.csv("23_pier_calls_755617_taxon_report.csv")
run9.raw <- read.csv("23_sewage_maremio_calls_755610_taxon_report.csv")
run10.raw <- read.csv("23_sewage_sign_calls_755611_taxon_report.csv")

# Apply consistent processing
r1 <- process_run(run1.raw, "pier_22")
r2 <- process_run(run2.raw, "espanola_22")
r3 <- process_run(run3.raw, "baquerizo_22")
r4 <- process_run(run4.raw, "mann_22")
r5 <- process_run(run5.raw, "carola_pipe_23")
r6 <- process_run(run6.raw, "mann_23")
r7 <- process_run(run7.raw, "marinos_23")
r8 <- process_run(run8.raw, "pier_23")
r9 <- process_run(run9.raw, "sewage_sw_23")
r10 <- process_run(run10.raw, "sewage_ne_23")

# Join all
all.runs <- reduce(
  list(r1, r2, r3, r4, r5, r6, r7, r8, r9, r10),
  ~ full_join(.x, .y, by = c("tax_id", "tax_level", "name", "category"))
)

# Remove names with "taxa" or "genus" or "uncultured"
all.runs <- all.runs %>%
  filter(!grepl("taxa|genus|uncultured", name))

# Select all taxa at Genus Level
all.runs.g <- all.runs %>%
  filter(tax_level == "2", category %in% c("eukaryota", "viruses", "bacteria", "archaea")) %>%
  select(name, pier_22, pier_23, mann_22, mann_23, marinos_23,
         carola_pipe_23, baquerizo_22, espanola_22, sewage_sw_23, sewage_ne_23)

# Set rownames, convert to matrix, handle NA
rownames(all.runs.g) <- all.runs.g$name
all.runs.g.mat <- all.runs.g %>%
  select(-name) %>%
  as.matrix()

all.runs.g.mat <- all.runs.g.mat
all.runs.g.mat[is.na(all.runs.g.mat)] <- 0.1

# Get top 50 taxa BEFORE log-transform
taxa_totals_g <- rowSums(all.runs.g.mat)
top50_taxa_g <- names(sort(taxa_totals_g, decreasing = TRUE))[1:50]

# Subset the matrix
all.runs.g.mat.top50 <- all.runs.g.mat[top50_taxa_g, ]

# Now log-transform
all.runs.g.mat.top50 <- log2(all.runs.g.mat.top50)

svg("czid_all_genus.svg", width = 5, height = 10) 

# Heatmap
pheatmap(
  all.runs.g.mat.top50,
  cluster_cols = TRUE,
  color = viridis(10),
  cellwidth = 7,
  cellheight = 7,
  fontsize = 5,
  fontsize_row = 7,
  fontsize_col = 7,
  treeheight_row = 20,
  treeheight_col = 20,
  main = "All Taxa \nGenus Level"
)

dev.off()

################################################################################

###Figure 3B###

# Load required libraries
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(showtext)
showtext_auto()  # Enables custom fonts
font_add("Arial", regular = "arial.ttf")  

# Load your data
df <- read.delim("amrfinder_combined_cleaned_R.tsv", header = TRUE, sep = "\t")

# Ensure uniqueness of Gene-Subclass mapping (important for clustering)
gene_class <- df %>%
  distinct(Gene, Class)

# Create a full list of samples and genes for presence/absence matrix
all_samples <- unique(df$Name)
all_genes <- unique(df$Gene)

# Count gene occurrences per sample
heatmap_data <- df %>%
  dplyr::count(Name, Gene) %>%
  tidyr::complete(Name = all_samples, Gene = all_genes, fill = list(n = 0)) %>%
  tidyr::pivot_wider(names_from = Name, values_from = n, values_fill = 0) %>%
  tibble::column_to_rownames("Gene")

# Reorder rows based on class
row_order <- gene_class %>%
  filter(Gene %in% rownames(heatmap_data)) %>%
  arrange(Class, Gene)

heatmap_data <- heatmap_data[row_order$Gene, , drop = FALSE]

# Optional: assign row annotation by class
row_annotation <- gene_class %>%
  filter(Gene %in% rownames(heatmap_data)) %>%
  column_to_rownames("Gene")

heatmap_data_binary <- (heatmap_data > 0) * 1

svg("amr_gene_presence.svg", width = 10, height = 30) 

# Generate the heatmap
x<- pheatmap(
  heatmap_data_binary,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  clustering_distance_cols = "binary",
  cellwidth = 5,
  cellheight = 5,
  annotation_row = row_annotation,
  color = c("white", "#005377"),
  legend_breaks = c(0, 1),
  legend_labels = c("Absent", "Present"),
  fontsize_row = 5,
  fontsize_col = 5,
  treeheight_col = 10,
  main = "AMR Gene Presence/Absence by Sample"
)

print(x)
dev.off()

################################################################################

###Figure3C###

library(ggplot2)
library(readr)

#data from seqkit_stats.txt

amr_data <- data.frame(
  Site = c("22 Playa Baquerizo", "22 Espanola", "22 Muelle de Pescadores", "22 Playa Mann", 
           "23 Muelle de Pescadores", "23 Playa Mann", "23 Punta Carola Pipe", "23 Laguna", 
           "23 Sewage SW", "23 Sewage NE"),
  AMR_per_10_Mbp = c(0.10, 0.00, 0.92, 0.00, 0.00, 0.18, 6.27, 0.42, 7.56, 8.05),
  year = c(2022, 2022, 2022, 2022, 2023, 2023, 2023, 2023, 2023, 2023),
  category = c("Ocean", "Ocean", "Ocean", "Ocean", "Ocean", "Ocean", "Outfall", "Outfall", "Sewage", "Sewage")
)

# Save plot
svg("amr_gene_density.svg", width = 10, height = 6) 

y <- ggplot(amr_data, aes(x = AMR_per_10_Mbp, y = reorder(Site, AMR_per_10_Mbp), fill = category)) +
  geom_col() +
  labs(
    title = "AMR Gene Density in Metagenome Samples",
    x = "AMR Genes per 10 Mb MAG",
    y = "Sampling Site"
  ) +
  scale_fill_manual(values = c("Ocean" = "#052f2f", "Outfall" = "#06a77d", "Sewage" = "#f1a208")) +
  theme_minimal(base_size = 14)

print(y)
dev.off()

################################################################################

###Figure3D###

# Load required libraries
library(pheatmap)
library(RColorBrewer)

# Read your matrix from CSV (from seqkit_stats.txt)
amr_matrix <- read.csv("amr_per_10mbp_matrix.csv", row.names = 1)

# Scale rows or log-transform if needed

amr_matrix <- log2(amr_matrix + .001)  # Use this if values vary widely

svg("amr_class_density.svg", width = 10, height = 5) 

# Create heatmap
pheatmap(
  mat = amr_matrix,
  color = viridis(20),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  border_color = "gray80",
  fontsize_row = 8,
  fontsize_col = 8,
  main = "AMR Gene Class Density per 10 Mb MAGs",
  angle_col = 45,
  cellwidth = 12,
  cellheight = 12,
  legend = TRUE,
  treeheight_row = 20,
  treeheight_col = 20,
  display_numbers = FALSE # change to TRUE to label values
)

dev.off()
