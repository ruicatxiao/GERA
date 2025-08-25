#!/usr/bin/env Rscript
# heatmap_generator.R
#
# This script reads in two AMR summary files (e.g., genome_amr_summary.tsv and plasmid_amr_summary.tsv),
library(pheatmap)
library(viridis)
library(grid)
library(tools)  

genome_file  <- "/data/ruicatxiao/gera_analysis/amr/genome_amr_summary.tsv"
plasmid_file <- "/data/ruicatxiao/gera_analysis/amr/plasmid_amr_summary.tsv"

palette_choice <- "magma"  

genome_data <- read.table(genome_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
plasmid_data <- read.table(plasmid_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)


genome_mat <- as.matrix(genome_data[,-1])
genome_mat <- apply(genome_mat, 2, as.numeric)

plasmid_mat <- as.matrix(plasmid_data[,-1])
plasmid_mat <- apply(plasmid_mat, 2, as.numeric)


global_min <- min(c(genome_mat, plasmid_mat))
global_max <- max(c(genome_mat, plasmid_mat))
breaks <- seq(global_min, global_max, length.out = 101)
if(palette_choice == "viridis"){
  colors <- viridis(100)
} else if(palette_choice == "magma"){
  colors <- magma(100)
} else {
  stop("Invalid palette_choice. Please choose either 'viridis' or 'magma'.")
}

create_heatmap <- function(mat, sample_ids, base_filename, title_text) {
  rownames(mat) <- sample_ids
  hm <- pheatmap(mat,
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 breaks = breaks,
                 color = colors,
                 angle_col = 315,
                 border_color = NA,
                 main = title_text,
                 silent = TRUE)
  g <- hm$gtable
  pdf_filename <- paste0(base_filename, ".pdf")
  pdf(pdf_filename, width = 8, height = 6)
  grid::grid.draw(g)
  dev.off()
  png_filename <- paste0(base_filename, ".png")
  png(png_filename, width = 800, height = 600)
  grid::grid.draw(g)
  dev.off()
  svg_filename <- paste0(base_filename, ".svg")
  svg(svg_filename, width = 8, height = 6)
  grid::grid.draw(g)
  dev.off()
  
  message("Heatmap saved as: ", pdf_filename, ", ", png_filename, ", and ", svg_filename)
}
base_genome  <- file_path_sans_ext(genome_file)
base_plasmid <- file_path_sans_ext(plasmid_file)
create_heatmap(genome_mat, genome_data[[1]], base_genome, "Genome AMR Heatmap")
create_heatmap(plasmid_mat, plasmid_data[[1]], base_plasmid, "Plasmid AMR Heatmap")
