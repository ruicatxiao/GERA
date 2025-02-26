#!/usr/bin/env Rscript
# heatmap_generator.R
#
# This script reads in two AMR summary files (e.g., genome_amr_summary.tsv and plasmid_amr_summary.tsv),
# computes a common color scale (global min and max) across both data sets, and then generates heatmaps.
# The heatmaps are saved to PDF, PNG, and SVG files (with the same base file name as the input files).
#
# Requirements:
#   - The first column of each file is the sample ID (non-numeric).
#   - The remaining columns (13 of them) are numeric counts for each AMR class.
#   - The heatmap will not perform any clustering (rows and columns are in the file order).
#   - Column names are rotated 45°.
#   - The user may choose either the "viridis" or "magma" palette from the viridis package.
#
# Load required libraries
library(pheatmap)
library(viridis)
library(grid)
library(tools)  # for file_path_sans_ext

# -------------------------------
# USER SETTINGS: Specify file names
# -------------------------------
# Change these to the paths of your two output summary files.
genome_file  <- "/data/ruicatxiao/gera_analysis/amr/genome_amr_summary.tsv"
plasmid_file <- "/data/ruicatxiao/gera_analysis/amr/plasmid_amr_summary.tsv"

# Palette option: choose "magma" or "viridis"
palette_choice <- "magma"  # change to "viridis" if desired

# -------------------------------
# READ THE DATA
# -------------------------------
# Read in each file assuming a tab-delimited format and a header row.
genome_data <- read.table(genome_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
plasmid_data <- read.table(plasmid_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# The first column is the sample ID; the remaining columns are numeric counts.
# Convert the numeric data to a matrix.
genome_mat <- as.matrix(genome_data[,-1])
genome_mat <- apply(genome_mat, 2, as.numeric)

plasmid_mat <- as.matrix(plasmid_data[,-1])
plasmid_mat <- apply(plasmid_mat, 2, as.numeric)

# -------------------------------
# DETERMINE COMMON COLOR SCALE
# -------------------------------
# Compute the overall minimum and maximum values across both matrices.
global_min <- min(c(genome_mat, plasmid_mat))
global_max <- max(c(genome_mat, plasmid_mat))

# Create a sequence of breaks for the color scale (here we use 101 break points for 100 color intervals).
breaks <- seq(global_min, global_max, length.out = 101)

# Choose the color palette from viridis.
if(palette_choice == "viridis"){
  colors <- viridis(100)
} else if(palette_choice == "magma"){
  colors <- magma(100)
} else {
  stop("Invalid palette_choice. Please choose either 'viridis' or 'magma'.")
}

# -------------------------------
# FUNCTION TO CREATE & SAVE HEATMAP
# -------------------------------
# This function creates a heatmap using pheatmap (with no clustering, column labels at 45°) and saves the
# plot in PDF, PNG, and SVG formats. The output file names are based on the input file's base name.
create_heatmap <- function(mat, sample_ids, base_filename, title_text) {
  # Set the row names of the matrix to the sample IDs.
  rownames(mat) <- sample_ids
  
  # Generate the heatmap without clustering (default order is retained)
  # and with column labels rotated 45°.
  hm <- pheatmap(mat,
                 cluster_rows = FALSE,
                 cluster_cols = FALSE,
                 breaks = breaks,
                 color = colors,
                 angle_col = 315,
                 border_color = NA,
                 main = title_text,
                 silent = TRUE)
  
  # Extract the ggplot/gtable object from pheatmap
  g <- hm$gtable
  
  # Create and save the heatmap in PDF format.
  pdf_filename <- paste0(base_filename, ".pdf")
  pdf(pdf_filename, width = 8, height = 6)
  grid::grid.draw(g)
  dev.off()
  
  # Create and save the heatmap in PNG format.
  png_filename <- paste0(base_filename, ".png")
  png(png_filename, width = 800, height = 600)
  grid::grid.draw(g)
  dev.off()
  
  # Create and save the heatmap in SVG format.
  svg_filename <- paste0(base_filename, ".svg")
  svg(svg_filename, width = 8, height = 6)
  grid::grid.draw(g)
  dev.off()
  
  message("Heatmap saved as: ", pdf_filename, ", ", png_filename, ", and ", svg_filename)
}

# -------------------------------
# CREATE THE HEATMAPS
# -------------------------------
# Use tools::file_path_sans_ext() to remove the .tsv extension from the file names for output.
base_genome  <- file_path_sans_ext(genome_file)
base_plasmid <- file_path_sans_ext(plasmid_file)

# Generate heatmap for the genome data.
create_heatmap(genome_mat, genome_data[[1]], base_genome, "Genome AMR Heatmap")

# Generate heatmap for the plasmid data.
create_heatmap(plasmid_mat, plasmid_data[[1]], base_plasmid, "Plasmid AMR Heatmap")
