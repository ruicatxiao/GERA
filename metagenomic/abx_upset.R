# Packages 
library(ggplot2) 
library(RColorBrewer) 
library(ComplexUpset) 

# ComplexUpset_1.3.3 RColorBrewer_1.1-3 ggplot2_3.5.1 

if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2_version <- packageVersion("ggplot2")
  
  # Check if version is newer than 3.5.1
  if (ggplot2_version > package_version("3.5.1")) {
    warning("This script was made with ggplot2 3.5.1. Your current version is newer (", 
            ggplot2_version, ") and may not work properly.")
  }
}

# reading in the data
FULL <- read.csv("kirby_bauer_data_amr.tsv", header = TRUE, sep = "\t")

#specifying names used for the antibiotic labels
antibiotics <-c("AmoxicillinClavulanicAcid","Azithromycin",
                "Cefoxitin","Ceftriaxone","Chloramphenicol","Gentamicin",
                "Ampicillin","Meropenem","Streptomycin","Tetracycline", 
                "TrimethoprimSulfamethoxazole", "Ciprofloxacin")

# constructing the plot with ggplot2 
svg("abx_upset.svg", width = 18, height = 6)

upset(
  FULL,
  antibiotics,
  name = "Resistance combinations of the twelve antibiotics",
  sort_intersections_by = "degree",  # sort by size
  sort_intersections = "ascending",              # let it sort them
  base_annotations = list(
    'Intersection size' = intersection_size(
      counts = FALSE,
      mapping = aes(fill = category)
    ) +
      geom_bar(stat = "identity", colour = "black") +
      scale_x_reverse() +                   # <-- reverse the axis
      ylab('Number of isolates') +
      labs(fill = "category") +
      scale_y_continuous()
  ),
  width_ratio = 0.1,
  stripes = 'white'
)

dev.off()
