library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(viridis)      
library(RColorBrewer) 


order_df   <- read_tsv("order_by_tree_v2.txt", col_types = cols())
sample_order <- order_df$sample




### 1. Heatmap for location.tsv
loc_df <- read_tsv("location_v2.tsv", col_types = cols())
loc_df <- reorder_df(loc_df)

# For a one‐column heatmap, create a dummy column (x = 1)
loc_df <- loc_df %>% mutate(dummy = 1)

custom_loc_palette <- c(
  "Fresh"  = "skyblue",   
  "Ocean"  = "dodgerblue",
  "Pipe"   = "orange",    
  "Sewage" = "sienna"     
)

ggplot(loc_df, aes(x = dummy, y = sample, fill = category)) +
  geom_tile(color="white") +
  scale_fill_manual(values = custom_loc_palette, na.value = "grey") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )


### 2. Heatmap for clermont_type.tsv
clermont_df <- read_tsv("clermont_type_v2.tsv", col_types = cols())
clermont_df <- reorder_df(clermont_df)
clermont_df <- clermont_df %>% mutate(dummy = 1)

ggplot(clermont_df, aes(x = dummy, y = sample, fill = claremont)) +
  geom_tile(color="white") +
  scale_fill_brewer(palette = "Set2", na.value = "grey") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())

# Define the exact column order you want
col_order <- c(
  "sample",
  "Gentamicin",
  "Streptomycin",
  "Amoxicillin/Clavulanic_Acid",
  "Ampicillin",
  "Cefoxitin",
  "Ceftriaxone",
  "Meropenem",
  "Ciprofloxacin",
  "Trimethoprim/Sulfamethoxazole",
  "Azithromycin",
  "Chloramphenicol",
  "Tetracycline"
)


reorder_df <- function(df) {
  df %>%
    filter(sample %in% sample_order) %>% 
    mutate(sample = factor(sample, levels = rev(sample_order))) %>%
    arrange(sample)
}

disk_df <- read_tsv("disk_v4.tsv", col_types = cols()) %>%
  select(all_of(col_order)) %>%
  reorder_df()

disk_long <- disk_df %>%
  pivot_longer(
    cols = -sample,
    names_to = "disk",
    values_to = "value"
  ) %>%
  mutate(
    disk = factor(disk, levels = col_order[-1])
  )

ggplot(disk_long, aes(x = disk, y = sample, fill = value)) +
  geom_tile(color="white") +
  scale_fill_manual(
    values   = c("R" = "#D61A46", "I" = "#FDDC22", "S" = "#236AB9"),
    na.value = "white"
  ) +
  theme_bw() +
  theme(
    axis.text.x    = element_text(angle = 45, hjust = 1)
  )

# Now plotting amr_combined_v4_summary.tsv


col_order_2 <- c("SAMPLEID","aadA1","aadA5","aph(3'')-Ib","aph(3')-Ia","aph(6)-Id","blaCTX-M-15","blaCTX-M-3","blaCTX-M-32","blaCTX-M-65","blaEC","blaTEM","blaTEM-1","blaTEM-176","cyaA_S352T","fosA3","uhpT_E350Q","erm(B)","mph(A)","floR","gyrA_D87N","gyrA_S83A","gyrA_S83L","parC_S57T","parC_S80I","parE_S458A","qnrB19","qnrS1","sat2","sul1","sul2","sul3","tet(A)","tet(B)","dfrA1","dfrA14","dfrA17","dfrA51","dfrA7")

amr_combined_df <- read_tsv("amr_combined_v4_summary.tsv", col_types = cols()) %>%
  rename(sample = SAMPLEID) %>%
  reorder_df()

amr_combined_long <- amr_combined_df %>%
  pivot_longer(
    cols = -sample,
    names_to = "antibiotic",
    values_to = "value"
  ) %>%
  mutate(
    antibiotic = factor(antibiotic, levels = col_order_2[-1]),
    value      = factor(value,      levels = c("G", "P", "GP"))
  )

custom_amr_pal <- c(
  "G"  = "#06A77D",
  "P"  = "#052F5F",
  "GP" = "#F1A208"
)

ggplot(amr_combined_long, aes(x = antibiotic, y = sample, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = custom_amr_pal, na.value = "grey") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x  = element_text(angle = 45, size = 18, hjust = 1)
  )