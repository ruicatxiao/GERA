library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(ggrepel)

# read data
df <- read.table("all_cord_freq.txt", header = FALSE, sep = "\t",
                 col.names = c("Sample","Gene","Position","AltFreq"),
                 stringsAsFactors = FALSE)

# 1. Bar chart of total SNP counts per sample
counts_df <- df %>% 
  group_by(Sample) %>% 
  summarise(SNP_Count = n())

ggplot(counts_df, aes(x = Sample, y = SNP_Count)) +
  geom_col() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Total SNP Counts per Sample",
       x     = "Sample",
       y     = "SNP Count")

# Saved to blatem1_snpcounts_per_sample

# 2. Violin plot of allele‐frequency distribution per sample
ggplot(df, aes(x = Sample, y = AltFreq)) +
  geom_violin(trim = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Distribution of SNP Allele Frequencies by Sample",
       x     = "Sample",
       y     = "Allele Frequency (%)")

# Saved to blatem1_snpfreq_per_sample

# 3. SNP counts per 10%‐gene‐body bins, colored by sample
df <- read.table("all_cord_freq.txt",
                 header    = FALSE,
                 sep       = "\t",
                 col.names = c("Sample","Gene","Position","AltFreq"),
                 stringsAsFactors = FALSE)

df <- df %>%
  mutate(Bin = cut(Position,
                   breaks = seq(0, 100, by = 10),
                   include.lowest = TRUE,
                   labels = paste0(seq(0, 90, 10), "-", seq(10, 100, 10))
  ))

bin_counts <- df %>%
  group_by(Bin, Sample) %>%
  summarise(SNP_Count = n(), .groups = "drop")

ggplot(bin_counts, aes(x = Bin, y = SNP_Count, fill = Sample)) +
  geom_col() +   # default is position = "stack"
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "SNP Counts per Gene‐Body Bin (stacked by Sample)",
       x     = "Gene Body (%)",
       y     = "SNP Count")

# Saved to blatem1_snpcounts_stackedsample

# 4 Kernal Density
ggplot(df, aes(x = Position, color = Sample)) +
  geom_density() +
  theme_bw() +
  labs(title = "Kernel Density of SNP Positions by Sample",
       x     = "Position on Gene (%)",
       y     = "Density")

# Saved to blatem1_kernal_density


# 5 Heatmap and clustering based on bin x sample similarity
heat_df <- df %>%
  group_by(Bin, Sample) %>%
  summarise(MeanFreq = mean(AltFreq), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = MeanFreq)

mat <- as.matrix(heat_df[,-1])
rownames(mat) <- heat_df$Bin

pheatmap(mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Mean Allele Frequency Heatmap\n(bin × sample)")

# Saved to blatem1_heatmap

# 6. PCA on per‐bin mean frequencies ───────────────────────────
# transpose so rows = samples, cols = bins
pca_data <- t(mat)
# replace NA with 0 (or you could impute)
pca_data[is.na(pca_data)] <- 0

p <- prcomp(pca_data, scale. = TRUE)
pca_df <- as.data.frame(p$x[,1:2])
pca_df$Sample <- rownames(pca_df)

ggplot(pca_df, aes(x = PC1, y = PC2, label = Sample)) +
  geom_point(size = 3) +
  geom_text_repel() +
  theme_bw() +
  labs(title = "PCA of Samples (on per‐bin allele frequencies)",
       x     = "PC1",
       y     = "PC2")

# saved to blatem1_pca

# 7 Hotspot
n_samples <- length(unique(df$Sample))

hotspot_df <- df %>%
  group_by(Position) %>%
  summarise(TotalSNPs   = n(),
            Samples     = n_distinct(Sample),
            MeanFreq    = mean(AltFreq),
            MaxFreq     = max(AltFreq),
            .groups     = "drop") %>%
  arrange(desc(TotalSNPs))

hotspots <- hotspot_df %>%
  filter(Samples == 16,   # seen in every sample
         MeanFreq >= 1)           # mean alt‐freq > 1%

ggplot(hotspots, aes(x = Position, y = TotalSNPs)) +
  geom_col() +
  theme_bw() +
  labs(title = "Hotspot SNP Positions\n(in all samples & mean freq > 1%)",
       x     = "Gene Position (%)",
       y     = "Total SNP")

# saved to blatem1_hotspot

ggplot(hotspots, aes(x = Position, y = MeanFreq)) +
     geom_point(color = "steelblue") +
     geom_smooth(method = "loess", formula = y ~ poly(x, 5), se = TRUE, linetype = "dashed", color = "grey") +
     theme_bw() + ylim (0,2) +
     labs(title = "Hotspot SNP Positions\n(in all samples & mean freq > 1%)",
          x     = "Gene Position (%)",
          y     = "Mean Alternative Allele Frequency (%)")
