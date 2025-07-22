library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(viridis)      
library(RColorBrewer)
library(cowplot)


order_df   <- read_tsv("order_by_tree.txt", col_types = cols())
sample_order <- order_df$sample

reorder_df <- function(df) {
  df %>%
    filter(sample %in% sample_order) %>% 
    mutate(sample = factor(sample, levels = rev(sample_order))) %>%
    arrange(sample)
}

### Tree rendered via ITOL online and saved as SVG
### Tree data below
(SW070919:0.0001339693,(SW07091:0.00005822084,TB06053:0.0001021984,SW070917:0.0004601076,MB06232:0.00009189716,PO06021:0.00004415348,(MB07152:0.0003191657,(SW071619:0.00003983498,SW071611:0.00003494364,SW07163:0.00003384095):0.0001121023[1]):0.0001101982[0.529801325],(SW071610:0.0003276469,(MB06021:0.00009676843,(PCP06101:0.0000411983,PCP06102:0.00003927395):0.0001029618[1]):0.0001017056[1]):0.00008378033[0.887417219]):0.0002081152[0.973509934],(MB06022:0.000212107,(PCP07083:0.0001942802,((PCP06173:0.000247636,(MB06233:0.0001396978,PCP06171:0.0001053928):0.0001736611[1]):0.0002498105[1],((PCP06021:0.0004863953,PO06101:0.0002820382):0.00026716[0.973509934],(TB06051:0.0004915913,(TB06052:0.0003300556,(PCP07014:0.0001609712,((PM06121:0.001588592,MN06021:0.0006454148):0.0002684115[1],(SW070912:0.0005330394,(SW070910:0.00004011681,SW071620:0.00003508365,PCP06238:0.00008006047):0.000274828[1]):0.000319134[1]):0.001249884[1]):0.0002236252[1]):0.0001225762[1]):0.0002873806[1],(MB06171:0.0003747471,(SW071614:0.00003686047,MB06173:0.00003936888,SW071618:0.0000359293,(SW07096:0.000139031,PCP06022:0.00004295396):0.0001017391[1],(PCP06235:0.0005238793,MB07081:0.00009357251,(MB06101:0.0002463126,MB06102:0.00003670895):0.00008505521[0.993377483]):0.0001007832[0.960264901],(MB07011:0.0001542031,(SW07164:0.0007176252,MB06234:0.0001400057):0.0001621799[1],(MB06231:0.0002167577,(MB06172:0.0001032897,SW07099:0.0001025036):0.00009766293[0.993377483]):0.00008880336[0.874172185]):0.00006891506[0.78807947]):0.0001057535[0.953642384]):0.0002244603[0.973509934]):0.000274164[0.973509934]):0.0002835989[1]):0.0001989811[0.973509934]):0.00008981592[0.973509934]);


# 0. Table plotting for MLST
mlst_df <- read_tsv("mlst.tsv", col_types = cols()) %>%
  rename(sample = sample, mlst = mlst) %>%
  reorder_df() %>%
  mutate(dummy = 1)

p0 <- ggplot(mlst_df, aes(x = dummy, y = sample)) +
  geom_tile(aes(fill = is.na(mlst)), color = "black") +
  scale_fill_manual(
    values = c(`TRUE` = "#EBEBEB", `FALSE` = "white"),
    guide  = FALSE
  ) +
  geom_text(aes(label = ifelse(is.na(mlst), "", mlst))) +
  scale_x_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.title   = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )


# 1. Heatmap for location.tsv
loc_df <- read_tsv("location.tsv", col_types = cols())
loc_df <- reorder_df(loc_df)
loc_df <- loc_df %>% mutate(dummy = 1)

custom_loc_palette <- c(
  "Fresh"  = "#F1A208",   
  "Ocean"  = "#052F5F",
  "Outfall"   = "#D5C67A",    
  "Sewage" = "#06A77D"     
)

p1 <- ggplot(loc_df, aes(x = dummy, y = sample, fill = category)) +
  geom_tile(width=2, color = "white") +
  scale_fill_manual(values = custom_loc_palette, na.value = "#EBEBEB") +
  theme_bw() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  )

# 2. Heatmap for clermont_type.tsv
clermont_df <- read_tsv("clermont_type.tsv", col_types = cols())
clermont_df <- reorder_df(clermont_df)
clermont_df <- clermont_df %>% mutate(dummy = 1)

custom_cler_palette <- c(
  "A"  = "#C6AED1",   
  "B1"  = "#CF9AAD",
  "B2"   = "#A3A00D",    
  "D" = "#BD457D",
  "E" = "#68457A",
  "G" = "#75B7B1"
)


p2 <- ggplot(clermont_df, aes(x = dummy, y = sample, fill = claremont)) +
  geom_tile(color="white") +
  scale_fill_manual(values = custom_cler_palette, na.value = "#EBEBEB") +
  theme_bw() +
  theme(axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1)
        )

# 3. Heatmap for Disk Assay
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

disk_df <- read_tsv("disk.tsv", col_types = cols()) %>%
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

p3 <- ggplot(disk_long, aes(x = disk, y = sample, fill = value)) +
  geom_tile(color="white") +
  scale_fill_manual(
    values   = c("R" = "#2B8CBE", "I" = "#A6BDDB", "S" = "#252525"),
    na.value = "#EBEBEB"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )

# 4. Combined Heatmap for RGI_CARD and AMRFinder+ results


col_order_2 <- c("SAMPLEID", "aadA1", "aadA5", "ANT(3'')-IIa", "aph(3'')-Ib", "aph(3')-Ia", "aph(6)-Id", "qacEdelta1", "qacJ", "qacL", "ampC", "CMY-132", "blaCTX-M-15", "blaCTX-M-3", "blaCTX-M-32", "blaCTX-M-65", "EC-13", "EC-14", "EC-15", "EC-18", "EC-8", "blaTEM", "blaTEM-1", "blaTEM-176", "TEM-244", "erm(B)", "emrR", "gyrA", "gyrA_D87N", "gyrA_S83A", "gyrA_S83L", "mdtH", "parC", "parC_S57T", "parC_S80I", "parE_S458A", "mph(A)", "Mrx", "acrE", "acrF", "acrS", "emrE", "evgA", "evgS", "gadW", "leuO", "mdtE", "mdtF", "mdtM", "mdtN", "mdtO", "mdtP", "ugd", "floR", "cyaA_S352T", "fosA3", "glpT", "uhpT", "uhpT_E350Q", "qnrB19", "QnrB5", "qnrS1", "sat2", "sul1", "sul2", "sul3", "tet(A)", "tet(B)", "tetR", "dfrA1", "dfrA14", "dfrA17", "dfrA51", "dfrA7")

amr_combined_df <- read_tsv("amr_rgi_combined_summary.tsv", col_types = cols()) %>%
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
  "G"  = "#F39B7F",
  "P"  = "#24325F",
  "GP" = "#008280"
)

p4 <- ggplot(amr_combined_long_2, aes(x = antibiotic, y = sample, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = custom_amr_pal, na.value = "#EBEBEB") +
  theme_bw() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x  = element_text(angle = 90, size = 18, hjust = 1)
  )


# Combined plot
plot_grid(p0,p1,p2,p3,p4, rel_widths = c(1,1.5,1,4,8), ncol=5)

