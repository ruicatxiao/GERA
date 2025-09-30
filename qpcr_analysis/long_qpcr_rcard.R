library(dplyr)
library(ggplot2)
library(forcats)
library(ggbreak)
library(scales)
library(RColorBrewer)
library(readxl)

# Library Versions
# readxl_1.4.3       scales_1.3.0       ggbreak_0.1.5      forcats_1.0.0     
# dplyr_1.1.4        ComplexUpset_1.3.3 RColorBrewer_1.1-3 ggplot2_3.5.1     


if (requireNamespace("ggplot2", quietly = TRUE)) {
  ggplot2_version <- packageVersion("ggplot2")
  
  # Check if version is newer than 3.5.1
  if (ggplot2_version > package_version("3.5.1")) {
    warning("This script was made with ggplot2 3.5.1. Your current version is newer (", 
            ggplot2_version, ") and may not work properly.")
  }
}

### Figure 2B ###

#LONGITUDINAL qPCR DATA BY SITE
KGV_qpcr_data <- read_excel("KGV_qpcr_data.xlsx")
KGV_qpcr_data$Date <- as.Date(KGV_qpcr_data$Date)
KGV_qpcr_data$Site <- fct_reorder(KGV_qpcr_data$Site, KGV_qpcr_data$Category, .fun = mean)
color_breaks <- c(0, 600, 800)
color_palette <- c("white", "#052f5f", "#052f5f")
values <- rescale(color_breaks, to = c(0,1))
qpcr_long_base <- ggplot(KGV_qpcr_data, aes(x=Date, y= Site, col=Human)) +
     geom_line(color="black") +
     geom_point(aes(fill= RFU_corrected_negatives), size = 5, shape = 21, stroke = 1, color="black") +
     scale_x_date(date_breaks = "1 week") +
    scale_fill_gradientn(colors = color_palette,
                         values = rescale(color_breaks, to = c(0,1)),
                         breaks = color_breaks)
print(qpcr_long_base)
qpcr_long_2 <- qpcr_long_base + 
  scale_x_break(c("2023-07-10", "2024-05-27")) +
  theme_bw() +
  theme(axis.text = element_text(size=15), 
        axis.title = element_text(size=0), 
        strip.text = element_text(size=11),
        legend.text = element_text(size=11),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
print(qpcr_long_2)
qpcr_long_3 <- qpcr_long_2 +
  theme(axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank())
print(qpcr_long_3)
svg("qpcr_long_3.svg", width = 24, height = 10)
print(qpcr_long_3)
dev.off()

### Figure S1 ###

#LONGITUDINAL RCard DATA BY SITE
KGV_rcard_data <- read_excel("KGV_rcard_data.xlsx")
KGV_rcard_data$Date <- as.Date(KGV_rcard_data$Date)
KGV_rcard_data$Site <- fct_reorder(KGV_rcard_data$Site, KGV_rcard_data$Category, .fun = mean)
rcolor_breaks <- c(0, 200, 5500)
color_palette <- c("white", "#052f5f", "#052f5f")
rvalues <- rescale(rcolor_breaks, to = c(0,1))

rcard_long_base <- ggplot(KGV_rcard_data, aes(x=Date, y= Site, col=Human)) +
  geom_line(color="black") +
  geom_point(aes(fill= Coliform), size = 5, shape = 21, stroke = 1, color="black") +
  scale_x_date(date_breaks = "1 week") +
  scale_fill_gradientn(colors = color_palette,
                       values = rescale(rcolor_breaks, to = c(0,1)),
                       breaks = rcolor_breaks)

print(rcard_long_base)
rcard_long_2 <- rcard_long_base + 
  scale_x_break(c("2023-07-10", "2024-05-27")) +
  theme_bw() +
  theme(axis.text = element_text(size=15), 
        axis.title = element_text(size=0), 
        strip.text = element_text(size=11),
        legend.text = element_text(size=11),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
print(rcard_long_2)
rcard_long_3 <- rcard_long_2 +
  theme(axis.text.x.top = element_blank(),
        axis.ticks.x.top = element_blank(),
        axis.line.x.top = element_blank())
print(rcard_long_3)

svg("rcard_long_3.svg", width = 12, height = 6)
print(rcard_long_3)
dev.off()