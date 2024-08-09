library(tidyverse)
library(Seurat)
source("../general_functions.R")

# For data download refer to: https://doi.org/10.1038/s42003-021-02810-x

ids <- c("prefilter_factorIDs", "harmony_factorIDs", "bbknn_factorIDs", "scanorama_factorIDs")
y_lab <- expression(paste(italic("Mustn1"), " expression"))
sub_title <- expression(paste("(McKellar ", italic("et al."), ")"))

### Plot for paper (all cells) ###
load("~/datasets/mckellar_scmuscle/scMuscle_mm10_slim_v1-1.RData")
# Make sure there is enough RAM available, takes 23.1GB
cols_pal_15 <- c("#636363FF", "#636363FF", "#636363FF", "#636363FF", "#636363FF",
                         "#636363FF", "#636363FF", "#3182BDFF", "#E6550DFF", "#636363FF",
                         "#31A354FF", "#636363FF", "#636363FF", "#756BB1FF", "#636363FF")
cols_pal_23 <- c("#636363FF", "#636363FF", "#636363FF", "#636363FF", "#636363FF",
                           "#636363FF", "#636363FF", "#636363FF", "#636363FF", "#636363FF",
                           "#3182BDFF", "#3182BDFF", "#3182BDFF", "#E6550DFF", "#31A354FF",
                           "#636363FF", "#636363FF", "#636363FF", "#636363FF", "#636363FF",
                           "#756BB1FF", "#636363FF", "#636363FF")

VlnPlot(scMuscle.slim.seurat, group.by = ids[2], features = c("Mustn1"), pt.size = 0, cols = cols_pal_23)  +
  default_theme() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("") + ylab(y_lab) + NoLegend() + ggtitle("Muscle scRNA-seq", subtitle=sub_title)

ggsave_fixed("Fig5_MuscleSingleCell.pdf", plot_height = 5, plot_width = 16, height = 12, width = 18)



### Plot for paper (myogenic lineage) ###
load("~/datasets/mckellar_scmuscle/myo_slim_seurat_v1-1.RData")
# Make sure there is enough RAM available, takes 4.2GB

VlnPlot(myo.slim.seurat, group.by = ids[4], features = c("Mustn1"), pt.size = 0, cols = c("#636363FF", "#636363FF", "#636363FF", "#636363FF", "#3182BDFF", "#E6550DFF", "#636363FF", "#636363FF", "#31A354FF"))  +
  default_theme() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("") + ylab(y_lab) + NoLegend() + ggtitle("Myogenic lineage scRNA-seq", subtitle=sub_title)

ggsave_fixed("Fig5_MyogenicSingleCell.pdf", plot_height = 5, plot_width = 7, height = 12, width = 9)


