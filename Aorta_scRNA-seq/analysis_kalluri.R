library(tidyverse)
library(Seurat)
source("../general_functions.R")

# For data download refer to DOI: 10.1161/CIRCULATIONAHA.118.038362
meta <- read_delim("META_DATA_Chow_12PCs_outfile.dms")
counts <- read_delim("Seurat_Chow_12PCs_outfile.mtx")

meta <- meta %>% dplyr::filter(!(NAME == "TYPE")) %>% column_to_rownames("NAME")
counts <- column_to_rownames(counts, "GENE")

kalluri_seurat <- CreateSeuratObject(counts, project = "kalluri", assay = "RNA", meta.data = meta)

kalluri_seurat <- SetIdent(kalluri_seurat, value = "Cluster")
kalluri_seurat <- SetIdent(kalluri_seurat, value = "Sub.Cluster")
levels(kalluri_seurat)

### Plot for paper ###
y_lab <- expression(paste(italic("Mustn1"), " expression"))
sub_title <- expression(paste("(Kalluri ", italic("et al."), ", SCP289)"))

VlnPlot(kalluri_seurat, features = c("Mustn1"), cols = c("#3182BDFF", "#E6550DFF", "#31A354FF", "#756BB1FF", "#636363FF", "#6BAED6FF"))  +
  default_theme() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("") + ylab(y_lab) + NoLegend() + ggtitle("Aorta scRNA-seq", subtitle=sub_title)

ggsave_fixed("Fig4_AortaSingleCell.pdf", plot_height = 5, plot_width = 10, height = 9, width = 13)
