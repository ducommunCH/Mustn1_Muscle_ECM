library(tidyverse)
library(pheatmap)
library(ggplotify)
source("../general_functions.R")

data_exp <- read_tsv("sample_expression.csv")

### Select the genes that contribute to term enrichment

# negative regulation of cell adhesion	GO:0007162
genes_adhesion <- str_to_title(c("CD74","H2-AA","THBS1","LAPTM5","TNC","ANXA1","NCKAP1L","IL4RA","MMP2","VSIG4","PTPRC","CD86","MMP14","FGL2","CYP1B1","SEMA3E","MYOC","SPI1","EPHB2"))
counts_adhesion <- data_exp %>% filter(SYMBOL %in% genes_adhesion) %>% 
  column_to_rownames("SYMBOL") %>% 
  dplyr::select(!(c("ENSEMBL", "ENTREZID", "GENENAME"))) %>%
  relocate(W1, W2, W3) %>% dplyr::rename("WT1" = W1, "WT2" = W2, "WT3" = W3, "HET1" = H1, "HET2" = H2, "HET3" = H3)

# regulation of vasculature development	GO:1901342
genes_vascular <- str_to_title(c("THBS1","SERPINF1","ECM1","EMILIN2","ITGB2","CCR2","C3AR1","PTGIS","CYBB","SFRP2","CMA1","ITGB3","DCN","C3","PF4","SPARC","CYP1B1","SEMA3E","CTSH"))
counts_vascular <- data_exp %>% filter(SYMBOL %in% genes_vascular) %>% 
  column_to_rownames("SYMBOL") %>% 
  dplyr::select(!(c("ENSEMBL", "ENTREZID", "GENENAME"))) %>%
  relocate(W1, W2, W3) %>% dplyr::rename("WT1" = W1, "WT2" = W2, "WT3" = W3, "HET1" = H1, "HET2" = H2, "HET3" = H3)

# collagen binding	GO:0005518
genes_collagen <- str_to_title(c("CTSS","THBS1","COL6A2","COL6A1","LOX","CTSK","PCOLCE","MRC2","PCOLCE2","ANTXR1","TGFBI","DPP4","DCN","SPARC","ITGA11","CRTAP","MAP1A"))
counts_collagen <- data_exp %>% filter(SYMBOL %in% genes_collagen) %>% 
  column_to_rownames("SYMBOL") %>% 
  dplyr::select(!(c("ENSEMBL", "ENTREZID", "GENENAME"))) %>%
  relocate(W1, W2, W3) %>% dplyr::rename("WT1" = W1, "WT2" = W2, "WT3" = W3, "HET1" = H1, "HET2" = H2, "HET3" = H3)

### Plot heatmaps ###
# Do separately for the 3 heatmaps
counts <- counts_adhesion
counts <- counts_vascular
counts <- counts_collagen

heatmap <- pheatmap::pheatmap(counts,
                              main = NA,
                              cluster_rows = T,
                              cluster_cols = F,
                              scale = "row",
                              clustering_distance_rows = "correlation",
                              # clustering_distance_cols = "euclidean",
                              show_rownames = T, show_colnames = T,
                              cellwidth = 5,
                              cellheight = 5,
                              treeheight_row = 5,
                              # treeheight_col = 20,
                              fontsize = 5,
                              border_color = NA,
                              gaps_col = 3,
                              # annotation_colors = pal_zfp697_mKO2,
                              # annotation_col = annot_cols,
                              legend = TRUE,
                              color = colorRampPalette(c("#3182BDFF", "white", "#E6550DFF"))(32)
)
as.ggplot(heatmap)

ggsave_fixed("Fig2_Heatmap_adhesion.pdf", plot_height = 5, plot_width = 5, height = 6, width = 6)
ggsave_fixed("Fig2_Heatmap_vascular.pdf", plot_height = 5, plot_width = 5, height = 6, width = 6)
ggsave_fixed("Fig2_Heatmap_collagen.pdf", plot_height = 5, plot_width = 5, height = 6, width = 6)

