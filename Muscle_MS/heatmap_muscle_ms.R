library(tidyverse)
library(pheatmap)
library(ggplotify)
source("../general_functions.R")

# Read in the expression values, use the imputed values
muscle_matrix <- read_tsv("SD_muscle_matrix_imputed_ttest.txt", comment = "#")

# Average the technical replicates
muscle_matrix <- muscle_matrix %>% mutate(Genes_trunc = str_extract(muscle_matrix$Genes, "[^;]+")) %>%
  #filter(`Student's T-test Significant KO_WT` == "+") %>%
  filter(!is.na(Genes)) %>%
  dplyr::rename_with(stringr::str_replace, 
                     pattern = "E:\\\\Mustang muscle\\\\TAduplicates\\\\20230322_EXPL4_MaMu_SA_M770_TA", replacement = "") %>%
  dplyr::rename_with(stringr::str_replace, pattern = ".raw", replacement = "") %>%
  mutate(KO1 = (rep1_KO_1 + rep2_KO_1 + rep3_KO_1)/3, .keep = "unused", .before = 1) %>%
  mutate(KO2 = (rep1_KO_2 + rep2_KO_2 + rep3_KO_2)/3, .keep = "unused", .before = 1) %>%
  mutate(KO3 = (rep1_KO_3 + rep2_KO_3 + rep3_KO_3)/3, .keep = "unused", .before = 1) %>%
  mutate(KO4 = (rep1_KO_4 + rep2_KO_4 + rep3_KO_4)/3, .keep = "unused", .before = 1) %>%
  mutate(WT1 = (rep1_WT_1 + rep2_WT_1 + rep3_WT_1)/3, .keep = "unused", .before = 1) %>%
  mutate(WT2 = (rep1_WT_2 + rep2_WT_2 + rep3_WT_2)/3, .keep = "unused", .before = 1) %>%
  mutate(WT3 = (rep1_WT_3 + rep2_WT_3 + rep3_WT_3)/3, .keep = "unused", .before = 1) %>%
  mutate(WT4 = (rep1_WT_4 + rep2_WT_4 + rep3_WT_4)/3, .keep = "unused", .before = 1)

### Select the genes that contribute to term enrichment
# extracellular matrix GO:0031012
genes <- str_to_title(c("TNXB","MYOC","ANXA1","FN1","BGN","FMOD","TGFBI","SERPINF1","COL6A2",
                        "COL6A1","COL12A1","OGN","POSTN","TNC","COL6A6","ASPN","PRELP","DPT","COMP","THBS4"))

counts <- muscle_matrix %>% filter(Genes_trunc %in% genes) %>% 
  column_to_rownames("Genes_trunc") %>% 
  dplyr::select(c("WT1", "WT2", "WT3", "WT4", "KO1", "KO2", "KO3", "KO4"))

### Plot heatmaps ###
heatmap <- pheatmap::pheatmap(counts,
                              main = NA,
                              cluster_rows = T,
                              cluster_cols = F,
                              scale = "row",
                              clustering_distance_rows = "correlation",
                              show_rownames = T, show_colnames = T,
                              cellwidth = 5,
                              cellheight = 5,
                              treeheight_row = 5,
                              fontsize = 5,
                              border_color = NA,
                              gaps_col = 4,
                              legend = TRUE,
                              color = colorRampPalette(c("#3182BDFF", "white", "#E6550DFF"))(32)
)
as.ggplot(heatmap)

ggsave_fixed("Fig6_Heatmap.pdf", plot_height = 5, plot_width = 5, height = 6, width = 6)

