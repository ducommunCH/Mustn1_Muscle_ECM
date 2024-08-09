library(tidyverse)
library(ggrepel)
source("../general_functions.R")

# Read in the data matrix, non-imputed values
aorta_matrix <- read_tsv("SD_aorta_matrix_non-imputed.txt", comment = "#")
write_excel_csv(aorta_matrix, file = "Supplementary_File_3_Aorta-MS.csv", col_names = TRUE)
aorta_matrix <- aorta_matrix %>% mutate(Genes_trunc = str_extract(aorta_matrix$Genes, "[^;]+"))

aorta_curve <- read_tsv("SD_aorta_curve_non-imputed.txt", comment = "#")
aorta_curve <- aorta_curve %>% dplyr::rename(`Student's T-test Difference KO_WT` = x, `-Log Student's T-test p-value KO_WT` = y)

# No significant differences
ggplot(aorta_matrix, mapping = aes(x=`Student's T-test Difference KO_WT`, y=`-Log Student's T-test p-value KO_WT`)) +
  geom_point(data = aorta_matrix, colour = "#878787", size = 0.1) +
  geom_vline(xintercept = 0, linetype="dotted") +
  geom_line(data = aorta_curve, colour = "#3182BD", size = 0.5) +
  default_theme() +
  ggtitle("Aorta proteomics", subtitle="(FDR = 0.05, s0 = 0.1)") +
  ylab("-log10(p-value)") +
  xlab("Difference (KO-WT)") +
  scale_y_continuous(expand = c(0, 0)) +
  ylim(0, 6) + xlim(-2,2)

ggsave_fixed("Fig4_AortaMS1.pdf", plot_height = 5, plot_width = 5, height = 8, width = 8)


# Read in the data matrix, imputed values
aorta_matrix_i <- read_tsv("SD_aorta_matrix_imputed.txt", comment = "#", col_types = "ddddddddcccddddddccccc")
write_excel_csv(aorta_matrix_i, file = "Supplementary_File_4_Aorta-MS_imputed.csv", col_names = TRUE)
aorta_matrix_i <- aorta_matrix_i %>% mutate(Genes_trunc = str_extract(aorta_matrix_i$Genes, "[^;]+"))
aorta_curve_i <- read_tsv("SD_aorta_curve_imputed.txt", comment = "#")
aorta_curve_i <- aorta_curve_i %>% dplyr::rename(`Student's T-test Difference KO_WT` = x, `-Log Student's T-test p-value KO_WT` = y)


# 1 significant protein (Mustn1), color it and label
ggplot(aorta_matrix_i, mapping = aes(x=`Student's T-test Difference KO_WT`, y=`-Log Student's T-test p-value KO_WT`)) +
  geom_point(data = filter(aorta_matrix_i, is.na(`Student's T-test Significant KO_WT`)), colour = "#878787", size = 0.1) +
  geom_point(data = filter(aorta_matrix_i, `Student's T-test Significant KO_WT` == "+"), colour = "#3182BD", size = 0.5) +
  geom_vline(xintercept = 0, linetype="dotted") +
  geom_line(data = aorta_curve_i, colour = "#3182BD", size = 0.5) +
  geom_text_repel(data = filter(aorta_matrix_i, Genes == "Mustn1"), aes(label = "Mustn1")) +
  default_theme() +
  ggtitle("Aorta proteomics", subtitle="(imputed, FDR = 0.05, s0 = 0.1)") +
  ylab("-log10(p-value)") +
  xlab("Difference (KO-WT)") +
  scale_y_continuous(expand = c(0, 0)) +
  ylim(0, 6) + xlim(-6,6)

ggsave_fixed("Fig4_AortaMS2.pdf", plot_height = 5, plot_width = 5, height = 8, width = 8)
