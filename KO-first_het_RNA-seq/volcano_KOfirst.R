library(tidyverse)
library(ggrepel)
source("../general_functions.R")

result_array_ids <- readRDS("result_array_ids.rds")
result_df <- result_array_ids$mustn1_hetz_vs_mustn1_wt

write_excel_csv(result_df, file = "Supplementary_File_2_RNA-seq_KO-first_het.csv", col_names = TRUE)

result_df <- result_df %>% drop_na() %>%
  mutate(significant = case_when(
    padj < 0.05 & log2FoldChange > 0 ~ "up",
    padj < 0.05 & log2FoldChange < 0 ~ "down",
    TRUE ~ "unchanged"))

sig_prot_up <- result_df %>% 
  filter(significant == "up") %>%
  dplyr::select(SYMBOL)
write_tsv(sig_prot_up, file = "sig_prot_up.txt", col_names = FALSE)

## Volcano plot ###
ggplot(result_df, mapping = aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(data = filter(result_df, padj>=0.05), colour = "grey80", size = 0.1) +
  geom_point(data = filter(result_df, padj<0.05, log2FoldChange<0), colour = "#3182BDFF", size = 0.5) +
  geom_point(data = filter(result_df, padj<0.05, log2FoldChange>0), colour = "#E6550DFF", size = 0.5) +
  ggtitle("Differential gene expression", subtitle = "(Soleus HET vs. WT)") +
  xlab("log2 Fold-Change") + 
  ylab("-log10 padj") +
  xlim(-1.1,1.1) + 
  ylim(0, 15) +
  geom_hline(yintercept = -log10(0.05), linetype="dotted") +
  geom_vline(xintercept = 0, linetype="dotted") +
  theme(text = element_text(family="Helvetica", color="black"),
        axis.text = element_text(size=6, colour = "black"),
        axis.ticks = element_line(size = 0.33, colour = "black"),
        plot.title = element_text(size = 7, hjust = 0.5),
        axis.title = element_text(size = 7),
        panel.background = element_rect(fill = "white"), 
        panel.border = element_rect(linetype = "solid", color = "black", size=0.75, fill = NA), 
        panel.grid.major = element_line(color = "grey93"), 
        panel.grid.minor = element_line(color = NA)) +
  geom_text_repel(data = filter(result_df, SYMBOL == "Mustn1"), parse = TRUE, aes(label = "italic('Mustn1')")) +
  annotate("text",
           x = Inf, y = Inf,
           hjust = 17, vjust = 1, size = 5,
           color = "#3182BDFF",
           label = paste(
             result_df %>%
               filter(significant == "down") %>%
               nrow()
           )
  ) +
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1, vjust = 1, size = 5,
           color = "#E6550DFF",
           label = paste(
             result_df %>%
               filter(significant == "up") %>%
               nrow()
           )
  ) +
  default_theme()

ggsave_fixed("Fig2_Hets_Volcano.pdf", plot_height = 5, plot_width = 5, height = 7, width = 7)


