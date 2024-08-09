library(tidyverse)
source("../general_functions.R")

load("1_unloading_reloading_deseq2.RData")
# For details of upstream analysis, please refer to GSE237099 and associated publication Correia et al.

## Volcano plot.
ggplot(Reloading_vs_Unloading, mapping = aes(x=log2FoldChange, y=-log10(padj))) +
  geom_point(data = filter(Reloading_vs_Unloading, padj>=0.05), colour = "grey80", size = 0.1) +
  geom_point(data = filter(Reloading_vs_Unloading, padj<0.05, abs(log2FoldChange)<=0.58), colour = "grey80", size = 0.1) +
  geom_point(data = filter(Reloading_vs_Unloading, padj<0.05, abs(log2FoldChange)>0.58), colour = "grey40", size = 0.2) +
  geom_point(data = filter(Reloading_vs_Unloading, mgi_symbol == "Mustn1"), colour = "#E6550D", size = 1) +
  geom_point(data = filter(Reloading_vs_Unloading, mgi_symbol == "Hspa1b"| mgi_symbol == "Atf3"| mgi_symbol == "Tnfrsf12a" | mgi_symbol == "Cryab" | mgi_symbol == "Dnaja4" | mgi_symbol == "Enah" | mgi_symbol == "Rgcc"), colour = "#3182BDFF", size = 1) +
  ggtitle("Reloading vs. Unloading") +
  xlab("log2 Fold-Change") + 
  ylab("-log10 padj") +
  xlim(-8,8) + 
  default_theme() +
  geom_hline(yintercept = -log10(0.05), linetype="dotted") +
  geom_vline(xintercept = -0.58, linetype="dotted") +
  geom_vline(xintercept = 0.58, linetype="dotted") +
  geom_text_repel(data = filter(Reloading_vs_Unloading, mgi_symbol == "Mustn1" | mgi_symbol == "Hspa1b"| mgi_symbol == "Atf3"| mgi_symbol == "Tnfrsf12a" | mgi_symbol == "Cryab" | mgi_symbol == "Dnaja4" | mgi_symbol == "Enah" | mgi_symbol == "Rgcc"), aes(label = mgi_symbol, fontface = "italic"))


ggsave_fixed("Fig1_Volcano_more.pdf", plot_height = 5, plot_width = 5, height = 7, width = 7)