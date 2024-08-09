library(eulerr)
library(ggplotify)
source("../general_functions.R")

# Make a Venn diagram with 220 RNA-seq up, 49 Mass-spec up. 15 in common.

list_common <- "Pcolce\nMcpt4\nMyoc\nCol6a1\nCma1\nPrelp\nDpp4\nTnc\nAnxa1\nDpysl3\nDpt\nBgn\nTgfbi\nCol6a2\nSerpinf1\n"


venn_2 <- euler(c("A" = 205, "B" = 34, "A&B" = 15))
p1 <- plot(venn_2, labels = list(labels = c("Transcriptomics", "Proteomics"), font = 1), fills = c("#3182BDFF", "#E6550DFF"), edges = list(lty = 1),
     quantities = list(type = "counts"))

p2 <- as.ggplot(p1) +
  ggtitle("Genes/proteins increased in Mustn1-HET/KO", subtitle = "Mustn1-HET (RNA-seq) or Mustn1-KO (MS) vs. WT") +
  theme(
  plot.title = element_text(size = 16, hjust = 0.5),
  plot.subtitle = element_text(size = 12, hjust = 0.5)) +
  annotate("text", x = 0, y = .5, hjust = 0, 
           label = list_common, size = 10/.pt,
           color = 'black')

ggsave_fixed("Fig6_Venn.pdf", plot = p2, plot_height = 8, plot_width = 10, height = 10, width = 14)


