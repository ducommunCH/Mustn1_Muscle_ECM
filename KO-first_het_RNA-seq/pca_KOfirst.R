library(tidyverse)
library(DESeq2)
source("../general_functions.R")

dds <- readRDS("dds.rds")

### PCA plot

# data transformation
dds_vst <- DESeq2::vst(dds, blind = FALSE)

pca_data <- DESeq2::plotPCA(dds_vst, intgroup = c("group"), returnData = TRUE)
pca_data <- mutate(pca_data, sample = c("HET_1", "HET_2", "HET_3", "WT_1", "WT_2", "WT_3"))
percent_var <- round(100 * attr(pca_data, "percentVar"))
segments <- pca_data %>%
  group_by(!!as.symbol("group")) %>%
  summarise(xend = mean(PC1), yend = mean(PC2))
pca_data <- merge(pca_data, segments, by = "group")
no_colors <- pca_data[, "group"] %>% unique() %>% length()


ggplot(pca_data, aes(PC1, PC2, color = group)) +
  geom_point(size = 3) +
  ggtitle("Principal component analysis", subtitle = "(Soleus HET vs. WT)") +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  geom_segment(aes(PC1, PC2, xend = xend, yend = yend)) +
  geom_point(data = segments, aes(x = xend, y = yend), size = 2) +
  geom_text_repel(data = pca_data, parse = TRUE, aes(label = sample)) +
  default_theme() +
  scale_color_manual(values=c("#3182BDFF", "#E6550DFF")) +
  theme(legend.position="none")

ggsave_fixed("Fig2_Hets_PCA.pdf", plot_height = 5, plot_width = 5, height = 8, width = 9)

