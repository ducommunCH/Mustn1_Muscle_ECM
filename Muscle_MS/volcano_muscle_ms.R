library(tidyverse)
library(ggrepel)
source("../general_functions.R")


muscle_matrix <- read_tsv("SD_muscle_matrix_non-imputed_ttest.txt", comment = "#")
muscle_matrix <- muscle_matrix %>% mutate(Genes_trunc = str_extract(muscle_matrix$Genes, "[^;]+"))

muscle_matrix$significant <- ""
muscle_matrix$significant[muscle_matrix$`Student's T-test Significant KO_WT` == "+" & muscle_matrix$Difference >= 0] <- "up"
muscle_matrix$significant[muscle_matrix$`Student's T-test Significant KO_WT` == "+" & muscle_matrix$Difference < 0] <- "down"

muscle_curve <- read_tsv("SD_muscle_curve_non-imputed.txt", comment = "#")
muscle_curve <- muscle_curve %>% dplyr::rename(Difference = x, `-Log(P-value)` = y)

# Extract list of significant proteins for pathway analysis
sig_prot_up <- muscle_matrix %>% 
  filter(`Student's T-test Significant KO_WT` == "+") %>% 
  filter(Difference > 0) %>%
  dplyr::select(Genes) %>% 
  filter(!is.na(Genes)) %>% 
  mutate(Genes = replace(Genes, Genes == "Mbnl1;Mbnl2", "Mbnl1"))
write_tsv(sig_prot_up, file = "sig_prot_up.txt", col_names = FALSE)

### Volcano plot for non-imputed values ###
muscle_matrix$label = ifelse(muscle_matrix$Genes == "Apcdd1" | muscle_matrix$Genes == "Tgfbi" | muscle_matrix$Genes == "Comp" | 
                               muscle_matrix$Genes == "Bdh1" | muscle_matrix$Genes == "Col12a1" | muscle_matrix$Genes == "Dpp4" |
                               muscle_matrix$Genes == "Pcolce",
                             muscle_matrix$Genes, "")

ggplot(muscle_matrix, mapping = aes(x=Difference, y=`-Log(P-value)`)) +
  geom_point(data = filter(muscle_matrix, significant == ""), colour = "#878787", size = 0.1) +
  geom_point(data = filter(muscle_matrix, significant == "up"), colour = "#E6550D", size = 0.5) +
  geom_point(data = filter(muscle_matrix, significant == "down"), colour = "#3182BD", size = 0.5) +
  geom_vline(xintercept = 0, linetype="dotted") +
  geom_line(data = muscle_curve, colour = "#878787", size = 0.5) +
  geom_text_repel(aes(label = label), 
                  max.overlaps = Inf, 
                  nudge_y = 0, 
                  nudge_x = 1.5, 
                  min.segment.length = 0) +
  default_theme() +
  ggtitle("Muscle proteomics", subtitle="(FDR = 0.05, s0 = 0.1)") +
  ylab("-log10(p-value)") +
  xlab("Difference (KO-WT)") +
  scale_y_continuous(expand = c(0, 0)) +
  ylim(0, 8) + xlim(-3, 3) +
  annotate("text",
           x = Inf, y = Inf,
           hjust = 17, vjust = 1, size = 5,
           color = "#3182BD",
           label = paste(
             muscle_matrix %>%
               filter(significant == "down") %>%
               nrow()
           )
  ) +
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1, vjust = 1, size = 5,
           color = "#E6550D",
           label = paste(
             muscle_matrix %>%
               filter(significant == "up") %>%
               nrow()
           )
  )
ggsave_fixed("Fig6_MSVolcano.pdf", plot_height = 7.5, plot_width = 7.5, height = 10, width = 9)


  