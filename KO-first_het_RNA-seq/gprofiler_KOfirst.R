#https://biit.cs.ut.ee/gplink/l/_fFeAsQFSv

library(tidyverse)
library(gprofiler2)
source("../general_functions.R")


# Select list of proteins increased in Mustn1-HET for enrichment analysis
sig_prot_up <- read_tsv(file = "sig_prot_up.txt", col_names = FALSE)

# Run enrichment analysis for GO terms and pathways (KEGG and REAC), exclude electronic annotation
gostres <- gost(query = list("Soleus muscle HET-WT" = sig_prot_up$X1), 
                organism = "mmusculus", ordered_query = FALSE, 
                multi_query = FALSE, significant = FALSE, exclude_iea = FALSE, 
                measure_underrepresentation = FALSE, evcodes = TRUE, 
                user_threshold = 0.05, correction_method = "g_SCS", 
                domain_scope = "annotated", custom_bg = NULL, 
                numeric_ns = "", sources = c("GO"), as_short_link = FALSE) #, "KEGG", "REAC"

### Manhattan plot of functional enrichment results ###

# Interactive plot
gostplot(gostres, capped = TRUE, interactive = TRUE, pal = c(`GO:MF` = "#3182BDFF", `GO:BP` = "#E6550DFF", `GO:CC` = "#31A354FF", KEGG =
                                                               "#756BB1FF", REAC = "#636363FF", WP = "#0099c6", TF = "#5574a6", MIRNA = "#22aa99", HPA =
                                                               "#6633cc", CORUM = "#66aa00", HP = "#990099"))
                                                             
# Static plot for publishing                                                                    
p <- gostplot(gostres, capped = FALSE, interactive = FALSE, pal = c(`GO:MF` = "#3182BDFF", `GO:BP` = "#E6550DFF", `GO:CC` = "#31A354FF", KEGG =
                                                                     "#756BB1FF", REAC = "#636363FF", WP = "#0099c6", TF = "#5574a6", MIRNA = "#22aa99", HPA =
                                                                     "#6633cc", CORUM = "#66aa00", HP = "#990099"))
                                                                   
# Plot with selected terms and table
pp <- publish_gostplot(p, highlight_terms = c("GO:0005201", "GO:0005518", "GO:0001968", "GO:0050840",
                                              "GO:0031012",
                                              "GO:0001568", "GO:0002376", "GO:0007155", "GO:0007162", "GO:0043062", "GO:0048870", "GO:1901342"
                                              ), 
                       width = NA, height = NA, filename = NULL)
#,"REAC:R-MMU-1474244", "REAC:R-MMU-1474290", "REAC:R-MMU-1474228"

ggsave_fixed("Fig2_GProfiler.pdf", plot_height = 18, plot_width = 15, height = 22, width = 22)

# Table only, can select columns
publish_gosttable(gostres, highlight_terms = c("GO:0005201", "GO:0005518", "GO:0001968", "GO:0050840",
                                               "GO:0031012",
                                               "GO:0001568", "GO:0002376", "GO:0007155", "GO:0007162", "GO:0043062", "GO:0048870", "GO:1901342"
                                               ), 
                  use_colors = FALSE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)
#,"REAC:R-MMU-1474244", "REAC:R-MMU-1474290", "REAC:R-MMU-1474228"
ggsave_fixed("Fig2_GProfiler_table.pdf", plot_height = 20, plot_width = 20, height = 22, width = 22)

