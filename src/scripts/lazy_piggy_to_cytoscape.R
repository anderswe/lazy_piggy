# Info --------------------------------------------------------------------
# lazy_piggy_to_cytoscape.R
# Anders E.
# Oct 5th, 2021


# Notes -------------------------------------------------------------------
# This script runs gProfiler on gCIS results from Lazy Piggy as:
# 1. combined hits (ranked gene lists)
# 2. separate (ranked gene lists)
# and preps both for visualization in Cytoscape

# Reimand / Bader lab paper at: https://www.nature.com/articles/s41596-018-0103-9

# GMT file for mouse GO:BP + GO:MF + REAC downloaded from gProfiler website on Oct 4, 2021
# at: /Volumes/ANDERS/data/cytoscape/inputs/gmts/20211004/gprofiler_full_mmusculus.name.gmt

# Libraries ---------------------------------------------------------------
library(HGNChelper)
library(magrittr)
library(janitor)
library(nichenetr)
library(gprofiler2)
library(tidyverse)
library(glue)
library(ggrepel)
library(clipr)
library(data.table)
`%notin%` <- Negate(`%in%`)


# Inputs ------------------------------------------------------------------
lookup <- readRDS("/Volumes/ANDERS/data/r_objs/mouse_human_lookup_table.rds") # just a conversion table from biomaRt

setwd("/Volumes/ANDERS/data/cytoscape/inputs")
gem_path <- glue("{getwd()}/gems")
plot_path <- "/Volumes/ANDERS/data/plots/gprofiler2" #path to output manhattan plots for prelim gprofiler results
# combined_gcis_path <- "/Volumes/ANDERS/data/Lazy Piggy/gCIS results/TAM+-/results_filtered_subclonal_high/gCIS.temp"

# Filtered
# These files are just filtered .tsv versions of the final output from src/scripts/gcis.R: sig_genes_tam_gcis_separate_no_recurrence_filter.rds
tam_pos_gcis_path <- "/Volumes/ANDERS/data/Lazy Piggy/20211208/Exp_group/TAM+/all_gcis_TAM+.tsv"  #"/Volumes/ANDERS/data/Lazy Piggy/gCIS results/PB_TAM+/results_filtered_subclonal_high/gCIS.temp"
tam_neg_gcis_path <- "/Volumes/ANDERS/data/Lazy Piggy/20211208/Exp_group/TAM-/all_gcis_TAM-.tsv"  #"/Volumes/ANDERS/data/Lazy Piggy/gCIS results/PB_TAM-/results_filtered_subclonal_high/gCIS.temp"

# Unfiltered
# tam_pos_gcis_path <- "/Volumes/ANDERS/data/Lazy Piggy/gCIS results/PB_TAM+/results_unfiltered_subclonal_high/gCIS.temp"
# tam_neg_gcis_path <- "/Volumes/ANDERS/data/Lazy Piggy/gCIS results/PB_TAM-/results_unfiltered_subclonal_high/gCIS.temp"
exclude <- c("En2", "Foxf2", "Sfi1") # false positives from gCIS
min_term_size <- 10 #pathway term size thresholds suggested from the paper linked in the Notes above
max_term_size <- 500
min_num_tumors <- 2 #threshold for including genes in g:Profiler run
bonferroni_cutoff <- 0.01 #think about using 0.1?? try diff cutoffs


# Functions ---------------------------------------------------------------
source("/Volumes/ANDERS/scripts/R/lazy_piggy/pre_cytoscape_functions.R")



# TAM+ vs TAM- ------------------------------------------------------------


# Read in
tam_pos <- fread(tam_pos_gcis_path, nThread = 8) %>% 
  clean_names() %>% 
  dplyr::rename(p_value = pvalue) %>% 
  as_tibble() %>% 
  dplyr::mutate(hum_gn = ifelse(gene %in% lookup$mgi,
                                plyr::mapvalues(gene, lookup$mgi, lookup$hgnc_new, warn_missing = F),
                                NA))

tam_neg <- fread(tam_neg_gcis_path, nThread = 8) %>% 
  clean_names() %>% 
  dplyr::rename(p_value = pvalue) %>% 
  as_tibble() %>% 
  mutate(hum_gn = ifelse(gene %in% lookup$mgi,
                      plyr::mapvalues(gene, lookup$mgi, lookup$hgnc_new, warn_missing = F),
                      NA))

# Filter genes
# NOTE TO SELF, MAYBE NEED TO CALCULATE BONFERRONI AFTER!!!! FILTERING BY OBSERVED SAMPLE COUNT!!!!!
tam_pos_genes <- tam_pos %>%
  filter(observed_sample_count >= min_num_tumors & gene %notin% exclude) %>%
  mutate(bonf = p.adjust(p_value, "bonferroni", length(p_value))) %>% 
  filter(bonf < bonferroni_cutoff) %>% #toggle
  arrange(desc(observed_sample_count), p_value) %>% #toggle
  # arrange(p_value) %>% #toggle
  pull(hum_gn) %>% na.omit() %>% as.character()
  
tam_neg_genes <- tam_neg %>%
  filter(observed_sample_count >= min_num_tumors & gene %notin% exclude) %>%
  mutate(bonf = p.adjust(p_value, "bonferroni", length(p_value))) %>% 
  filter(bonf < bonferroni_cutoff) %>% #toggle
  arrange(desc(observed_sample_count), p_value) %>% #toggle
  # arrange(p_value) %>% #toggle
  pull(hum_gn) %>% na.omit() %>% as.character()
  


# Run gprofiler2 
gres2 <- gost(list("TAM+" = tam_pos_genes, "TAM-" = tam_neg_genes),
             organism = "hsapiens",
             multi_query = F, #multi_query = F is counter-intuitive. needs to be F for easy output to enrichmentmap/cytoscape. Doesn't actually change analysis in gprofiler.
             ordered_query = T,
             significant = T,
             exclude_iea = T,
             sources = c("GO:BP","GO:MF", "GO:CC","REAC", "KEGG", "CORUM", "WP"),
             evcodes = T)

# Filter gprofiler2 results
res2 <- gres2$result %>% filter(term_size >= min_term_size & term_size <= max_term_size)
gres2$result <- res2

# temp
print(gostplot(gres2, capped = T, interactive = T))

# Plot
pdf(glue("{plot_path}/lazy_piggy_tam_pos_and_neg.pdf"), height = 10, width = 10)
print(gostplot(gres2, capped = T, interactive = F)+
        geom_text_repel(label = gres2$result$term_name))
dev.off()


# Format gprofiler results as gems
gpro_gems <- map(res2$query %>% unique(), \(x){
  res2 %>% filter(query == x) %>% gprofiler_to_gem()
})
names(gpro_gems) <- res2$query %>% unique()

# REACTOME digression
# copy paste these into reactome here:
# https://reactome.org/PathwayBrowser/#/TOOL=AT
# then export as gem
write_clip(tam_pos_genes)
write_clip(tam_neg_genes)
reac_gem_pos <- reactome_to_gem("/Volumes/ANDERS/data/cytoscape/inputs/reactome/mm_no_bonf_1_tumor_min/result_tam_pos.csv") 
reac_gem_neg <- reactome_to_gem("/Volumes/ANDERS/data/cytoscape/inputs/reactome/mm_no_bonf_1_tumor_min/result_tam_neg.csv") 

# merge gems
gem_pos <- rbind(gpro_gems$`TAM+`, reac_gem_pos)
gem_neg <- rbind(gpro_gems$`TAM-`, reac_gem_neg)

# Export gems
export_as_gem(gem_pos, "pos")
export_as_gem(gem_neg, "neg")

# Export .rnks
export_rnk(tam_pos, "TAM+")
export_rnk(tam_neg, "TAM-")

# Export gmt
# reac_gmt <- read.GMT("/Volumes/ANDERS/data/cytoscape/inputs/gmts/20211012/hs_reac_only/gprofiler_full_hsapiens.name.gmt")





# # TAM+ vs TAM- CONVERTED TO HUMAN --------------------------------------------------------------------
# 
# # Read in
# tam_pos <- read_gcis(tam_pos_gcis_path)
# tam_neg <- read_gcis(tam_neg_gcis_path)
# 
# # Export .rnks
# export_rnk(tam_pos, "TAM+")
# export_rnk(tam_neg, "TAM-")
# 
# # Filter genes and pull human
# tam_pos_genes <- tam_pos %>%
#   filter(num_tumors >= min_num_tumors & gene %notin% exclude) %>%
#   filter(!is.na(hum_gn)) %>% #new
#   arrange(p_value) %>%
#   pull(hum_gn)
# tam_neg_genes <- tam_neg %>%
#   filter(num_tumors >= min_num_tumors & gene %notin% exclude) %>%
#   filter(!is.na(hum_gn)) %>% #new
#   arrange(p_value) %>%
#   pull(hum_gn)
# 
# # Run gprofiler2 
# gres2 <- gost(list("TAM+" = tam_pos_genes, "TAM-" = tam_neg_genes),
#               organism = "hsapiens",
#               multi_query = F, #multi_query = F is counter-intuitive. needs to be F for easy output to enrichmentmap/cytoscape. Doesn't actually change analysis in gprofiler.
#               ordered_query = T,
#               significant = T,
#               exclude_iea = T,
#               sources = c("GO:BP","GO:MF", "GO:CC", "KEGG", "CORUM", "WP"),
#               evcodes = T)
# 
# # Filter gprofiler2 results
# res2 <- gres2$result %>% filter(term_size >= min_term_size & term_size <= max_term_size)
# gres2$result <- res2
# 
# # Plot
# pdf(glue("{plot_path}/lazy_piggy_tam_pos_and_neg.pdf"), height = 10, width = 10)
# print(gostplot(gres2, capped = T, interactive = F)+
#         geom_text_repel(label = gres2$result$term_name))
# dev.off()
# 
# 
# # Send to cytoscape
# map(res2$query %>% unique(), \(x){
#   res2 %>% filter(query == x) %>% export_as_gem()
# })



# # Combined TAM+/- genes -----------------------------------------------------------
# 
# # Read in
# gcis <- read_gcis(combined_gcis_path)
# 
# # Filter genes
# genes <- gcis %>%
#   filter(num_tumors >= min_num_tumors & gene %notin% exclude) %>%
#   arrange(p_value) %>%
#   pull(gene)
# 
# # Export .rnk
# export_rnk(gcis, "TAM+-")
# 
# # Run gprofiler2
# gres <- gost(list("TAM+-" = genes),
#              organism = "mmusculus",
#              ordered_query = T,
#              significant = T,
#              exclude_iea = T,
#              sources = c("GO:BP","GO:MF", "GO:CC", "REAC", "KEGG", "CORUM", "WP"),
#              evcodes = T)
# 
# # Filter gprofiler2 results
# res <- gres$result %>% filter(term_size >= min_term_size & term_size <= max_term_size)
# gres$result <- res
# 
# # Plot
# pdf(glue("{plot_path}/lazy_piggy_combined_tam.pdf"), height = 5, width = 10)
# print(gostplot(gres, capped = T, interactive = F)+
#         geom_text_repel(label = gres$result$term_name))
# dev.off()
# 
# # Send to cytoscape
# export_as_gem(res)

