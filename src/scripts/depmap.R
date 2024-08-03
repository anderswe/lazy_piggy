# Info --------------------------------------------------------------------
#
# Anders E.
# Taylor lab
# May 7, 2024
# 
# goal here is to plot depmap dependency scores of LP screen hits
# 


# Libraries ---------------------------------------------------------------
library(magrittr)
library(ggplot2)
source("src/scripts/utils.R")
source("src/scripts/configs.R")



# Dirs --------------------------------------------------------------------

outdir <- glue::glue("{repo_dir}/outs/")


# Import ------------------------------------------------------------------


# lp screen hits
genes <- readRDS(glue::glue("{awe_local_dir}20220524/sig_genes_tam_gcis_separate_no_recurrence_filter.rds")) %>% 
  dplyr::filter(gene %notin% c("Sfi1", "En2", "Foxf2", "Pisd-ps3", "Pisd-ps1")) %>% 
  dplyr::pull(gene) %>% 
  unique() %>%
  gprofiler2::gconvert() %>% 
  dplyr::pull(name) %>% 
  grep("MIR", ., value = T, invert = T) %>% 
  c("KCNB2")

# depmap dependency scores
deps <- data.table::fread("src/data/depmap_23q4/CRISPRGeneEffect.csv") %>% # downloaded 2024-05-07 from: https://depmap.org/portal/download/all/
  dplyr::rename(cell_line = V1)

# mb line achilles accession codes
mb_lines <- data.table::fread("src/data/depmap_23q4/mb.txt") %>% names() # corresponding to: "D283MED;D341Med;D384;D425;D458;D556;DAOY;ONS76;SUMB002;UW228"



# clean -------------------------------------------------------------------

new_names <- colnames(deps) %>% stringr::str_split(" ") %>% lapply(\(x) x[[1]]) %>% unlist()

df <- deps %>% 
  magrittr::set_colnames(new_names) %>% 
  dplyr::select(cell_line, !!genes) %>% 
  tidyr::pivot_longer(-cell_line, names_to = "gene", values_to = "score") %>% 
  dplyr::mutate(mb_status = ifelse(cell_line %in% mb_lines, "mb", "other")) # seems to be missing 3 of the MB lines


pdf(glue::glue("{outdir}/depmap.pdf"), w = 6, h = 4)
ggplot(df, aes(gene, score, fill = mb_status, color = gene, shape = mb_status))+
  geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.6)+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45))+
  scale_shape_manual(values = c(2, 16), name = "", labels = c("MB", "Other"))+
  scale_color_manual(values = c(tableau20, tableau20))+
  labs(x = "", y = "DepMap Gene Effect Score")+
  theme(axis.text.x = element_text(hjust = 1, vjust = 1), legend.position = c(0.1,0.8))+
  guides(color = "none", fill = "none")+
  scale_y_reverse()+
  geom_hline(yintercept = -1)
dev.off()


pdf(glue::glue("{outdir}/depmap_boxplots.pdf"), w = 6, h = 4)
ggplot(df, aes(gene, score, fill = mb_status, color = gene, shape = mb_status))+
  # geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.6)+
  geom_boxplot()+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45))+
  scale_shape_manual(values = c(2, 16), name = "", labels = c("MB", "Other"))+
  scale_color_manual(values = c(tableau20, tableau20))+
  scale_fill_manual(values = c(tableau20, tableau20))+
  labs(x = "", y = "DepMap Gene Effect Score")+
  theme(axis.text.x = element_text(hjust = 1, vjust = 1), legend.position = c(0.1,0.8))+
  guides(color = "none", fill = "none")+
  scale_y_reverse()
dev.off()

# agnostic mb vs other comparison -----------------------------------------

diffs <- deps %>% 
  magrittr::set_colnames(new_names) %>% 
  tidyr::pivot_longer(-cell_line, names_to = "gene", values_to = "score") %>% 
  dplyr::mutate(mb_status = ifelse(cell_line %in% mb_lines, "mb", "other")) %>% 
  dplyr::group_by(mb_status, gene) %>% 
  dplyr::summarise(med = median(score)) %>% 
  tidyr::pivot_wider(names_from = "mb_status", values_from = "med") %>% 
  dplyr::mutate(diff = other - mb) %>% 
  dplyr::arrange(desc(diff))


# session info ------------------------------------------------------------
sessionInfo()



