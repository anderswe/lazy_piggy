# Info --------------------------------------------------------------------

# Figure 2 reproducibility code
# Dec 1, 2022
# Anders E.
# Taylor lab


# Libraries -------------------------------------------------------------------------
library(readxl)
library(tidyverse)
library(lubridate)
library(survival)
library(janitor)
library(glue)
library(survminer)
library(magrittr)
library(DESeq2)
library(ggrepel)
library(cowplot)
library(viridis)
library(tidytable)
library(HGNChelper)
library(ComplexHeatmap)

source("src/scripts/configs.R")
source("src/scripts/utils.R")




# Panel A -----------------------------------------------------------------
# OS by genotype

# inventory sheet V3
df <- read_excel(glue("{proj_dir}Sample_Information/Lazy Piggy Project Mouse Inventory V3.xlsx"), sheet = 1, skip = 22) %>%
  clean_names() %>% 
  filter(!is.na(genotype))

#ptc ctrl mice
ptc <- read_excel(glue("{proj_dir}/Survival/SB_Survival.xlsx")) %>% 
  clean_names() %>% 
  select(days, ptc) %>% 
  filter(!is.na(ptc))

#full genotype LP-ptc mice
lp_ptc <- read_tsv(glue("{proj_dir}/Lazy_Piggy_Analysis/SPLINK/LP_sample_metadata.txt")) %>% clean_names() %>% 
  filter(experiment_group != "CTRL")

# LP mice no ptc
lp_no_ptc <- df %>% filter(!is.na(month_mm_12) & #filter out samples without a sacrifice date
                             is.na(ptc) & #filter out ptc+/-
                             grepl("LP-", genotype) & #needs LP transposon
                             nestin_cre == 1 & #needs nestin:cre
                             n_luc_sb100 == 1 & #needs nestin:luc-SB100
                             r26_pber == 1) %>% #needs r26-pb-er
  mutate(birthday = glue("{month_mm_8}-{day_dd_9}-{year_yyyy_10}") %>% as_date(format = '%m-%d-%Y'),
         sacday = glue("{month_mm_12}-{day_dd_13}-{year_yyyy_14}") %>% as_date(format = '%m-%d-%Y'),
         survival_days = difftime(sacday, birthday, units = "days") %>% as.numeric(),
         group = "lp_no_ptc")

# rbind
sdf <- rbind(ptc %>% dplyr::rename(survival_days = days) %>% mutate(group = "ptc") %>% select(-ptc),
             lp_no_ptc %>% select(survival_days, group),
             lp_ptc %>% dplyr::rename(survival_days = total_surival) %>% mutate(group = "lp_ptc") %>% select(survival_days, group)) %>% 
  dplyr::mutate(group = factor(group, levels = c("lp_no_ptc", "ptc", "lp_ptc")))


# create survival objects
surv <- Surv(time = sdf$survival_days, event = rep(1, length(sdf$survival_days)))
fit <- survfit(surv ~ group, data = sdf)

# Prepare KM plot
caption <- ifelse(surv_pvalue(fit, sdf)$pval < 0.0001,
                  "Log-rank\np < .0001",
                  glue("p = {surv_pvalue(fit, sdf)$pval %>% signif(3)}"))
gg <- ggsurvplot(fit, data = sdf, conf.int = F, xlab = "Time (Days)", pval = FALSE,
                 xlim=c(0,350), risk.table = T, #max(df$total_surival)+5
                 legend.title = " ", 
                 # legend = "right",
                 tables.theme = theme_cleantable(),
                 break.time.by = 50,
                 risk.table.fontsize = 4,
                 risk.table.height = 0.20, 
                 surv.median.line = c("h"),
                 legend.labs = c("Nestin:Cre+/-; N-Luc-SB100+/-; R26-PB-ER+/-; LP+/-", "Ptc+/-", "Nestin:Cre+/-; N-Luc-SB100+/-; R26-PB-ER+/-; LP+/-; Ptc+/-"),
                 censor.shape = c("|"), censor.size = 6, 
                 axes.offset = T, tables.y.text = FALSE,
                 palette = as.character(c(pal["blue"], "black", pal["red"]))
)+
  guides(colour = guide_legend(nrow = 3))
gg$plot <- gg$plot+
  annotate(geom = "text", x = 0, y = 0, label = caption, hjust = 0, vjust = 0, size = 4)+
  theme(legend.justification=c(0,0), legend.position = "top") #legend.position = c(0.02, 0.1)
gg$table <- gg$table + 
  theme(plot.title = element_text(hjust = 0, size = 12))
gg

# Export plot
pdf(glue("{repo_dir}outs/survival_lp_by_genotype_20220113.pdf"),height = 4.5, width = 6)
print(gg, newpage = FALSE)
dev.off()


# Panel B -----------------------------------------------------------------
# lazy piggy oncoplot

# import
sig_genes <- readRDS(glue("{awe_local_dir}20220524/sig_genes_tam_gcis_separate_no_recurrence_filter.rds")) # TODO: confirm
metadata <- glue("{awe_local_dir}20220524/metadata_20220524.csv") %>% 
  data.table::fread() %>% 
  select.(-barcode_seq, -tnp_end, -V1) %>% 
  mutate.(donor = (donor %>% str_replace("^7", "chr7") %>% str_replace("^10", "chr10")),
          comment = case_when(
            mouse_id == "3-27-15M" ~ "Biological Control (BlBPcre+; EZH2)",
            mouse_id == "2-16-17T" ~ "LP Triple Control",
            mouse_id == "2-26-16M" ~ "Biological LP Control",
            mouse_id %in% c("11-11-15U", "2-24-17T", "3-13-17T") ~ "Biological Control",
            grepl("4-7-15", mouse_id) ~ "To omit. Unclear genotype + NA tam status. Marked as a CTRL in RNA metadata.",
            TRUE ~ as.character(comment))) %>% 
  unique()

# identify false positives
`%notin%` <- Negate(`%in%`)
false_pos <- c("Sfi1", "En2", "Foxf2", "Pisd-ps3", "Pisd-ps1")

# import and remove false pos
sig_gcis <- sig_genes %>% filter.(!(gene %in% false_pos))

# get dims
genes <- unique(sig_gcis$gene)
row_num <- length(genes)
tamp_samples <- metadata %>% dplyr::filter(tam_status == "TAM+") %>% pull(mouse_id) %>% unique() %>% paste0(.,"_PT")
tamn_samples <- metadata %>% dplyr::filter(tam_status == "TAM-") %>% pull(mouse_id) %>% unique() %>% paste0(.,"_PT")
samples <- sig_gcis$observed_samples %>% paste0(collapse = ";") %>% str_split(";") %>% unlist() %>% unique()
# all_samples <- unique(c(tamp_samples, tamn_samples))
col_num <- length(samples)

# create empty matrix list
mtx <- matrix(0,row_num,col_num)
colnames(mtx) <- samples
rownames(mtx) <- genes
mtx_list <- list(mtx, mtx) %>% set_names(c("tamp", "tamn"))


# assign values to tamp matrix
for(i in genes){ 
  avail_samples <- sig_gcis %>%
    filter(gene == i) %>%
    pull(observed_samples) %>% 
    str_split(";") %>% 
    unlist() %>% 
    unique() %>% 
    .[. %in% tamp_samples]
  mtx_list[["tamp"]][i, avail_samples] <- 1
}
# assign values to tamn matrix
for(i in genes){ 
  avail_samples <- sig_gcis %>%
    filter(gene == i) %>%
    pull(observed_samples) %>% 
    str_split(";") %>% 
    unlist() %>% 
    unique() %>% 
    .[. %in% tamn_samples]
  mtx_list[["tamn"]][i, avail_samples] <- 1
}

# rename Sox2ot to Sox2, 4831426I19Rik to Syne3
mtx_list %<>% purrr::map(\(x){ 
  rownames(x) <- ifelse(rownames(x) == "Sox2ot", "Sox2", rownames(x))
  rownames(x) <- ifelse(rownames(x) == "4831426I19Rik", "Syne3", rownames(x))
  return(x)})


# oncoprint
col = c(tamp = as.character(pal['red']), tamn = as.character(pal['blue']))
op <- ComplexHeatmap::oncoPrint(mtx_list,
                alter_fun = list(
                  background = alter_graphic("rect", fill = "#CCCCCC"),
                  tamp = alter_graphic("rect", fill = col["tamp"]),
                  tamn = alter_graphic("rect", fill = col["tamn"])),
                col = col,
                show_pct = TRUE,
                pct_side = "right",
                row_names_side = "left",
                right_annotation = NULL,
                row_names_gp = gpar(fontsize = 10, col = dplyr::case_when(
                  rownames(mtx_list[[1]]) %in% c("Kcnb1", "Kcnh2") ~ pal[3], #Kv hits
                  rownames(mtx_list[[1]]) %in% c("Sox2", "Smarca4", "Dync1h1") ~ pal[5], #known SHH alterations / stem cell compartment marker cf. Morrissy, Vanner
                  TRUE ~ "black")),
                pct_gp = gpar(fontsize = 10),
                show_heatmap_legend = FALSE)
                
pdf(glue("{repo_dir}outs/gcis_oncoprint.pdf"), h = 10, w = 7)
draw(op)
dev.off()


# Panel C&E -----------------------------------------------------------------
# done in Cytoscape + EnrichmentMap

# Panel D -----------------------------------------------------------------
# RNAseq volcano

# read in
res_filt <- readRDS(glue("{repo_dir}/outs/rna_deseq_res_filt.rds"))

# cleaning volcano df
dict <- gprofiler2::gconvert(rownames(res_filt), "mmusculus", "MGI")
p_adj_cutoff <- 0.05
lfc_cutoff <- 1
lp_volc <- res_filt %>%  
  mutate(symbol = ifelse(rownames(.) %in% dict$input, plyr::mapvalues(rownames(.),dict$input,dict$name), NA)) %>% 
  filter(!is.na(symbol) & !is.na(padj)) %>% 
  mutate(colour = ifelse(padj > p_adj_cutoff | abs(log2FoldChange) < lfc_cutoff, "not_sig",
                         ifelse(log2FoldChange > 0, "up", "down")),
         label = ifelse(symbol %in% (filter(., colour != "not_sig") %>% slice_min(order_by = padj, n = 10) %>% pull(symbol)) | (grepl("Kcn", symbol) & colour == "up") | colour == "down",symbol,""))

# Set plot boundaries
xbound <- max(abs(lp_volc$log2FoldChange)) + 0.1

# palette
mypal <- c(as.character(pal['red']), "grey", as.character(pal['blue']))

# manually remove Scd2 for visualization purposes
lp_volc %<>% filter(symbol != "Scd2")

# ggplot
pdf(glue("{repo_dir}/outs/volcano_deseq2_lazy_piggy_rnaseq_tam_pos_vs_neg.pdf"), h = 5, w = 5)
print(ggplot(lp_volc, aes(log2FoldChange, -log10(padj), label = label, colour = colour))+
        geom_point()+
        geom_label_repel(size = ifelse(grepl("Kcn", lp_volc$symbol), 3.5 ,3.5),
                        segment.size = 0.2, box.padding = 0.5,
                        color = ifelse(grepl("Kcn", lp_volc$symbol), pal[3],"black"),
                        segment.color = "black", max.overlaps = Inf)+ 
        scale_color_manual(values = alpha(c(mypal), 0.9),
                           breaks = c("up", "not_sig", "down"),
                           labels = c("TAM+", "Not Sig", "TAM-"))+
        theme(panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              plot.title = element_text(size=10),
              legend.key = element_blank(),
              legend.title = element_blank(),
              legend.position = "none", #c(.95, .95)
              legend.justification = c("right", "top"),
              legend.margin = margin(rep(1,4)))+
        xlab(expression('log'[2]*'FoldChange'))+
        ylab(expression('-log'[10]*'P'[adj]))+
        xlim(-xbound, xbound)+
        guides(colour = guide_legend(override.aes = list(size=3)))
)
dev.off()


# Panel F -----------------------------------------------------------------
# MB vs normal brain volcanos

# import
dds <- readRDS("src/objs/dds_all_subgroup.RDS") #MAGIC cohort DESeq2 results from Hendrikse et al. 2022 Nature

# K channel list
# from: https://www.genenames.org/data/genegroup/#!/group/183
k_channels <- read_tsv("src/objs/group-183.tsv") %>% 
  clean_names() %>% 
  pull(approved_symbol)

# loop by subgroup
for(subgroup in c("Group4", "SHH", "Group3", "WNT")){
  
  
  # lfcShrink
  resLFC <- lfcShrink(dds, coef=glue("RNA_seq_subgroup_allMB_k4_{subgroup}_vs_normal"), type="apeglm")
  
  # plot setup
  p_adj_cutoff <- 1e-3
  mypal <- c(as.character(pal['red']), "grey", as.character(pal['blue']))
  dds_k <- resLFC %>% as.data.frame() %>% 
    mutate(symbol = str_split_fixed(rownames(.), "___", 2) %>% as.data.frame() %>% pull(V1)) %>% 
    filter(symbol %in% k_channels) %>% 
    mutate(label = ifelse(symbol %in% (slice_min(., order_by = padj, n = 35) %>% pull(symbol)),symbol,""),
           colour = ifelse(padj > p_adj_cutoff | abs(log2FoldChange) < 1, "not_sig",
                           ifelse(log2FoldChange > 0, "up", "down")))
  
  # Set plot boundaries
  xbound <- max(abs(dds_k$log2FoldChange)) + 0.25
  
  # ggplot
  outdir <- glue("{repo_dir}/outs")
  pdf(glue("{outdir}/volcano_deseq2_{subgroup}_vs_normal_brain_k_channels.pdf"), 
      height = 5, width = 5)
  print(ggplot(dds_k, aes(log2FoldChange, -log10(padj), label = label, colour = colour))+
          geom_point()+
          geom_text_repel(size = 3.5, segment.size = 0.2, box.padding = 0.3, color = ifelse(dds_k$symbol == "KCNB2", "red" ,"black"), segment.color = "black")+
          scale_color_manual(values = alpha(c(mypal), 0.9),
                             breaks = c("up", "not_sig", "down"),
                             labels = c(glue("{subgroup} MB"), "Not Sig", "Ctrl brain"))+
          theme(panel.background = element_blank(),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                plot.title = element_text(size=12, hjust = 0.5),
                legend.position = "none")+
          xlab(expression('log'[2]*'(FoldChange)'))+
          ylab(expression('-log'[10]*'(P-adjusted)'))+
          labs(title = glue("{str_replace(subgroup,'Group', 'Group ')} MB vs Normal Cerebellum"))  +
          xlim(-xbound, xbound)+
          guides(colour = guide_legend(override.aes = list(size=3)))
  )
  dev.off()
  
  # save DE results
  write_xlsx(dds_k, glue("{outdir}/de_results_{subgroup}.xlsx"))
}



# Panel G -----------------------------------------------------------------
# done in bioRender



# Bonus: kcnb1 / kcnh2 RNAseq boxplots by insertion status -------------------------------------------------------------------

# outdir
outdir <- glue::glue("{repo_dir}/outs")
rna_dir <- glue("{proj_dir}/Lazy_Piggy_Analysis/RNAseq/Expression")

# insertions
input_file <- glue::glue("{awe_local_dir}20211201/all_annotations_20211201.annot.gz")
insertions <- data.table::fread(input_file) %>% 
  tidyr::separate(sample, c("patient", "tissue", "library"), sep = "_")

# rna count matrix
LP_count_matrix <- readRDS(glue::glue("{outdir}/lp_counts_clean.rds")) %>% as.matrix() %>% DGEobj.utils::convertCounts(unit = "CPM", log = TRUE)
LP_metadata <- read.table(file = glue("{rna_dir}/LP_sample_metadata.txt"), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t")


# loop
genes <- c("Kcnb1", "Kcnh2")
mypal <- c(as.character(pal['red']), as.character(pal['blue']))
for(gene in genes){
  
  # insertion samples
  insertion_samples <- insertions %>% dplyr::filter(gene == !!gene) %>% dplyr::pull(patient) %>% unique()
  
  # count number of samples with insertions at gene ... with >1 supporting read
  insertions %>% dplyr::filter(gene == !!gene & count > 1) %>% dplyr::pull(patient) %>% unique() %>% length() %>% print()
  
  # and of those...how many have support from a PBL/PBR library? ... and are not from a CTRL?
  insertions %>% dplyr::filter(gene == !!gene & count > 1 & grepl("PB", library) & tissue == "PT") %>% dplyr::pull(patient) %>% unique() %>% length() %>% print()
  
  # and how many non control sample are there with PBL/PBR libraries?
  insertions %>% dplyr::filter(tissue == "PT" & grepl("PB", library)) %>% dplyr::pull(patient) %>% unique() %>% length()
  
  # gene key
  gene_key <- gene %>% gprofiler2::gconvert(organism = "mmusculus", target = "ENSG") %>% dplyr::pull(target)
  
  # subset matrix
  df <- LP_count_matrix[gene_key,] %>% 
    tibble::enframe() %>% 
    dplyr::rename(sample = name, counts = value) %>% 
    dplyr::left_join(LP_metadata %>% dplyr::select(sample = Mouse, tam = TamStatus), by = "sample") %>% 
    dplyr::mutate(insertion_status = ifelse(sample %in% insertion_samples, "Insertion", "NoInsertion"),
                  insertion_tam = glue::glue("{insertion_status}_{tam}"))
  
  pdf(glue::glue("{outdir}/boxplot_expression_{gene}_no_pvalues.pdf"), w = 4, h = 3.5)
  print(ggplot(df, aes(insertion_status, counts, color = insertion_status))+
          geom_violin(width=0.6)+
          geom_boxplot(width=0.1, color="grey", alpha=0.5, outlier.shape = NA)+
          geom_jitter(width=0.1, alpha = 0.5)+
          # stat_compare_means()+
          theme_classic()+
          ylab('Expression (log CPM)')+
          theme(legend.position = "none", axis.title.x = element_blank())+
          scale_color_manual(values = mypal))
  dev.off()
  
}




# bonus: volcano gcis shh vs normal cb ------------------------------------

# import
dds <- readRDS("src/objs/dds_all_subgroup.RDS") #MAGIC cohort DESeq2 results from Hendrikse et al. 2022 Nature

# gcis genes
false_pos <- c("Sfi1", "En2", "Foxf2", "Pisd-ps3", "Pisd-ps1")
sig_genes <- readRDS(glue("{awe_local_dir}20220524/sig_genes_tam_gcis_separate_no_recurrence_filter.rds")) %>% 
  dplyr::filter(gene %notin% false_pos) %>% 
  dplyr::pull(gene) %>% 
  gprofiler2::gorth(source_organism = "mmusculus", target_organism = "hsapiens") %>% 
  dplyr::pull(ortholog_name)


# loop by subgroup
for(subgroup in c("Group4", "SHH", "Group3", "WNT")){
  
  
  # lfcShrink
  resLFC <- lfcShrink(dds, coef=glue("RNA_seq_subgroup_allMB_k4_{subgroup}_vs_normal"), type="apeglm")
  
  # plot setup
  p_adj_cutoff <- 1e-3
  mypal <- c(as.character(pal['red']), "grey", as.character(pal['blue']))
  dds_g <- resLFC %>% as.data.frame() %>% 
    mutate(symbol = str_split_fixed(rownames(.), "___", 2) %>% as.data.frame() %>% pull(V1)) %>% 
    filter(symbol %in% sig_genes) %>% 
    mutate(label = ifelse(symbol %in% (slice_min(., order_by = padj, n = 25) %>% pull(symbol) %>% c(., "SMARCA4", "KCNH2", "DYNC1H1", "KCNB1")) ,symbol,""),
           colour = ifelse(padj > p_adj_cutoff | abs(log2FoldChange) < 1, "not_sig",
                           ifelse(log2FoldChange > 0, "up", "down")))
  
  # Set plot boundaries
  xbound <- max(abs(dds_g$log2FoldChange)) + 0.25
  
  # ggplot
  outdir <- glue("{repo_dir}/outs")
  pdf(glue("{outdir}/volcano_deseq2_{subgroup}_vs_normal_brain_gcis_genes.pdf"), 
      height = 5, width = 5)
  print(ggplot(dds_g, aes(log2FoldChange, -log10(padj), label = label, colour = colour))+
          geom_point()+
          geom_text_repel(size = 2.75, segment.size = 0.2, box.padding = 0.35, color = ifelse(dds_g$symbol == "KCNB2", "red" ,"black"), segment.color = "black")+
          scale_color_manual(values = alpha(c(mypal), 0.9),
                             breaks = c("up", "not_sig", "down"),
                             labels = c(glue("{subgroup} MB"), "Not Sig", "Ctrl brain"))+
          theme(panel.background = element_blank(),
                panel.border = element_rect(colour = "black", fill=NA, size=1),
                plot.title = element_text(size=12, hjust = 0.5),
                legend.position = "none")+
          xlab(expression('log'[2]*'(FoldChange)'))+
          ylab(expression('-log'[10]*'(P-adjusted)'))+
          labs(title = glue("{str_replace(subgroup,'Group', 'Group ')} MB vs Normal Cerebellum"))  +
          xlim(-xbound, xbound)+
          guides(colour = guide_legend(override.aes = list(size=3)))
  )
  dev.off()
  
  # save DE results
  write_xlsx(dds_g, glue("{outdir}/de_results_gcis_genes_{subgroup}.xlsx"))
}





# session info ------------------------------------------------------------
sessionInfo()


