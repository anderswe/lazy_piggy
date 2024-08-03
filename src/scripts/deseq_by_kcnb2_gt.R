# Info --------------------------------------------------------------------
#
# Anders E.
# Taylor lab
# Nov 10, 2023
# 
# goal here is to run DESeq2 on Math1-Cre; SmoM2; Kcnb2 +/+ vs Kcnb2 null samples
# 


# Libraries ---------------------------------------------------------------
library(DESeq2)
library(gprofiler2)
library(ggrepel)
library(magrittr)
library(cowplot)
library(viridis)
library(tidyverse)
library(HGNChelper)
library(glue)
library(janitor)
library(ggsci)
library(writexl)
source("src/scripts/utils.R")
source("src/scripts/configs.R")

# palette
mypal <- pal[c("red", "grey", "blue")] %>% magrittr::set_names(c("Knockout", "not_sig", "Control"))

# Inputs ------------------------------------------------------------------

outdir <- glue::glue("{repo_dir}/outs/")
counts <- readRDS(glue("{repo_dir}/src/data/star_counts/counts.rds"))
md <- glue::glue("{repo_dir}/src/data/fastqs/KCNB2_merge/metadata.csv") %>% 
  data.table::fread() %>% 
  dplyr::mutate(sample = sample %>% stringr::str_replace("Mspos", "MS_pos")) %>% 
  dplyr::filter(condition != "Piezo2_ko") %>% 
  dplyr::mutate(gt = ifelse(grepl("ko", condition), "ko", "wt") %>% factor(levels = c("wt", "ko"))) %>% 
  dplyr::select(sample, gt) %>% 
  dplyr::arrange(factor(sample, levels = colnames(counts))) %>% 
  tibble::column_to_rownames("sample")
    



# Cleaning ------------------------------------------------------------------


# keep only genes with non-zero expression in at least 2 samples per group
ko_samples <- md %>% dplyr::filter(gt == "ko") %>% rownames()
wt_samples <- md %>% dplyr::filter(gt == "wt") %>% rownames()
na_mtx <- counts #duplicate count matrix
na_mtx[na_mtx == 0] <- NA #convert 0s to NAs in the duplicate
na_counts <- data.frame(gene = rownames(counts), ko = "", wt = "") #create an empty df to populate with counts of NAs
na_counts$ko <- apply(na_mtx,1,function(x) sum(is.na(x[ko_samples]))) #count em
na_counts$wt <- apply(na_mtx,1,function(x) sum(is.na(x[wt_samples]))) # ""
na_counts$ko_non_zero <- length(ko_samples) - na_counts$ko #get number of non-zero samples per group
na_counts$wt_non_zero <- length(wt_samples) - na_counts$wt # ""
genes_keep <- na_counts %>% filter(ko_non_zero >= 2 & wt_non_zero >= 2) %>% pull(gene) #get gene names
counts <- counts[genes_keep, ]


# Run DESeq2 --------------------------------------------------------------

# deseq2
dds <- DESeqDataSetFromMatrix(countData = counts, colData = md, design = ~ gt)
dds <- DESeq(dds)
res <- lfcShrink(dds, coef = "gt_ko_vs_wt", type = "apeglm")

# get gene symbol dict
dict <- res %>% 
  rownames() %>% 
  gprofiler2::gconvert(organism = "mmusculus", target = "MGI") %>% 
  dplyr::select(input, target)

# convert to symbols and filter out NA
res_filt <- res %>% 
  as.data.frame() %>% 
  dplyr::mutate(gene = plyr::mapvalues(rownames(.), dict$input, dict$target, warn_missing = FALSE)) %>% 
  dplyr::filter(!is.na(padj) & !grepl("ENSMUS", gene))


# save for plotting
saveRDS(res_filt, glue("{repo_dir}/outs/rna_deseq_res_filt.rds"))
data.table::fwrite(res_filt, glue::glue("{outdir}/deseq_kcnb2_ko_vs_wt.csv"))


# TPM calculation ---------------------------------------------------------

# Imports
ref <- readr::read_tsv("https://data.broadinstitute.org/Trinity/CTAT/cnv/mouse_gencode.GRCm38.p6.vM25.basic.annotation.by_gene_name.infercnv_positions",
                       col_names = c("symbol", "chr", "start", "end")) %>% 
  dplyr::mutate(length = abs(end - start),
                kb_length = length / 1000,
                ensmus = ifelse(symbol %in% dict$target,
                                plyr::mapvalues(symbol, dict$target, dict$input, warn_missing = FALSE),
                                NA))

# Intersect and order
intsx <- intersect(rownames(counts), ref$ensmus)
ref_filt <- ref %>% dplyr::filter(ensmus %in% intsx)
counts_filt <- counts %>% 
  dplyr::filter(rownames(.) %in% intsx) %>% 
  dplyr::arrange(factor(rownames(.), levels = ref_filt$ensmus))


# convert to TPM for visualization:
# 1. Divide the read counts by the length of each exon in kilobases. This gives you reads per kilobase (RPK).
# 2. Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
# 3. Divide the RPK values by the “per million” scaling factor. This gives you TPM.
rpk <- (counts_filt / ref_filt$kb_length) %>% as.matrix()
per_million <- colSums(rpk) / 1000000
tpm <- rpk %*% diag(per_million^-1)
colnames(tpm) <- colnames(counts_filt)

ggvln <- tpm %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("ensmus") %>% 
  dplyr::mutate(gene = plyr::mapvalues(ensmus, dict$input, dict$target, warn_missing = FALSE)) %>% 
  dplyr::select(-ensmus) %>% 
  tidyr::pivot_longer(-gene, names_to = "sample", values_to = "tpm") %>% 
  dplyr::mutate(gt = plyr::mapvalues(sample, rownames(md), as.character(md$gt), warn_missing = FALSE))

# violin plot
gene <- "Kcnb2"
pdf(glue::glue("{outdir}/{gene}_tpm_violinplot.pdf"), w = 4, h  = 4)
ggplot(ggvln %>% dplyr::filter(gene == !!gene), aes(gt, tpm, colour = gt))+
  geom_violin(width=0.6)+
  geom_boxplot(width=0.1, color="grey", alpha=0.2)+
  geom_jitter(width=0.1, alpha = 0.8)+
  ggpubr::stat_compare_means()+
  theme_classic()+
  ylab(glue::glue("{gene} expression (TPM)"))+
  theme(legend.position = "none",axis.title.x = element_blank())+
  scale_color_manual(labels = c("Kcnb2 KO", "Kcnb2 WT"), values = as.character(pal[c("red", "blue")]))+
  scale_x_discrete(labels = c("Kcnb2 KO", "Kcnb2 WT"))
dev.off()


# Volcano plots -----------------------------------------------------------


lfc_thresh <- 1
padj_thresh <- 0.05
ggdf <- res_filt %>%
  mutate(colour = ifelse(padj > padj_thresh | abs(log2FoldChange) < lfc_thresh, "not_sig",
                         ifelse(log2FoldChange > 0, "Knockout", "Control")),
         symbol = make.unique(gene),
         label = ifelse(colour == "not_sig", NA, symbol)) %>% 
  dplyr::filter(!grepl("^Gm|Rik$", symbol))# filter out unannotated transcripts for visualization

xbound <- max(abs(ggdf$log2FoldChange)) + 0.25




# ggplot
highlights <- c("Kcnb2", "Kcnc4", "Kcnt1", "Kcnab3")
pdf(glue("{outdir}/volcano_kcnb2_ko_vs_wt.pdf"), h = 10, w = 11)
print(ggplot(ggdf, aes(log2FoldChange, -log10(padj), label = label, colour = colour))+
        geom_point()+
        geom_text_repel(size = ifelse(ggdf$label %in% highlights, 6, 5),
                        segment.size = 0.2, box.padding = 0.5,
                        color = ifelse(ggdf$label %in% highlights, "red", "black"),
                        segment.color = "black", max.overlaps = Inf)+ 
        scale_color_manual(values = alpha(c(mypal), 0.9),
                           breaks = c("Knockout", "Control"))+
        theme(panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              plot.title = element_text(size=10),
              legend.key = element_blank(),
              legend.title = element_blank(),
              # legend.position = c(.1, .95),
              legend.justification = c("left", "top"),
              legend.margin = margin(rep(1,4)))+
        xlab(expression('log'[2]*'(FoldChange)'))+
        ylab(expression('-log'[10]*'P'['adj']))+
        xlim(-xbound, xbound)+
        guides(colour = guide_legend(override.aes = list(size=3)))
)
dev.off()



# rtk subset volcano ------------------------------------------------------

rtks <- purrr::map(c("321", "1095", "1096"), \(x){
  glue::glue("{repo_dir}/src/objs/group-{x}.csv") %>% # downloaded from e.g. https://www.genenames.org/data/genegroup/#!/group/183
    data.table::fread(skip = 1) %>% 
    return()}) %>% 
  data.table::rbindlist() %>% 
  dplyr::pull(`Approved symbol`) %>% 
  gprofiler2::gorth() %>% 
  dplyr::pull(ortholog_name)


lfc_thresh_rtks <- 0.1
p_thresh_rtks <- 0.1
pdf(glue::glue("{outdir}/volcano_rtks_kcnb2_ko_vs_wt.pdf"), w = 6, h = 5)
res_filt %>%
  mutate(colour = ifelse(pvalue > p_thresh_rtks | abs(log2FoldChange) < lfc_thresh_rtks, "not_sig",
                         ifelse(log2FoldChange > 0, "Knockout", "Control"))) %>% 
  dplyr::filter(gene %in% rtks) %>%
  ggplot(aes(log2FoldChange, -log10(pvalue), label = gene, colour = colour))+
        geom_point()+
        geom_text_repel(size = 3.5,
                        segment.size = 0.2, box.padding = 0.5,
                        color = "black",
                        segment.color = "black", max.overlaps = Inf)+ 
        scale_color_manual(values = alpha(c(mypal), 0.9),
                           breaks = c("Knockout", "Control"))+
        theme(panel.background = element_blank(),
              panel.border = element_rect(colour = "black", fill=NA, size=1),
              plot.title = element_text(size=10),
              legend.key = element_blank(),
              legend.title = element_blank(),
              # legend.position = c(.1, .95),
              legend.justification = c("left", "top"),
              legend.margin = margin(rep(1,4)))+
        xlab(expression('log'[2]*'(FoldChange)'))+
        ylab(expression('-log'[10]*'P'))+
        xlim(-2, 2)+
        guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()



# endocytosis subset volcano ----------------------------------------------

# using gene set: KEGG_ENDOCYTOSIS

# reread
res_filt <- readRDS(glue::glue("{repo_dir}/outs/rna_deseq_res_filt.rds"))

# genes
endo_genes <- glue::glue("{repo_dir}/src/objs/KEGG_ENDOCYTOSIS.v2023.2.Hs.grp") %>% 
  data.table::fread(skip = 1) %>% 
  dplyr::pull() %>% 
  gprofiler2::gorth() %>% 
  dplyr::pull(ortholog_name) %>% 
  unique()

lfc_thresh_endo <- 0.1
p_thresh_endo <- 0.1
pdf(glue::glue("{outdir}/volcano_endo_kcnb2_ko_vs_wt.pdf"), w = 6, h = 5)
res_filt %>%
  mutate(colour = ifelse(pvalue > p_thresh_endo | abs(log2FoldChange) < lfc_thresh_endo, "not_sig",
                         ifelse(log2FoldChange > 0, "Knockout", "Control")),
         label = ifelse(colour != "not_sig", gene, "")) %>% 
  dplyr::filter(gene %in% endo_genes) %>%
  ggplot(aes(log2FoldChange, -log10(pvalue), label = label, colour = colour))+
  geom_point()+
  geom_text_repel(size = 2.5,
                  segment.size = 0.2, box.padding = 0.3,
                  color = "black",
                  segment.color = "black", max.overlaps = Inf)+ 
  scale_color_manual(values = alpha(c(mypal), 0.9),
                     breaks = c("Knockout", "Control"))+
  theme(panel.background = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(size=10),
        legend.key = element_blank(),
        legend.title = element_blank(),
        # legend.position = c(.1, .95),
        legend.justification = c("left", "top"),
        legend.margin = margin(rep(1,4)))+
  xlab(expression('log'[2]*'(FoldChange)'))+
  ylab(expression('-log'[10]*'P'))+
  xlim(-2, 2)+
  guides(colour = guide_legend(override.aes = list(size=3)))
dev.off()



# go enrichment -----------------------------------------------------------

go_ko <- ggdf %>% 
  dplyr::filter(padj < padj_thresh & log2FoldChange > lfc_thresh) %>% 
  dplyr::pull(gene) %>% 
  gprofiler2::gost(ordered_query = TRUE, exclude_iea = TRUE) %>% 
  .[["result"]] %>% 
  dplyr::select(term_name, p_value:intersection_size, source)

go_wt <- ggdf %>% 
  dplyr::filter(padj < padj_thresh & log2FoldChange < lfc_thresh) %>% 
  dplyr::pull(gene) %>% 
  gprofiler2::gost(ordered_query = TRUE, exclude_iea = TRUE) %>% 
  .[["result"]] %>% 
  dplyr::select(term_name, p_value:intersection_size, source)

pdf(glue::glue("{outdir}/go_enrichment_kcnb2_ko.pdf"), w = 7, h = 8)
go_ko %>% 
  dplyr::filter(source != "TF") %>% 
  ggplot(aes(x = forcats::fct_reorder(term_name, -p_value), y = -log10(p_value), fill = -p_value))+
  geom_col()+
  coord_flip()+
  viridis::scale_fill_viridis(option = "rocket")+
  theme_classic()+
  theme(legend.position = "none", axis.text.y = element_text(size = 12, color = "black"))+
  xlab("")
dev.off()

pdf(glue::glue("{outdir}/go_enrichment_kcnb2_wt.pdf"), w = 7, h = 6)
go_wt %>% 
  dplyr::filter(source != "TF") %>% 
  ggplot(aes(x = forcats::fct_reorder(term_name, -p_value), y = -log10(p_value), fill = -p_value))+
  geom_col()+
  coord_flip()+
  viridis::scale_fill_viridis(option = "mako")+
  theme_classic()+
  theme(legend.position = "none")+
  xlab("")
dev.off()



