# Info --------------------------------------------------------------------
#
#
# Anders E.
# Taylor lab
# Jan 31, 2023
# modified from Patryk's original script


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


# Inputs ------------------------------------------------------------------
rna_dir <- glue("{proj_dir}/Lazy_Piggy_Analysis/RNAseq/Expression")
outdir <- glue("{repo_dir}/outs/")
LP_count_matrix <- read.table(file = glue("{rna_dir}/gene_matrix.txt"), header = TRUE, row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
LP_metadata <- read.table(file = glue("{rna_dir}/LP_sample_metadata.txt"), header = TRUE, stringsAsFactors = FALSE, check.names = FALSE, sep = "\t")

# Cleaning ------------------------------------------------------------------

# Remove missing sample from the matrix
LP_count_matrix <- LP_count_matrix[,!colnames(LP_count_matrix) %in% c("4-7-15KWC")]

# Only include metadata samples present in the matrix
LP_metadata <- LP_metadata[LP_metadata$Mouse %in% colnames(LP_count_matrix),]

# Sort the metadata and the matrix to the same order
LP_count_matrix <- LP_count_matrix[,order(colnames(LP_count_matrix))]
LP_metadata <- LP_metadata[order(LP_metadata$Mouse),]

# Remove symbols in the factor names
LP_metadata[LP_metadata$TamStatus == "TAM+",]$TamStatus="TAMPLUS"
LP_metadata[LP_metadata$TamStatus == "TAM-",]$TamStatus="TAMMINUS"

# convert design variables to factors with proper level order
LP_metadata$TamStatus=factor(LP_metadata$TamStatus, levels = c("TAMMINUS", "TAMPLUS"))
LP_metadata$Sex=factor(LP_metadata$Sex, levels = c("M", "F"))

# keep only genes with non-zero expression in at least 2 samples per group
pos_samples <- LP_metadata %>% filter(TamStatus == "TAMPLUS") %>% pull(Mouse) #get names of TAM+ samples
neg_samples <- LP_metadata %>% filter(TamStatus == "TAMMINUS") %>% pull(Mouse) #get names of TAM- samples
na_mtx <- LP_count_matrix #duplicate count matrix
na_mtx[na_mtx == 0] <- NA #convert 0s to NAs in the duplicate
na_counts <- data.frame(gene = rownames(LP_count_matrix), pos = "", neg = "") #create an empty df to populate with counts of NAs
na_counts$pos <- apply(na_mtx,1,function(x) sum(is.na(x[pos_samples]))) #count em
na_counts$neg <- apply(na_mtx,1,function(x) sum(is.na(x[neg_samples]))) # ""
na_counts$pos_non_zero <- length(pos_samples) - na_counts$pos #get number of non-zero samples per group
na_counts$neg_non_zero <- length(neg_samples) - na_counts$neg # ""
genes_keep <- na_counts %>% filter(pos_non_zero >= 2 & neg_non_zero >= 2) %>% pull(gene) #get gene names
LP_count_matrix <- LP_count_matrix[genes_keep, ]



# Run DESeq2 --------------------------------------------------------------

# deseq2
dds <- DESeqDataSetFromMatrix(countData = LP_count_matrix, colData = LP_metadata, design = ~ Sex + TamStatus)
dds <- DESeq(dds)
res <- lfcShrink(dds, coef = "TamStatus_TAMPLUS_vs_TAMMINUS", type = "apeglm")
res_filt <- res %>% as.data.frame() %>% filter(!is.na(pvalue))

# save for plotting
saveRDS(res_filt, glue("{repo_dir}/outs/rna_deseq_res_filt.rds"))




# PCA prep ------------------------------------------------

#Transform the data
vsd <- varianceStabilizingTransformation(dds, blind = FALSE)

# outdir
output_dir <- glue("{repo_dir}/outs")

# save
saveRDS(vsd, glue("{repo_dir}/outs/vsd_rna_seq.rds"))

# PCA plotting in figure_s2.R

# UMAP calculation ------------------------------------------------

normalized_counts <- assay(vsd) %>% t()
umap <- umap::umap(normalized_counts)

umap_df <- umap$layout %>% 
  as.data.frame() %>% 
  dplyr::rename(UMAP1 = 1, UMAP2 = 2) %>% 
  as_tibble(rownames = "Mouse") %>%
  left_join(LP_metadata, by = "Mouse")


# UMAP plotting ------------------------------------------------
umap_by_tam <- ggplot(umap_df, aes(UMAP1, UMAP2))+
  geom_point(aes(fill=TamStatus),colour="black",pch=21, size=3)+
  theme_classic()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.margin=margin(0,10,0,0),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.position = c(0.9,0.9))+
  scale_fill_manual(values = c(mypal[1], mypal[3]),
                    breaks=c("TAMPLUS", "TAMMINUS"),
                    labels = c("TAM+", "TAM-"))

umap_by_sex <- ggplot(umap_df, aes(UMAP1, UMAP2))+
  geom_point(aes(fill=Sex),colour="black",pch=21, size=3)+
  theme_classic()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.margin=margin(0,10,0,0),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.position = c(0.85,0.9))+
  scale_fill_manual(values = c("#F26463", "#A7D1E8"),
                    breaks=c("F", "M"),
                    labels = c("Female", "Male"))

umap_by_tam_days <- ggplot(umap_df, aes(UMAP1, UMAP2))+
  geom_point(aes(fill=DaysOnTamoxifin),colour="black",pch=21, size=3)+
  theme_classic()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        legend.key = element_blank(),
        legend.margin=margin(0,10,0,5),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.position = c(0.75,0.9),
        legend.direction = "horizontal")+
  scale_fill_viridis(name = "Tamoxifen\n(days)")

pdf(glue("{output_dir}/umaps_lazy_piggy_rna_seq.pdf"), height = 15, width = 5)
plot_grid(umap_by_tam, umap_by_tam_days, umap_by_sex, ncol = 1)
dev.off()





