# Info --------------------------------------------------------------------

# Figure S4 reproducibility code
# Aug 28, 2023
# Anders E.
# Taylor lab


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(readxl)
library(janitor)
library(glue)
library(magrittr)
library(Seurat)
library(Nebulosa)
library(BiocFileCache)
library(GEOquery)
library(harmony)
library(data.table)
library(destiny)
library(viridis)
library(job)
library(cowplot)
source("src/scripts/utils.R")
source("src/scripts/configs.R")

# outdir
outdir <- glue::glue("{repo_dir}outs")


# Panel B (Riemondy) -----------------------------------------------------------------

# this object was made by downloading harmonized umap embeddings, metadata,
# and the expression matrix for this object from:
# UCSC Cell Browser: https://d33sxa6bpqwi51.cloudfront.net/?ds=human

so <- readRDS(glue::glue("{riemondy_dir}/SHH.rds")) # so = seurat object
md <- read_excel(glue::glue("{riemondy_dir}/sampleInfo_suppTable_1.xlsx"), skip = 1) %>% janitor::clean_names()

# phase
so %<>% Seurat::CellCycleScoring(s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes, set.ident = F)

# featureplots for manuscript
fp_list_riem <- purrr::map(c("SOX2", "DCX", "KCNB2"), \(x){
  FeaturePlot(so, x, order = T, cols = rev(RColorBrewer::brewer.pal(11, "RdBu")))+NoAxes()+guides(colour = guide_colorbar(ticks = FALSE))
})
phase_pal <- c("#2B83BA", "#FDAE61", "#D7191C")
phase_plot_riem <- DimPlot(so, group.by = "Phase", cols = phase_pal)+NoAxes()+ggtitle("")

pdf("{outdir}/riemondy_featureplots.pdf", h = 4, w = 4)
phase_plot_riem + fp_list_riem[[1]] + fp_list_riem[[2]] + fp_list_riem[[3]] + patchwork::plot_layout(guides = "collect")
dev.off()



# Panel B (Hovestadt) -----------------------------------------------------------------


setwd(northcott_dir)

#accession code
acc <- "GSE119926" #ref: Hovestadt V, Smith KS, Bihannic L, Filbin MG et al. Resolving medulloblastoma cellular architecture by single-cell genomics. Nature 2019 Aug;572(7767):74-79

# get the series matrix
gse <- getGEO(acc)[[1]]
md <- pData(gse) %>% clean_names() %>% filter(grepl("SHH", characteristics_ch1_2))
keys <- md$title %>% paste0(collapse = "|")#get accessions for SHH samples  --> went and downloaded supplementary from GEO, and stored in getwd()

# import all
raw <- map(list.files(path = getwd(), full.names = T, pattern = keys), \(x){
  df <- read_tsv(x, skip = 1, col_names = F) %>% column_to_rownames("X1")
  colnames(df) <- read_tsv(x, col_names = F)[1,] %>% t() %>% as.character()
  return(df)
})

# specs
dims <- 10 #picked from ElbowPlot(so) 
res <- 0.2

# convert to seurat objects
sos <- map(raw, \(x){
  
  #NOTE: the input values are actually TPM!!!
  CreateSeuratObject(counts = x, assay = "RNA") %>% 
    # SCTransform(verbose = FALSE) %>%
    FindVariableFeatures() %>% 
    ScaleData(verbose = F) %>% 
    RunPCA(verbose = F) %>% 
    RunUMAP(dims = 1:dims, verbose = F) %>%
    FindNeighbors(dims = 1:dims, verbose = F) %>%
    FindClusters(resolution = res, verbose = F) %>% 
    return()
}) %>% setNames(c("MUV41", "SJ454", "SJ577", "RCMB18", "RCMB24")) #hard-coding here. confirm correct.

# set orig.ident labels
sos$MUV41$orig.ident <- "MUV41"
sos$SJ454$orig.ident <- "SJ454"
sos$SJ577$orig.ident <- "SJ577"
sos$RCMB18$orig.ident <- "RCMB18"
sos$RCMB24$orig.ident <- "RCMB24"

# merge all
all <- merge(sos[[1]], sos[2:5]) %>%
  # SCTransform(assay = "RNA", verbose = F) %>%
  FindVariableFeatures(verbose = F) %>% 
  ScaleData(verbose = F) %>% 
  RunPCA(verbose = F) %>%
  RunHarmony(group.by.vars = "orig.ident",
             reduction = "pca",
             assay.use = "RNA", #SCT
             theta = 2) %>%
  FindNeighbors(dims = 1:dims, verbose = FALSE, reduction = "harmony") %>%
  FindClusters(resolution = res, verbose = FALSE) %>%
  RunUMAP(dims = 1:dims, verbose = FALSE, reduction = "harmony")

all %<>% CellCycleScoring(s.features = cc.genes.updated.2019$s.genes, g2m.features = cc.genes.updated.2019$g2m.genes)


# export for manuscript
fp_list <- purrr::map(c("SOX2", "DCX", "KCNB2"), \(x){
  FeaturePlot(all, x, order = T, cols = rev(RColorBrewer::brewer.pal(11, "RdBu")))+NoAxes()+guides(colour = guide_colorbar(ticks = FALSE))
  })
phase_pal <- c("#2B83BA", "#FDAE61", "#D7191C")
phase_plot <- DimPlot(all, group.by = "Phase", cols = phase_pal)+NoAxes()+ggtitle("")

pdf(glue::glue("{outdir}/northcott_featureplots.pdf"), h = 4, w = 4)
phase_plot + fp_list[[1]] + fp_list[[2]] + fp_list[[3]] + patchwork::plot_layout(guides = "collect")
dev.off()

