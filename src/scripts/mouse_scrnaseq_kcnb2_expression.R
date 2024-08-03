# Info --------------------------------------------------------------------
#
# Anders E.
# Taylor lab
# 14 May 2024
# 
# goal here is to visualize Kcnb2 expression in mouse Math1Cre; SmoM2 tumours
# 


# Libraries ---------------------------------------------------------------
# salloc -N 1 -c 1 --mem 48G -t 3:00:00
# srun --pty bash
# cd /hpf/largeprojects/mdtaylor/aerickson/data/clones/lazy_piggy
# module load R/4.2.1
# R
library(magrittr)
library(Seurat)
library(ggplot2)
library(rhdf5)
source("src/scripts/utils.R")


# Import ---------------------------------------------------------------
# read in Siyi objects
siyi_objs <- list.files("src/data/scrnaseq/GSE197402_RAW", full.names = T)

siyi_ol <- purrr::map(siyi_objs, \(x){ 
  
  h5 <- rhdf5::h5read(x, name = "matrix")
  
  counts <- Matrix::sparseMatrix(
    dims = h5$shape,
    i = as.numeric(h5$indices),
    p = as.numeric(h5$indptr),
    x = as.numeric(h5$data),
    index1 = FALSE) # thank you: https://gist.github.com/slowkow/d3c4b77c9bf2a75f6dad4843d7d3aefc
  
  bc <- h5$barcodes
  gn <- h5$features$name %>% make.names(unique = T)
  
  colnames(counts) <- bc
  rownames(counts) <- gn
  
  CreateAssayObject(counts) %>% 
  CreateSeuratObject() %>% 
    # SCTransform(method = "glmGamPoi", verbose = F) %>%
    # RunPCA(verbose = F) %>% 
    # RunUMAP(dims = 1:20, verbose = F) %>%
    # FindNeighbors(dims = 1:20, verbose = F) %>%
    # FindClusters(resolution = 0.2, verbose = F) %>% 
    return()
  
})
siyi_ol[[1]]$sample_id <- "P21_1"
siyi_ol[[2]]$sample_id <- "P21_2"

# read in Jerry object
jerry_path <- "src/data/scrnaseq/mskp7wt_x6"

mtx <- Matrix::readMM(list.files(jerry_path, pattern = "mtx", full.names = T))
bc <- data.table::fread(list.files(jerry_path, pattern = "barcodes", full.names = T), header = F) %>% dplyr::pull()
gn <- data.table::fread(list.files(jerry_path, pattern = "features", full.names = T), header = F) %>% dplyr::pull(V2) %>% make.names(unique = T)

rownames(mtx) <- gn
colnames(mtx) <- bc

jo <- CreateSeuratObject(mtx)# %>% 
  # SCTransform(method = "glmGamPoi", verbose = F) %>%
  # RunPCA(verbose = F) %>% 
  # RunUMAP(dims = 1:20, verbose = F) %>%
  # FindNeighbors(dims = 1:20, verbose = F) %>%
  # FindClusters(resolution = 0.2, verbose = F)
jo$sample_id <- "P7"



# Integrate ---------------------------------------------------------------

so <- merge(jo, siyi_ol) %>% 
  # SCTransform(assay = "RNA", method = "glmGamPoi", verbose = T) %>%
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(verbose = T) %>% 
  harmony::RunHarmony(group.by.vars = "sample_id",
             reduction = "pca",
             assay.use = "SCT") %>% 
  FindNeighbors(dims = 1:20, verbose = T, reduction = "harmony") %>%
  FindClusters(resolution = 0.2, verbose = T) %>%
  RunUMAP(dims = 1:20, verbose = T, reduction = "harmony") %>% 
  JoinLayers()

# cell cycle
orth <- data.table::fread("src/data/orthologs/new_ortho.csv")
s_genes <- cc.genes.updated.2019$s.genes %>% plyr::mapvalues(orth$gene_name, orth$mouse_gene_name, warn_missing = FALSE) %>% unique()
g2m_genes <- cc.genes.updated.2019$g2m.genes %>% plyr::mapvalues(orth$gene_name, orth$mouse_gene_name, warn_missing = FALSE) %>% unique()

so %<>% CellCycleScoring(s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)


# quick save
saveRDS(so, glue::glue("src/data/scrnaseq/integrated.rds"))


# Initial Plots ---------------------------------------------------------------

p1 <- DimPlot(so, group.by = c("seurat_clusters", "sample_id"), raster = T, label = T) & NoAxes() & ggtitle("")
p2 <- FeaturePlot(so, "Kcnb2", raster = T) & NoAxes() & NoLegend()

pdf("outs/mouse_math1cre_smom2_kcnb2_expression.pdf", w = 15, h = 5)
(p1 | p2) + patchwork::plot_layout(widths = c(2,1))
dev.off()


# Subset to just tumour cells ---------------------------------------------------------------


ss <- so %>% 
  subset(subset = seurat_clusters %in% c(0:6, 11)) %>% 
  NormalizeData(normalization.method = "LogNormalize", scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>% 
  RunPCA(verbose = T) %>% 
  harmony::RunHarmony(group.by.vars = "sample_id",
             reduction = "pca",
             assay.use = "SCT") %>% 
  FindNeighbors(dims = 1:20, verbose = T, reduction = "harmony") %>%
  FindClusters(resolution = 0.2, verbose = T) %>%
  RunUMAP(dims = 1:20, verbose = T, reduction = "harmony") %>% 
  JoinLayers()

saveRDS(ss, glue::glue("src/data/scrnaseq/integrated_subset.rds"))

phase_plot <- DimPlot(ss, group.by = "Phase", cols = phase_pal)+NoAxes()+ggtitle("")
fp_list <- FeaturePlot(ss, c("Sox2", "Dcx", "Kcnb2"), order = T, ncol = 3) & NoAxes() & NoLegend() #, cols = rev(RColorBrewer::brewer.pal(11, "RdBu"))

pdf("outs/mouse_math1cre_smom2_kcnb2_expression_tumor_subset.pdf", w = 8, h = 2)
(phase_plot | fp_list) + patchwork::plot_layout(widths = c(1, 3), guides = "collect")
dev.off()
