# Info --------------------------------------------------------------------

# Figure S2 reproducibility code
# Feb 4, 2023
# Anders E.
# Taylor lab


# Libraries ---------------------------------------------------------------
library(tidyverse)
library(readxl)
library(janitor)
library(glue)
library(survival)
library(survminer)
library(magrittr)
library(data.table)
library(tidytable)
library(gggenes)
source("src/scripts/utils.R")
source("src/scripts/configs.R")

# outdir
outdir <- glue::glue("{repo_dir}outs")



# Panels D & E ------------------------------------------------------------
# insertion plots 


# to add


# Panel F -----------------------------------------------------------------
# PCA reductions for bulk RNAseq from lazy piggy tumours

# palette
mypal <- c(as.character(pal['red']), "grey", as.character(pal['blue']))

# imports
vsd <- readRDS(glue("{repo_dir}/outs/vsd_rna_seq.rds"))

# get labels
labels <- plotPCA(vsd, intgroup = c("Sex"))$labels

# by tam status
pca_by_tam <- ggplot(plotPCA(vsd, intgroup = c("TamStatus"), returnData = TRUE), aes(x=PC1,y=PC2 )) +
  geom_point(aes(fill=TamStatus),colour="black",pch=21, size=3) +
  scale_fill_manual(values=alpha(c(mypal[1], mypal[3]), 1), breaks=c("TAMPLUS", "TAMMINUS"), labels = c("Tamoxifen+", "Tamoxifen-")) +
  xlab(labels$x) +
  ylab(labels$y) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.margin=margin(0,10,0,0),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.position = c(0.15,0.1))

# by sex
pca_by_sex <- ggplot(plotPCA(vsd, intgroup = c("Sex"), returnData = TRUE), aes(x=PC1,y=PC2 )) +
  geom_point(aes(fill=Sex),colour="black",pch=21, size=3) +
  scale_fill_manual(values = c("#F26463", "#A7D1E8"), breaks=c("F", "M"), labels = c("Female", "Male")) +
  xlab(labels$x) +
  ylab(labels$y) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        legend.key = element_blank(),
        legend.title = element_blank(),
        legend.margin=margin(0,10,0,0),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.position = c(0.15,0.1))


# by tam days
pca_by_tam_days <- ggplot(plotPCA(vsd, intgroup = c("Sex"), returnData = TRUE) %>% left_join((LP_metadata %>% dplyr::select(name = Mouse, tam_days = DaysOnTamoxifin)), by = "name"),
                          aes(x=PC1,y=PC2 )) +
  geom_point(aes(fill=tam_days),colour="black",pch=21, size=3) +
  scale_fill_viridis(name = "Tamoxifen\n(days)")+ 
  xlab(labels$x) +
  ylab(labels$y) +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(colour = "black", size=1),
        legend.margin=margin(0,10,0,5),
        legend.box.margin=margin(-10,-10,-10,-10),
        legend.direction = "horizontal",
        legend.position = c(0.3,0.1))

# print
pdf(glue("{output_dir}/pca_lazy_piggy_rna_seq.pdf"), height = 15, width = 5)
plot_grid(pca_by_tam,pca_by_tam_days,pca_by_sex,ncol = 1)
dev.off()


