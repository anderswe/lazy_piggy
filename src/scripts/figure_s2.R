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
library(Rsubread)
library(biomaRt)
library(gprofiler2)
source("src/scripts/utils.R")
source("src/scripts/configs.R")

# outdir
outdir <- glue::glue("{repo_dir}outs")

# Panel B ------------------------------------------------------------
# produced in src/scripts/overlap_map.R


# Panels D & E ------------------------------------------------------------
# insertion plots 

# get gene ref from biomaRt
ref <- getInBuiltAnnotation("mm9") %>% clean_names()

# convert entrez id to gene symbol
mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
dict <- getBM(attributes = c("mgi_symbol", "entrezgene_id"), filters = "entrezgene_id", values = ref$gene_id, bmHeader = T, mart = mart) %>% clean_names()
ref$symbol <- plyr::mapvalues(ref$gene_id, dict$ncbi_gene_formerly_entrezgene_id, dict$mgi_symbol, warn_missing = F)


# Inputs
genes <- c("Kcnh2", "Kcnb1", "Kcnip3")
buffer <- 1000

# filtered insertions
insertion_path <- glue("{awe_local_dir}20211201/all_annotations_20211201.annot.gz") # for reviewers: this is the same as for_reviewers/processed_data/wgs_lp_tumours/all_annotations_20211201.csv.gz

# import insertion locations
pb <- read_csv(insertion_path, col_names = T)

# loop
for(i in genes){

      # subset
      insertions <- pb %>% filter(gene == i)
      y_val <- length(insertions$loc) + 1
      insertions %<>% mutate(start = loc+buffer,
            end = loc-buffer,
            y_coord = seq(2,y_val, (y_val-2)/length(gene_location))[-1],
            orientation = ifelse(gene_orientation == "SAME", 1, 0)) #counterintuitive
      
      # get coords for whole gene
      gene_coords <- ref %>% filter(symbol == i) %>%
      mutate(orientation = ifelse(strand == "-", 0, 1),
            start = min(start),
            end = max(end)) %>% 
      dplyr::rename(gene = symbol) %>% .[1,]

      # get coords for individual transcripts
      exon_coords <- ref %>% filter(symbol == i) %>%
      dplyr::rename(from = start, to = end, gene = symbol)

      # number exons in order, based on transcription orientation
      if(exon_coords$strand[1] == "-"){
            exon_coords %<>% arrange(from) %>% 
                  mutate(exon_number = length(exon_coords$from):1)
      }else{
            exon_coords %<>% arrange(from) %>% 
                  mutate(exon_number = 1:length(exon_coords$from))
      }

      # final cleaning on exon_coords
      exon_coords %<>% mutate(orientation = ifelse(strand == "-", 0, 1),
            exon = glue("Ex{exon_number}"),
            start = gene_coords$start,
            end = gene_coords$end)

      # plot
      gg <- ggplot(gene_coords, aes(xmin = start, xmax = end, y = 1, fill = gene, forward = orientation, label = gene)) + 
            geom_gene_arrow(fill = "white")+
            geom_gene_arrow(data=insertions, aes(xmin = start, xmax = end, y = y_coord, forward = orientation, fill = gene_orientation),
                              arrowhead_width = grid::unit(0.01, "npc"),#grid::unit(abs(min(insertions$loc) - max(insertions$loc)) / 57000, "mm"),#grid::unit(1.5, "mm"),
                              arrowhead_height = grid::unit((y_val - 2)/10, "mm"),#grid::unit(2, "mm"),
                              arrow_body_height = grid::unit((y_val - 2)/20, "mm"))+#grid::unit(1, "mm"))+
            geom_subgene_arrow(data = exon_coords,
                              aes(xmin = start, xmax = end, y = 1, fill = gene,
                                    xsubmin = from, xsubmax = to), color="black", alpha=.7)+
            geom_gene_label()+ #for central Kcnb1 label
            geom_text(data=insertions, aes(x=loc, label = "|"), nudge_y = 1, size = 1.9)+ #bars for exact insertion sites
            geom_text(data=exon_coords %>% mutate(spot = (from + to)/2) %>% filter(exon_number %in% c(max(exon_number), min(exon_number))), aes(x=spot, label = exon), nudge_y = -1, size = 2)+ #for exon labels
            theme_bw() +
            scale_fill_manual(values = c("#7FC97F", "#F0027F", "#386CB0"))+
            theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  panel.background = element_rect(colour = "black", size=1),
                  legend.key = element_blank(),
                  legend.title = element_blank(),
                  legend.position = "none",
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(),
                  axis.title.x=element_text(size=8))+
            ylab("")+
            xlab(glue("Chromosome {gene_coords$chr %>% str_replace('chr','') %>% unique()}"))+
            ylim(c(0, y_val))
      
      pdf(glue("{outdir}/insertion_plot_{i}.pdf"), width = 5, height = 4)
      print(gg)
      dev.off()

}



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
pdf(glue("{outdir}/pca_lazy_piggy_rna_seq.pdf"), height = 15, width = 5)
plot_grid(pca_by_tam,pca_by_tam_days,pca_by_sex,ncol = 1)
dev.off()


