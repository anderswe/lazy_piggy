# Info --------------------------------------------------------------------
#
# Anders E.
# Taylor lab
# Jan 31, 2023
# modified from Patryk's original overlap map script


# Notes -------------------------------------------------------------------
# here we will plot jaccard similarity matrices for IR / JX / PB libraries.

# Libraries ---------------------------------------------------------------
suppressMessages(library(tidyverse))
suppressMessages(library(fastverse))
suppressMessages(library(tidytable))
suppressMessages(library(foreach))
suppressMessages(library(dtplyr))
suppressMessages(library(glue))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(gespeR))
source("src/scripts/utils.R")
source("src/scripts/configs.R")


# Inputs -------------------------------------------------------------------
input_file <- glue("outs/subclonal_insertions.csv.gz")
annot_path <- glue("{awe_local_dir}Annotated_files_20170802")
metadata_file <- glue("{awe_local_dir}20220524/metadata_20220524.csv")
libs <- c("IR", "JX", "PB")


# Import -------------------------------------------------------------

# load insertions
insertions <- data.table::fread(input_file) %>% 
  mutate.(lib = dplyr::case_when(grepl("IR", orientation) ~ "IR", grepl("PB", orientation) ~ "PB", grepl("JX", orientation) ~ "JX"))


# load metadata
metadata <- fread(metadata_file, nThread = 8) %>% 
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


# Clean -------------------------------------------------------------------


# split libraries
ovlp <- map(libs, \(x) {insertions %>% filter.(lib == x) %>% mutate.(insertion_map = glue("{chr}_{loc}_{gene_orientation}"))}) %>% 
  setNames(libs)


mice <- metadata %>%
  filter.(grepl("TAM", tam_status)) %>% 
  arrange.(tam_status) %>% 
  pull.(mouse_id) %>% 
  unique() %>% 
  .[. %in% unique(insertions$patient)]

# calculate jaccard scores
ir_jx_mtx <- jaccard_matrix(ovlp$IR, ovlp$JX, mice) #takes about a minute for each
ir_pb_mtx <- jaccard_matrix(ovlp$IR, ovlp$PB, mice)
pb_jx_mtx <- jaccard_matrix(ovlp$PB, ovlp$JX, mice)


# prep to plot
label <- "orig_insertions_no_filters"
selection <- "ir_jx_mtx"
mat <- ir_jx_mtx %>% as.matrix()
lp_cols <- pal[c("red", "blue")] %>% as.character()

# reorder matrix
order <- diag(mat) %>%
  as.data.frame() %>% 
  dplyr::rename(jaccard = 1) %>%
  rownames_to_column(var = "mouse_id") %>% 
  left_join.(metadata, by = "mouse_id") %>% 
  arrange.(tam_status, desc(jaccard)) %>% 
  pull.(mouse_id)

hm <- ComplexHeatmap::Heatmap(mat,
                              name = glue("Jaccard Index"),
                              col = viridis::inferno(100), #c("white", "red3"),
                              show_row_names = FALSE,
                              show_column_names = FALSE,
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              row_order = order,
                              column_order = order,
                              row_split = c(rep("A", (metadata %>% filter.(mouse_id %in% mice & tam_status == "TAM+") %>% nrow())),
                                            rep("B", (metadata %>% filter.(mouse_id %in% mice & tam_status == "TAM-") %>% nrow()))),
                              column_split = c(rep("A", (metadata %>% filter.(mouse_id %in% mice & tam_status == "TAM+") %>% nrow())),
                                               rep("B", (metadata %>% filter.(mouse_id %in% mice & tam_status == "TAM-") %>% nrow()))),
                              border = TRUE,
                              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = lp_cols),
                                                                                  labels = c("Tamoxifen+", "Tamoxifen-"), 
                                                                                  labels_gp = gpar(col = "white", fontsize = 10))),
                              left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = lp_cols),
                                                                               labels = c("Tamoxifen+", "Tamoxifen-"), 
                                                                               labels_gp = gpar(col = "white", fontsize = 10))),
                              row_title = "IR Libraries",
                              column_title = "JX Libraries",
                              heatmap_legend_param = list(direction = "vertical",
                                                          # at = c(0, 0.5, 1),
                                                          legend_height = unit(4, "cm"),
                                                          title_position = "leftcenter-rot"))
draw(hm)

pdf(glue("outs/heatmap_jaccard_{selection}_{label}.pdf"), w = 10, h = 10)
draw(hm)
dev.off()

# violin plot
pdf(glue("outs/violin_jaccard_{selection}_{label}.pdf"), w = 3, h = 3)
diag(mat) %>%
  as.data.frame() %>% 
  dplyr::rename(jaccard = 1) %>%
  rownames_to_column(var = "mouse_id") %>% 
  left_join.(metadata, by = "mouse_id") %>% 
  mutate.(tam_status = (ifelse(tam_status == "TAM+", "Tamoxifen+", "Tamoxifen-") %>% factor(levels = c("Tamoxifen+", "Tamoxifen-")))) %>% 
  ggplot(aes(tam_status, jaccard, color = tam_status))+
  geom_violin(width=0.8)+
  geom_boxplot(width=0.2, color="black", alpha = 0.2, outlier.shape = NA)+
  geom_jitter(width=0.1, alpha = 0.5)+
  ggpubr::stat_compare_means(label.x.npc = "center", method = "t.test")+
  theme_classic()+
  scale_color_manual(values = lp_cols)+
  ylab('Jaccard Index (IR vs. JX)')+
  theme(legend.position = "none",axis.title.x = element_blank())
dev.off()


# rbo overlap map ---------------------------------------------------------
# BIG_LETTERS <- c(LETTERS,
#                  do.call("paste0",CJ(LETTERS,LETTERS)),
#                  do.call("paste0",CJ(LETTERS,LETTERS,LETTERS)),
#                  do.call("paste0",CJ(LETTERS,LETTERS,LETTERS,LETTERS)))
# test <- "01-07-16M"
# x <- ovlp$IR %>%
#   filter.(patient == test) %>% 
#   arrange.(desc(count)) %>% 
#   pull.(insertion_map) %>% 
#   unique()
# y <- ovlp$JX %>%
#   filter.(patient == test) %>% 
#   arrange.(desc(count)) %>% 
#   pull.(insertion_map) %>% 
#   unique()
# 
# # probably need to recode insertions to random numerical values and name with SAME ordered names
# dict <- data.frame(insertions = union(x,y),
#                    index = 1:length(union(x,y)))
# a <- plyr::mapvalues(x, dict$insertions, dict$index, warn_missing = FALSE) %>% as.numeric()
# b <- plyr::mapvalues(y, dict$insertions, dict$index, warn_missing = FALSE) %>% as.numeric()
# names(a) <- BIG_LETTERS[1:length(a)]
# names(b) <- BIG_LETTERS[1:length(b)]
# rbo(a, b, p = 0.9)

# abandoning this since rbo values for multiple test mice are all 0...suspect this is related to list length...

# percentile rank overlap -------------------------------------------------

# presumably heavily vulnerable/susceptible to uncertainty as size of intersection between samples decreases
percentile_similarity_matrix <- function(in_dt, ref_dt, mice){
  require(foreach)
  
  # initialize matrix
  mtx <- as.data.frame(matrix(NA,length(mice),length(mice)))
  rownames(mtx) <- mice
  colnames(mtx) <- mice
  
  # nested loop
  foreach(i = mice) %do% {
    foreach(j = mice) %do% {
      mtx[i,j] <- in_dt %>% 
        filter.(patient == i) %>% 
        mutate.(rank = frank(count),
                percentile = rank / length(.$insertion_map)) %>% 
        filter.(insertion_map %in% (ref_dt %>% filter.(patient == j) %>% pull.(insertion_map))) %>% 
        pull.(percentile) %>% 
        sum() %>% 
        `/`(nrow(ref_dt %>% filter.(patient == j)))
    }
  }
  
  return(mtx)
}

# calculate percentile similarity scores
ir_jx_mtx <- percentile_similarity_matrix(ovlp$IR, ovlp$JX, mice) #takes ~10 minutes for each
ir_pb_mtx <- percentile_similarity_matrix(ovlp$IR, ovlp$PB, mice)
pb_jx_mtx <- percentile_similarity_matrix(ovlp$PB, ovlp$JX, mice)

# map version
# using an inverted method: calculating avg percentile of ref_df in in_df (based only on overlapping insertions)
mice <- metadata %>%
  filter.(grepl("TAM", tam_status)) %>% 
  arrange.(tam_status) %>% 
  pull.(mouse_id) %>% 
  unique() %>% 
  .[. %in% unique(insertions$patient)]

percentile_one_vs_many <- function(in_df, ref_df){
  # NB: this crashes for ir/pb and pb/jx (or ir/jx if using subclonal)
  # with error:
  # Error in `dplyr::bind_cols()`:
  # ! Can't recycle `..1` (size 37) to match `..2` (size 57).
  # presumably because there is zero overlap between some sample pairs --> summarise.() step will collapse to a shorter table that map_dfc then tries to cbind
  # can find a manual fix to input zeros
  ref_df %>%
    filter.(patient %in% mice) %>% 
    add_count.(patient, name = "n_insertions") %>% 
    mutate.(rank = frank(count),
            percentile = rank / n_insertions,
            .by = patient) %>%
    filter.(insertion_map %in% in_df$insertion_map) %>% 
    unique() %>%
    summarise.(mean_pctile = mean(percentile), .by = patient) %>% 
    as_tibble() %>% 
    column_to_rownames("patient") %>%
    dplyr::rename(!!unique(in_df$patient) := mean_pctile) %>% 
    return()
  
}
future::plan("multisession")
ir <- ovlp$IR
jx <- ovlp$JX
pb <- ovlp$PB
ir_jx_mtx <- furrr::future_map_dfc(mice, \(y){ ir %>% filter.(patient == y) %>% percentile_one_vs_many(jx) }) %>% .[mice,mice]
ir_pb_mtx <- furrr::future_map_dfc(mice, \(y){ ir %>% filter.(patient == y) %>% percentile_one_vs_many(pb) }) %>% .[mice,mice]
pb_jx_mtx <- furrr::future_map_dfc(mice, \(y){ pb %>% filter.(patient == y) %>% percentile_one_vs_many(jx) }) %>% .[mice,mice]


# prep to plot
label <- "subclonal"
selection <- "ir_jx_mtx"
mat <- ir_jx_mtx %>% as.matrix()
lp_cols <- c("#E773A0", "#526A8F")

# reorder matrix --> pick up here...why renaming???
order <- diag(mat) %>%
  as.data.frame() %>% 
  dplyr::rename(pctile_ovlp = 1) %>%
  rownames_to_column(var = "mouse_id") %>% 
  left_join.(metadata, by = "mouse_id") %>% 
  arrange.(tam_status, desc(pctile_ovlp)) %>% 
  pull.(mouse_id)

hm <- ComplexHeatmap::Heatmap((mat*100),
                              name = glue("Average Percentile in Overlap"),
                              col = viridis::inferno(100), #c("white", "red3"),
                              show_row_names = FALSE,
                              show_column_names = FALSE,
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              row_order = order,
                              column_order = order,
                              row_split = c(rep("A", (metadata %>% filter.(mouse_id %in% mice & tam_status == "TAM+") %>% nrow())),
                                            rep("B", (metadata %>% filter.(mouse_id %in% mice & tam_status == "TAM-") %>% nrow()))),
                              column_split = c(rep("A", (metadata %>% filter.(mouse_id %in% mice & tam_status == "TAM+") %>% nrow())),
                                               rep("B", (metadata %>% filter.(mouse_id %in% mice & tam_status == "TAM-") %>% nrow()))),
                              border = TRUE,
                              top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = lp_cols),
                                                                                  labels = c("TAM+", "TAM-"), 
                                                                                  labels_gp = gpar(col = "white", fontsize = 10))),
                              left_annotation = rowAnnotation(foo = anno_block(gp = gpar(fill = lp_cols),
                                                                               labels = c("TAM+", "TAM-"), 
                                                                               labels_gp = gpar(col = "white", fontsize = 10))),
                              row_title = "IR Libraries",
                              column_title = "JX Libraries",
                              heatmap_legend_param = list(direction = "vertical",
                                                          # at = c(0, 0.5, 1),
                                                          legend_height = unit(4, "cm"),
                                                          title_position = "leftcenter-rot"))
draw(hm)

pdf(glue("{output_dir}/heatmap_percentile_overlap_{selection}_{label}.pdf"), w = 10, h = 10)
draw(hm)
dev.off()

# violin plot
pdf(glue("{output_dir}/violin_percentile_overlap_{selection}_{label}.pdf"), w = 3, h = 3)
diag(mat*100) %>%
  as.data.frame() %>% 
  dplyr::rename(jaccard = 1) %>%
  rownames_to_column(var = "mouse_id") %>% 
  left_join.(metadata, by = "mouse_id") %>% 
  mutate.(tam_status = factor(tam_status, levels = c("TAM+", "TAM-"))) %>% 
  ggplot(aes(tam_status, jaccard, color = tam_status))+
  geom_violin(width=0.8)+
  geom_boxplot(width=0.2, color="black", alpha = 0.2, outlier.shape = NA)+
  geom_jitter(width=0.1, alpha = 0.5)+
  ggpubr::stat_compare_means(label.x.npc = "center", method = "t.test")+
  theme_classic()+
  scale_color_manual(values = lp_cols)+
  ylab('Avg. Percentile in Overlap (IR vs. JX)')+
  theme(legend.position = "none",axis.title.x = element_blank())
dev.off()
