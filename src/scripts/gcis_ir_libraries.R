# Info --------------------------------------------------------------------
#
# Anders E.
# Taylor lab
# Jan 31, 2023
# modified from Patryk's original gCIS script


# Notes -------------------------------------------------------------------
# here we will to establish PBR and PBL
# cutoffs separately first, THEN filter insertions, THEN combine left and right, THEN calculate + plot jaccard similarity matrices.
# and will do the same for IR and JX libraries


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
input_file <- glue("{awe_local_dir}20211201/all_annotations_20211201.annot.gz") # for reviewers: this is the same as for_reviewers/processed_data/wgs_lp_tumours/all_annotations_20211201.csv.gz
output_dir <- glue("{repo_dir}outs") 
metadata_file <- glue("{awe_local_dir}20220524/metadata_20220524.csv")
gene_annot_input <- glue("{awe_local_dir}20211104/gene_annot_gCIS_15000.txt")
lib_type <- "restriction"
min_rec <- 3
window_size <- 5
pvalue_correction <- "bonferroni"
qvalue_threshold <- 0.1
recurrence_filter <- FALSE 

#Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)


# load metadata -----------------------------------------------------------
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



# Read processing ---------------------------------------------------------

## load 
print("Loading insertions")
insertions <- fread(input_file, nThread = 8) %>% #read in
  rename.(Sample_key = sample) %>% 
  mutate.(strand = "+", # strand unavailable
          ULP = NA) %>% 
  separate.(Sample_key,
            into = c("patient", "sample", "orientation"),
            sep = "_",
            remove = FALSE)


## Calculate clusters
# If an insertion is found with {window_size} bp of another insertion in the same library, both insertions are but within the same cluster.
print(paste("Clustering Insertions within ", window_size, " bp", sep=""))
insertions %<>% 
  mutate.(cluster=NA) %>%
  arrange.(Sample_key, chr, loc) %>% 
  lazy_dt() %>% 
  group_by(Sample_key) %>% 
  group_modify(cluster_insertions_group) %>% 
  as.data.table()

## Collapse the clusters to pick the most likely insertion
# i.e. assign aggregated cluster counts to insertion site with most counts...within each cluster
print("Collapsing Insertions...")
insertions %<>% 
  mutate.(count_sum = sum(count), .by = c(Sample_key, cluster)) %>% 
  slice_max.(count, n = 1, .by = c(Sample_key, cluster)) %>% 
  arrange.(Sample_key, chr, loc) %>% 
  select.(-count) %>% 
  rename.(count = count_sum)

## Add tam status info
insertions %<>% mutate.(tam_status = plyr::mapvalues(patient, metadata$mouse_id, metadata$tam_status, warn_missing = FALSE))

## save clustered insertions
print("Saving the Data...")
fwrite(insertions, file = glue("{output_dir}/all_insertions_clustered.csv.gz"), nThread = 8)
# insertions <- fread(file = glue("{output_dir}/all_insertions_clustered.csv.gz"), nThread = 8)

# Prelim filters ----------------------------------------------------------

## Non-standard chromosomes
print("Removing insertions in non-standard chromosomes...")
standard_chr <- c(glue("chr{1:22}"), "chrX", "chrY") 
insertions %<>% filter.(chr %in% standard_chr)

## Donor Filter
print("Removing insertions in donor chromosomes...")
donor_md <- metadata %>% select.(mouse_id, donor) %>% rename.(patient = mouse_id) %>% unique()
insertions %<>% left_join.(donor_md, by = "patient") %>% filter.(chr != donor) 

## Remove insertions with single read support
insertions %<>% filter.(count > 1)

## Recurrent/contaminant insertions + salvage
if(recurrence_filter == TRUE){
  print("Removing insertions found in multiple mice...")
  
  insertions %<>% mutate.(insertion_key = paste(chr,loc,strand, sep="_"))
  
  #Blacklist of recurrent insertions (recurrent insertions within the same mouse are preserved)
  keys <- insertions %>% 
    select.(patient, insertion_key) %>%
    unique() %>% 
    select.(insertion_key) %>%
    pull.() 
  reccurent_insertions <- unique(keys[duplicated(keys)])
  
  #Salvage Insertions
  salvaged_insertions <- insertions %>% 
    mutate.(insertion_key=paste(chr,loc,strand, sep="_")) %>% 
    filter.(insertion_key %in% reccurent_insertions) %>% 
    slice_max.(count, n = 1, .by = insertion_key) %>% 
    select.(-insertion_key)
  
  
  #Remove blacklisted insertions
  insertions %<>% 
    filter.(!(insertion_key %in% reccurent_insertions)) %>% 
    select.(-insertion_key)
  
  #Add back salvaged insertions 
  print(paste("Salvaged ", dim(salvaged_insertions)[1], " insertions...", sep=""))
  insertions <- unique(bind_rows(insertions, salvaged_insertions))
}

## Remove blacklisted insertions i.e. insertions in CTRL mice libraries
print("Removing insertions found in control mice...")
blacklist_samples <- "PTCH"

blacklist_keys <- insertions %>%
  filter.(patient %in% blacklist_samples) %>% 
  mutate.(key = glue("{chr}_{loc}_{orientation}")) %>% 
  pull.(key) %>% 
  unique()

insertions %<>% 
  mutate.(insertion_key = glue("{chr}_{loc}_{orientation}")) %>% 
  filter.(!(insertion_key %in% blacklist_keys)) %>% 
  select.(-insertion_key)


## filter out non experimental samples
ctrl_list <- metadata %>% filter.(tam_status == "CTRL") %>% pull.(mouse_id) %>% unique()
insertions %<>% filter.(!(patient %in% ctrl_list))


# calculate gcis thresholds ---------------------------------------------------------

## Top 1 percent read count filter & Top 0.1 percent read sum filter

print("Top 1 percent read count filter and Top 0.1 percent read sum filter...")
insertions %<>% 
  mutate.(top_1per_max_insertion = max(count)*0.01,
          sum_reads_0.1_per = sum(count)*0.001,
          .by = Sample_key)



## 95% percentile of negative binomial
print("95% percentile of negative binomial filter......")
insertions %<>% 
  lazy_dt() %>% 
  group_by(Sample_key) %>% 
  group_modify(neg_binom_thresh_groups) %>% 
  as.data.table()
insertions$neg_binom_thresh[insertions$neg_binom_thresh == 0] <- NA



# order thresholds by stringency
insertions %<>% 
  rowwise() %>% 
  mutate(clonal_thresh = max(c(top_1per_max_insertion, sum_reads_0.1_per, neg_binom_thresh), na.rm = T), #notice the call to mutate() not mutate.() --> rowwise() doesn't play nice with {dtplyr} or {tidytable}
         subclonal_thresh = median(c(top_1per_max_insertion, sum_reads_0.1_per, neg_binom_thresh), na.rm = T),
         low_thresh = min(c(top_1per_max_insertion, sum_reads_0.1_per, neg_binom_thresh), na.rm = T)) %>% 
  ungroup() %>% 
  as.data.table()


## save clustered insertions with filter thresholds 
fwrite(insertions, file = glue("{output_dir}/all_insertions_clustered_with_thresholds.csv.gz"), nThread = 8)
# insertions <- fread(file = glue("{output_dir}/all_insertions_clustered_with_thresholds.csv.gz"), nThread = 8)

# plot insertion counts by tam status
ggdf <- fread(file = glue("{output_dir}/all_insertions_clustered_with_thresholds.csv.gz"), nThread = 8) %>%
  mutate.(lib = str_sub(orientation, 1L, 2L)) %>% 
  slice_max.(count, n = 1, .by = c(patient, lib, chr, loc, gene_orientation)) %>% 
  summarise.(sum_counts = sum(count), .by = c(patient, tam_status))

vln_pal <- c("TAM+" = "red", "TAM-" = "blue")
vln <- ggplot(ggdf, aes(tam_status, sum_counts, color = tam_status))+
  geom_violin(width=0.6)+
  geom_jitter(width=0.1, alpha = 0.5)+
  geom_boxplot(width=0.1, color="black", outlier.shape = NA)+
  ggpubr::stat_compare_means()+
  theme_classic()+
  scale_color_manual(values = vln_pal)+
  theme(legend.position = "none")+
  xlab("")+
  ylab("Num. SB insertion counts per sample")

pdf(glue::glue("{output_dir}/counts_by_sample.pdf"), w = 5, h = 5)
print(vln)
dev.off()

# filter to clonal/subclonal/lowest insertions ---------------------------------------------

# first, save all combinations
insertions %>% filter.(count >= clonal_thresh) %>% fwrite(file = glue("{output_dir}/clonal_insertions.csv.gz"), nThread = 8)
insertions %>% filter.(count >= subclonal_thresh) %>% fwrite(file = glue("{output_dir}/subclonal_insertions.csv.gz"), nThread = 8)
insertions %>% filter.(count >= low_thresh) %>% fwrite(file = glue("{output_dir}/low_insertions.csv.gz"), nThread = 8)

# now, filter to subclonal
insertions %<>% filter.(count >= subclonal_thresh)

# or read in directly if starting here:
# insertions <- fread(file = glue("{output_dir}/subclonal_insertions.csv.gz"), nThread = 8)

# merge left and right ----------------------------------------------------
insertions %<>%
  mutate.(lib = str_sub(orientation, 1L, 2L)) %>% 
  slice_max.(count, n = 1, .by = c(patient, lib, chr, loc, gene_orientation))



# gCIS prep --------------------------------------------------------------------

## Load the gene annotation and probabilities
gene_annot <- fread(gene_annot_input)

# now we need to split by library
hold <- insertions

## Parse the insertions for gCIS
insertions %<>% 
  filter.(lib == "IR") %>% # tweak here
  mutate.(name = paste(patient, sample, sep="_"))


## Count Insertions in Each Gene
cl <- parallel::makeForkCluster(4)
doParallel::registerDoParallel(cl)

fe <- foreach(gene = gene_annot$gene,
              g_chr = gene_annot$chr, 
              g_start = gene_annot$start_buff,
              g_end = gene_annot$end_buff) %dopar% {
                insertions %>% 
                  filter.(chr == g_chr &
                            loc >= g_start &
                            loc <= g_end) %>% 
                  mutate.(gene = gene)
              }
foreach::registerDoSEQ()
gene_annot %<>% mutate.(insertions = fe)


## String of All Samples with Insertion
gene_annot$observed_samples = mapply(function(x){paste(x$name, collapse =";")}, x = gene_annot$insertions)

## Count Observed Insertions in Each Gene
gene_annot$observed_sample_count = mapply(function(x){length(x$name)}, x = gene_annot$insertions)

## Count insertions per gene, by tam status 
tam_pos <- metadata %>% filter.(tam_status == "TAM+") %>% pull.(mouse_id) %>% unique() 
tam_neg <- metadata %>% filter.(tam_status == "TAM-") %>% pull.(mouse_id) %>% unique() 
gene_annot$observed_tam_pos_insertion_count = mapply(function(x){x %>% filter.(patient %in% tam_pos) %>% pull.(name) %>% length()}, x = gene_annot$insertions)
gene_annot$observed_tam_neg_insertion_count = mapply(function(x){x %>% filter.(patient %in% tam_neg) %>% pull.(name) %>% length()}, x = gene_annot$insertions)

## Count mice per gene
gene_annot$observed_patient_count = mapply(function(x){length(unique(x$patient))}, x = gene_annot$insertions)

## Count numbers of tam+ and tam- patients
gene_annot$tam_pos_patients = mapply(function(x){x %>% filter.(patient %in% tam_pos) %>% pull.(name) %>% unique() %>% length()}, x = gene_annot$insertions)
gene_annot$tam_neg_patients = mapply(function(x){x %>% filter.(patient %in% tam_neg) %>% pull.(name) %>% unique() %>% length()}, x = gene_annot$insertions)

## Calculate Expected Insertions in Each Gene
gene_annot$expected_sample_count = gene_annot$TA_insertion_prob*nrow(unique(insertions))


## Calculate Chi-Squared Test Statistic and Pvalue
n <- length(unique(insertions$name))
n_tam_pos <- length(tam_pos)
n_tam_neg <- length(tam_neg)
gene_annot %<>% mutate.(chi_test_stat = ((observed_sample_count - expected_sample_count)^2 / (expected_sample_count*(1-(expected_sample_count / n)))),
                        chi_test_stat_tam_pos = ((observed_tam_pos_insertion_count - expected_sample_count)^2 / (expected_sample_count*(1-(expected_sample_count / n)))),
                        chi_test_stat_tam_neg = ((observed_tam_neg_insertion_count - expected_sample_count)^2 / (expected_sample_count*(1-(expected_sample_count / n)))))

#Chi-Squared P-value
gene_annot$Pvalue = pchisq(gene_annot$chi_test_stat, 1, lower.tail = FALSE)
gene_annot$Pvalue_tam_pos = pchisq(gene_annot$chi_test_stat_tam_pos, 1, lower.tail = FALSE)
gene_annot$Pvalue_tam_neg = pchisq(gene_annot$chi_test_stat_tam_neg, 1, lower.tail = FALSE)

## Correct the P-value
gene_annot$Pvalue_adj = p.adjust(gene_annot$Pvalue, method = pvalue_correction)
gene_annot$Pvalue_adj_tam_pos = p.adjust(gene_annot$Pvalue_tam_pos, method = pvalue_correction)
gene_annot$Pvalue_adj_tam_neg = p.adjust(gene_annot$Pvalue_tam_neg, method = pvalue_correction)

## quick save
saveRDS(gene_annot, glue::glue("{output_dir}/final_gene_annot_ir_libraries.rds"))
# gene_annot <- readRDS(glue::glue("{output_dir}/final_gene_annot_ir_libraries.rds"))

## Filter to Find Significant Genes
print(paste("Filtering the genes using qvalue threshold: ", qvalue_threshold, " and observed sample count min: ",min_rec,sep=""))  
qvalue_threshold <- 0.1

sig_genes <- gene_annot %>% #PICKUP filter by PVAL not sep
  filter.((Pvalue_adj < qvalue_threshold & observed_patient_count >= min_rec)) %>% 
  # filter.((Pvalue_adj_tam_pos < qvalue_threshold & tam_pos_patients >= min_rec) |
  #           (Pvalue_adj_tam_neg < qvalue_threshold & tam_neg_patients >= min_rec)) %>%
  arrange.(Pvalue_adj_tam_pos, Pvalue_adj_tam_neg)

# save
saveRDS(sig_genes, glue("{output_dir}/sig_genes_tam_gcis_separate_no_recurrence_filter_ir_libraries.rds"))



# Panel X -----------------------------------------------------------------
# SB oncoplot

# import
sig_genes <- readRDS(glue("{output_dir}/sig_genes_tam_gcis_separate_no_recurrence_filter_ir_libraries.rds"))
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
                                pct_gp = gpar(fontsize = 8),
                                show_heatmap_legend = FALSE)

pdf(glue("{output_dir}/gcis_oncoprint_ir_libraries.pdf"), h = 26, w = 6)
draw(op)
dev.off()


# go enrichment
go <- sig_genes %>% 
  dplyr::filter(!(gene %in% false_pos)) %>% 
  dplyr::pull(gene) %>% 
  gprofiler2::gost(ordered_query = TRUE, exclude_iea = TRUE) %>% 
  .[["result"]] %>% 
  dplyr::select(term_name, p_value:intersection_size, source)

pdf(glue::glue("{output_dir}/go_enrichment_ir_oncoplot.pdf"), h = 4, w = 6)
ggplot(go %>% dplyr::filter(source != "TF"), aes(x = forcats::fct_reorder(term_name, -p_value), y = -log10(p_value), fill = -p_value))+
  geom_col()+
  coord_flip()+
  viridis::scale_fill_viridis(option = "mako")+
  theme_classic()+
  theme(legend.position = "none")+
  xlab("")
dev.off()


