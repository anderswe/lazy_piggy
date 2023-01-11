## =============================================
##    INPUT
## =============================================

input_file <- "/Volumes/ANDERS/data/Lazy Piggy/20211201/all_annotations_20211201.annot" 
output_dir <- "/Users/anders/Library/CloudStorage/OneDrive-UniversityofToronto/PhD/Projects/Sleeping_Beauty/Lazy Piggy/20220524_overlap_maps" 
metadata_file <- "/Volumes/ANDERS/data/Lazy Piggy/20211201/metadata_lp_tam_status.txt"
gene_annot_input="/Volumes/ANDERS/data/Lazy Piggy/20211104/gene_annot_gCIS_15000.txt"
lib_type <- "restriction"
clonality_threshold=0.05
min_rec=3
pvalue_correction="bonferroni"
qvalue_threshold=0.1
lib_type <- "restriction"
suffix <- "LP"
recurrence_filter <- TRUE

#Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## =============================================
##  Notes
## =============================================
# the difference here is we will try to establish PBR and PBL
# cutoffs separately first, THEN filter insertions, THEN combine, THEN perform gCIS.

## =============================================
##  Libraries
## =============================================

suppressMessages(library(tidyverse))
suppressMessages(library(fastverse))
suppressMessages(library(glue))
suppressMessages(library(ComplexHeatmap))

## =============================================
##  Functions
## =============================================
source("/Volumes/ANDERS/scripts/R/lazy_piggy/whole_pipeline_functions_20211208.R")

## =============================================
##  READ PROCESSING
## =============================================
  
  #------------------------------------
  # Load the Insertions
  #------------------------------------
    print("Loading insertions")

    insertions <- fread(input_file, nThread = 8) %>% #read in
      rename(Sample_key = sample) %>% 
      mutate(strand = "+", #Patryk suggested set all to "+" here otherwise need to go back to alignment step to get strand info. haven't looked into how much of a difference this will make.
             ULP = NA, #same as with strand
             patient = str_split_fixed(Sample_key, "_", 3) %>% as.data.frame() %>% pull(V1), #split Sample_key into three variables for compatibility with the rest of the script
             sample = str_split_fixed(Sample_key, "_", 3) %>% as.data.frame() %>% pull(V2),
             orientation = str_split_fixed(Sample_key, "_", 3) %>% as.data.frame() %>% pull(V3)) %>% 
      relocate(Sample_key) %>%  #chuck this up front
      filter(grepl("PB", orientation) | grepl("Primary|PTCH", patient)) #filter to PB library samples, plus controls for blacklist

  #-------------------------------------------
  # Insertion clustering
  #-------------------------------------------
  # If an insertion is found with 5 bp of another insertion (window_size) in the same library, both insertions are but within the same cluster. Each cluster is then collapsed to a single insertion.
    window_size=5
    print(paste("Clustering Insertions within ", window_size, " bp", sep=""))
    insertions %<>% mutate(cluster=NA) %>% group_by(Sample_key, strand) %>% arrange(chr, loc, .by_group = TRUE) %>% group_modify(~ cluster_insertions(.x, window_size)) %>% ungroup()  
    
    #Collapse the clusters to pick the most likely insertion  
    print("Collapsing Insertions...")
    insertions %<>% group_by(patient, sample, orientation, Sample_key, cluster, strand) %>% group_modify(~ collapse_clust(.x, type=lib_type)) %>% ungroup()
    
  #-------------------------------------------
  # Estimate the clonality of each insertion
  #-------------------------------------------
    if (lib_type == "shear"){
      print("Calculating Insertion Clonality...")
      insertions %<>% mutate(clonality=NA) %>% group_by(Sample_key) %>% group_modify(~ calc_clonality(.x)) %>% ungroup() 
    }
    
  #-------------------------------------------
  # Save clustered insertions
  #-------------------------------------------
    print("Saving the Data...")
    
    fwrite(insertions, file = glue("{output_dir}/all_insertions_clustered.csv.gz"), nThread = 8)
    # insertions <- fread(file = glue("{output_dir}/all_insertions_clustered.csv.gz"), nThread = 8)
    

# read metadata -----------------------------------------------------------
metadata <- fread(metadata_file, nThread = 8)
# manual modification here to the metadata so that all controls get included in the eventual cohort split
metadata %<>% mutate(tnp_end = case_when( # won't actually make a difference when creating blacklist, which will just merge these all again.
  grepl("L-P", Mouse_ID) ~ "PBL",
  grepl("R-P", Mouse_ID) ~ "PBR",
  TRUE ~ as.character(tnp_end) 
)) %>% as.data.frame()    

# split cohorts -----------------------------------------------------------
# here we are splitting into PBR:TAM+, PBR:TAM-, PBR:CTRL, and same three for PBL.
    
    print("Spliting insertions into groups using using the 'Exp_' columns...")
    #Create a new folder for each experimental group and cohort
    exp_groups = colnames(metadata)[grepl(colnames(metadata), pattern = "^Exp_")]
    
    for (name in exp_groups){
      
      print(paste("...",name, sep=""))
      for (cohort in (metadata[,name] %>% unique() %>% na.omit() %>% as.character())){
        
        for(leftright in (metadata$tnp_end %>% unique() %>% na.omit() %>% as.character())){
        print(paste("......",cohort, sep=""))
        cohort_dir <- glue("{output_dir}/{name}/{cohort}/{leftright}")
        dir.create(cohort_dir, recursive = T)
        
        metadata_subset = metadata[metadata[,name] == cohort,] %>% filter(!is.na(Mouse_ID) & tnp_end == leftright) %>% mutate(full_name = paste(Mouse_ID , Sample_ID , sep="_"))
        
        cohort_samples = metadata_subset[metadata_subset[,name] == cohort,] %>% dplyr::select(full_name) %>% pull()
        
        #Subset the gCIS file and place into the folder
        if(grepl("TAM", cohort)){
          insertions %>% mutate(full_name = paste(patient, sample, sep="_")) %>% 
            dplyr::filter(full_name %in% cohort_samples & grepl(leftright, Sample_key)) %>%  #select the correct sample and correct orientation only
            dplyr::select(-full_name) %>% write.table(., file = glue("{cohort_dir}/clustered_insertions_{cohort}_{leftright}.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
        }else{
          insertions %>% mutate(full_name = paste(patient, sample, sep="_")) %>% 
            dplyr::filter(full_name %in% cohort_samples & ((grepl("Primary", patient) & (str_sub(orientation,3,3) == str_sub(leftright,3,3))) | orientation == leftright)) %>%  #get 1. Primary IRR and Primary JXR and PTCH PBR only, during the CTRL loop with leftright == PBR, for example.
            dplyr::select(-full_name) %>% write.table(., file = glue("{cohort_dir}/clustered_insertions_{cohort}_{leftright}.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
          
        } 
        
        } #end orientation loop: PBR and PBL
        
      }#end cohort loop: TAM+/- and CTRL
    }#end group loop: tam status
    
# question - will it make a difference to split cohorts here into tam+ and tam-? or should I do all tam, then gcis, then retro-attribute to either tam+/-
    # well, thresholds will be calculated on a per-sample basis, so splitting there won't make a difference. chi-squared test in gcis is agnostic to total number/count
    # of insertions
# question - should it have made a difference at all by splitting PBR and PBL insertions? I guess yes. for example if PBL has fewer insertions/less read depth or whatever
    # if PBL is biased AGAINST. then combining all insertions leaves those PBL insertions (that for some reason don't have a matched PBR) in a pool of all insertions
    # that is predominated by PBR insertions, and for which the thresholds are determined heavily by PBR read counts, etc. so yes, this is a good approach.

    
# blacklist creation ------------------------------------------------------
tam_samples <- metadata %>% 
      filter(grepl("TAM", Exp_group) & !duplicated(Mouse_ID)) %>% 
      select(Mouse_ID, Sample_ID)
    
    pre_blacklist <- map(list.files(path = glue("{output_dir}/Exp_group/CTRL"), full.names = T, recursive = T),
                         \(x){ read_tsv(x) }) %>% rbindlist() %>% 
      mutate(insertion = glue("{chr}_{loc}_{strand}"),
             Sample_ID = "PT") %>% 
      dplyr::select(Sample_ID, insertion)
    
    blacklist_insertion <- inner_join(tam_samples, pre_blacklist, by = "Sample_ID") %>% 
      dplyr::rename(patient = Mouse_ID, sample = Sample_ID, insertion_key = insertion) %>% 
      write_tsv(glue("{output_dir}/blacklist.tsv"))
    blacklist_insertion_loc <- glue("{output_dir}/blacklist.tsv")
    

# gcis thresholds ---------------------------------------------------------
# quick stopcheck
if (!(lib_type %in% c("shear", "restriction"))){print ("Error: please select library type: \"shear\" or \"restriction\""); exit()}

#Load the insertions, start loop
rm(insertions)
tam_folders <- glue("{output_dir}/Exp_group") %>% list.files(pattern = "TAM", full.names = T)
insertion_files <- map(tam_folders, \(x){ list.files(path = x, full.names = T, recursive = T, pattern = "clustered") }) %>% unlist()

for(y in insertion_files){
  
  # organize working dir
  current_dir <- dirname(y)
  setwd(current_dir)
  
  # read in
  insertions=read.table(y, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  if (lib_type == "restriction"){
    print("Assuming it is a RESTRICTION type library preparation...")
    
    #------------------------------------
    # Top 1 percent read count filter
    #------------------------------------
    print("Top 1 percent read count filter...")
    insertions %<>% group_by(patient, sample) %>% group_modify(~ top_1per_max_insertion(.x)) %>% ungroup()
    
    #------------------------------------
    # Top 0.1 percent read sum
    #------------------------------------
    print("Top 0.1 percent read sum filter...")
    insertions %<>% group_by(patient, sample) %>% group_modify(~ sum_reads_0.1_per(.x)) %>% ungroup()
    
    #------------------------------------
    # 95% percentile of negative binomial
    #------------------------------------
    print("95% percentile of negative binomial filter......")
    insertions %<>% group_by(patient, sample) %>% group_modify(~ neg_binom_thresh(.x)) %>% ungroup()
    
    #------------------------------------
    # Save thresholds for later inspection/curiosity
    #------------------------------------
    print("Saving threshold info......")
    insertions %>% dplyr::select(Sample_key, top_1per_max_insertion, sum_reads_0.1_per, neg_binom_thresh) %>% unique() %>%
      fwrite(file = glue("{current_dir}/thresholds.csv.gz"), nThread = 8)
    
    #------------------------------------
    # Pick 2nd most stringent of the flags for the clonality filter #EDIT 2021-11-10: using "subclonal" filter i.e. 2nd most stringent threshold per insertion
    #------------------------------------
    print("Picking 2nd most stringent filter......")
    insertions %<>% rowwise() %>% mutate(clonality_threshold = max(top_1per_max_insertion, sum_reads_0.1_per, neg_binom_thresh),
                                         subclonality_threshold = median(c(top_1per_max_insertion, sum_reads_0.1_per, neg_binom_thresh)),
                                         lowest_threshold = min(top_1per_max_insertion, sum_reads_0.1_per, neg_binom_thresh)) %>%
      dplyr::filter(count >= subclonality_threshold) %>% #select stringency here!!!
      select(-top_1per_max_insertion, -sum_reads_0.1_per, -neg_binom_thresh, -clonality_threshold, -subclonality_threshold, -lowest_threshold)
    
    #------------------------------------
    # Filter insertions with a single read support
    #------------------------------------
    insertions = insertions[insertions$count > 1,]
    
  } else if (lib_type == "shear"){
    
    print("Assuming it is a SHEAR type library preparation...")
    #------------------------------------
    # ULP Flag
    #------------------------------------
    #Filter out insertions with a single ULP support
    print("Filtering out ULP = 1 insertions...")
    insertions = insertions[insertions$ULP > 1,]
    
    #------------------------------------
    # Clonality Filter
    #------------------------------------  
    #  Filter out non clonal insertions using the provided clonality threshold
    print(paste("Filtering out insertions with less then ", clonality_threshold, " clonality...", sep=""))
    insertions = insertions[insertions$clonality > clonality_threshold,]
  }
  
  #------------------------------------
  # Save Brett-filtered insertions
  #------------------------------------
  print("Saving Brett-filtered insertions......")
  insertions %>% fwrite(file = glue("{current_dir}/brett_filtered_insertions_{basename(y) %>% str_replace('.tsv','') %>% str_replace('clustered_insertions_','')}_subclonal.csv.gz"), nThread = 8)
  
} #end y loop for separate PBR/PBL thresholds



# merge IRR and IRL,  proceed to gCIS -------------------------------------

for(y in c("TAM+", "TAM-", "all_tam")){
  
  #Load the gene annotation and probabilities
  gene_annot=read.table(gene_annot_input, stringsAsFactors = FALSE, header = TRUE)
  
  print("Merging IRL and IRR Insertions...")
  current_dir <- glue("{output_dir}/Exp_group/{y}")
  dir.create(current_dir, showWarnings = F)
  setwd(current_dir)
  if(y == "all_tam"){
    brett_paths <- list.files("..", pattern = "brett", full.names = T, recursive = T) %>% grep(pattern = "\\+|-", ., value = T)
  }else{
    brett_paths <- list.files(path = current_dir, recursive = T, pattern = "brett", full.names = T) 
  }
   
  insertions <- map(brett_paths, \(x){ fread(x, nThread = 8) }) %>% rbindlist()
  
  # merge left and right
  insertions %<>%
    group_by(patient, sample, chr, loc, strand) %>% slice_max(count, n = 1) %>%
    ungroup() %>%
    select(-orientation, -Sample_key, -cluster)
  
  #------------------------------------
  # Standard Chromosome Flag
  #------------------------------------
  #Filter out insertions in non-standard chr
  print("Removing insertions in non-standared chromsomes...")
  standard_chr=c(paste("chr",1:22, sep=""), "chrX", "chrY")
  
  insertions = insertions[insertions$chr %in% standard_chr,]
  
  #------------------------------------
  # Donor Filter
  #------------------------------------ 
  print("Removing insertions in donor chromosomes...")
  insertions %<>% dplyr::left_join(unique(metadata[,c("Mouse_ID", "Donor")]), by = c("patient" = "Mouse_ID")) %>% rowwise() %>% dplyr::filter(!grepl(paste(chr,"$", sep=""), Donor)) %>% dplyr::filter(!grepl(paste(chr,",", sep=""), Donor))
  
  #------------------------------------
  # Recurrence Filter
  #------------------------------------
  if (recurrence_filter == TRUE){
    print("Removing insertions found in multiple mice...")
    #Blacklist of recurrent insertions (recurrent insertions within the same mouse a preserved)
    keys=insertions %>% mutate(insertion_key=paste(chr,loc,strand, sep="_")) %>% dplyr::select(patient, insertion_key) %>% unique() %>% dplyr::select(insertion_key) %>% pull() 
    reccurent_insertions=unique(keys[duplicated(keys)])
    
    #Salvage Insertions
    salvage <-function(tbl, type){
      
      if (type == "shear"){
        keep_patient = tbl %>% arrange(desc(ULP), desc(count)) %>%  head(., n=1) %>% select(patient) %>% pull()
        return(tbl %>% dplyr::filter(patient %in% keep_patient))
      } else if (type == "restriction"){
        keep_patient = tbl %>% arrange(desc(count)) %>% head(., n=1) %>% select(patient) %>% pull()
        return(tbl %>% dplyr::filter(patient %in% keep_patient))
        
      }
    }
    
    salvaged_insertions = insertions %>% mutate(insertion_key=paste(chr,loc,strand, sep="_")) %>% rowwise() %>% 
      dplyr::filter(insertion_key %in% reccurent_insertions) %>% group_by(insertion_key) %>% 
      group_modify(~ salvage(.x, lib_type)) %>% ungroup() %>% dplyr::select(-insertion_key)
    
    
    #Remove blacklisted insertions
    insertions %<>% mutate(insertion_key=paste(chr,loc,strand, sep="_")) %>% dplyr::filter(!(insertion_key %in% reccurent_insertions)) %>% 
      dplyr::select(-insertion_key)
    
    #Add back salvaged insertions 
    print(paste("Salvaged ", dim(salvaged_insertions)[1], " insertions...", sep=""))
    insertions = unique(bind_rows(insertions, salvaged_insertions))
  }
  
  #------------------------------------
  # Remove supplied blacklisted insertions i.e. insertions in CTRL mice
  #------------------------------------
  if (!is.null(blacklist_insertion_loc)){
    
    print("Removing insertions found in the supplied blacklist...")
    #Extra step for some customized filtering. These could be the insertions found in controls. Must specify the mouse name and insertion key to remove. 
    blacklist_keys = blacklist_insertion %>% mutate(full_key = paste(patient, sample, insertion_key, sep="_")) %>% dplyr::select(full_key) %>% pull()
    
    insertion_keys = insertions %>% mutate(full_key = paste(patient, sample, chr, loc, strand, sep="_")) %>% dplyr::select(full_key) %>% pull()
    
    insertions = insertions[!(insertion_keys %in% blacklist_keys),]
  }
  #------------------------------------
  # Parse the insertions for gCIS
  #------------------------------------
  print("Parsing insertions for output...")
  if (lib_type == "shear"){
    insertions %<>% mutate(name=paste(patient, sample, sep="_")) %>% dplyr::select(name, patient, sample, gene, gene_location, functional_effect, chr, loc, count, gene_orientation, strand, ULP, clonality)  
  } else if (lib_type == "restriction"){
    insertions %<>% mutate(name=paste(patient, sample, sep="_")) %>% dplyr::select(name, patient, sample, gene, gene_location, functional_effect, chr, loc, count, gene_orientation, strand)  
  }
  
  print("Saving pre-gCIS insertions......")
  insertions %>% fwrite(file = glue("{current_dir}/pre_gcis_insertions_{y}_subclonal.csv.gz"), nThread = 8)
  
  #------------------------------------
  # Count Insertions in Each Gene
  #------------------------------------
  print("Intersecting insertions with the gene bounds...")
  chr_split=list()
  for (chr in unique(insertions$chr)){
    print(paste("...",chr, sep=""))
    
    #Subset the search space by chromosome
    gene_annot_subset = gene_annot[gene_annot$chr == chr,]
    insertions_subset = insertions[insertions$chr == chr,]
    
    #Find Insertions that intersect each gene bound
    gene_annot_subset$insertions = mapply(function(gene, insertions, g_chr, g_start, g_end)
    {insertions %>% filter(loc >= g_start & loc <= g_end) %>% mutate(gene=gene) %>% as.data.frame()}, 
    gene=gene_annot_subset$gene, insertions=list(insertions_subset), 
    g_chr=gene_annot_subset$chr, g_start=gene_annot_subset$start_buff, g_end=gene_annot_subset$end_buff, SIMPLIFY = FALSE)
    
    chr_split[[chr]]=gene_annot_subset
    
  }
  
  #Merge all the chromosomes back together
  gene_annot = as.data.frame(rbindlist(chr_split))
  
  #-----------------------------------------
  # String of All Samples with Insertion
  #-----------------------------------------
  
  gene_annot$observed_samples = mapply(function(x){paste(x$name, collapse =";")}, x = gene_annot$insertions)
  
  #-----------------------------------------
  # Count Observed Insertions in Each Gene
  #-----------------------------------------
  print("Counting observed insertions in each gene...")  
  #Observed samples count
  gene_annot$observed_sample_count = mapply(function(x){length(unique(x$name))}, x = gene_annot$insertions)
  
  #-----------------------------------------
  # Count Number of Unique Patients
  #----------------------------------------- 
  print("Counting number of unique patients...")   
  #Observed samples count
  gene_annot$observed_patient_count = mapply(function(x){length(unique(x$patient))}, x = gene_annot$insertions)
  
  #---------------------------------------------
  # Calculate Expected Insertions in Each Gene
  #---------------------------------------------
  print("Calculating expected number of insertions in each gene...")    
  #Expected samples count
  gene_annot$expected_sample_count=gene_annot$TA_insertion_prob*dim(unique(insertions))[1]
  
  #----------------------------------------------------
  # Calculate Chi-Squared Test Statistic and Pvalue
  #----------------------------------------------------
  print("Calculating Chi-Squared test statistic...")   
  #Chi-Squared Test statistic
  n=length(unique(insertions$name))
  gene_annot = gene_annot %>% mutate(chi_test_stat = ((observed_sample_count - expected_sample_count)^2 / (expected_sample_count*(1-(expected_sample_count/n)))))
  
  #Chi-Squared P-value
  gene_annot$Pvalue=pchisq(	gene_annot$chi_test_stat, 1, lower.tail = FALSE)
  
  #-------------------------------------------------
  # Correct the P-value
  #------------------------------------------------- 
  print(paste("Correcting the pvalue using ", pvalue_correction, sep="")) 
  #Corrected P-value
  gene_annot$Pvalue_adj=p.adjust(gene_annot$Pvalue, method = pvalue_correction)
  
  #-------------------------------------------------
  # Filter to Find Significant Genes
  #-------------------------------------------------  
  print(paste("Filtering the genes using qvalue threshold: ", qvalue_threshold, " and observed sample count min: ",min_rec,sep=""))  
  #Filter the significant genes
  sig_genes = gene_annot %>% filter(Pvalue_adj < qvalue_threshold & observed_patient_count >= min_rec) %>% arrange(Pvalue_adj)
  
  saveRDS(sig_genes, file = glue("{current_dir}/sig_genes_{y}.rds"))
  print(paste("Detected  ", dim(sig_genes)[1], " significant genes", sep="")) 
  
  
  ## =============================================
  ##    OUTPUT PARSING
  ## =============================================
  print("Outputing significant genes and insertions...") 
  
  #-------------------------------------------------
  # Output all the Genes
  #-------------------------------------------------   
  
  gene_annot %>% dplyr::select(gene, chr, start, end, strand, observed_sample_count, observed_patient_count, Pvalue, Pvalue_adj,observed_samples) %>%
    arrange(Pvalue_adj) %>%
    write_tsv(glue("{current_dir}/all_gcis_{y}.tsv"))
  
  #-------------------------------------------------
  # Output the Significant Genes
  #-------------------------------------------------  
  
  sig_genes %>% dplyr::select(gene, chr, start, end, strand, observed_sample_count, observed_patient_count, Pvalue, Pvalue_adj,observed_samples) %>%
    write_tsv(glue("{current_dir}/sig_gcis_{y}.tsv"))
  
  #-------------------------------------------------
  # Output the Significant Insertions
  #-------------------------------------------------   
  
  #Ensure the insertions within gene window (+/ set buffer) have a gene name
  sig_genes$insertions=mapply(function(insertions, gene){insertions$gene = gene; return(insertions)}, sig_genes$insertions, sig_genes$gene, SIMPLIFY = FALSE)
  
  if (lib_type == "shear"){
    sig_insertions = as.data.frame(rbindlist(sig_genes$insertions)) %>% dplyr::select(name, gene, gene_location, functional_effect, gene_orientation, chr, loc, strand, count, ULP,clonality) 
  } else if (lib_type == "restriction"){
    sig_insertions = as.data.frame(rbindlist(sig_genes$insertions)) %>% dplyr::select(name, gene, gene_location, functional_effect, gene_orientation, chr, loc, strand, count) 
  }
  
  write_tsv(sig_insertions, glue("{current_dir}/sig_insertions_{y}.tsv"))
  
  #-------------------------------------------------
  # Output all Insertions, stratified by tam status
  #-------------------------------------------------
  # tamp_samples <- metadata %>% filter(Exp_group == "TAM+") %>% pull(Mouse_ID) %>% unique() %>% paste0(.,"_PT")
  # tamn_samples <- metadata %>% filter(Exp_group == "TAM-") %>% pull(Mouse_ID) %>% unique() %>% paste0(.,"_PT")
  # 
  # 
  # gene_annot %>% filter(grepl(paste0(tamp_samples, collapse = "|"), observed_samples)) %>%
  #   dplyr::select(-insertions) %>%
  #   fwrite(glue("{current_dir}/gcis_tam_neg.csv"), nThread = 8)
  # gene_annot %>% filter(grepl(paste0(tamn_samples, collapse = "|"), observed_samples)) %>%
  #   dplyr::select(-insertions) %>% 
  #   fwrite(glue("{current_dir}/gcis_tam_pos.csv"), nThread = 8)
  
  
  #-------------------------------------------------
  # Output the Bed Files (View in IGV)
  #-------------------------------------------------   
  print("Outputing bedtracks to view in IGV (red = + strand, blue = - strand) ...")  
  convert_to_bed <- function(insertions, type)
  {
    V1=insertions$chr
    V2=insertions$loc
    V3=insertions$loc
    
    if (type == "shear"){
      V4=paste(insertions$name, insertions$gene, insertions$strand, round(insertions$clonality,3), sep="_")
    } else if (type == "restriction"){
      V4=paste(insertions$name, insertions$gene, insertions$strand, insertions$count, sep="_")
    }
    
    V5="."
    V6="."
    V7="."
    V8=ifelse(insertions$strand  == "SAME", "255,0,0", "0,0,255")
    return(data.frame(V1,V2,V3,V4,V5,V6,V7,V8))
  }
  
  write.table(convert_to_bed(sig_insertions, lib_type), file = glue("{current_dir}/sig_insertions_{y}.bed"), quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)    
  write.table(convert_to_bed(insertions, lib_type), file = glue("{current_dir}/all_insertions_{y}.bed"), quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)    
  
  
}



# oncoprint ---------------------------------------------------------------
# identify false positives
`%notin%` <- Negate(`%in%`)
false_pos <- c("Sfi1", "En2", "Foxf2")

# import and remove false pos
# tamp_sig_gcis <- read_tsv(list.files(path = output_dir, recursive = T, full.names = T, pattern = "sig_gcis_TAM\\+")) %>% #be careful with regex
#   filter(gene %notin% false_pos)
# tamn_sig_gcis <- read_tsv(list.files(path = output_dir, recursive = T, full.names = T, pattern = "sig_gcis_TAM-")) %>% 
#   filter(gene %notin% false_pos)
sig_gcis <- read_tsv(list.files(path = output_dir, recursive = T, full.names = T, pattern = "sig_gcis_all_tam")) %>% #be careful with regex
  filter(gene %notin% false_pos)

# get dims
genes <- unique(sig_gcis$gene)
row_num <- length(genes)
tamp_samples <- metadata %>% filter(Exp_group == "TAM+") %>% pull(Mouse_ID) %>% unique() %>% paste0(.,"_PT")
tamn_samples <- metadata %>% filter(Exp_group == "TAM-") %>% pull(Mouse_ID) %>% unique() %>% paste0(.,"_PT")
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
col = c(tamp = "#E773A0", tamn = "#526A8F") #("#f56262","#1ba2ce")
op <- oncoPrint(mtx_list,
          # graphics
          alter_fun = list(
            background = alter_graphic("rect", fill = "#CCCCCC"),
            tamp = alter_graphic("rect", fill = col["tamp"]),
            tamn = alter_graphic("rect", fill = col["tamn"])),
          # colors
          col = col,
          # flip gene names / percentages
          pct_side = "right", row_names_side = "left",
          # remove barplot
          right_annotation = NULL,
          # adjust rowname / percent font sizes
          row_names_gp = gpar(fontsize = 8),
          pct_gp = gpar(fontsize = 8),
          # relabel legend
          heatmap_legend_param = list(title = NULL, labels = c("TAM+", "TAM-"))
            )
plot_height <- 6
plot_width <- 5
pdf(glue("/Users/Anders/OneDrive - University of Toronto/Lazy Piggy KCNB2/Bioinformatics/Anders_MDT_20211214/gcis_oncoprint.pdf"), h = plot_height, w = plot_width)
draw(op)
dev.off()
pdf(glue("{output_dir}/gcis_oncoprint.pdf"), h = plot_height, w = plot_width)
draw(op)
dev.off()



