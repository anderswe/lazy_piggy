# Info --------------------------------------------------------------------

# Nov 10, 2023
# Anders E.
# Taylor lab

# this script writes (and submits)
# jobs 

# Setup -------------------------------------------------------------------
# module load R/4.2.1
# R
library(glue)
library(magrittr)
source("src/scripts/configs.R")

# samples
parent <- xh_dir
samples <- glue::glue("{xh_dir}/metadata.csv") %>% 
    data.table::fread() %>% 
    dplyr::pull(sample) %>% 
    stringr::str_replace("Mspos", "MS_pos")

# directory mgmt
scripts_dir <- glue::glue("{xh_dir}/scripts")
logs_dir <- glue::glue("{xh_dir}/logs")
dir.create(scripts_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(logs_dir, showWarnings = FALSE, recursive = TRUE)


# Create sh scripts ------------------------------------------------------
for(sample in samples){

    # get R1 and R2
    r1_dir <- list.files(xh_dir, full.names = TRUE) %>% grep(sample, ., value = T) %>% grep("md5|fastqc", ., value = T, invert = T) %>% grep("R1", ., value = T)
    r2_dir <- list.files(xh_dir, full.names = TRUE) %>% grep(sample, ., value = T) %>% grep("md5|fastqc", ., value = T, invert = T) %>% grep("R2", ., value = T)
    

    # Generate the shell script content using glue
    shell_script_content <- glue::glue('
    #!/bin/bash
    #SBATCH -J {sample}
    #SBATCH --output {logs_dir}/star_align_{sample}_outs.out
    #SBATCH -e {logs_dir}/star_align_{sample}_errors.out
    #SBATCH --mem 118G
    #SBATCH --time 2-00:00:00

    # check
    echo "Processing sample: {sample}"
    cd {parent};
    echo $PWD;

    # prep env and dirs
    module load star/2.5.4b;
    mkdir -p $PWD/aligned;

    # get reads
    read1={r1_dir};
    read2={r2_dir};

    # run
    STAR --runMode alignReads \\\\\n
    --runThreadN 4 \\\\\n
    --genomeDir {mdt_dir}/Genome_file/STAR_GRCm38/ \\\\\n
    --outSAMtype BAM SortedByCoordinate \\\\\n
    --quantMode TranscriptomeSAM GeneCounts \\\\\n
    --readFilesCommand zcat \\\\\n
    --readFilesIn $read1 $read2 \\\\\n
    --outFileNamePrefix $PWD/aligned/{sample} \\\\\n
    --outFilterType BySJout \\\\\n
    --outFilterMultimapNmax 20 \\\\\n
    --alignSJoverhangMin 8 \\\\\n
    --alignMatesGapMax 200000 \\\\\n
    --alignIntronMax 200000 \\\\\n
    --chimSegmentReadGapMax parameter 3 \\\\\n
    --alignSJDBoverhangMin 10 \\\\\n
    --alignSJstitchMismatchNmax 5 -1 5 5 \\\\\n
    --outSAMmultNmax 20 \\\\\n
    --chimSegmentMin 12 \\\\\n
    --chimJunctionOverhangMin 12 \\\\\n
    --twopassMode Basic \\\\\n
    --outSAMattrRGline ID:{sample} PU:Illumina SM:{sample}
    ')

    # Write the shell script to a file
    shell_script_filename <- glue::glue("{scripts_dir}/star_align_{sample}.sh")
    writeLines(shell_script_content, shell_script_filename)

}


# Submit shell scripts ---------------------------------------------------

for(sample in samples){ system(glue::glue("sbatch {scripts_dir}/star_align_{sample}.sh")) }

system("squeue --me")


# Concatenate counts -----------------------------------------------------

# read in all unstranded counts, output a cbound dataframe for all samples
counts <- purrr::map(list.files(pattern = "ReadsPerGene", path = glue::glue("{xh_dir}/aligned"), full.names = T), function(x){
  
  lab <- basename(x) %>% stringr::str_replace("ReadsPerGene.out.tab", "")
  
  readr::read_tsv(x, col_names = F, skip = 4) %>% 
    dplyr::select(gene = X1, !!rlang::quo_name(lab) := X2) %>% 
    return()
  
}) %>%
  purrr::reduce(dplyr::left_join, by = "gene") %>% 
  tibble::column_to_rownames("gene")

# subset out Piezo2 KO samples
to_keep <- glue::glue("{xh_dir}/metadata.csv") %>% 
    data.table::fread() %>% 
    dplyr::mutate(sample = sample %>% stringr::str_replace("Mspos", "MS_pos")) %>% 
    dplyr::filter(condition != "Piezo2_ko") %>% 
    dplyr::pull(sample)

final_counts <- counts[,to_keep]
saveRDS(final_counts, glue::glue("{xh_dir}/counts.rds"))

