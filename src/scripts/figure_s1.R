# Info --------------------------------------------------------------------

# Figure S1 reproducibility code
# Feb 1, 2023
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
source("src/scripts/utils.R")
source("src/scripts/configs.R")

# outdir
outdir <- glue::glue("{repo_dir}outs")



# Panel D -----------------------------------------------------------------
# qPCR copy number determination for lazy piggy founder mouselines



# import and clean
path <- glue::glue("{proj_dir}/from_kevin/LP_copy_number/Copy number determination Dec 19, 2013.xls") #yeah, 2013
df <- read_excel(path, skip = 7) %>% 
  janitor::clean_names() %>% 
  dplyr::select(-average_ct, -x6, -x7, -copy, -ct_2, -x10, -x11, -x12) %>% 
  tidyr::unite("sample_key",
               sample_name, target_name,
               remove = FALSE) %>% 
  .[1:50,] %>% 
  mutate(label = case_when(
    grepl("h2", sample_name) ~ "water",
    grepl("48|54|56|129|137|139|143", sample_name) ~ "founder",
    TRUE ~ "standard"),
    input_copies = ifelse(label == "standard",
                          as.numeric(sample_name),
                          NA))

# fit lm
std_fit <- lm(ct ~ input_copies, data = df) %>% broom::tidy()
model_label <- glue::glue("y = {std_fit$estimate[2] %>% signif(3)}x + {std_fit$estimate[1] %>% signif(3)}")

df %<>% mutate(predicted_copies = ifelse(label == "founder",
                                         (ct - std_fit$estimate[1]) / std_fit$estimate[2],
                                         NA))

# plot
pdf(glue("{outdir}/lp_founder_copy_number_determination.pdf"), h = 3, w = 4)
ggplot(filter(df, label == "standard"), aes(input_copies, ct))+
  geom_smooth(method='lm', alpha = 0.3, color = "grey")+
  geom_jitter(color="black", size=2, alpha=0.5, width = 0.5)+
  geom_jitter(data = filter(df, label == "founder"), aes(predicted_copies, ct, fill = sample_name), pch=21, size = 2)+
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(9,'Spectral'))(7),
                    name = "Founder")+
  theme_classic()+
  ylab(expression(C[T]))+
  xlab("Copy Number (Standard, Black; Founder, Colour)")+
  annotate(geom = "text", x = 0, y = (df %>% filter(label == "standard") %>% pull(ct) %>% min()), label = model_label, hjust = 0, vjust = 0, size = 4)
dev.off()          






# Panel H -----------------------------------------------------------------
# lp survival by tam status


# Import
df <- read_tsv(glue::glue("{proj_dir}/Lazy_Piggy_Analysis/SPLINK/LP_sample_metadata.txt")) %>% clean_names() %>% 
  filter(experiment_group != "CTRL")
df$experiment_group %<>% factor(levels = c("TAM+", "TAM-"))

# Palette
mypal <- pal[c("red", "blue")] %>% as.character()

# Survival object
surv <- Surv(time = df$total_surival, event = rep(1, length(df$total_surival)))
fit <- survfit(surv ~ experiment_group, data = df)

# Prepare KM plot
caption <- glue("p = {surv_pvalue(fit, df)$pval %>% signif(3)}")
gg <- ggsurvplot(fit, data = df, conf.int = F, xlab = "Time (Days)", pval = FALSE,
                 xlim=c(0,200), risk.table = T, #max(df$total_surival)+5
                 legend.title = " ", 
                 # legend = "right",
                 tables.theme = theme_cleantable(), # I really liked the look of the clean table theme, but it was unclear which timepoints the at-risk table values were associated with.
                 break.time.by = 50,
                 risk.table.fontsize = 4,
                 risk.table.height = 0.20, 
                 # risk.table.title = "test",
                 # tables.col = "strata",
                 # risk.table.y.text = TRUE,
                 # risk.table.y.text.col = FALSE,
                 # risk.table.pos = "in",
                 surv.median.line = c("h"),
                 legend.labs = c("Tamoxifen +", "Tamoxifen -"),
                 censor.shape = c("|"), censor.size = 6, 
                 axes.offset = T, tables.y.text = FALSE,
                 palette = mypal)+
  guides(colour = guide_legend(nrow = 2))
gg$plot <- gg$plot+
  annotate(geom = "text", x = 0, y = 0, label = caption, hjust = 0, vjust = 0, size = 4)+
  theme(legend.justification=c(0.4, 0.1), legend.position = c(0.75,0.75)) 
gg$table <- gg$table + 
  theme(plot.title = element_text(hjust = 0, size = 12))
gg

# Export plot
pdf(glue::glue("{outdir}/survival_lp_by_tam_status.pdf"),height = 4, width = 4)
print(gg, newpage = FALSE)
dev.off()


# Survival by tam status, subgrouped by donor chromosome
for(i in unique(df$donor)){

  df_chr <- df %>% dplyr::filter(donor == i)
  surv <- Surv(time = df_chr$total_surival, event = rep(1, length(df_chr$total_surival)))
  fit <- survfit(surv ~ experiment_group, data = df_chr)
  copies <- ifelse(i == "7", "1000", "600")
  caption <- glue::glue("Donor chromosome {i}\n({copies}+ LP copies)\np = {surv_pvalue(fit, df_chr)$pval %>% signif(3)}")
  
  gg <- ggsurvplot(fit, data = df_chr, conf.int = F, xlab = "Time (Days)", pval = FALSE,
                   xlim=c(0,200), risk.table = T, #max(df$total_surival)+5
                   legend.title = " ", 
                   # legend = "right",
                   tables.theme = theme_cleantable(), # I really liked the look of the clean table theme, but it was unclear which timepoints the at-risk table values were associated with.
                   break.time.by = 50,
                   risk.table.fontsize = 4,
                   risk.table.height = 0.20, 
                   # risk.table.title = "test",
                   # tables.col = "strata",
                   # risk.table.y.text = TRUE,
                   # risk.table.y.text.col = FALSE,
                   # risk.table.pos = "in",
                   surv.median.line = c("h"),
                   legend.labs = c("Tamoxifen +", "Tamoxifen -"),
                   censor.shape = c("|"), censor.size = 6, 
                   axes.offset = T, tables.y.text = FALSE,
                   palette = mypal)+
    guides(colour = guide_legend(nrow = 2))
  gg$plot <- gg$plot+
    annotate(geom = "text", x = 0, y = 0, label = caption, hjust = 0, vjust = 0, size = 4)+
    theme(legend.justification=c(0.4, 0.1), legend.position = c(0.75,0.75)) 
  gg$table <- gg$table + 
    theme(plot.title = element_text(hjust = 0, size = 12))
  gg
  
  pdf(glue::glue("{outdir}/survival_lp_by_tam_status_donor_chr{i}.pdf"),height = 4, width = 4.5)
  print(gg, newpage = FALSE)
  dev.off()
  
  
}



# survival by donor chromosome
# Survival object
surv <- Surv(time = df$total_surival, event = rep(1, length(df$total_surival)))
fit <- survfit(surv ~ donor, data = df)

# Prepare KM plot
caption <- glue("Log-rank\np = {surv_pvalue(fit, df)$pval %>% signif(3)}")
gg <- ggsurvplot(fit, data = df, conf.int = F, xlab = "Time (Days)", pval = FALSE,
                 xlim=c(0,200), risk.table = T, #max(df$total_surival)+5
                 legend.title = " ", 
                 # legend = "right",
                 tables.theme = theme_cleantable(), # I really liked the look of the clean table theme, but it was unclear which timepoints the at-risk table values were associated with.
                 break.time.by = 50,
                 risk.table.fontsize = 4,
                 risk.table.height = 0.20, 
                 # risk.table.title = "test",
                 # tables.col = "strata",
                 # risk.table.y.text = TRUE,
                 # risk.table.y.text.col = FALSE,
                 # risk.table.pos = "in",
                 surv.median.line = c("h"),
                 legend.labs = c("Chr10 Donor", "Chr7 Donor"),
                 censor.shape = c("|"), censor.size = 6, 
                 axes.offset = T, tables.y.text = FALSE,
                 palette = mypal)+
  guides(colour = guide_legend(nrow = 2))
gg$plot <- gg$plot+
  annotate(geom = "text", x = 0, y = 0, label = caption, hjust = 0, vjust = 0, size = 4)+
  theme(legend.justification=c(0.4, 0.1), legend.position = c(0.75,0.75)) 
gg$table <- gg$table + 
  theme(plot.title = element_text(hjust = 0, size = 12))
gg

# Export plot
pdf(glue::glue("{outdir}/survival_lp_by_donor.pdf"),height = 4, width = 4)
print(gg, newpage = FALSE)
dev.off()

# survival by donor, subgrouped by tam status
for(i in unique(df$experiment_group)){
  
  df_chr <- df %>% dplyr::filter(experiment_group == i) %>% dplyr::mutate(donor = ifelse(donor == "7", "chr7 (1000+)",  "chr10 (600+)"))
  surv <- Surv(time = df_chr$total_surival, event = rep(1, length(df_chr$total_surival)))
  fit <- survfit(surv ~ donor, data = df_chr)
  caption <- glue::glue("p = {surv_pvalue(fit, df_chr)$pval %>% signif(3)}")
  
  gg <- ggsurvplot(fit, data = df_chr, conf.int = F, xlab = "Time (Days)", pval = FALSE,
                   xlim=c(0,200), risk.table = T, #max(df$total_surival)+5
                   legend.title = " ", 
                   # legend = "right",
                   tables.theme = theme_cleantable(), 
                   break.time.by = 50,
                   risk.table.fontsize = 4,
                   risk.table.height = 0.20, 
                   # risk.table.title = "test",
                   # tables.col = "strata",
                   # risk.table.y.text = TRUE,
                   # risk.table.y.text.col = FALSE,
                   # risk.table.pos = "in",
                   surv.median.line = c("h"),
                   # legend.labs = c("Tamoxifen +", "Tamoxifen -"),
                   censor.shape = c("|"), censor.size = 6, 
                   axes.offset = T, tables.y.text = FALSE,
                   palette = mypal)+
    guides(colour = guide_legend(nrow = 2))
  gg$plot <- gg$plot+
    annotate(geom = "text", x = 0, y = 0, label = caption, hjust = 0, vjust = 0, size = 4)+
    theme(legend.justification=c(0.4, 0.1), legend.position = c(0.75,0.75)) 
  gg$table <- gg$table + 
    theme(plot.title = element_text(hjust = 0, size = 12))
  gg
  
  pdf(glue::glue("{outdir}/survival_lp_by_donor_{i}.pdf"),height = 4, width = 4.5)
  print(gg, newpage = FALSE)
  dev.off()
  
  
}


# strat by donor AND tam status
dtdf <- df %>% dplyr::mutate(donor_tam = glue::glue("{donor}_{experiment_group}"))
surv <- Surv(time = dtdf$total_surival, event = rep(1, length(dtdf$total_surival)))
fit <- survfit(surv ~ donor_tam, data = dtdf)

# Prepare KM plot
caption <- glue("Log-rank\np = {surv_pvalue(fit, dtdf)$pval %>% signif(3)}")
gg <- ggsurvplot(fit, data = dtdf, conf.int = F, xlab = "Time (Days)", pval = FALSE,
                 xlim=c(0,200), risk.table = T, #max(df$total_surival)+5
                 legend.title = " ", 
                 # legend = "right",
                 tables.theme = theme_cleantable(), 
                 break.time.by = 50,
                 risk.table.fontsize = 4,
                 risk.table.height = 0.20, 
                 surv.median.line = c("h"),
                 censor.shape = c("|"), censor.size = 6, 
                 axes.offset = T, tables.y.text = FALSE)+
  guides(colour = guide_legend(nrow = 2))
gg$plot <- gg$plot+
  annotate(geom = "text", x = 0, y = 0, label = caption, hjust = 0, vjust = 0, size = 4)+
  theme(legend.justification=c(0.4, 0.1), legend.position = c(0.7,0.75)) 
gg$table <- gg$table + 
  theme(plot.title = element_text(hjust = 0, size = 12))
gg

# Export plot
pdf(glue::glue("{outdir}/survival_lp_by_donor_tam.pdf"),height = 4, width = 4)
print(gg, newpage = FALSE)
dev.off()


# bonus -------------------------------------------------------------------

# survival by: quad vs quint gt
xldf <- readxl::read_excel(glue::glue("{proj_dir}Sample_Information/Lazy Piggy Project Mouse Inventory V3.xlsx"), skip = 22) %>% 
  janitor::clean_names() %>% 
  dplyr::select(-control_exp, -x2) %>% 
  dplyr::filter(sac_n_y == "Y") %>% 
  dplyr::mutate(gt = tolower(genotype)) %>%
  dplyr::mutate(gt_simple = dplyr::case_when(
    nestin_cre == 1 & n_luc_sb100 == 1 & ptc == 1 & r26_pber == 1 & (lp_137_chr10 == 1 | lp_129_chr7 == 1) ~ "quint",
    nestin_cre == 1 & n_luc_sb100 == 1 & ptc == 1 & r26_pber == 1 ~ "quad",
    TRUE ~ "other"
  ))

# now clean up genotypes, etc etc lubridate to calc OS









