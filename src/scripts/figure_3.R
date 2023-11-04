# Info --------------------------------------------------------------------

# Figure 3 reproducibility code
# Nov 4, 2023
# Anders E.
# Taylor lab

# analyses originally performed on: 2021-10-20 and 2022-01-18
# cleaned and posted to github on 2023-11-04


# Libraries -------------------------------------------------------------------------
library(readxl)
library(tidyverse)
library(lubridate)
library(survival)
library(janitor)
library(glue)
library(survminer)
library(magrittr)
library(DESeq2)
library(ggrepel)
library(cowplot)
library(viridis)
library(tidytable)
library(HGNChelper)
library(ComplexHeatmap)

source("src/scripts/configs.R")
source("src/scripts/utils.R")


# Panel A -----------------------------------------------------------------
# Non-tumour OS by Kcnb2 genotype

path <- glue::glue("{surv_data_dir}/MSK nontumor combo survival.xlsx")

# Data import and cleaning
df <- read_excel(path, sheet = 1)
names(df) <- c("days_survival", "double_pos", "het", "double_neg")
df %<>% mutate(id = 1:length(days_survival)) %>% 
  pivot_longer(cols = grep("double|het", names(.), value = T), names_to = "genotype") %>% 
  filter(!is.na(value)) %>% 
  mutate(genotype = factor(genotype, levels = c("double_neg", "het", "double_pos"))) %>% 
  arrange(genotype)

# Survival object
surv <- Surv(time = df$days_survival, event = df$value)
fit <- survfit(surv ~ genotype, data = df)

# Prepare KM plot
# https://github.com/kassambara/survminer/issues/246
mypal <- c("#FB6467FF", "#B7E4F9FF", "#24325FFF")
caption <- glue("p = {surv_pvalue(fit, df)$pval %>% signif(3)}")
gg <- ggsurvplot(fit, data = df, conf.int = F, xlab = "Time (Days)", pval = FALSE,
                 xlim=c(0,300), risk.table = T,
                 legend.title = " ", 
                 tables.theme = theme_cleantable(), # I really liked the look of the clean table theme, but it was unclear which timepoints the at-risk table values were associated with.
                 # break.time.by = 10,
                 risk.table.fontsize = 4,
                 risk.table.height = 0.20, 
                 # risk.table.title = "test",
                 # tables.col = "strata",
                 # risk.table.y.text = TRUE,
                 # risk.table.y.text.col = FALSE,
                 # risk.table.pos = "in",
                 surv.median.line = c("h"),
                 # legend.labs = c(expression(Kcnb2^{"+/+"}), "Kcnb2+/-", "Kcnb2-/-"),
                 censor = F,
                 # censor.shape = c("|"), censor.size = 6,
                 axes.offset = T, tables.y.text = FALSE,
                 palette = mypal)+
  guides(colour = guide_legend(nrow = 3))
gg$plot <- gg$plot+
  # annotate(geom = "text", x = 0, y = 1.2, label = caption, hjust = 0, vjust = 0, size = 4)+
  labs(tag = caption)+
  theme(
    # legend.position = c(0.15,0.3),
        # legend.box.background = element_rect(colour = "black"),
    plot.tag.position = c(0.2, 0.9),
    plot.tag = element_text(hjust = 0, vjust = 0, size = 12),
        legend.title = element_blank(),
        legend.text.align = 0)+ #legend.justification=c(0,0),
  scale_color_manual(labels = c(expression(Kcnb2^{"-/-"}), expression(Kcnb2^{"+/-"}), expression(Kcnb2^{"+/+"})), values = mypal)
gg$table <- gg$table + 
  theme(plot.title = element_text(hjust = 0, size = 12))
gg

# Export plot
outdir <- glue("{repo_dir}/outs")
pdf(glue::glue("{outdir}/survival_non_tumor_kcnb2_20220118.pdf"), height = 4.5, width = 4.5)
print(gg, newpage = FALSE)
dev.off()

# Pairwise comparisons
pairwise_survdiff(Surv(days_survival, value) ~ genotype, data = df)



# Panel G -----------------------------------------------------------------
# Math1-Cre SmoM2 MB OS by Kcnb2 genotype

# Data import and cleaning
df <- read_excel(glue::glue("{surv_data_dir}/MSK survival.xlsx")) %>% 
  clean_names() %>% 
  mutate(id = 1:length(days_survival)) %>% 
  pivot_longer(cols = grep("msk", names(.), value = T), names_to = "genotype") %>% 
  filter(!is.na(value)) %>% 
  mutate(value = NULL,
         genotype = str_replace(genotype, "msk_3", "double_neg") %>% str_replace("msk_2", "het") %>% str_replace("msk", "double_pos") %>% factor(levels = c("double_pos", "het", "double_neg"))) %>% 
  arrange(genotype)
  
# Survival object
surv <- Surv(time = df$days_survival, event = rep(1, length(df$days_survival)))
fit <- survfit(surv ~ genotype, data = df)

# Prepare KM plot
# https://github.com/kassambara/survminer/issues/246
mypal <- pal_rickandmorty(palette = "schwifty")(12)[3:5]
caption <- glue("p = {surv_pvalue(fit, df)$pval %>% signif(3)}")
gg <- ggsurvplot(fit, data = df, conf.int = F, xlab = "Time (Days)", pval = FALSE,
           xlim=c(0,max(df$days_survival)+5), risk.table = T,
           legend.title = " ", 
           tables.theme = theme_cleantable(), # I really liked the look of the clean table theme, but it was unclear which timepoints the at-risk table values were associated with.
           break.time.by = 10,
           risk.table.fontsize = 4,
           risk.table.height = 0.20, 
           # risk.table.title = "test",
           # tables.col = "strata",
           # risk.table.y.text = TRUE,
           # risk.table.y.text.col = FALSE,
           # risk.table.pos = "in",
           surv.median.line = c("h"),
           legend.labs = c("Kcnb2+/+", "Kcnb2+/-", "Kcnb2-/-"),
           censor.shape = c("|"), censor.size = 6, 
           axes.offset = T, tables.y.text = FALSE,
           palette = mypal)
gg$plot <- gg$plot + annotate(geom = "text", x = 0, y = 0, label = caption, hjust = 0, vjust = 0, size = 4)
gg$table <- gg$table + 
  theme(plot.title = element_text(hjust = 0, size = 12))
gg

# Export plot
outdir <- glue("{repo_dir}/outs")
pdf(glue::glue("{outdir}/survival_msk_20211020.pdf"), height = 5, width = 5)
print(gg, newpage = FALSE)
dev.off()

# Pairwise comparisons, if needed
pairwise_survdiff(Surv(days_survival, rep(1, length(days_survival))) ~ genotype, data = df)
