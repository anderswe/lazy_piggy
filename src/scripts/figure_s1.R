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






# Panel G -----------------------------------------------------------------
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


