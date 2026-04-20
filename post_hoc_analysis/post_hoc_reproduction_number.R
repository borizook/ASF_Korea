library(tidyverse)
library(ggsci)
scen_num = ""
setwd("~/save/Mechanistic_modelling/Marryland/Model_comparison/Candidates_comparisons")
folder_scen = paste0("Model_B",scen_num)

source("/save_home/jslim/Mechanistic_modelling/Marryland/OOD_Senario_results/Model_3_1/code/function_create_wd.R")
create_wd(folder_scen)

setwd(scen_folder_wd)
source(("/save_home/jslim/Mechanistic_modelling/Marryland/code/Model_4_3/init_model_wb.R"))

if (!any(list.files(file.path(scen_folder_wd,"/init_data" )) %in% "initdata.Rdata")){
  
  init.model(end_week = end_week, 
             dt = 1/10, 
             n.breaks = 4, 
             cell.adj.order = "first")
  
  setwd(paste0(file.path(scen_folder_wd,"/init_data" )))
  save(list = ls(), file = "initdata.Rdata", compress = TRUE)
}


setwd(paste0(file.path(scen_folder_wd,"/init_data" )))

load("initdata.Rdata", .GlobalEnv)

library(bayestestR)


# priors for EasyABC ------------------------------------------------------
priors <- list(
  c("unif", 0, 1),       # transmission rate
  c("unif", 1/4.66, 1/4.66), # sigma - detection rate for wild boar cells (Iu_d to Id)
  c("unif", 1/5.88, 1/5.88), # gamma - sanitation rate for wild boar cells (Id to S)
  c("unif", 0.073, 0.073), # eta - recovery rate (Iu_u to S)
  c("unif", 0, 1),   # 9.  detectibility
  c("unif", 0, 1)    # 9.  relative susceptibility/infectibility of cell with zero forest density (should be between 0 - 1)
  
)


# essential_pars ---
#  sigma - detection rate for wild boar cells (Iu_d to Id)
#  gamma - sanitation rate for wild boar cells (Id to S)
#  eta - Recovery Rate from Iu to S

pars_estimate = c("transmission rate", 
                  "detectibility",
                  "rel_sus_infect")
pars = c(0.09, 1/4.66, 1/5.88, 0.073, 0.5, 1)
pars_d = matrix(ncol = length(pars), pars )
colnames(pars_d) = c("transmission rate", 
                     "detection rate",
                     "sanitation rate", 
                     "recovery rate",
                     "detectibility", 
                     "rel_sus_infect")


results_ABC <-  ABC_Lenormand <-readRDS(paste0(file.path(scen_folder_wd,"/Posterior distribution/results_ABC.rds" )))

post_ABC<- results_ABC$param


colnames(post_ABC)  = pars_estimate

priors[pars_estimate] -> priors_esti



recov_p = 1/pars_d[4]

sani_p = 1/pars_d[3]







detectability = post_ABC[,2]

# Probability of infection from forested to forested
prob_f2f <- (1 - exp(-post_ABC[,1]))

# Probability of infection from forested to forested
prob_f2uf <- (1 - exp(-(post_ABC[,1] * post_ABC[,3])))

# Probability of infection from forested to forested
prob_uf2uf <- (1 - exp(-(post_ABC[,1] * post_ABC[,3] *post_ABC[,3] )))

# Probability of infection from forested to forested
prob_uf2f <- (1 - exp(-(post_ABC[,1]  *post_ABC[,3] )))
# 



## reproduction number
f2f_detected <- (1 - exp(-post_ABC[,1]* sani_p)) * 6 
f2f_undetected <-(1 - exp(-post_ABC[,1]* recov_p)) * 6 
uf2f_detected <- (1 - exp(-post_ABC[,1]* sani_p* post_ABC[,2])) * 6 
uf2f_undetected <- (1 - exp(-post_ABC[,1]* recov_p* post_ABC[,2])) * 6 
f2uf_detected <- (1 - exp(-post_ABC[,1]* sani_p* post_ABC[,2])) * 6 
f2uf_undetected <- (1 - exp(-post_ABC[,1]* recov_p* post_ABC[,2])) * 6 
uf2uf_detected <- (1 - exp(-post_ABC[,1]* sani_p* post_ABC[,2]* post_ABC[,2])) * 6 
uf2uf_undetected <- (1 - exp(-post_ABC[,1]* recov_p* post_ABC[,2]* post_ABC[,2])) * 6 




tibble(f2f_detected, f2f_undetected, uf2f_detected, uf2f_undetected, f2uf_detected, f2uf_undetected, uf2uf_detected, uf2uf_undetected) %>%
  pivot_longer(names_to ="type", values_to = "prob", cols = everything()) %>%
  mutate(detect = case_when(
    grepl(pattern = "_detected", x =type) ~ "detected", 
    grepl(pattern = "_undetected", x =type) ~ "undetected", 
  ), 
  from_forest = case_when(
    grepl(pattern = "\\bf2", x = type) ~ "from_forest",
    grepl(pattern = "uf2", x = type) ~ "from_unforest",
  ), 
  to_forest = case_when(
    grepl(pattern = "2f", x = type) ~ "to_forest",
    grepl(pattern = "2uf", x = type) ~ "to_unforest",
    
    
  )) -> dtf_r

dtf_r %>%
  
  mutate(f_type = factor(type, 
                         levels = c("f2f_detected", "f2uf_detected","uf2f_detected", "uf2uf_detected", 
                                    "f2f_undetected", "f2uf_undetected","uf2f_undetected", "uf2uf_undetected" ))) -> dtf_r

dtf_r %>%
  mutate(f_type2 = factor(type, 
                          levels = c("f2f_detected","f2f_undetected", "f2uf_detected", "f2uf_undetected","uf2f_detected", "uf2f_undetected", "uf2uf_detected",  "uf2uf_undetected"
                          ))) -> dtf_r
library(bbplot)

dtf_r %>% group_by(type) %>%
  summarise(hdi_95 = hdi(prob))


tibble(f2f_detected, f2f_undetected, uf2f_detected, uf2f_undetected, f2uf_detected, f2uf_undetected, uf2uf_detected, uf2uf_undetected) %>%
  apply(MARGIN = 2, FUN = hdi)

dtf_r %>% ggplot() +
  geom_violin(aes(y = prob, x = detect ,group = detect, col = detect, fill = detect)) +
  geom_boxplot(aes(y = prob, x = detect,group = detect), width = 0.1) +
  scale_y_continuous(limits = c(0, 5)) +
  scale_color_manual(values = pal_simpsons("springfield")(16)[c(13,5)])  +
  scale_fill_manual(values = pal_simpsons("springfield")(16)[c(13,5)]) +
  bbplot::bbc_style()



dtf_r %>% ggplot() +
  geom_violin(aes(y = prob, x = c(from_forest  ) ,group = from_forest , col = from_forest , fill = from_forest )) +
  geom_boxplot(aes(y = prob, x = c(from_forest  ),group = from_forest ), width = 0.1) +
  scale_y_continuous(limits = c(0, 5)) +
  scale_color_manual(values = pal_simpsons("springfield")(16)[c(13,5)])  +
  scale_fill_manual(values = pal_simpsons("springfield")(16)[c(13,5)]) +
  bbplot::bbc_style()



dtf_r %>% ggplot() +
  geom_violin(aes(y = prob, x = c(to_forest  ) ,group = to_forest , col = to_forest , fill = to_forest )) +
  geom_boxplot(aes(y = prob, x = c(to_forest  ),group = to_forest ), width = 0.1) +
  scale_y_continuous(limits = c(0, 5)) +
  scale_color_manual(values = pal_simpsons("springfield")(16)[c(13,5)])  +
  scale_fill_manual(values = pal_simpsons("springfield")(16)[c(13,5)]) +
  bbplot::bbc_style()



library(scales)

show_col(pal_simpsons("springfield")(16))


setwd("/save_home/jslim/Mechanistic_modelling/Marryland/Model_comparison/Candidates_comparisons/Model_B/post_hoc")




dtf_r %>%
  filter(!(f_type %in% c("uf2f_detected", "uf2f_undetected"))) %>%
  group_by(f_type) %>%
  summarise(lo = quantile(prob, prob = 0.025), 
            med = quantile(prob, prob = 0.5),
            up = quantile(prob, prob = 0.975)
  )


dtf_r %>%
  filter(!(f_type %in% c("uf2f_detected", "uf2f_undetected"))) %>%
  ggplot() +
  geom_violin(aes(y = prob, x = f_type,group = f_type, col = f_type, fill = f_type)) +
  geom_boxplot(aes(y = prob, x = f_type,group = f_type), width = 0.1) +
  scale_y_continuous(limits = c(0, 6),breaks = seq(0,6, by =1), "Reproduction number") +
  scale_x_discrete(labels = c("F2F", "F2UF/UF2F", "UF2UF", "F2F", "F2UF/UF2F", "UF2UF")) +
  geom_hline(yintercept = 6, size = 0.7) +
  geom_hline(yintercept = 0, size = 0.7) +
  scale_color_manual(values = pal_simpsons("springfield")(16)[c(13,5, 1, 7, 2, 10)])  +
  scale_fill_manual(values = pal_simpsons("springfield")(16)[c(13,5, 1, 7, 2, 10)]) +
  bbplot::bbc_style() +
  theme(legend.position = "none", 
        axis.title.y = element_text(size = 18), 
        axis.title.x = element_blank()) -> cond_R
cond_R
print(cond_R)



library(ggplot2)
library(dplyr)
library(scales)
pal_lancet <- c(
  
  # 🔴 Red side
  "#B2001D",   # deep Lancet red
  "#E56B6F",  # medium crimson
  
  "#F4A3A1",  # light rose red
  # 🔵 Blue side
  "#01497C",  # deep Lancet blue
  "#2A6EBB",  # medium blue
  "#74B3E7"  # light blue
  
)


# Plot ---------------------------------------------------------------------
dtf_r %>%
  filter(!(f_type %in% c("uf2f_detected", "uf2f_undetected"))) %>%
  ggplot(aes(x = f_type, y = prob, group = f_type,
             colour = f_type, fill = f_type)) +
  geom_violin(trim = FALSE, alpha = 0.6, linewidth = 0.6) +
  geom_boxplot(width = 0.15, outlier.shape = NA,
               fill = "white", colour = "#444444", linewidth = 0.5) +
  geom_hline(yintercept = 1, colour = "#999999", linewidth = 0.4, linetype = "dashed") +
  scale_y_continuous(
    name = "Reproduction number (R)",
    limits = c(0, 6),
    breaks = seq(0, 6, 1)
  ) +
  scale_x_discrete(
    labels = c("F2F", "F2UF/UF2F", "UF2UF", "F2F", "F2UF/UF2F", "UF2UF")
  ) +
  scale_color_manual(values = pal_lancet) +
  scale_fill_manual(values = pal_lancet) +
  theme_minimal() +
  theme(
    # text / titles
    plot.title = element_text(size = 18, face = "bold", colour = "#111111"),
    plot.subtitle = element_text(size = 13, colour = "#555555", margin = margin(b = 6)),
    plot.caption = element_text(size = 9, colour = "#777777", margin = margin(t = 8)),
    axis.title.y = element_text(size = 14, margin = margin(r = 8)),
    axis.title.x = element_blank(),
    axis.text = element_text(size = 12, colour = "#222222"),
    # grid
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#D9D9D9", linewidth = 0.4),
    # legend
    legend.position = "none",
    # margins
    plot.margin = margin(10, 16, 10, 10)
  ) -> cond_R2;cond_R2
