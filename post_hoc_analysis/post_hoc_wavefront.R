
library(units)
library(tidyr)
library(ggsci)




library(tidyverse)
#saveRDS(l_ss, file = "wavefront_data.RDS")
l_ss<-readRDS( file = "wavefront_data.RDS")

l_ss %>% lapply(FUN = function(x){
  x %>% group_by(week) %>%
    arrange(type) %>%
    mutate(diff = -diff(cumu_max)) -> wavefront_data_dif
  return(wavefront_data_dif)
}) %>% bind_rows(.id = "sim") -> dtf_ss


dtf_ss %>% group_by(type, week) %>%
  summarise(lo_dist = quantile(cumu_max/1000, prob = 0.025), 
            med_dist = median(cumu_max/1000),
            up_dist = quantile(cumu_max/1000, prob = 0.975), 
            lo_diff =  quantile(diff/1000, prob = 0.025),
            med_diff = median(diff/1000),
            up_diff = quantile(diff/1000, prob = 0.975),
  ) %>%
  ungroup() -> cl_dtf_ss

cl_dtf_ss$med_diff < 6

cl_dtf_ss %>%
  ggplot(aes(x = week,col = type, fill = type)) +
  geom_step(aes( y =med_dist), alpha = 1) + 
  scale_color_manual(labels = c("True", "Apparent"),values = c("grey80", "black"))  +
  bbplot::bbc_style() 

cl_dtf_ss %>%
  ggplot(aes(x = week,col = type, fill = type)) +
  geom_step(aes( y =med_dist), alpha = 1) + 
  geom_ribbon(aes(ymin = lo_dist, ymax = up_dist), alpha = 0.2) +
  scale_color_manual(labels = c("True", "Apparent"),values = pal_frontiers()(10)[c(1,5)])  +
  scale_fill_manual(labels = c("True", "Apparent"),values = pal_frontiers()(10)[c(1,5)]) + 
  
  bbplot::bbc_style() 



setwd("~/save/Mechanistic_modelling/Marryland/Model_comparison/Candidates_comparisons/Model_B/Trajectories")

png(units = "in", width = 22, height = 8, filename = "longitudinal_distance.png",res = 600)
font <- "Helvetica"

cl_dtf_ss %>%
  mutate(week = ymd("2019-09-01") + weeks(week - 1)) %>%
  ggplot(aes(x = week,col = type, fill = type)) +
  geom_step(aes( y =med_dist), alpha = 1) + 
  geom_step(aes(y = up_dist), alpha = 0.2, linetype = "dashed") +
  geom_step(aes(y = lo_dist), alpha = 0.2, linetype = "dashed") + 
  scale_color_manual(labels = c("True", "Apparent"),values = pal_frontiers()(10)[c(1,5)])  +
  scale_fill_manual(labels = c("True", "Apparent"),values = pal_frontiers()(10)[c(1,5)]) + 
  
  bbplot::bbc_style() +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "3 month") +
  scale_y_continuous("Longitudial distance", labels = scales::label_number(scale = 1, suffix = " km")) +
  ggplot2::theme(plot.title = ggplot2::element_text(family = font, 
                                                    size = 28, face = "bold", color = "#222222"), 
                 plot.subtitle = ggplot2::element_text(family = font, size = 22, margin = ggplot2::margin(9, 0, 9, 0)), plot.caption = ggplot2::element_blank(), 
                 legend.position = "top", legend.text.align = 0, legend.background = ggplot2::element_blank(), 
                 legend.title = ggplot2::element_blank(), 
                 legend.key = ggplot2::element_blank(), 
                 legend.text = ggplot2::element_text(family = font, size = 18, 
                                                     color = "#222222"), 
                 axis.title.y = ggplot2::element_text( size = 18),
                 axis.title.x = ggplot2::element_blank(), 
                 
                 axis.text = ggplot2::element_text(family = font, size = 18, 
                                                   color = "#222222"), 
                 axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5,b = 10)), axis.ticks = ggplot2::element_blank(), 
                 axis.line = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(), 
                 panel.grid.major.y = ggplot2::element_line(color = "#cbcbcb"), 
                 panel.grid.major.x = ggplot2::element_blank(), 
                 panel.background = ggplot2::element_blank(), 
                 strip.background = ggplot2::element_rect(fill = "white"), 
                 strip.text = ggplot2::element_text(size = 22, hjust = 0))


dev.off()


cl_dtf_ss

cl_dtf_ss$up_diff

cl_dtf_ss %>% dplyr::select(week, lo_diff, med_diff, up_diff) %>%
  mutate(week = ymd("2019-09-01") + weeks(week - 1)) %>%
  pivot_longer(names_to = "type", values_to = "value", col = -week) %>% 
  mutate(type = factor(type, levels = c("up_diff", "med_diff", "lo_diff")))%>%
  ggplot(aes(x = week, y = value, group = type, col = type)) +
  geom_step(direction = "vh") +
  #geom_hline(yintercept = 25, linetype = "dashed", color = "red", linewidth = 1) +
  scale_color_manual(values = c( "grey70", "black",  "grey70"), 
                     labels = c("97.5%", "50%","2.5%" ) )+ 
  scale_x_date(date_labels = "%Y-%m", date_breaks = "3 month") +
  scale_y_continuous("Longitudial distance", labels = scales::label_number(scale = 1, suffix = " km")) +
  ggplot2::theme(plot.title = ggplot2::element_text(family = font, 
                                                    size = 28, face = "bold", color = "#222222"), 
                 plot.subtitle = ggplot2::element_text(family = font, size = 22, margin = ggplot2::margin(9, 0, 9, 0)), plot.caption = ggplot2::element_blank(), 
                 legend.position = "top", legend.text.align = 0, legend.background = ggplot2::element_blank(), 
                 legend.title = ggplot2::element_blank(), 
                 legend.key = ggplot2::element_blank(), 
                 legend.text = ggplot2::element_text(family = font, size = 18, 
                                                     color = "#222222"), 
                 axis.title.y = ggplot2::element_text( size = 18),
                 axis.title.x = ggplot2::element_blank(), 
                 
                 axis.text = ggplot2::element_text(family = font, size = 18, 
                                                   color = "#222222"), 
                 axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5,b = 10)), axis.ticks = ggplot2::element_blank(), 
                 axis.line = ggplot2::element_blank(),
                 panel.grid.minor = ggplot2::element_blank(), 
                 panel.grid.major.y = ggplot2::element_line(color = "#cbcbcb"), 
                 panel.grid.major.x = ggplot2::element_blank(), 
                 panel.background = ggplot2::element_blank(), 
                 strip.background = ggplot2::element_rect(fill = "white"), 
                 strip.text = ggplot2::element_text(size = 22, hjust = 0))
