
whole_hex %>% st_as_sf %>% 
  mutate(r_med =r_med_sp_undetected_whole )%>% 
  tm_shape() +
  tm_fill("r_med", title = "Reproduction number",
          palette = "-RdYlBu",
          breaks = c(0, 0.5, 1, 2, 3, 4), 
          alpha = 1, 
          legend.reverse = TRUE
  ) +
  tm_shape(NK_map) +
  tm_fill(fill = "grey60")+
  tm_layout(  legend.position = c("right","bottom"), 
              legend.outside = TRUE,
              legend.frame = FALSE,
              frame = FALSE)  ->spatial_R
print(spatial_R)


pal_lancet_div <- c(
  "#023858", "#045A8D", "#74A9CF", "#F7F7F7",  # 파랑→흰색
  "#FEE0D2", "#FB6A4A", "#CB181D"               # 흰색→붉은색
)


whole_hex %>% st_as_sf %>% 
  mutate(r_over_1 =r_sp_undetected_whole_over_1 )%>% 
  tm_shape() +
  tm_fill("r_over_1", title = "Probability of R>1",
          palette = pal_lancet_div,
          style="cont",
          alpha = 1 , 
          legend.na.show = FALSE,
          legend.reverse = TRUE)    +
  tm_shape(NK_map) +
  tm_fill(fill = "grey60")+
  tm_legend(na.show = FALSE, reverse = FALSE) +
  
  tm_layout(legend.position = c("right","bottom"), 
            legend.outside = FALSE,
            frame = FALSE) +
  tm_compass(position = c("right", "top")) +
  tm_scale_bar(position = c("right", "top"), 
               breaks  = c(0, 50, 100))                  -> spatial_R_over_1
# tm_shape(whole_hex) +
# tm_borders() -> spatial_R_over_1
print(spatial_R_over_1)

setwd("/save_home/jslim/Mechanistic_modelling/Marryland/Model_comparison/Candidates_comparisons/Model_B/post_hoc")



pdf(file = "two_R_plots.pdf",   # The directory you want to save the file in
    width = 14, # The width of the plot in inches
    height = 6) 
tmap_arrange(spatial_R, spatial_R_over_1, widths = c(.5, .5), ncol = 2)
dev.off()

