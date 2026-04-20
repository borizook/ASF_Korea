risk_ob_plot <-
  tm_shape(tm_kor, bbox = hex_bbox) +
  tm_polygons(col = "grey80", lwd = 1) +
  
  tm_shape(sf_hex) +
  tm_fill(
    "risk_ob",
    title   = "Prob. of ASF detection",
    palette = c("#FFF5F0", "#FCAe91", "#FB6A4A", "#DE2D26", "#A50F15"),
    style   = "cont",
    alpha   = 1,
    legend.reverse = TRUE      # ✅ 여기 추가!
  ) +
  
  tm_shape(NK_map) +
  tm_fill(fill = "grey60") +
  
  tm_shape(uni_sf_hex) +
  tm_borders(col = "black") +
  
  tm_shape(ob_sf_hex) +
  tm_borders() +
  
  tm_layout(
    frame = FALSE,
    legend.frame = FALSE,
    legend.position = c("right", "bottom")
  ) +
  
  tm_compass(position = c("right", "top"), size = 3) + 
  
  tm_scalebar(position = c("right", "top"), breaks = c(0, 50, 100), text.size = 0.7) 

risk_ob_plot


risk_inci_plot <-
  tm_shape(tm_kor, bbox = hex_bbox) +
  tm_polygons(col = "grey80", lwd = 1) +
  
  tm_shape(sf_hex) +
  tm_fill(
    "risk_inci",
    title   = "Risk of ASF occurrence",
    palette = c("#E0F3F8", "#74A9CF", "#0570B0", "#023858"),  # Lancet blue 계열
    style   = "cont",
    breaks  = seq(0, 1, by = 0.2),
    alpha   = 1,
    legend.reverse = TRUE        # 🔁 범례 순서 반전 (0 아래, 1 위)
  ) +
  
  tm_shape(NK_map) +
  tm_fill(fill = "grey60") +
  
  tm_shape(uni_sf_hex) +
  tm_borders(col = "black") +
  
  tm_layout(
    frame = FALSE,
    legend.frame = FALSE,       # ✅ 범례 테두리 제거
    legend.bg.color = NA,       # ✅ 범례 배경 투명
    legend.position = c("right", "bottom")
  ) +
  
  tm_compass(position = c("right", "top"), size = 3) + 
  
  tm_scalebar(position = c("right", "top"), breaks = c(0, 50, 100), text.size = 0.7) 

risk_inci_plot

setwd("/save_home/jslim/Mechanistic_modelling/Marryland/Model_comparison/Candidates_comparisons/Model_B/post_hoc")

pdf(file = "two_spatial_plots.pdf",   # The directory you want to save the file in
    width = 18, # The width of the plot in inches
    height = 6) 
tmap_arrange(risk_ob_plot, risk_inci_plot, widths = c(.5, .5), ncol = 2)
dev.off()




