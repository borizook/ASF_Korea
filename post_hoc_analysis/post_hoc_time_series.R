library(ggplot2)
library(dplyr)
setwd("/save_home/jslim/Mechanistic_modelling/Marryland/Model_comparison/Candidates_comparisons/Model_B/post_hoc")

pdf(file = "figure 4.pdf",   # The directory you want to save the file in
    width = 12, # The width of the plot in inches
    height = 6) 


# 색 정의 (Lancet 톤)
blue95 <- "#A6CEE3"; blue50 <- "#5BA3E6"; blueLn <- "#01497C"  # undetected
red95  <- "#F4A3A1"; red50  <- "#E56B6F"; redLn  <- "#B2001D"  # detected
obs_black <- "#222222"

# y축 최대값을 이전 블루 그래프와 동일하게
y_max <- y_max_blue  # 예: 120

p <- ggplot() +
  # --------- RIBBONS: 상태=색(fill), 구간=투명도(alpha) ---------
# Detected (Red)
stat_summary(data = dt_trajec[[1]],
             aes(x = ob_week, y = n,
                 fill  = "Detected",
                 alpha = "95%"),
             geom = "ribbon", fun.data = mean_cl_quantile_95,
             colour = NA) +
  stat_summary(data = dt_trajec[[1]],
               aes(x = ob_week, y = n,
                   fill  = "Detected",
                   alpha = "50%"),
               geom = "ribbon", fun.data = mean_cl_quantile_50,
               colour = NA) +
  # Undetected (Blue)
  stat_summary(data = dt_trajec[[6]],
               aes(x = ob_week, y = n,
                   fill  = "Undetected",
                   alpha = "95%"),
               geom = "ribbon", fun.data = mean_cl_quantile_95,
               colour = NA) +
  stat_summary(data = dt_trajec[[6]],
               aes(x = ob_week, y = n,
                   fill  = "Undetected",
                   alpha = "50%"),
               geom = "ribbon", fun.data = mean_cl_quantile_50,
               colour = NA) +
  
  # ------------------------ MEDIAN LINES ------------------------
stat_summary(data = dt_trajec[[6]],
             aes(x = ob_week, y = n, colour = "Undetected median"),
             geom = "line", fun = median, linewidth = 1.1) +
  stat_summary(data = dt_trajec[[1]],
               aes(x = ob_week, y = n, colour = "Detected median"),
               geom = "line", fun = median, linewidth = 1.1) +
  
  # ------------------------ OBSERVED LINE -----------------------
geom_line(data = sim_same_par_traj_reported,
          aes(x = min_Id_week, y = n, colour = "Observed"),
          linewidth = 0.9) +
  
  # ---------------------------- SCALES --------------------------
# fill: Status만 표시 (색으로 상태 구분)
scale_fill_manual(
  name   = "Status",
  values = c("Undetected" = blue50, "Detected" = red50)
) +
  # alpha: Interval만 표시 (투명도로 구간 구분)
  scale_alpha_manual(
    name   = "Interval",
    values = c("95%" = 0.25, "50%" = 0.45)
  ) +
  # colour: 라인 항목 (Observed 먼저, 그다음 median들)
  scale_color_manual(
    name   = NULL,
    breaks = c("Observed", "Undetected median", "Detected median"),
    values = c("Observed" = obs_black,
               "Undetected median" = blueLn,
               "Detected median"   = redLn)
  ) +
  
  scale_x_continuous(limits = c(x_min, x_max + 1)) +
  scale_y_continuous(limits = c(0, y_max),
                     breaks = seq(0, y_max, by = 20),
                     name = "Number of detected hexagons") +
  
  guides(
    fill   = guide_legend(order = 1, override.aes = list(alpha = 1, colour = NA)),
    alpha  = guide_legend(order = 2),
    colour = guide_legend(order = 3)
  ) +
  
  labs(x = "Week") +
  theme_minimal(base_family = font) +
  theme(
    axis.title.y  = element_text(size = 16, margin = margin(r = 8)),
    axis.title.x  = element_text(size = 15, margin = margin(t = 8)),
    axis.text     = element_text(size = 13, colour = "#222"),
    axis.ticks    = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#D9D9D9", linewidth = 0.4),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.text = element_text(size = 12, colour = "#222"),
    legend.key  = element_blank(),
    plot.margin = margin(10, 16, 10, 10)
  )

p
dev.off()

