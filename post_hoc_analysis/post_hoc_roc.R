
## Calculating the risk of ASF reports
dt_trajec[[3]] %>%
  group_by(unit) %>%
  summarise(numb = sum(value ==3),
            total_simu = times_permute,
            risk = numb/total_simu, 
            
  ) %>% dplyr::select(risk) -> risk_ob_hexa





risk_ob_hexa %>% 
  mutate(num = 1:n(), 
         binary_ob = case_when(num %in% unit_ID_ob ~ 1, 
                               !num %in% unit_ID_ob ~ 0, )
  ) %>%
  ggplot() +
  geom_violin(aes(x = factor(binary_ob), y = risk, col = factor(binary_ob), fill = factor(binary_ob))) +
  bbc_style() +
  theme(plot.title = element_text(size = 15)) +
  scale_x_discrete(labels = c("Undetected", "Detected")) +
  ggtitle("Predicted risk of detecting ASF-positive hexagon (Model fitting)") +
  scale_y_continuous(limits = c(0, 1)) 

unit_ID_ob<- readRDS("/save_home/jslim/Mechanistic_modelling/Marryland/data/data/unit_ID_ob.rds")

risk_ob_hexa %>% 
  mutate(num = 1:n(), 
         binary_ob = case_when(num %in% unit_ID_ob ~ 1, 
                               !num %in% unit_ID_ob ~ 0, )
  ) -> ssss



library(pROC)

roc(ssss$binary_ob, ssss$risk)-> roc_risk


plot(roc_risk)


