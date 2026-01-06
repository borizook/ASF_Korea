
library(doParallel)
library(foreach)
library(tidyverse)
library(bbplot)
library(ggplot2)
library(EasyABC)
library(tidyverse)
###################################################################
# 0. Define the folder where the results will be saved
###################################################################
rm(list = ls())



###################################################################
###################################################################
# Specify scenario number, number of particles for ABC-APMC, 
# number of cores to be used, and the study period

scen_num = "" # scenario number
n.cores = n_cores = 30 # number of cores
end_week = 155 # Study period
nb_simul = 1000 # number of particles

###################################################################
###################################################################

setwd("/Results/")
folder_scen = paste0("Model_I",scen_num)

source("/Model/function_create_wd.R")
create_wd(folder_scen)

setwd(scen_folder_wd)

# 
# folder_prior_pred = "Prior predictive"
# if (folder_prior_pred %in% dir(getwd()) == FALSE) {
#   dir.create(file.path(getwd(), folder_prior_pred))}
# 


folder_posterior_points = "Posterior points"
if (folder_posterior_points %in% dir(getwd()) == FALSE) {
  dir.create(file.path(getwd(), folder_posterior_points))}


###################################################################
# 1. Simulate the model with set of parameters
###################################################################


# Saving the data needed to initialize the model

source(("/Model/init_model_wb.R"))

if (!any(list.files(file.path(scen_folder_wd,"/init_data" )) %in% "initdata.Rdata")){
  
  init.model(end_week = end_week, 
             dt = 1/10, 
             n.breaks = 4, 
             cell.adj.order = "first")
  
  setwd(paste0(file.path(scen_folder_wd,"/init_data" )))
  save(list = ls(), file = "initdata.Rdata", compress = TRUE)
}







# priors for EasyABC ------------------------------------------------------
priors <- list(
  c("unif", 0, 1),       # amplitude of beta
  c("unif", 5/2, 5/2),       # the first peak time in a time unit of this model, 
  c("unif", 0, 1),       # minimum of beta
  c("unif", 2, 2),       # number of cycle per year
  c("unif", 1/4.66, 1/4.66), # sigma - detection rate for wild boar cells (Iu_d to Id)
  c("unif", 1/5.88, 1/5.88), # gamma - sanitation rate for wild boar cells (Id to S)
  c("unif", 0.073, 0.073), # eta - recovery rate (Iu_u to S)
  c("unif", 0, 1),   # 9.  detectibility
  c("unif", 0, 1),   # 9.  rel_suscep (should be between 0 - 1)
  c("unif", 0, 1)    # 9.  rel_infect (should be between 0 - 1)  
)


# essential_pars ---
#  sigma - detection rate for wild boar cells (Iu_d to Id)
#  gamma - sanitation rate for wild boar cells (Id to S)
#  eta - Recovery Rate from Iu to S

pars_estimate = c("amplitutde of beta", 
                  "minimum of beta",
                  "detectibility",
                  "rel_suscep", 
                  "rel_infect")
pars = c(0.07, 5/2, 0.37, 2 ,1/4.66, 1/5.88, 0.073, 0.2, 0.5, 0.5)
pars_d = matrix(ncol = length(pars), pars )
colnames(pars_d) = c("amplitutde of beta", 
                     "h_shift of beta", 
                     "minimum of beta",
                     "number of cycle", 
                     "detection rate",
                     "sanitation rate", 
                     "recovery rate",
                     "detectibility", 
                     "rel_suscep", 
                     "rel_infect")


model <- function(x){
  # I0 <- c(sample(I0, 1), sample(unitData$unit.id[unitData$county == "Tulcea" & unitData$dp == 0 & unitData$edge == 1], 1))
  
  source(("/Model/cal_sum_stats.sim_ver3_extended.R"))
  source(("/Model/Simulate_model_wb_general.R"))
  source(("/Model/gen.matsum.sim.R"))
  load("initdata.Rdata", .GlobalEnv)
  
  
  library(tidyverse)
  
  set.seed(x[1])
  SimulateModel(end_week = end_week, 
                N = N, 
                adjacency = adjacency,
                I0 = I0, 
                seasonality = c(TRUE),
                seasonality_beta = x[(2:5)],
                essential_pars = x[c(c(6:8))],
                detect_pars  = x[9], 
                susc_infect_pars = x[c(10,11)], 
                unitData = unitData, 
                dt = 1/10,
                trans.mode = "dens", 
                forest_unit = "decimal", 
                detect = "constant", 
                susc_infect = "different constant")
  
  
  matsum <- gen.matsum.sim(matr.states, I0)
  sum.stats <- calc.sum.stats.sim(matsum, n.breaks = 4,
                                  sum_stats_output = c("incidence", "forest", "dist", "weekly", "ele"),
                                  sum_stats_spatial = c("user_defined"),
                                  ob_indicator = 3)
  return(sum.stats)
  rm(list = ls())
}








if (
  !any(
    list.files(
      file.path(scen_folder_wd,"/Posterior distribution/" )
    ) %in% "results_ABC.rds")){
  
  
  
  sum.stats.obs = readRDS("/data/obs.soi-nobuff.rds")
  
  # Run APMC parameter estimation -------------------------------------------
  setwd(paste0(file.path(scen_folder_wd,"/init_data" )))
  set.seed(1024)
  
  ABC_Lenormand <- ABC_sequential(method="Lenormand", model=model, prior=priors,
                                  nb_simul=nb_simul, summary_stat_target=sum.stats.obs,
                                  p_acc_min=0.05, verbose = TRUE, use_seed = TRUE, n_cluster = n.cores)
  
  
  
  saveRDS(ABC_Lenormand, file = paste0(file.path(scen_folder_wd,"/Posterior distribution/results_ABC.rds" )))
}

ABC_Lenormand <- readRDS(paste0(file.path(scen_folder_wd,"/Posterior distribution/results_ABC.rds" )))



post_param <- ABC_Lenormand$param %>% apply(MARGIN = 2, FUN = median)



if (!any(list.files(file.path(scen_folder_wd, folder_posterior_points)) %in% "posterior_points.RDS")){
  
  
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  
  doParallel::registerDoParallel(cores = my.cluster)
  foreach::getDoParRegistered()
  foreach(i = 1:(nb_simul/2), .combine = rbind, .packages = c("tidyverse") , .export = "scen_folder_wd") %dopar% {
    
    source(("/model/cal_sum_stats.sim_ver3_extended.R"))
    source(("/model/Simulate_model_wb_general.R"))
    source(("/model/gen.matsum.sim.R"))
    setwd(paste0(file.path(scen_folder_wd,"/init_data" )))
    
    load("initdata.Rdata", .GlobalEnv)
    
    
    library(tidyverse)
    
    set.seed(i)
    
    pars <-c(
      post_param[1],
      runif(1, 5/2, 5/2),
      post_param[2],
      runif(1, 2, 2),
      runif(1, 1/4.66, 1/4.66),
      runif(1,  1/5.88,  1/5.88),
      runif(1, 0.073, 0.073),
      post_param[3],
      post_param[4], 
      post_param[5]
      
    )
    names(pars) = c("amplitutde of beta", 
                       "h_shift of beta", 
                       "minimum of beta",
                       "number of cycle", 
                       "detection rate",
                       "sanitation rate", 
                       "recovery rate",
                       "detectibility", 
                       "rel_suscep", 
                       "rel_infect")
    
    
    
    
    set.seed(i)
    
    SimulateModel(end_week = end_week, 
                  N = N, 
                  adjacency = adjacency,
                  I0 = I0, 
                  seasonality = c(TRUE),
                  seasonality_beta = pars[1:4],
                  essential_pars = pars[5:7],
                  detect_pars  = pars[c(8)], 
                  susc_infect_pars = pars[c(9, 10)], 
                  unitData = unitData, 
                  dt = 1/10,
                  trans.mode = "dens", 
                  forest_unit = "decimal", 
                  detect = "constant", 
                  susc_infect = "different constant")
    
    matsum <- gen.matsum.sim(matr.states, I0)
    sum.stats <- calc.sum.stats.sim(matsum, n.breaks = 4,
                                    sum_stats_output = c("incidence", "forest", "dist", "weekly", "ele"),
                                    sum_stats_spatial = c("user_defined"),
                                    ob_indicator = 3)
    return(list(pars,matr.states, sum.stats))
    
  } -> post_points
  
  parallel::stopCluster(cl = my.cluster)
  
  
  saveRDS(post_points, paste0(file.path(scen_folder_wd, folder_posterior_points),"/posterior_points.RDS"))
}



library(tidyverse)
library(ggplot2)
library(bbplot)
###################################################################
# 3. Simulated the trajectories by using the estimated parameters
###################################################################
names(priors) = colnames(pars_d)


results_ABC <-  ABC_Lenormand <-readRDS(paste0(file.path(scen_folder_wd,"/Posterior distribution/results_ABC.rds" )))

post_ABC<- results_ABC$param


colnames(post_ABC)  = pars_estimate

priors[pars_estimate] -> priors_esti

lapply(seq_len(ncol(post_ABC)), FUN = function(i){ post_ABC[,i]}) %>% 
  setNames(pars_estimate) -> list_post_ABC


source("/post_hoc_analysis/overlap_area_function.R")



lapply(seq_len(length(list_post_ABC)), FUN = function(x){
  
  return(overlap_area(distri_prior = priors_esti[[x]][c(2,3)] %>% as.numeric ,posterior = post_ABC[,x], priors_names  = pars_estimate[[x]]))
  
  
}) -> l_overlap_result


post_ABC %>% as.data.frame %>%
  pivot_longer(cols =everything(),
               names_to = "var", 
               values_to = "post") %>%
  group_by(var) %>%
  summarise(low_95 = quantile(post, probs = 0.025), 
            median = quantile(post, probs = 0.5),
            high_95 = quantile(post, probs = 0.975)) %>%
  arrange(factor(var, levels = pars_estimate )) -> post_summary




setwd(paste0(scen_folder_wd, "/Posterior distribution"))
write_csv(post_summary, "post_summary.csv")
post_summary %>% t() %>% as.data.frame -> format_post_summary
write_csv(format_post_summary, "format_post_summary.csv")



setwd(paste0(scen_folder_wd, "/Posterior distribution/overlap"))
for (i in seq_len(length(l_overlap_result))){
  jpeg(paste0("overlap_",pars_estimate[i],".jpg"), width = 600, height = "600", quality = 100)
  print(l_overlap_result[[i]][[2]])
  dev.off()
}



setwd(paste0(scen_folder_wd, "/Posterior distribution/overlap"))
for (i in seq_len(length(l_overlap_result))){
  pdf(file = paste0("overlap_",pars_estimate[i],".pdf"),   # The directory you want to save the file in
      width = 6, # The width of the plot in inches
      height = 6) 
  print(l_overlap_result[[i]][[2]])
  dev.off()
}


setwd(scen_folder_wd)
# list_temporal_trajec = list(NA)

if (!any(list.files(paste0(scen_folder_wd,"/Trajectories/" )) %in% "simulated_trajectories.RDS")){
  
  library(doParallel)
  library(foreach)
  
  
  
  my.cluster <- parallel::makeCluster(
    n.cores, 
    type = "PSOCK"
  )
  
  doParallel::registerDoParallel(cores = my.cluster)
  foreach::getDoParRegistered()
  traj_simul = nb_simul/2
  setwd(file.path(scen_folder_wd,"/init_data" ))
  foreach(i = 1:traj_simul, .combine = rbind, .packages = c("tidyverse"), .export = c("post_ABC")) %dopar% {
    
    
    source(("/save_home/jslim/Mechanistic_modelling/Marryland/code/Model_4_3/cal_sum_stats.sim_ver3_extended.R"))
    source(("/save_home/jslim/Mechanistic_modelling/Marryland/code/Model_4_3/Simulate_model_wb_general.R"))    
    source(("/save_home/jslim/Mechanistic_modelling/Marryland/code/Model_4_3/gen.matsum.sim.R"))
    load("initdata.Rdata", .GlobalEnv)
    which(colnames(pars_d) %in% pars_estimate) -> pl_estimated
    pars[pl_estimated] = post_ABC[i, seq_along(pars_estimate)]
    
    set.seed(i)
    
    SimulateModel(end_week = end_week, 
                  N = N, 
                  adjacency = adjacency,
                  I0 = I0, 
                  seasonality = c(TRUE),
                  seasonality_beta = pars[(1:4)],
                  essential_pars = pars[c(5:7)],
                  detect_pars  = pars[8], 
                  susc_infect_pars = pars[c(9,10)], 
                  unitData = unitData, 
                  dt = 1/10,
                  trans.mode = "dens", 
                  forest_unit = "decimal", 
                  detect = "constant", 
                  susc_infect = "different constant")
    
    
    
    ## incidence of time series
    
    
    matr.states %>% gen.matsum.sim(I0)    %>% 
      filter(value == 3) %>%
      arrange(week) %>%
      summarise(ob_week = min(week)) %>%
      pull ->ob_start
    
    
    
    matr.states %>% gen.matsum.sim(I0 =  I0) %>% 
      filter( value %in% c(1,2)) %>%
      mutate(ob_week = week - ob_start + 1 ) %>% 
      group_by(ob_week) %>%
      summarise(n =n()) %>% 
      complete(ob_week = min(ob_week):c(end_week), fill = list(n=0)) %>%
      arrange(ob_week)-> sim_same_par_traj_inci_total
    
    
    
    matr.states %>% gen.matsum.sim(I0)    %>% 
      filter(value == 3) %>%
      arrange(week) %>%
      summarise(ob_week = min(week)) %>%
      pull ->ob_start
    
    
    
    matr.states %>% gen.matsum.sim(I0 =  I0) %>% 
      filter( value %in% c(1)) %>%
      mutate(ob_week = week - ob_start + 1 ) %>% 
      group_by(ob_week) %>%
      summarise(n =n()) -> temp_l
    if(nrow(temp_l)>0) {
      
      temp_l %>% complete(ob_week = min(ob_week):c(end_week), fill = list(n=0)) %>%
        arrange(ob_week)-> sim_same_par_traj_inci_undetected
      
      
    } else {
      
      tibble(ob_week = 1:155, n = 0) -> sim_same_par_traj_inci_undetected
      
    }
    
    
    matr.states %>% gen.matsum.sim(I0 =  I0) %>% 
      filter( value %in% c(2)) %>%
      mutate(ob_week = week - ob_start + 1 ) %>% 
      group_by(ob_week) %>%
      summarise(n =n()) %>% 
      complete(ob_week = min(ob_week):c(end_week), fill = list(n=0)) %>%
      arrange(ob_week)-> sim_same_par_traj_inci_to_be_detected
    
    
    ## reported time series
    
    matr.states %>% gen.matsum.sim(I0) %>% summarise(max_week = max(week)) %>% pull -> max_week
    
    
    
    matr.states %>% gen.matsum.sim(I0) %>%
      filter(value == 3) %>%
      arrange(week) %>%
      summarise(ob_week = min(week)) %>%
      pull ->ob_start
    
    
    
    matr.states %>% gen.matsum.sim(I0 =  I0) %>% 
      filter( value == 3) %>%
      mutate(ob_week = week - ob_start + 1 ) %>% 
      group_by(ob_week) %>%
      summarise(n =n()) %>% 
      filter(ob_week <=end_week) %>%
      complete(ob_week = min(ob_week):end_week, fill = list(n=0)) %>%
      arrange(ob_week)-> sim_same_par_traj_reported
    
    
    
    
    
    ## spatital_reported
    
    
    gen.matsum.sim(matr.states, I0) %>%
      filter(value == 3) %>%
      dplyr::select(unit, value) %>%
      complete(unit = 1:N, fill = list(value = 0)) %>% 
      distinct(unit, .keep_all = TRUE) %>%
      arrange(unit)-> sim_same_par_traj_spatial_reported
    
    
    
    ## spatital_incidence
    
    gen.matsum.sim(matr.states, I0) %>%
      filter( value %in% c(1,2)) %>%
      dplyr::select(unit, value) %>%
      complete(unit = 1:N, fill = list(value = 0)) %>%
      distinct(unit, .keep_all = TRUE) %>%
      arrange(unit)-> sim_same_par_traj_spatial_inci
    
    
    
    ## The number of reinfection
    
    matr.states %>% gen.matsum.sim(I0 =  I0) %>% 
      filter(value %in% 1:2) %>% 
      group_by(unit) %>% summarise(n = n()) %>% 
      mutate(num_reinfect = n - 1) -> table_num_reinfect
    
    
    
    return(list(sim_same_par_traj_reported,
                sim_same_par_traj_inci_total,
                sim_same_par_traj_spatial_reported, sim_same_par_traj_spatial_inci, ob_start,
                sim_same_par_traj_inci_undetected,sim_same_par_traj_inci_to_be_detected, table_num_reinfect, matr.states))
    
    
  } -> list_trajec
  parallel::stopCluster(cl = my.cluster)
  
  saveRDS(list_trajec, paste0(scen_folder_wd, "/Trajectories/", "simulated_trajectories.RDS"))
  
}
# 
# list_temporal_trajec %>% 
#   unlist(recursive = FALSE) %>% 
#   do.call(what = rbind.data.frame) %>% 
#   group_by(ob_week) %>% 
#   summarise(median = median(n), 
#             lower = quantile(n, 0.025), 
#             higher = quantile(n, 0.975)) -> dtf_trajectories

# library(bbplot)
# devtools::install_github('bbc/bbplot')
# ggplot(data = dtf_trajectories) +
#   geom_line(aes(x = ob_week, y = median), linewidth = 2) +
#   geom_line(aes(x = ob_week, y = lower), col = "grey50", linewidth = 1, linetype = 2) +
#   geom_line(aes(x = ob_week, y = higher) , col = "grey50", linewidth = 1, linetype = 2) + 
#   bbc_style()
# 




###################################################################
# 4. Plotting the posterior distribution of the estimated parameters
###################################################################

gc()

## ggplot for histogram and 95% interval

rm(list=ls()[!ls() %in% c("l_overlap_result","list_post_ABC","pars_estimate","n_cores" ,"n.cores", "unobserved_day", "end_week","priors", "scen_folder_wd", "pars_d")])


results_ABC <-   readRDS(file.path(scen_folder_wd,"/Posterior distribution/results_ABC.rds" ))

# infect_non_forest = results_ABC$param[,2]
# sus_non_forest = results_ABC$param[,3]
# detect_rate = results_ABC$param[,4]
library(bayesplot)
library(cowplot)


post_ABC<- results_ABC$param



colnames(post_ABC)  = pars_estimate

priors[pars_estimate] -> priors_esti

lapply(seq_len(ncol(post_ABC)), FUN = function(i){ post_ABC[,i]}) %>% 
  setNames(pars_estimate) -> list_post_ABC






post_plot_hist_dens = function(mat_post, par_name){
  
  mat_post %>% quantile(probs = c(0.025, 0.5, 0.975)) -> cri_var
  
  density_plot<<- matrix(ncol= 1, mat_post) %>% as.data.frame  %>% ggplot(data = .) +
    geom_density(aes(x = V1), fill = "deepskyblue3", col = "white",size = 0.3) + 
    geom_segment(aes(x = cri_var[[1]], y = 0, xend = cri_var[[3]], yend = 0), col = "black", size = 0.5) +
    geom_point(aes(x = cri_var[[2]], y= 0),size =4)+
    
    annotate(geom="text", 
             x=cri_var[[1]], 
             y=0.05, 
             label=round(cri_var[[1]], 2),
             size = 5,
             color="black", 
             fontface = "bold") +
    annotate(geom="text", 
             x=cri_var[[3]], 
             y=0.05, 
             label=round(cri_var[[3]], 2),
             size = 5,
             color="black", 
             fontface = "bold") +
    annotate(geom="text", 
             x=cri_var[[2]], 
             y=0.05, 
             label=round(cri_var[[2]], 2),
             size = 5,
             color="black", 
             fontface = "bold") +
    bbc_style() 
  
  hist_plot <<-matrix(ncol= 1, mat_post) %>% as.data.frame  %>% ggplot(data = .) +
    geom_histogram(aes(x = V1, y = after_stat(density)), fill = "deepskyblue3",col = "white", size = 0.3) + 
    geom_segment(aes(x = cri_var[[1]], y = 0, xend = cri_var[[3]], yend = 0), col = "black", size = 0.5) +
    geom_point(aes(x = cri_var[[2]], y= 0),size =4) +
    
    annotate(geom="text", 
             x=cri_var[[1]], 
             y=0.05, 
             label=round(cri_var[[1]], 2),
             size = 5,
             color="black", 
             fontface = "bold") +
    annotate(geom="text", 
             x=cri_var[[3]], 
             y=0.05, 
             label=round(cri_var[[3]], 2),
             size = 5,
             color="black", 
             fontface = "bold") +
    annotate(geom="text", 
             x=cri_var[[2]], 
             y=0.05, 
             label=round(cri_var[[2]], 2),
             size = 5,
             color="black", 
             fontface = "bold") +
    bbc_style() 
  
  var_plot<<- plot_grid(hist_plot, density_plot )
  
  
}
setwd(paste0(file.path(scen_folder_wd, "/Posterior distribution"), "/dist_plot"))

for(i in seq_len(length(list_post_ABC))){
  jpeg(paste0("post_",pars_estimate[i],".jpg"), width = 1200, height = "600", quality = 100)
  post_plot_hist_dens(mat_post = list_post_ABC[[i]], pars_estimate[[i]])
  title <- ggdraw() + draw_label(pars_estimate[[i]], fontface='bold', size = 25)
  print(plot_grid(title,var_plot, ncol=1, rel_heights=c(0.1, 1))) # rel_heights values control title margins
  dev.off()
}


for (i in seq_len(length(list_post_ABC))){
  pdf(file = paste0("post_",pars_estimate[i],".pdf"),   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 4) 
  post_plot_hist_dens(mat_post = list_post_ABC[[i]], pars_estimate[[i]])
  title <- ggdraw() + draw_label(pars_estimate[[i]], fontface='bold', size = 25)
  print(plot_grid(title,var_plot, ncol=1, rel_heights=c(0.1, 1))) # rel_heights values control title margins
  dev.off()
}






post_plot_hist_dens = function(mat_post, par_name){
  
  mat_post %>% quantile(probs = c(0.025, 0.5, 0.975)) -> cri_var
  
  density_plot<<- matrix(ncol= 1, mat_post) %>% as.data.frame  %>% ggplot(data = .) +
    geom_density(aes(x = V1), fill = "deepskyblue3", col = "white",size = 0.3) + 
    geom_segment(aes(x = cri_var[[1]], y = 0, xend = cri_var[[3]], yend = 0), col = "black", size = 0.5) +
    geom_point(aes(x = cri_var[[2]], y= 0),size =4)+
    annotate(geom="text", 
             x=cri_var[[1]], 
             y=0.05, 
             label=round(cri_var[[1]], 2),
             size = 5,
             color="black", 
             fontface = "bold") +
    annotate(geom="text", 
             x=cri_var[[3]], 
             y=0.05, 
             label=round(cri_var[[3]], 2),
             size = 5,
             color="black", 
             fontface = "bold") +
    annotate(geom="text", 
             x=cri_var[[2]], 
             y=0.05, 
             label=round(cri_var[[2]], 2),
             size = 5,
             color="black", 
             fontface = "bold") +
    bbc_style() 
  
  hist_plot <<-matrix(ncol= 1, mat_post) %>% as.data.frame  %>% ggplot(data = .) +
    geom_histogram(aes(x = V1, y = after_stat(density)), fill = "deepskyblue3",col = "white", size = 0.3) + 
    geom_segment(aes(x = cri_var[[1]], y = 0, xend = cri_var[[3]], yend = 0), col = "black", size = 0.5) +
    geom_point(aes(x = cri_var[[2]], y= 0),size =4) +
    annotate(geom="text", 
             x=cri_var[[1]], 
             y=0.05, 
             label=round(cri_var[[1]], 2),
             size = 5,
             color="black", 
             fontface = "bold") +
    annotate(geom="text", 
             x=cri_var[[3]], 
             y=0.05, 
             label=round(cri_var[[3]], 2),
             size = 5,
             color="black", 
             fontface = "bold") +
    annotate(geom="text", 
             x=cri_var[[2]], 
             y=0.05, 
             label=round(cri_var[[2]], 2),
             size = 5,
             color="black", 
             fontface = "bold") +
    bbc_style() 
  
  var_plot<<- plot_grid(hist_plot, density_plot )
  
  
}
setwd(paste0(file.path(scen_folder_wd, "/Posterior distribution"), "/dist_plot"))



for(i in seq_len(length(list_post_ABC))){
  jpeg(paste0("post_",pars_estimate[i],".jpg"), width = 1200, height = "600", quality = 100)
  post_plot_hist_dens(mat_post = list_post_ABC[[i]], pars_estimate[[i]])
  title <- ggdraw() + draw_label(pars_estimate[[i]], fontface='bold', size = 25)
  print(plot_grid(title,var_plot, ncol=1, rel_heights=c(0.1, 1)))# rel_heights values control title margins
  dev.off()
}


for (i in seq_len(length(list_post_ABC))){
  pdf(file = paste0("post_",pars_estimate[i],".pdf"),   # The directory you want to save the file in
      width = 12, # The width of the plot in inches
      height = 4) 
  post_plot_hist_dens(mat_post = list_post_ABC[[i]], pars_estimate[[i]])
  title <- ggdraw() + draw_label(pars_estimate[[i]], fontface='bold', size = 25)
  print(plot_grid(title,var_plot, ncol=1, rel_heights=c(0.1, 1))) # rel_heights values control title margins
  dev.off()
}






amp <- post_ABC[,c(1)]

h_shift <- rep(pars_d[2],nrow(post_ABC))

minimum <- post_ABC[,c(2)]

num_cyc_year <- rep(pars_d[4],nrow(post_ABC))


real_time = 1:155



beta_par<-list(amp, h_shift, minimum, num_cyc_year)

names(beta_par) = c("amp", "h_shift","minimum", "num_cyc_year")



beta_par %>% bind_cols %>% apply(MARGIN = 1, FUN = function(X){

  amp <- X[1]
  base_line <- X[3]
  num_cyc_year <- X[4]
  h_shift <- X[2]

  x = 1:155
  return(amp * sin(num_cyc_year * 2 * pi * (x - h_shift)/52) + base_line + amp)



}) %>% t %>% as_tibble -> beta_traj_post



names(beta_traj_post) = 1:155

beta_traj_post %>% mutate(sim = 1:n()) %>% pivot_longer(cols = `1`:`155`, names_to = "name") %>%
  mutate(name_num = as.numeric(name),
         sim_char = as.character(sim))->long_beta_traj



library(bbplot)


mean_cl_quantile_95 <- function(x, q = c(0.025, 0.975), na.rm = TRUE){
  dat <- data.frame(y = mean(x, na.rm = na.rm),
                    ymin = quantile(x, probs = q[1], na.rm = na.rm),
                    ymax = quantile(x, probs = q[2], na.rm = na.rm))
  return(dat)
}


mean_cl_quantile_50 <- function(x, q = c(0.25, 0.75), na.rm = TRUE){
  dat <- data.frame(y = mean(x, na.rm = na.rm),
                    ymin = quantile(x, probs = q[1], na.rm = na.rm),
                    ymax = quantile(x, probs = q[2], na.rm = na.rm))
  return(dat)
}

font <- "Helvetica"
library(bbplot)
library(scales)


long_beta_traj %>% mutate(date = weeks(long_beta_traj$name_num) + ymd("2019-09-01")) -> long_beta_traj_n


long_beta_traj_n %>%
  ggplot(data = .) +
  # geom_line(alpha = 0.05, aes(x = date, y = value, group = sim)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_quantile_50, aes(x = date,y = value, fill = "50%"), alpha = 0.3) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_quantile_95, aes(x = date,y = value, fill = "95%"), alpha = 0.15) +
  stat_summary(geom = "line", aes(x = date,y = value), col = "red", alpha = 1, size = 1) +
  scale_fill_manual("",values = c("50%" = "red", "95%" = "red")) +
  scale_x_date(date_breaks = "2 months", date_labels = "%B") +
  ggplot2::theme(

    plot.title = ggplot2::element_text(family = font,
                                       size = 28, face = "bold", color = "#222222"),
    plot.subtitle = ggplot2::element_text(family = font, size = 22, margin = ggplot2::margin(9, 0, 9, 0)),
    plot.caption = ggplot2::element_blank(),
    legend.position = "top",
    legend.text = ggplot2::element_text(family = font, size = 18),
    axis.title = element_text(family = font, size = 20),
    axis.text = ggplot2::element_text(family = font, size = 18,
                                      color = "#222222"),
    axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5,
                                                                 b = 10), size = 10),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(color = "#cbcbcb"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "white"),
    strip.text = ggplot2::element_text(size = 22, hjust = 0)) +
  labs(x = "Time", y = "Transmission rate") -> beta_seas
setwd(paste0(file.path(scen_folder_wd, "/Posterior distribution"), "/dist_plot"))

pdf(file = paste0("post_beta_season.pdf"),   # The directory you want to save the file in
    width = 16, # The width of the plot in inches
    height = 4)
plot(beta_seas)
dev.off()



long_beta_traj_n %>%
  ggplot(data = .) +
  geom_line(alpha = 0.05, aes(x = date, y = value, group = sim)) +
  stat_summary(geom = "line", aes(x = date,y = value), col = "red", alpha = 1, size = 1) +
  scale_fill_manual("",values = c("50%" = "red", "95%" = "red")) +
  scale_x_date(date_breaks = "2 months", date_labels = "%B") +
  ggplot2::theme(

    plot.title = ggplot2::element_text(family = font,
                                       size = 28, face = "bold", color = "#222222"),
    plot.subtitle = ggplot2::element_text(family = font, size = 22, margin = ggplot2::margin(9, 0, 9, 0)),
    plot.caption = ggplot2::element_blank(),
    legend.position = "top",
    legend.text = ggplot2::element_text(family = font, size = 18),
    axis.title = element_text(family = font, size = 20),
    axis.text = ggplot2::element_text(family = font, size = 18,
                                      color = "#222222"),
    axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5,
                                                                 b = 10), size = 10),
    axis.ticks = ggplot2::element_blank(),
    axis.line = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.grid.major.y = ggplot2::element_line(color = "#cbcbcb"),
    panel.grid.major.x = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    strip.background = ggplot2::element_rect(fill = "white"),
    strip.text = ggplot2::element_text(size = 22, hjust = 0)) +
  labs(x = "Time", y = "Transmission rate") -> beta_seas_2



setwd(paste0(file.path(scen_folder_wd, "/Posterior distribution"), "/dist_plot"))

pdf(file = paste0("post_beta_season_2.pdf"),   # The directory you want to save the file in
    width = 16, # The width of the plot in inches
    height = 4)
plot(beta_seas_2)
dev.off()







setwd(scen_folder_wd)

###################################################################
# 5. Plotting the trajectories based on the estimated parameters
###################################################################
rm(list=ls()[!ls() %in% c("l_overlap_result","pars_estimate","n_cores" ,"n.cores", "unobserved_day", "end_week","priors", "scen_folder_wd", "pars_d")])

setwd(paste0(file.path(scen_folder_wd,"/init_data" )))
load("initdata.Rdata", .GlobalEnv)

# Prepare the simulated trajectories using the posterior distributions.
# list(sim_same_par_traj_reported,sim_same_par_traj_inci, sim_same_par_traj_spatial_reported, sim_same_par_traj_spatial_inci))
# 
# return(list(sim_same_par_traj_reported,
#             sim_same_par_traj_inci_total,
#             sim_same_par_traj_spatial_reported, sim_same_par_traj_spatial_inci, ob_start,
#             sim_same_par_traj_inci_undetected,sim_same_par_traj_inci_to_be_detected))
# 

list_trajec = readRDS(paste0(scen_folder_wd, "/Trajectories/", "simulated_trajectories.RDS"))

list_trajec %>% apply(MARGIN = 1, FUN= function(x){
  
  x[[3]]
  
})->list_ob_hex


list_trajec %>% apply(MARGIN = 1, FUN= function(x){
  
  x[[4]]
  
})->list_inci_hex

list_trajec %>% apply(MARGIN = 2, FUN = function(x){
  
  return(do.call(what = rbind.data.frame, x))
  
}) -> dt_trajec


nrow(dt_trajec[[1]])/(end_week) -> times_permute

dt_trajec[[2]] %>% dplyr::summarise(min_week = min(ob_week)) %>% pull -> min_week

# dt_trajec[[3]] %>% dplyr::mutate(week = rep(rep(min_week:end_week, each = N), time = times_permute))-> dt_trajec[[3]]
dt_trajec[[3]]$posterior = rep(1:times_permute, each = N )

# dt_trajec[[4]] %>% dplyr::mutate(week = rep(rep(min_week:end_week, each = N), time = times_permute))-> dt_trajec[[4]]
dt_trajec[[4]]$posterior = rep(1:times_permute, each = N )

dt_trajec[[5]][[1]] %>% quantile(probs = c(0.025, 0.5, 0.975)) %>% as.data.frame -> dist_first_delay

names(dist_first_delay) = "The delayed time period from the primary case"


setwd(paste0(file.path(scen_folder_wd,"/Posterior distribution/")))

write_csv(dist_first_delay, "delayed_time_primary_case.csv")


names(dt_trajec) = c("sim_same_par_traj_reported",
                     "sim_same_par_traj_inci_total", 
                     "sim_same_par_traj_spatial_reported", 
                     "sim_same_par_traj_spatial_inci", 
                     "sim_same_par_traj_inci_undetected",
                     "sim_same_par_traj_inci_to_be_detected")



mean_cl_quantile_95 <- function(x, q = c(0.025, 0.975), na.rm = TRUE){
  dat <- data.frame(y = mean(x, na.rm = na.rm),
                    ymin = quantile(x, probs = q[1], na.rm = na.rm),
                    ymax = quantile(x, probs = q[2], na.rm = na.rm))
  return(dat)
}


mean_cl_quantile_50 <- function(x, q = c(0.25, 0.75), na.rm = TRUE){
  dat <- data.frame(y = mean(x, na.rm = na.rm),
                    ymin = quantile(x, probs = q[1], na.rm = na.rm),
                    ymax = quantile(x, probs = q[2], na.rm = na.rm))
  return(dat)
}




# Prepare the trajectory that was used to fit the model.


source(("/model/gen.matsum.sim.R"))





source(("/model/gen.matsum.R"))


## reported time series from fitted data
sim_same_par_traj_reported = readRDS("/save_home/jslim/Mechanistic_modelling/Marryland/data/data/sim_same_par_traj_reported.rds")
# Plotting the trajectories that were simulated by using posterior distribution and was used for model fitting.

font <- "Helvetica"

dt_trajec[[1]] %>% 
  pull(ob_week) %>% min ->x_min

dt_trajec[[1]] %>% 
  pull(ob_week) %>% max ->x_max

setwd((file.path(scen_folder_wd, "/Trajectories/Temporal/")))



jpeg(paste0("trajectory_ob_temporal.jpg"), width = 1000, height = "600", quality = 100)
ggplot(data = dt_trajec[[1]], aes(x = ob_week)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_quantile_95, aes(y = n, fill = "95%"), alpha = 0.1) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_quantile_50, aes(y = n, fill = "50%"), alpha = 0.3) +
  stat_summary(geom = "line", fun = median, aes(y = n, color = "Median"), alpha = 1, size = 1) +
  scale_fill_manual("",values = c("50%" = "red", "95%" = "red")) +
  scale_alpha_manual("",values = c(0.1, 0.3, 1)) +
  geom_line(data = sim_same_par_traj_reported, aes(x = min_Id_week, y = n, color = "Simulated Data")) +
  scale_color_manual("",values = c("Median" = "red", "Simulated Data" = "black")) +
  scale_x_continuous(limits = c(x_min, x_max + 1)) +
  guides(
    fill = guide_legend(title = "", override.aes = list(alpha = c(0.3, 0.1), size = 4)),
    color = guide_legend(title = NULL)) +
  ggplot2::theme(
    
    plot.title = ggplot2::element_text(family = font, 
                                       size = 28, face = "bold", color = "#222222"), 
    plot.subtitle = ggplot2::element_text(family = font, size = 22, margin = ggplot2::margin(9, 0, 9, 0)), 
    plot.caption = ggplot2::element_blank(), 
    legend.position = "top", 
    legend.text = ggplot2::element_text(family = font, size = 18), 
    axis.title = element_text(family = font, size = 20),
    axis.text = ggplot2::element_text(family = font, size = 18, 
                                      color = "#222222"), 
    axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5, 
                                                                 b = 10)), 
    axis.ticks = ggplot2::element_blank(), 
    axis.line = ggplot2::element_blank(), 
    panel.grid.minor = ggplot2::element_blank(), 
    panel.grid.major.y = ggplot2::element_line(color = "#cbcbcb"), 
    panel.grid.major.x = ggplot2::element_blank(), 
    panel.background = ggplot2::element_blank(), 
    strip.background = ggplot2::element_rect(fill = "white"), 
    strip.text = ggplot2::element_text(size = 22, hjust = 0)) +
  labs(x = "Week", y = "Number of detected hexagons") -> trj_ob;trj_ob;
print(trj_ob)
dev.off()





setwd((file.path(scen_folder_wd, "/Trajectories/Temporal/")))
pdf(file = paste0("trajectory_ob_temporal.pdf"),   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 6) 
print(trj_ob)
dev.off()






setwd((file.path(scen_folder_wd, "/Trajectories/Temporal/")))
jpeg(paste0("trajectory_ob_temporal.jpg"), width = 1000, height = "600", quality = 100)
print(trj_ob)
dev.off()



setwd((file.path(scen_folder_wd, "/Trajectories/Temporal/")))

ggplot(data = dt_trajec[[2]], aes(x = ob_week)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_quantile_95, aes(y = n, fill = "95%"), alpha = 0.1) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_quantile_50, aes(y = n, fill = "50%"), alpha = 0.3) +
  stat_summary(geom = "line", fun = median, aes(y = n, color = "Median"), alpha = 1, size = 1) +
  scale_fill_manual("",values = c("50%" = "blue", "95%" = "blue")) +
  scale_alpha_manual("",values = c(0.1, 0.3, 1)) +
  #geom_line(data = sim_same_par_traj_inci, aes(x = ob_week, y = n, color = "Simulated Data")) +
  scale_color_manual("",values = c("Median" = "blue", "Simulated Data" = "black")) +
  scale_x_continuous(limits = c(x_min, x_max + 1)) +
  guides(
    fill = guide_legend(title = "", override.aes = list(alpha = c(0.3, 0.1), size = 4)),
    color = guide_legend(title = NULL)) +
  ggplot2::theme(
    
    plot.title = ggplot2::element_text(family = font, 
                                       size = 28, face = "bold", color = "#222222"), 
    plot.subtitle = ggplot2::element_text(family = font, size = 22, margin = ggplot2::margin(9, 0, 9, 0)), 
    plot.caption = ggplot2::element_blank(), 
    legend.position = "top", 
    legend.text = ggplot2::element_text(family = font, size = 18), 
    axis.title = element_text(family = font, size = 18),
    axis.text = ggplot2::element_text(family = font, size = 18, 
                                      color = "#222222"), 
    axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5, 
                                                                 b = 10)), 
    axis.ticks = ggplot2::element_blank(), 
    axis.line = ggplot2::element_blank(), 
    panel.grid.minor = ggplot2::element_blank(), 
    panel.grid.major.y = ggplot2::element_line(color = "#cbcbcb"), 
    panel.grid.major.x = ggplot2::element_blank(), 
    panel.background = ggplot2::element_blank(), 
    strip.background = ggplot2::element_rect(fill = "white"), 
    strip.text = ggplot2::element_text(size = 22, hjust = 0)) +
  labs(x = "Week", y = "Number of newly ASF-positive hexagon") -> trj_inci;trj_inci





setwd((file.path(scen_folder_wd, "/Trajectories/Temporal/")))
jpeg(paste0("trajectory_inci_temporal.jpg"), width = 1000, height = "600", quality = 100)
trj_inci  
dev.off()



setwd((file.path(scen_folder_wd, "/Trajectories/Temporal/")))
pdf(file = paste0("trajectory_inci_temporal.pdf"),   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 6) 
trj_inci  
dev.off()





ggplot(data = dt_trajec[[6]], aes(x = ob_week)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_quantile_95, aes(y = n, fill = "95%"), alpha = 0.1) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_quantile_50, aes(y = n, fill = "50%"), alpha = 0.3) +
  stat_summary(geom = "line", fun = median, aes(y = n, color = "Median"), alpha = 1, size = 1) +
  scale_fill_manual("",values = c("50%" = "blue", "95%" = "blue")) +
  scale_alpha_manual("",values = c(0.1, 0.3, 1)) +
  #geom_line(data = sim_same_par_traj_inci, aes(x = ob_week, y = n, color = "Simulated Data")) +
  scale_color_manual("",values = c("Median" = "blue", "Simulated Data" = "black")) +
  scale_x_continuous(limits = c(x_min, x_max + 1)) +
  guides(
    fill = guide_legend(title = "", override.aes = list(alpha = c(0.3, 0.1), size = 4)),
    color = guide_legend(title = NULL)) +
  ggplot2::theme(
    
    plot.title = ggplot2::element_text(family = font, 
                                       size = 28, face = "bold", color = "#222222"), 
    plot.subtitle = ggplot2::element_text(family = font, size = 22, margin = ggplot2::margin(9, 0, 9, 0)), 
    plot.caption = ggplot2::element_blank(), 
    legend.position = "top", 
    legend.text = ggplot2::element_text(family = font, size = 18), 
    axis.title = element_text(family = font, size = 18),
    axis.text = ggplot2::element_text(family = font, size = 18, 
                                      color = "#222222"), 
    axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5, 
                                                                 b = 10)), 
    axis.ticks = ggplot2::element_blank(), 
    axis.line = ggplot2::element_blank(), 
    panel.grid.minor = ggplot2::element_blank(), 
    panel.grid.major.y = ggplot2::element_line(color = "#cbcbcb"), 
    panel.grid.major.x = ggplot2::element_blank(), 
    panel.background = ggplot2::element_blank(), 
    strip.background = ggplot2::element_rect(fill = "white"), 
    strip.text = ggplot2::element_text(size = 22, hjust = 0)) +
  labs(x = "Week", y = "Number of newly ASF-positive hexagon (undetected)") -> trj_inci_undetected;trj_inci_undetected


setwd((file.path(scen_folder_wd, "/Trajectories/Temporal/")))
jpeg(paste0("trajectory_inci_undetected_temporal.jpg"), width = 1000, height = "600", quality = 100)
trj_inci_undetected  
dev.off()



setwd((file.path(scen_folder_wd, "/Trajectories/Temporal/")))
pdf(file = paste0("trajectory_inci_undetected_temporal.pdf"),   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 6) 
trj_inci_undetected  
dev.off()




ggplot(data = dt_trajec[[7]], aes(x = ob_week)) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_quantile_95, aes(y = n, fill = "95%"), alpha = 0.1) +
  stat_summary(geom = "ribbon", fun.data = mean_cl_quantile_50, aes(y = n, fill = "50%"), alpha = 0.3) +
  stat_summary(geom = "line", fun = median, aes(y = n, color = "Median"), alpha = 1, size = 1) +
  scale_fill_manual("",values = c("50%" = "blue", "95%" = "blue")) +
  scale_alpha_manual("",values = c(0.1, 0.3, 1)) +
  #geom_line(data = sim_same_par_traj_inci, aes(x = ob_week, y = n, color = "Simulated Data")) +
  scale_color_manual("",values = c("Median" = "blue", "Simulated Data" = "black")) +
  scale_x_continuous(limits = c(x_min, x_max + 1)) +
  guides(
    fill = guide_legend(title = "", override.aes = list(alpha = c(0.3, 0.1), size = 4)),
    color = guide_legend(title = NULL)) +
  ggplot2::theme(
    
    plot.title = ggplot2::element_text(family = font, 
                                       size = 28, face = "bold", color = "#222222"), 
    plot.subtitle = ggplot2::element_text(family = font, size = 22, margin = ggplot2::margin(9, 0, 9, 0)), 
    plot.caption = ggplot2::element_blank(), 
    legend.position = "top", 
    legend.text = ggplot2::element_text(family = font, size = 18), 
    axis.title = element_text(family = font, size = 18),
    axis.text = ggplot2::element_text(family = font, size = 18, 
                                      color = "#222222"), 
    axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5, 
                                                                 b = 10)), 
    axis.ticks = ggplot2::element_blank(), 
    axis.line = ggplot2::element_blank(), 
    panel.grid.minor = ggplot2::element_blank(), 
    panel.grid.major.y = ggplot2::element_line(color = "#cbcbcb"), 
    panel.grid.major.x = ggplot2::element_blank(), 
    panel.background = ggplot2::element_blank(), 
    strip.background = ggplot2::element_rect(fill = "white"), 
    strip.text = ggplot2::element_text(size = 22, hjust = 0)) +
  labs(x = "Week", y = "Number of newly ASF-positive hexagon (to be detected)") -> trj_inci_to_be_detected;trj_inci_to_be_detected


setwd((file.path(scen_folder_wd, "/Trajectories/Temporal/")))
jpeg(paste0("trajectory_inci_to_be_detected_temporal.jpg"), width = 1000, height = "600", quality = 100)
trj_inci_to_be_detected  
dev.off()



setwd((file.path(scen_folder_wd, "/Trajectories/Temporal/")))
pdf(file = paste0("trajectory_inci_to_be_detected_temporal.pdf"),   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 6) 
trj_inci_to_be_detected  
dev.off()





list_trajec %>% apply(MARGIN = 1, FUN = function(x){
  
  x[[8]] %>% 
    dplyr::select(unit, n) %>%
    complete(unit = 1:N, fill = list(n = 0)) %>% 
    dplyr::select(n)
  
  
  
}) %>%
  bind_cols() %>% 
  mutate(unit = 1:N) %>% 
  relocate(unit) %>% 
  apply(MARGIN = 1, FUN = function(x){
    
    
    x[-1] %>% quantile(probs = c(0.025, 0.5, 0.975))
    
    
  }) -> dist_num_reinfect_per_unit

dist_num_reinfect_per_unit %>% t %>% as_tibble %>% dplyr::select(`50%`) %>% pull(`50%`) %>% hist

dist_num_reinfect_per_unit %>% t %>% as_tibble %>%
  geom_violin(aes(y = "50%", x = ""), fill = "tomato") +
  geom_boxplot(aes(y = "50%", x = ""), width = 0.2) +
  scale_y_continuous(breaks = seq(0, 30, by = 10), limits = c(0, 30))




list_trajec %>% apply(MARGIN = 1, FUN = function(x){
  
  x[[8]] %>% summarise(num = n(), 
                       n_reinfect = sum(num_reinfect >= 1), 
                       perc = n_reinfect/num) %>% return
  
  
  
  
}) %>% bind_rows(.id = "sim") -> dt_reinfect




dt_reinfect



font <- "Helvetica"



dt_reinfect %>% ggplot() +
  geom_violin(aes(y = num, x = ""), fill = "tomato") +
  geom_boxplot(aes(y = num, x = ""), width = 0.2) +
  scale_y_continuous(breaks = seq(0, 1500, by = 500), limits = c(0, 1500)) +
  ggplot2::theme(plot.title = ggplot2::element_text(family = font, 
                                                    size = 12, face = "bold", color = "#222222"), plot.subtitle = ggplot2::element_text(family = font, 
                                                                                                                                        size = 22, margin = ggplot2::margin(9, 0, 9, 0)), plot.caption = ggplot2::element_blank(), 
                 legend.position = "top", legend.text.align = 0, legend.background = ggplot2::element_blank(), 
                 legend.title = ggplot2::element_blank(), legend.key = ggplot2::element_blank(), 
                 legend.text = ggplot2::element_text(family = font, size = 18, 
                                                     color = "#222222"), axis.title = ggplot2::element_blank(), 
                 axis.text = ggplot2::element_text(family = font, size = 18, 
                                                   color = "#222222"), axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5, 
                                                                                                                                    b = 10)), axis.ticks = ggplot2::element_blank(), 
                 axis.line = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), 
                 panel.grid.major.y = ggplot2::element_line(color = "#cbcbcb"), 
                 panel.grid.major.x = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                 strip.background = ggplot2::element_rect(fill = "white"), 
                 strip.text = ggplot2::element_text(size = 22, hjust = 0)) +
  ggtitle("The number of detected hexagons") -> num_infect;num_infect





dt_reinfect %>% ggplot() +
  geom_violin(aes(y = perc, x = ""), fill = "cornflowerblue") +
  geom_boxplot(aes(y = perc, x = ""), width = 0.2) +
  scale_y_continuous(breaks = seq(0,1,by = 0.2), limits = c(0, 1)) +
  ggplot2::theme(plot.title = ggplot2::element_text(family = font, 
                                                    size = 12, face = "bold", color = "#222222"), plot.subtitle = ggplot2::element_text(family = font, 
                                                                                                                                        size = 22, margin = ggplot2::margin(9, 0, 9, 0)), plot.caption = ggplot2::element_blank(), 
                 legend.position = "top", legend.text.align = 0, legend.background = ggplot2::element_blank(), 
                 legend.title = ggplot2::element_blank(), legend.key = ggplot2::element_blank(), 
                 legend.text = ggplot2::element_text(family = font, size = 18, 
                                                     color = "#222222"), axis.title = ggplot2::element_blank(), 
                 axis.text = ggplot2::element_text(family = font, size = 18, 
                                                   color = "#222222"), axis.text.x = ggplot2::element_text(margin = ggplot2::margin(5, 
                                                                                                                                    b = 10)), axis.ticks = ggplot2::element_blank(), 
                 axis.line = ggplot2::element_blank(), panel.grid.minor = ggplot2::element_blank(), 
                 panel.grid.major.y = ggplot2::element_line(color = "#cbcbcb"), 
                 panel.grid.major.x = ggplot2::element_blank(), panel.background = ggplot2::element_blank(), 
                 strip.background = ggplot2::element_rect(fill = "white"), 
                 strip.text = ggplot2::element_text(size = 22, hjust = 0)) +
  ggtitle("The proportion of the affected hexagons with two or more infected periods") -> prop_reinfect;prop_reinfect


setwd((file.path(scen_folder_wd, "/Trajectories/")))
pdf(file = paste0("proprotions_re_detected.pdf"),   # The directory you want to save the file in
    width = 8, # The width of the plot in inches
    height = 6) 
print(prop_reinfect)
dev.off()

setwd((file.path(scen_folder_wd, "/Trajectories/")))
jpeg(paste0("proprotions_re_detected.jpg"), width = 600, height = "600", quality = 100)
print(prop_reinfect)
dev.off()


























# Spatial risk and agreement
{
  
  # Spatial risk
  
  
  
  dt_trajec[[3]] %>%
    group_by(unit) %>% 
    summarise(occurrence = sum(value)) %>%
    mutate(binary = case_when(occurrence != 0 ~ 1,
                              occurrence ==0 ~ 0)) -> spatial_extent_ob
  
  dt_trajec[[3]] %>% 
    group_by(posterior) %>%
    filter(value == 3) %>% 
    dplyr::select(unit, posterior) %>% 
    group_map(.f = ~.x[,1] %>% pull) -> id_unit_report
  
  
  
  
  # Spatial risk
  dt_trajec[[4]] %>%
    group_by(unit) %>% 
    summarise(occurrence = sum(value)) %>%
    mutate(binary = case_when(occurrence != 0 ~ 1,
                              occurrence ==0 ~ 0)) -> spatial_extent_inci
  
  dt_trajec[[4]] %>% 
    group_by(posterior) %>%
    filter(value %in% c(1,2)) %>% 
    dplyr::select(unit, posterior) %>% 
    group_map(.f = ~.x[,1] %>% pull)  -> id_unit_inci
  
  
  
  
  unit_ID_ob<- readRDS("/data/unit_ID_ob.rds")
  
  
  
  ## Function to calculate spatial agreement
  sp_agree = function(id_unit_report, unit_ID_ob){
    
    ele_inter = intersect(id_unit_report, unit_ID_ob)
    ele_uni = union(id_unit_report, unit_ID_ob)
    sp_agree = length(ele_inter) / length(ele_uni)
    return(sp_agree)
    
  }
  
  
  ## Applying the function to whole element of the list
  lapply(X = id_unit_report, FUN = function(x, unit_ID_ob){
    
    return(sp_agree(id_unit_report = x, unit_ID_ob))
    
  }, unit_ID_ob = unit_ID_ob) %>% 
    unlist() -> sp_agree_dist
  
  quantile(sp_agree_dist, prob = c( 0.025, 0.5,0.975)) -> sp_agree_ob
  
  
  
}
setwd((file.path(scen_folder_wd, "/Diagnostic")))
rbind(sp_agree_ob ) %>% as.data.frame %>% rownames_to_column %>% write_csv("spatial_agreement.csv")

## For plotting, load hexagon file

hex <- readRDS("/data/hexagon.RDS")
library(sf)
hex %>% st_as_sf ->sf_hex




library(tmap)

# 

# Calculating the risk of ASF reports

## Calculating the risk of ASF reports
dt_trajec[[3]] %>%
  group_by(unit) %>%
  summarise(numb = sum(value ==3),
            total_simu = times_permute,
            risk = numb/total_simu, 
            
  ) %>% dplyr::select(risk) -> risk_ob_hexa


sf_hex$risk_ob = risk_ob_hexa[[1]]
sf_hex[I0[[1]],"risk_ob"] = 1

sf_hex[unit_ID_ob,] ->ob_sf_hex



setwd((file.path(scen_folder_wd, "/Trajectories/Spatial/")))
jpeg(paste0("Risk_ob_hex.jpeg"), width = 800, height = "800", quality = 100)
tm_shape(sf_hex) +
  tm_fill("risk_ob") + 
  tm_shape(ob_sf_hex) +
  tm_borders()
dev.off()


pdf(file = paste0("Risk_ob_hex.pdf"),   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 6) 
tm_shape(sf_hex) +
  tm_fill("risk_ob") + 
  tm_shape(ob_sf_hex) +
  tm_borders()
dev.off()





## Calculating the risk of ASF occurrence
dt_trajec[[4]] %>%
  group_by(unit) %>%
  summarise(numb = sum(value %in% c(1,2)),
            total_simu = times_permute,
            risk = numb/total_simu, 
            
  ) %>% dplyr::select(risk) -> risk_inci_hexa


sf_hex$risk_inci = risk_inci_hexa[[1]]
sf_hex[I0[[1]],"risk_inci"] = 1




setwd((file.path(scen_folder_wd, "/Trajectories/Spatial/")))
jpeg(paste0("Risk_inci_hex.jpeg"), width = 800, height = "800", quality = 100)

tm_shape(sf_hex) +
  tm_fill("risk_inci") 
dev.off()
pdf(file = paste0("Risk_inci_hex.pdf"),   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 6) 

tm_shape(sf_hex) +
  tm_fill("risk_inci") 

dev.off()

sf_hex %>% dplyr::select(risk_ob, risk_inci) -> risk_sf_hex

setwd((file.path(scen_folder_wd, "/Trajectories/Spatial/")))
st_write(layer = "risk_sf_hex", obj = risk_sf_hex, dsn = ".", layer_options = "OVERWRITE=true", driver = "ESRI Shapefile", update = TRUE)


setwd((file.path(scen_folder_wd, "/Trajectories/Spatial/")))
st_write(layer = "extent_ob_hex", obj = ob_sf_hex, dsn = ".", layer_options = "OVERWRITE=true", driver = "ESRI Shapefile",update = TRUE)

library(bbplot)





setwd((file.path(scen_folder_wd, "/Diagnostic")))
jpeg(paste0("violin_predict_VS_data_obhex.jpeg"), width = 800, height = "800", quality = 100)
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
dev.off()


setwd((file.path(scen_folder_wd, "/Diagnostic")))
jpeg(paste0("violin_predict_VS_data_incihex.jpeg"), width = 800, height = "800", quality = 100)
risk_inci_hexa %>% 
  mutate(num = 1:n(), 
         binary_inci = case_when(num %in% unit_ID_inci ~ 1, 
                                 !num %in% unit_ID_inci ~ 0, )
  ) %>%
  ggplot() +
  geom_violin(aes(x = factor(binary_inci), y = risk, col = factor(binary_inci), fill = factor(binary_inci))) +
  bbc_style() +
  theme(plot.title = element_text(size = 15)) +
  scale_x_discrete(labels = c("Un-occurred", "Occured")) +
  ggtitle("Predicted risk of ASF occurrence in a hexagon (for underlying diz. dynamics)") +
  scale_y_continuous(limits = c(0, 1)) 
dev.off()


pdf(file = paste0("violin_predict_VS_data_obhex.pdf"),   # The directory you want to save the file in
    width = 10, # The width of the plot in inches
    height = 6) 
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
dev.off()


# 
# setwd((file.path(scen_folder_wd, "/Diagnostic")))
# jpeg(paste0("violin_predict_VS_data_incihex.jpeg"), width = 800, height = "800", quality = 100)
# risk_inci_hexa %>% 
#   mutate(num = 1:n(), 
#          binary_inci = case_when(num %in% unit_ID_inci ~ 1, 
#                                  !num %in% unit_ID_inci ~ 0, )
#   ) %>%
#   ggplot() +
#   geom_violin(aes(x = factor(binary_inci), y = risk, col = factor(binary_inci), fill = factor(binary_inci))) +
#   bbc_style() +
#   theme(plot.title = element_text(size = 15)) +
#   scale_x_discrete(labels = c("Un-occurred", "Occured")) +
#   ggtitle("Predicted risk of ASF occurrence in a hexagon (for underlying diz. dynamics)") +
#   scale_y_continuous(limits = c(0, 1)) 
# dev.off()

# 
# pdf(file = paste0("violin_predict_VS_data_incihex.pdf"),   # The directory you want to save the file in
#     width = 10, # The width of the plot in inches
#     height = 6) 
# risk_inci_hexa %>% 
#   mutate(num = 1:n(), 
#          binary_inci = case_when(num %in% unit_ID_inci ~ 1, 
#                                  !num %in% unit_ID_inci ~ 0, )
#   ) %>%
#   ggplot() +
#   geom_violin(aes(x = factor(binary_inci), y = risk, col = factor(binary_inci), fill = factor(binary_inci))) +
#   bbc_style() +
#   theme(plot.title = element_text(size = 15)) +
#   scale_x_discrete(labels = c("Un-occurred", "Occured")) +
#   ggtitle("Predicted risk of ASF occurrence in a hexagon (for underlying diz. dynamics)") +
#   scale_y_continuous(limits = c(0, 1)) 
# dev.off()




## Function to calculate spatial sensitivity and specificity
spatial_sens_spec = function(trajectory_ID_inci, data_ID_inci){
  data_ID_inci = unique(data_ID_inci)
  pos = length(data_ID_inci)
  v_true_pos = intersect(data_ID_inci, trajectory_ID_inci)
  true_pos = length(v_true_pos)
  sens_n = true_pos/pos
  
  
  v_negative = setdiff(1:N, data_ID_inci)
  neg_n = length(v_negative)
  
  v_true_neg = setdiff(v_negative, trajectory_ID_inci)
  true_neg_n = length(v_true_neg)
  
  spec_n= true_neg_n/neg_n
  
  as.data.frame(cbind(sens_n, spec_n)) -> results
  
  return( results)
  
}



lapply(id_unit_report, FUN = spatial_sens_spec, data_ID_inci = unit_ID_ob) %>%
  do.call(what = rbind.data.frame) %>%
  apply(MARGIN = 2, FUN = quantile, probs = c(0.025, 0.5, 0.975)) %>%
  t() %>% as.data.frame %>% rownames_to_column()-> sens_spec_summary




setwd((file.path(scen_folder_wd, "/Diagnostic")))
write_csv(sens_spec_summary, "sens_spec_summary_dtf.csv")


risk_ob_hexa %>% 
  mutate(num = 1:n(), 
         binary = case_when(num %in% unit_ID_ob ~ 1, 
                            !num %in% unit_ID_ob ~ 0, )
  ) %>% group_by(binary) %>% 
  summarise(lo_95 = quantile(risk, probs = 0.025),
            med = median(risk), 
            up_95 = quantile(risk, probs = 0.975)) -> predicted_risk_ob


risk_inci_hexa %>% 
  mutate(num = 1:n(), 
         binary = case_when(num %in% unit_ID_inci ~ 1, 
                            !num %in% unit_ID_inci ~ 0, )
  ) %>% group_by(binary) %>% 
  summarise(lo_95 = quantile(risk, probs = 0.025),
            med = median(risk), 
            up_95 = quantile(risk, probs = 0.975)) -> predicted_risk_inci


bind_rows(predicted_risk_ob, predicted_risk_inci) %>% bind_cols(id = data.frame(id = rep(c("ob", "inci"), each = 2)) ) -> dtf_predicted_risk

setwd((file.path(scen_folder_wd, "/Diagnostic")))
write_csv(dtf_predicted_risk, "predicted_risk_dtf.csv")


cbind(paste0(round(sens_spec_summary, 2)[1,3] , "(",round(sens_spec_summary, 2)[1,2], " - ", round(sens_spec_summary, 2)[1,4],")"),
      paste0(round(sens_spec_summary, 2)[2,3] , "(",round(sens_spec_summary, 2)[2,2], " - ", round(sens_spec_summary, 2)[2,4],")")) %>% as.data.frame -> format_sens_spec_summary


setwd((file.path(scen_folder_wd, "/Diagnostic")))
write_csv(sens_spec_summary, "sens_spec_summary.csv")


