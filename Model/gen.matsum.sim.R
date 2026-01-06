gen.matsum.sim <- function(matr.states, I0) {
  library(tidyverse)
  I0_seed = I0
  
  matsum <-
    reshape2::melt(matr.states, varnames = c("unit", "timestep")) %>%
    arrange(unit, timestep) %>%
    group_by(unit) %>%
    filter(value != lag(value, order_by = timestep)) %>% # get only dt of state change
    ungroup() %>%
    add_row(unit = I0_seed[[1]],
            timestep = ((I0_seed[[2]]-1) / dt) + 1,
            value = 2) %>%  # add initial infected to assume infection by WB
    arrange(unit, timestep) %>%
    mutate(week = ceiling(timestep * dt)) %>%
    dplyr::select(-timestep) %>%
    relocate(unit, week, value)
  
  
  matsum_undected_incidence <-
    reshape2::melt(matr.states, varnames = c("unit", "timestep")) %>%
    arrange(unit, timestep) %>%
    group_by(unit) %>%
    filter( value == 1 & lag(value, order_by = timestep) == 0  ) %>% # get only dt of state change
    # add_row(unit = I0_seed[[1]],
    #         timestep = I0_seed[[2]] / dt,
    #         value = 2) %>%  # add initial infected to assume infection by WB
    arrange(unit, timestep) %>%
    mutate(week = ceiling(timestep * dt)) %>%
    dplyr::select(-timestep) %>%
    relocate(unit, week, value) %>%
    arrange(week, unit)
  
  
  reshape2::melt(matr.states, varnames = c("unit", "timestep")) %>%
    arrange(timestep, unit ) %>% summarise(max_tp = max(timestep)) %>% pull->max_tp

  matsum_undected_num <-
    reshape2::melt(matr.states, varnames = c("unit", "timestep")) %>%
    arrange(timestep, unit ) %>%
    group_by(timestep) %>%
    filter(value == 1) %>%
    summarise(n = n()) %>%
    mutate(week_raw = timestep * dt) %>%
    filter(week_raw == ceiling(week_raw)) %>%
    complete(fill = list(n = 0), week_raw = 1:ceiling(max_tp * dt)) %>%
    dplyr::rename(week = week_raw) %>%
    dplyr::select(week,n)
    
  
  matsum_dected_num <-
    reshape2::melt(matr.states, varnames = c("unit", "timestep")) %>%
    arrange(timestep, unit ) %>%
    group_by(timestep) %>%
    filter(value == 3) %>%
    summarise(n = n()) %>%
    mutate(week_raw = timestep * dt) %>%
    filter(week_raw == ceiling(week_raw)) %>%
    complete(fill = list(n = 0), week_raw = 1:ceiling(max_tp * dt)) %>%
    dplyr::rename(week = week_raw) %>%
    dplyr::select(week,n)
    
  sim.timesteps = ncol(matr.states)
  sim.n.weeks = ceiling(sim.timesteps * dt)
  
  LIST.sim = list(matsum = matsum, 
                  matsum_undected_incidence = matsum_undected_incidence,
                  sim.timesteps = sim.timesteps, 
                  sim.n.weeks = sim.n.weeks, 
                  matsum_undected_num = matsum_undected_num, 
                  matsum_detected_num = matsum_dected_num)
  list2env(LIST.sim, .GlobalEnv)
  return(matsum)
}
