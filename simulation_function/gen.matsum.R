gen.matsum<- function(matr.states, I0) {
  I0_seed = I0
  
  matsum <-
    reshape2::melt(matr.states, varnames = c("unit", "week")) %>%
    arrange(unit, week) %>%
    group_by(unit) %>%
    filter(value != lag(value, order_by = week)) %>% # get only dt of state change
    ungroup() %>%
    add_row(unit = I0_seed[[1]],
            week = c(I0_seed[[2]]-1),
            value = 4) %>%  # add initial infected to assume infection by WB
    arrange(unit, week) %>%
    relocate(unit, week, value)
  
  
  matsum_undected_incidence <-
    reshape2::melt(matr.states, varnames = c("unit", "week")) %>%
    arrange(unit, week) %>%
    group_by(unit) %>%
    filter( lag(value, order_by = week) == 0 & value %in% c(2,3) | lag(value, order_by = week) %in% c(2,3) & value == 0) %>% # get only dt of state change
    filter(lead(value, order_by = week) == 0 & value %in% c(2,3) | lag(value, order_by = week) %in% c(2,3) & value == 0 ) %>%
    # add_row(unit = I0_seed[[1]],
    #         timestep = I0_seed[[2]] / dt,
    #         value = 2) %>%  # add initial infected to assume infection by WB
    arrange(unit, week) %>%
    relocate(unit, week, value)
  

  sim.timesteps = ncol(matr.states)
  sim.n.weeks = ceiling(sim.timesteps * dt)
  
  LIST.sim = list(matsum = matsum, 
                  matsum_undected_incidence = matsum_undected_incidence,
                  sim.timesteps = sim.timesteps, 
                  sim.n.weeks = sim.n.weeks)
  list2env(LIST.sim, .GlobalEnv)
  return(matsum)
}
