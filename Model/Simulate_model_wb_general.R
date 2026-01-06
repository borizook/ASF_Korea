SimulateModel<-function(end_week, 
                        N,
                        adjacency,
                        I0, 
                        essential_pars,
                        detect_pars, 
                        susc_infect_pars ,
                        unitData, 
                        dt, 
                        trans.mode = c("dens", "freq"), 
                        seasonality = c(FALSE, TRUE),
                        seasonality_beta,
                        forest_unit = c("percent", "decimal"), 
                        detect = c("constant", "varying"), 
                        susc_infect = c("single constant", 
                                        "different constant", 
                                        "single varying", 
                                        "different varying")
)
{
  # R function to run an SIuIdS model assuming p% of detection of unit and associated recovery rate.
  #
  #
  # Inputs ----
  #   I0 - unitID(s) of initial infected(s)
  #   pars - vector of estimated parameter values
  # essential_pars ---
  #  sigma - detection rate for wild boar cells (Iu_d to Id)
  #  gamma - recovery rate for wild boar cells (Id to S)
  #  eta - Rate from Iu to S
  # detect_pars ---
  #  When detect = "constant"
  #     constant value - detectibility ranging from 0 to 1 (when detect = "constant" )
  #  When detect = "varying"
  #     Maximum value - maximum detectibility ranging from 0 to 1 
  #     Slope - decreaing rate of detectibilty with the function of elevation 
  #
  # susc_infect_pars --- phi
  #  When susc_infect = "single constant"
  #     Phi - relative susceptibility and infectibility
  #  When susc_infect = "different constant"
  #     Phi - 
  #     Psi - 
  #  When susc_infect = "single varying"
  #  When susc_infect = "different varying"
  
  #  
  #   unitData - data frame with unit ID, forest percent cover per cell, host, cell, user defined county, number of adjacency cells,
  #   dt - time step (in fraction of week)
  #   trans.mode - set whether transmission mode from wb cells to wb cells is density-dependent or frequency-dependent
  #
  # Output ----
  #   matr.states - matrix of epi unit status per time step
  
  
  # Paramaters ---- depending on the model structures you choose (detect, and susc_infect argument)
  #   beta - transmission rate from cell to cell   #a x sin(2pi x (t-b)/55) + a + c
  #      beta[1]: amplitutde; 
  #      beta[2]: the first peak time in a time unit of this model, 
  #      beta[3]: minimum beta
  #      beta[4]: the number of cycle per year; e.g., when two peaks in a year, beta[4] should be 2
  #   phi_nonforested - relative susceptibility of cell with zero forest density (should be between 0 - 1)
  #   eta - Rate from Iu to S
  #   sens - surveillance sensitivity (Probability that Iu will be detected)
  # The number of parameters
  #   detect
  #     "constant" -> one / "varying" -> two
  #   susc_infect
  #     "single constant" -> one / "varying constant" -> two / "single varying" -> one / "different varying " -> two
  # 
  # essential = c("beta","sigma", "gamma", "eta")
  # matr.states number codes:
  # 0 = state S
  # 1 = State Iu_u
  # 2 = State Iu_d
  # 3 = State Id
  
  # Time steps  --------------------------------------------------------------
  # Set the time step
  
  {
    essential = c("sigma", "gamma", "eta")
    
    n_pars = length(essential)
    detect = match.arg(detect)
    susc_infect = match.arg(susc_infect)

    if(seasonality == FALSE){
      
      if(length(seasonality_beta) != 1){
        
        stop("Only single parameter for transmission rate should be provided")
        
        
      }
      
      n_pars <- n_pars + 1
      
    } else {
      
      
      
      if(length(seasonality_beta) != 4){
        
        stop("Four parameters for transmission rate should be provided")
        
        
      }
      
      
      n_pars <- n_pars + 4
      
    }
    
    
    
    
    
    
    if(detect == "constant"){
      
      
      if(length(detect_pars) != 1 ){
        
        
        stop("Only single parameter for detectibility should be provided")
        
        
        
      }
      
      
      n_pars <- n_pars + 1
      
    } else {
      
      if(length(detect_pars) != 2 ){
        
        
        stop("Two parameters for detectibility should be provided")
        
        
        
      }
      
      n_pars <- n_pars + 2
      
    }
    
    
    
    if(susc_infect == "single constant"){
      
      
      
      if(length(susc_infect_pars) != 1 ){
        
        
        stop("Only single parameter for susceptibility/infectibility should be provided")
        
        
        
      }
      
      
      n_pars <- n_pars + 1
      
    }
    if(susc_infect == "different constant"){
      
      
      
      if(length(susc_infect_pars) != 2 ){
        
        
        stop("Two parameters for susceptibility/infectibility should be provided")
        
        
        
      }
      
      n_pars <- n_pars + 2
      
    }
    if(susc_infect == "single varying"){
      
      
      if(length(susc_infect_pars) != 1 ){
        
        
        stop("Only single parameter for susceptibility/infectibility should be provided")
        
        
        
      }
      n_pars <- n_pars + 1
      
    }
    
    if(susc_infect == "different varying"){
      
      if(length(susc_infect_pars) != 2 ){
        
        
        stop("Two parameters for susceptibility/infectibility should be provided")
        
        
        
      }
      
      n_pars <- n_pars + 2
      
    }
    # 
    # 
    # if(redetection == FALSE){
    #   
    #   if(!is.null(redetection_par)){
    #     
    #     stop("Should be NULL for redetection parameter")
    #     
    #     
    #   }
    #   
    #   n_pars <- n_pars
    #   
    # } else {
    #   
    #   
    #   
    #   if(length(redetection_par) != 1){
    #     
    #     stop("Single parameter for re-detection should be provided")
    #     
    #     
    #   }
    #   
    #   
    #   n_pars <- n_pars + 1
    #   
    # }
    
    
  }
  
  dt=dt
  
  # Extract parameters  -----------------------------------------------------
  ## estimated pars -----
  beta = seasonality_beta
  sigma = essential_pars[1]
  gamma = essential_pars[2]
  eta = essential_pars[3]
  
  
  
  if(detect == "constant"){
    max_sens = detect_pars[1]
    adj_sens = rep(max_sens, N)
    
    
  } else {
    
    max_sens = detect_pars[1]
    slope_sens = detect_pars[2]
    # std of elevation ranging 0 ~ 1
    ele_std <- (unitData$ele - min(unitData$ele)) / (max(unitData$ele) - min(unitData$ele))
    sens_vari = 1 - slope_sens * ele_std
    adj_sens  = max_sens * sens_vari 
  }
  
  
  
  
  
  if(susc_infect == "single constant"){
    
    
    phi_nonforested = susc_infect_pars[1] # domain of this parameter should be between 0 and 1
    # Modifying the unit of forest coverage
    {
      forest_unit <- match.arg(arg = forest_unit)
      
      if (forest_unit == "percent") {
        
        perc_for  = (unitData$for_cov / 100)
        
        
      }
      
      if (forest_unit == "decimal") {
        
        perc_for  = unitData$for_cov
        
      }
    }
    
    phi <- rep(1, N)
    
    phi[perc_for<0.8] = phi_nonforested
    
    
  }
  
  
  if(susc_infect == "different constant"){
    
    
    phi_nonforested = susc_infect_pars[1] # domain of this parameter should be between 0 and 1
    psi_nonforested = susc_infect_pars[2] # domain of this parameter should be between 0 and 1    
    
    # Modifying the unit of forest coverage
    {
      forest_unit <- match.arg(arg = forest_unit)
      
      if (forest_unit == "percent") {
        
        perc_for  = (unitData$for_cov / 100)
        
        
      }
      
      if (forest_unit == "decimal") {
        
        perc_for  = unitData$for_cov
        
      }
    }
    
    phi <- rep(1, N)
    psi <- rep(1, N)
    
    phi[perc_for<0.8] = phi_nonforested
    psi[perc_for<0.8] = psi_nonforested
    
    
    
  }
  if(susc_infect == "single varying"){
    
    
    phi_nonforested = susc_infect_pars[1] # domain of this parameter should be between 0 and 1
    
    
    # Modifying the unit of forest coverage
    {
      forest_unit <- match.arg(arg = forest_unit)
      
      if (forest_unit == "percent") {
        
        perc_for  = (unitData$for_cov / 100)
        
        
      }
      
      if (forest_unit == "decimal") {
        
        perc_for  = unitData$for_cov
        
      }
    }
    
    # modifying the forest coverage into function of linear assumption
    phi =(1 - phi_nonforested) * perc_for + phi_nonforested
    
  }
  
  if(susc_infect == "different varying"){
    phi_nonforested = susc_infect_pars[1] # domain of this parameter should be between 0 and 1
    psi_nonforested = susc_infect_pars[2] # domain of this parameter should be between 0 and 1
    
    # Modifying the unit of forest coverage
    {
      forest_unit <- match.arg(arg = forest_unit)
      
      if (forest_unit == "percent") {
        
        perc_for  = (unitData$for_cov / 100)
        
        
      }
      
      if (forest_unit == "decimal") {
        
        perc_for  = unitData$for_cov
        
      }
    }
    
    phi =(1 - phi_nonforested) * perc_for + phi_nonforested
    psi =(1 - psi_nonforested) * perc_for + psi_nonforested
    
    
  }
  # 
  # 
  # if(redetection == TRUE) {
  #   
  #   rel_inf_reinfect = redetection_par
  #   
  # }
  # 
  
  
  
  tStart <- 0 
  tStop <- end_week + (1/sigma)
  n.timesteps <- (tStop/dt)-(tStart/dt)
  n.weeks <- n.timesteps*dt
  
  
  # Set the times
  tVals = seq(tStart+dt,tStop-dt,dt)
  
  
  # modifying the forest coverage into function of linear assumption
  
  # modifying the elevation to exponential distribution function
  # sens_vari =  sens * exp(-sens * ele_std)
  
  # Set vectors of unit-specific parameters 
  sigma <- matrix(sigma,N)
  eta <- matrix(eta, N)
  gamma <- matrix(gamma,N)
  sens_surv <- matrix(adj_sens)
  
  # Create object for in-model clock
  time_R=matrix(-1,N)
  
  # Initialize the populations -------------------------------------------
  # Host states
  
  I00 <- I0[I0[,2] == 1,1]
  matr.states <- matrix(0, N, length(tVals)+1)
  matr.states[I00,1] <- 2
  
  # Matrix of units within surveillance zone
  # within_sz <- matrix(0,N,length(tVals)+1)
  
  # vector of units to remove from surveillance
  # surv.remove <- NULL  
  
  time_step=1
  
  # matrix of lambdas for each timestep  
  lambda.table = matrix(0,N,ncol = 1) # col1=wb2wb
  colnames(lambda.table) = c( "wb2wb") 
  
  # For loop sim -------------------------------------------------------------
  for(stp in tVals){
    # reset lambdas
    lambda.table[,] = 0
    
    
    
    #lambda.table[1,1] = 0.5
    # set current model time
    time_step=time_step+1
    
    # copy data from previous time step
    matr.states[,time_step]=matr.states[,time_step-1]

    
    ## Seasonality of transmission rate
    if(length(beta) > 1 ) {
      
      real_time <- ceiling(time_step * dt)
      
      
      amp <- beta[1]
      minimum <- beta[3]
      num_cyc_year <- beta[4]
      
      h_shift <- beta[2]
        
      beta <- amp * sin(num_cyc_year * 2 * pi * (real_time - h_shift)/52 ) + minimum + amp
      
    }
    
    
    # introduction of seed cases
    pseed_case <-  ((c(I0[,2] - 1) / dt) + 1) %in% time_step
    
    
    if(any(pseed_case)){
      
      matr.states[I0[pseed_case, 1],time_step] <- 2
      
    }
    ## Generate unit state vectors for current time step --------------------
    # matr.states number codes:
    # 0 = state S
    # 1 = State Iu_u
    # 2 = State Iu_d
    # 3 = State Id
    
    # logical vector of all epi units that are infectious (detected or undetected)
    is_infectious_unit <- (matr.states[, time_step] %in% 1:3)
    
    # logical vector of all units in state_Iu and state_Id
    state_Iu_u <- (matr.states[, time_step] == 1)
    state_Iu_d <- (matr.states[, time_step] == 2)
    state_Id <- (matr.states[, time_step] == 3)
    
    ## vector of susceptible units (formerly state_S)
    is_unit_Susc = matr.states[,time_step]==0
    
    # vector of susceptible units with relative susceptibility 
    # 
    # if(redetection == TRUE) {
    #   re_inf = rep(1, N)
    #   
    #   if(time_step == 2){
    #     ind_re_inf = (matr.states[,1:(time_step-1)] %in% 1:2) > 0
    #     
    #   } else {
    #     
    #     ind_re_inf = matr.states[,1:(time_step-1)] %>% 
    #       apply(X = . , MARGIN = 1, FUN = function(x){
    #         
    #         sum(x %in% c(1:2)) > 0
    #       } 
    #       
    #       )
    #   }
    #   
    #   re_inf[ind_re_inf] = rel_inf_reinfect
    #   RelSusc_unit_Susc = phi*is_unit_Susc * re_inf
    #   
    # } else {
    #   
    #   RelSusc_unit_Susc = phi*is_unit_Susc 
    #   
    # }
    
    RelSusc_unit_Susc = phi*is_unit_Susc 
    
    if(grepl(susc_infect, pattern = "single")){
      
      
      psi <- phi    
      
    }
    
    
    
    
    # Calculate FoIs ----------------------------------------------------------
    if (any(is_infectious_unit)){
      indexes=adjacency[is_infectious_unit[adjacency[,1]] & is_unit_Susc[adjacency[,2]],]
      
      
      trans.mode <- match.arg(arg = trans.mode)
      
      if(is.matrix(indexes)){
        if(trans.mode == "freq"){
          if (any(duplicated(indexes[,2]))){
            tmp =
              RelSusc_unit_Susc[indexes[,2]] *
              (psi[indexes[,1]] *
                 beta / 
                 unitData$n_adj[indexes[,2]])
            tmp=by(tmp,indexes[,2],sum)
            lambda.table[sort(unique(indexes[,2])),1] = as.numeric(tmp)
          } else {
            tmp = 
              RelSusc_unit_Susc[indexes[,2]] *
              (psi[indexes[,1]] *
                 beta / 
                 unitData$n_adj[indexes[,2]])
            lambda.table[sort(unique(indexes[,2])),1] = as.numeric(tmp)
            
          }
        }
        
        if(trans.mode == "dens"){
          if (any(duplicated(indexes[,2]))){
            tmp =
              RelSusc_unit_Susc[indexes[,2]] *
              (psi[indexes[,1]] *
                 beta)
            tmp=by(tmp,indexes[,2],sum)
            lambda.table[sort(unique(indexes[,2])),1] = as.numeric(tmp)
          } else {
            tmp = 
              RelSusc_unit_Susc[indexes[,2]] *
              (psi[indexes[,1]] *
                 beta)
            lambda.table[sort(unique(indexes[,2])),1] = as.numeric(tmp)
          }
        }
      } else {
        if(trans.mode == "freq"){
          tmp = 
            RelSusc_unit_Susc[indexes[2]] *
            (psi[indexes[1]] *
               beta / 
               unitData$n_adj[indexes[2]])
          lambda.table[sort(unique(indexes[2])),1] = as.numeric(tmp)
        }
        if(trans.mode == "dens"){
          if (any(duplicated(indexes[2]))){
            tmp =
              RelSusc_unit_Susc[indexes[2]] *
              (psi[indexes[1]] *
                 beta)
            tmp=by(tmp,indexes[2],sum)
            lambda.table[sort(unique(indexes[2])),1] = as.numeric(tmp)
          } else {
            tmp = 
              RelSusc_unit_Susc[indexes[2]] *
              (psi[indexes[1]] *
                 beta)
            lambda.table[sort(unique(indexes[2])),1] = as.numeric(tmp)
          }
        }
      }
      
    }
    
    # Determine state transitions -------------------------------------------------------------------------
    ## S (0) to Iu (1) transition -------------------------------------------
    ## Define whether or not an infection occurs:
    if(any(is_unit_Susc)){
      # vector of units to potentially transition
      may.trans <- which(rowSums(lambda.table)!=0)
      
      if(length(may.trans) > 0){
        prob = rep(0, times = length(may.trans))
        prob = 1-(exp(-lambda.table[may.trans,]*dt))
        
        bi_infect <- (rbinom(n = length(may.trans),size = 1,prob=prob) ==1)
      } 
      
      if(any(bi_infect)){
        
        ind_infect <- may.trans[bi_infect]
        
        matr.states[ind_infect,time_step] <- 1 + rbinom(length(ind_infect), 1, prob = sens_surv[ind_infect])
        
      }
      
    }
    
    
    
    
    
    
    
    
    ### Iu to Id or S transition ------------------------------------------
      
      ## Iu_u (1) to S (0) transition -------------------------------------------
    if(any(state_Iu_u)){
      
      probs_recov <- (1-exp(-eta*dt))[state_Iu_u]
      
      matr.states[state_Iu_u,time_step] = matr.states[state_Iu_u,time_step] +
        (- matr.states[state_Iu_u,time_step]) * rbinom(n = sum(state_Iu_u),size = 1,prob=probs_recov)
      
    }
    
    
    
    
    
    
      ## Iu_d (2) to Id (3) transition -------------------------------------------
    if(any(state_Iu_d)){
      
      probs_det <- (1-exp(-sigma*dt))[state_Iu_d]
      
      matr.states[state_Iu_d,time_step] = matr.states[state_Iu_d,time_step] +
        (3 - matr.states[state_Iu_d,time_step]) * rbinom(n = sum(state_Iu_d),size = 1,prob=probs_det)
      
      
    }
    
    
    
    ## Id (3) to S (0) transition ------------------------------------------
    # Probability of Id to S transition
    # gamma
    if(any(state_Id)){
      prob_recov <- (1-exp(-gamma*dt))[state_Id]
      
      # removing NA
      
      # Sampling random values (0: Susceptible)
      # Instead of using if statement where the process is dependent on the random number from rbinom with probability of "prob_recov", 
      # This function used rbinom and simple calculation
      matr.states[state_Id,time_step]=matr.states[state_Id,time_step] +
        (-matr.states[state_Id,time_step]) * rbinom(sum(state_Id),1,prob=prob_recov)
    }
    # Close for-loop ----------------------------------------------------------
  }
  
  # Close function ----------------------------------------------------------------------------------------
  matr.states <<- matr.states
}






