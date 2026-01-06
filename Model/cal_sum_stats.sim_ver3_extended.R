calc.sum.stats.sim <- function(matsum, 
                               n.breaks, 
                               sum_stats_output = c("incidence", "forest", "dist", "weekly" ,"ele"), 
                               sum_stats_spatial = c("county", "user_defined"), 
                               ob_indicator = 3) {
  
  if(0 == matsum %>% 
     filter(value == ob_indicator) %>%
     nrow) {
    
    matsum -> ob.matsum
    
  } else {
    matsum %>% 
      filter(value == 3) %>%
      arrange(week) %>%
      summarise(ob_week = min(week)) %>%
      pull ->ob_start
    
    matsum %>% filter(week >= ob_start &
                        week <= end_week) -> ob.matsum
    
  }
  sum.stats.ob.init <- ob.matsum %>%
    filter(value == ob_indicator) %>%
    mutate(div = cut(
      week,
      breaks = (0:n.breaks) * ceiling(end_week / n.breaks),
      right = TRUE)
    ) %>%
    mutate(div.n = as.integer(div))
  
  # sum.stats <- sum.stats %>%
  #   left_join(unitData, by = c("unit" = "unit.id")) %>%
  #   # dplyr::select(unit, div, div.n, county, epi_unit) %>%
  #   # group_by(county, epi_unit, div.n) %>%
  #   dplyr::select(unit, div, div.n, county, host) %>%
  #   group_by(county, host, div.n) %>%
  #   count(name="inc") %>%
  #   ungroup() %>%
  #   complete(
  #     county = unique(unitData$county),
  #     host = unique(unitData$host),
  #     # epi_unit = unique(unitData$epi_unit),
  #     # div = factor(levels(sum.stats$div)),
  #     div.n = 1:n.breaks,
  #     fill = list(inc = 0)
  #   ) %>%
  #   mutate(div.n = factor(div.n)) %>%
  #   # arrange(epi_unit, county, div.n) %>%
  #   arrange(host, county, div.n) %>%
  #   group_by(host, county) %>%
  #   # group_by(epi_unit, county) %>%
  #   mutate(cum.inc = cumsum(inc)) %>%
  #   mutate(
  #     sum.stat = case_when(
  #       host == "domestic pig" ~ inc,  # to set if sum stat is inc or cum.inc for each epi_unit/host
  #       host == "wild boar" ~ inc
  #     )
  #   ) %>% 
  #   # mutate(
  #   #   sum.stat = case_when(
  #   #     epi_unit == "village" | epi_unit == "industrial" ~ inc,
  #   #     epi_unit == "cell" ~ inc
  #   #   )
  #   # ) %>%
  #   ungroup()
  
  tot_sum.stats = NA
  
  if(c("incidence") %in% sum_stats_output){
    
    if(sum_stats_spatial == "county"){
      
      
      sum.stats.ob <- sum.stats.ob.init %>%
        left_join(unitData, by = c("unit" = "unit.id")) %>%
        # dplyr::select(unit, div, div.n, county, epi_unit) %>%
        # group_by(county, epi_unit, div.n) %>%
        #dplyr::select(unit, div.n, county) %>%
        dplyr::select(div.n, county) %>%
        group_by_all() %>%
        dplyr::count(name="inc") %>%
        ungroup() %>%
        dplyr::rename(region = county) %>%
        complete(
          region = unique(unitData %>% dplyr::select(county)) %>% pull,
          # epi_unit = unique(unitData$epi_unit),
          # div = factor(levels(sum.stats$div)),
          div.n = 1:n.breaks,
          fill = list(inc = 0)
        ) %>%
        mutate(div.n = factor(div.n)) %>%
        # arrange(epi_unit, county, div.n) %>%
        arrange(region, div.n) %>%
        # group_by(epi_unit, county) %>%
        # mutate(cum.inc = cumsum(inc)) %>%
        mutate(
          sum.stat = inc
        ) %>% 
        # mutate(
        #   sum.stat = case_when(
        #     epi_unit == "village" | epi_unit == "industrial" ~ inc,
        #     epi_unit == "cell" ~ inc
        #   )
        # ) %>%
        ungroup()
      
      
      
      
      sum.stats <- sum.stats.ob
      
      overall <- sum.stats %>%
        group_by(div.n) %>% 
        dplyr::summarise(dplyr::across(where(is.numeric), sum), .groups="drop") %>% 
        mutate(region="All") %>% 
        relocate(region, .before=div.n)
      
      
      sum.stats_incidence <- rbind(overall, sum.stats)
      
      tot_sum.stats <- sum.stats_incidence$sum.stat
      
      
    } 
    if(sum_stats_spatial == "user_defined") {
      sum.stats.ob <- sum.stats.ob.init %>%
        left_join(unitData, by = c("unit" = "unit.id")) %>%
        # dplyr::select(unit, div, div.n, county, epi_unit) %>%
        # group_by(county, epi_unit, div.n) %>%
        #dplyr::select(unit, div.n, county) %>%
        dplyr::select(div.n, user_defined) %>%
        group_by_all() %>%
        dplyr::count(name="inc") %>%
        ungroup() %>%
        dplyr::rename(region = user_defined) %>%
        complete(
          region = unique(unitData %>% dplyr::select(user_defined)) %>% pull,
          # epi_unit = unique(unitData$epi_unit),
          # div = factor(levels(sum.stats$div)),
          div.n = 1:n.breaks,
          fill = list(inc = 0)
        ) %>%
        mutate(div.n = factor(div.n)) %>%
        # arrange(epi_unit, county, div.n) %>%
        arrange(region, div.n) %>%
        # group_by(epi_unit, county) %>%
        # mutate(cum.inc = cumsum(inc)) %>%
        mutate(
          sum.stat = inc
        ) %>% 
        # mutate(
        #   sum.stat = case_when(
        #     epi_unit == "village" | epi_unit == "industrial" ~ inc,
        #     epi_unit == "cell" ~ inc
        #   )
        # ) %>%
        ungroup()
      
      
      
      
      sum.stats <- sum.stats.ob
      
      overall <- sum.stats %>%
        group_by(div.n) %>% 
        dplyr::summarise(dplyr::across(where(is.numeric), sum), .groups="drop") %>% 
        mutate(region="All") %>% 
        relocate(region, .before=div.n)
      
      
      sum.stats_incidence <- rbind(overall, sum.stats)
      
      tot_sum.stats <- sum.stats_incidence$sum.stat
      
      
      
      
      
      
      
      
      
      
      
    }
    
    
    
    
    
    
    
  }
  
  
  
  
  if(c("forest") %in% sum_stats_output){
    
    
    if(sum_stats_spatial == "county"){
      # 
      # sum_forest =  sum.stats.ob.init %>% 
      #   left_join(unitData, by = c("unit" = "unit.id")) %>%
      #   dplyr::select(div.n, county, for_cov) %>% 
      #   dplyr::rename(region = county) %>%
      #   complete(div.n = 1:n.breaks,
      #            region = unique(unitData %>% dplyr::select(county)) %>% pull,
      #            fill = list(for_cov = 0)) %>%
      #   group_by(div.n, region) %>%
      #   dplyr::summarise(sum_for_cov = sum(for_cov)) %>% 
      #   arrange(region, div.n)
      
      sum_forest_global =  sum.stats.ob.init %>% 
        left_join(unitData, by = c("unit" = "unit.id")) %>%
        dplyr::select(div.n, for_cov) %>% 
        complete(div.n = 1:n.breaks,
                 fill = list(for_cov = 0)) %>%
        group_by(div.n) %>%
        dplyr::summarise(sum_for_cov = sum(for_cov)) %>%
        mutate(region = "All", .before = sum_for_cov )
      
      sum.stats_forest <-rbind(sum_forest_global )
      
      
    } 
    if(sum_stats_spatial == "user_defined"){
      
      
      # 
      # sum_forest =  sum.stats.ob.init %>% 
      #   left_join(unitData, by = c("unit" = "unit.id")) %>%
      #   dplyr::select(div.n, user_defined, for_cov) %>% 
      #   dplyr::rename(region = user_defined) %>%
      #   complete(div.n = 1:n.breaks,
      #            region = unique(unitData %>% dplyr::select(user_defined)) %>% pull,
      #            fill = list(for_cov = 0)) %>%
      #   group_by(div.n, region) %>%
      #   dplyr::summarise(sum_for_cov = sum(for_cov)) %>% 
      #   arrange(region, div.n)
      
      sum_forest_global =  sum.stats.ob.init %>% 
        left_join(unitData, by = c("unit" = "unit.id")) %>%
        dplyr::select(div.n, for_cov) %>% 
        complete(div.n = 1:n.breaks,
                 fill = list(for_cov = 0)) %>%
        group_by(div.n) %>%
        dplyr::summarise(sum_for_cov = mean(for_cov)) %>%
        mutate(region = "All", .before = sum_for_cov )
      
      sum.stats_forest <-rbind(sum_forest_global )
      
      
      
      
    }
    
  }
  
  
  
  ob.matsum %>%
    filter(value == 3) %>%
    group_by(week) %>%
    summarise(newly_ob = n()) %>%
    complete(week = 1:end_week, fill = list(newly_ob = 0) ) %>% pull(newly_ob) -> v_newly_ob
  
  
  
  dist_m <- readRDS("/save_home/jslim//Mechanistic_modelling/Marryland/data/data/dist_m.rds")
  
  sum.stats.ob.init %>%
    left_join(unitData, by = c("unit" = "unit.id")) %>% 
    dplyr::select(unit, div.n, user_defined) %>%
    complete(div.n = 1:n.breaks, user_defined = unique(unitData %>% dplyr::select(user_defined)) %>% pull) %>%
    group_split(user_defined, div.n) %>%
    lapply(FUN = function(x){
      
      x %>% pull(unit) -> index
      
      dist_m[index,index] %>% sum / (2 * choose(k = 2, n = length(index))) -> avg_dist
      
      x$user_defined %>% unique ->region_ind
      x$div.n %>% unique ->div_ind
      
      return(c(div_ind, region_ind, avg_dist))
      
    }) %>% do.call(what = rbind.data.frame) -> dtf_dist
  
  names(dtf_dist) = c( "div.n", "region","avg_dist")
  
  dtf_dist %>% mutate(avg_dist = case_when(is.na(avg_dist) == TRUE ~ 0, 
                                           TRUE ~ avg_dist)) -> avg_dist_local
  
  avg_dist_local %>% group_by(div.n) %>% summarise(avg_dist= sum(avg_dist)) %>%
    mutate(region = "All") %>% 
    relocate(region, .before=avg_dist) -> avg_dist_global
  
  
  sum.stats_avg_dist<- rbind(avg_dist_global,avg_dist_local)
  
  
  
  
  
  sum_ele_global =  sum.stats.ob.init %>% 
    left_join(unitData, by = c("unit" = "unit.id")) %>%
    dplyr::select(div.n, ele) %>% 
    complete(div.n = 1:n.breaks,
             fill = list(ele = 0)) %>%
    group_by(div.n) %>%
    dplyr::summarise(avg_ele = mean(ele)) %>%
    mutate(region = "All", .before = avg_ele )
  
  sum.stats_ele <-rbind(sum_ele_global )
  
  
  # 
  # 
  # sum.stats.ob.init %>% 
  #   group_by(unit) %>%
  #   summarise(n = n()) %>% 
  #   group_by(n) %>% 
  #   summarise(n = n())-> n_reports
  # c(sum(n_reports[1,1]), sum(n_reports[-1,]))
  
  
  
  
  if(c("dist") %in% sum_stats_output){
    
    
    tot_sum.stats <- c(tot_sum.stats, sum.stats_avg_dist$avg_dist)  
    
  } 
  
  if(c("forest") %in% sum_stats_output){
    
    
    tot_sum.stats <- c(tot_sum.stats,sum.stats_forest$sum_for_cov )  
    
    
  } 
  
  
  if(c("weekly") %in% sum_stats_output){
    
    
    tot_sum.stats <- c(tot_sum.stats,v_newly_ob )  
    
    
  } 
  
  if(c("ele") %in% sum_stats_output){
    
    
    tot_sum.stats <- c(tot_sum.stats,sum.stats_ele$avg_ele )  
    
    
  } 
  
  
  return(tot_sum.stats)
  
  
  
}
