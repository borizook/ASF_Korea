init.model <- function(end_week,
                       dt = 1/10,
                       n.breaks,
                       cell.adj.order=c("first","second")){
  ## Set model specifics ------------
  #### buffer or no buffer landscape
  mod.buff = "-nobuff"
  
  
  setwd("/save_home/jslim/Mechanistic_modelling/Marryland/data/data")
  
  I0_dtf <-read.csv("I0_dtf.csv")
  
  ## Derive unitData values -------------------------------------------------
  herds <- read_csv(paste0("herds",mod.buff,".csv"))
  unitData <- dplyr::select(herds, unit.id, dp, ind, cell_id, for_cov, n_vil, n_ind, edge, county, user_defined , epi_unit, host, ele)
  
  # Total number of epi units
  N <- nrow(unitData)
  # Total number of herds and their ids
  # dp.ids <- unitData$unit.id[unitData$dp==1]
  # n.herds <- length(dp.ids)
  # # Total number of industrial sites and their ids
  # ind.ids <- unitData$unit.id[unitData$ind==1]
  # n.ind <- length(ind.ids)
  # # Total number of vils and their ids
  # vil.ids <- unitData$unit.id[unitData$dp==1 & unitData$ind==0]
  # n.vil <- length(vil.ids)
  # # Total number of boar cells and their ids
  cell.ids <- unitData$unit.id[unitData$dp==0]
  n.cells <- length(cell.ids)
  # Set forest coverage threshold indicator (1 if forest coverage >= forThresh, else 0)
  #unitData$for_thresh=as.numeric((unitData$for_cov>=forThresh))
  
  
  ## Set model timing parameters -------------------------------------------------------------
  #  tStart - estimated number of weeks prior to first detection to serve as infection start point for simulation
  #  tStop - stop week of simulation (in weeks)
  #  dt - time step (in fraction of week)
  
  ## Build relationship matrices  ----------------------------------------------------
  # Distance matrix between herds
  # setwd("D:/onedrive/data/paper/2. Working/[2023]ASF_modelling/2. Korea_model/prep_data")
  # herdDist <- readRDS(paste0("herd.dist",mod.buff,".rds"))
  # colnames(herdDist) <- 1:ncol(herdDist)
  # 
  # Adjacency matrix of cells without names
  cell.adj.order=match.arg(arg=cell.adj.order)
  if (cell.adj.order=="first"){
    cell.adj=1
  }
  if (cell.adj.order=="second"){
    cell.adj=2
  }
  setwd("/save_home/jslim/Mechanistic_modelling/Marryland/data/data")
  cellAdj <- readRDS(paste0("cell.adj.", cell.adj, mod.buff,".rds"))
  rownames(cellAdj) <- seq((1), N)
  colnames(cellAdj) <- seq((1), N)
  # 
  # # Herd-cell container matrix
  # ### set cell-herd infection order
  # cell.herd.order=match.arg(arg=cell.herd.order)
  # if (cell.herd.order=="zero"){
  #   cell.herd=0
  # }
  # if (cell.herd.order=="first"){
  #   cell.herd=1
  # }
  # setwd("D:/onedrive/data/paper/2. Working/[2023]ASF_modelling/2. Korea_model/prep_data")
  # cellHerd <- readRDS(paste0("cell.herd.", cell.herd, mod.buff, ".rds"))
  # colSums(cellHerd)[1:10]
  # cellHerd[1:5,1:5]
  # cellHerd["2090",1:5]
  # to fix sim-model line 348, 363 (anywhere with unitData$cell_id)
  
  # unitData column for number of adjacent cells to each cell, calculated from cellAdj
  unitData$n_adj <- NA
  
  a=matrix(0,nrow(unitData),nrow(unitData))
  for (i in as.numeric(colnames(cellAdj))){
    a[as.numeric(colnames(cellAdj)[1]):as.numeric(tail(colnames(cellAdj),1)),i] =
      cellAdj[,as.character(unitData$cell_id[i])]
  }
  cellAdj<-a
  
  unitData$n_adj <- colSums(cellAdj,na.rm=TRUE)
  
  # matrix of adjacency pairs
  adjacency <- which(cellAdj==1, arr.ind = TRUE)
  
  # List of herds in each herd's surveillance zone
  # identify herds within each unit's sz
  # herdSZ <- herdDist < sz 
  # # create list of herds with herds within sz
  # herdSZ <- apply(herdSZ, 1, function(x) names(which(x)))
  # # set list values as numerics
  # herdSZ <- lapply(herdSZ, function(x) as.numeric(x))
  # 
  ### number of herds within a herds' maxDist if using frequency-dependent transmission
  # herds.in.maxdist <- rowSums(herdDist < maxDist)
  # unitData$herds_in_maxdist <- 0
  # #unitData$herds_in_maxdist[1:n.herds] <- herds.in.maxdist
  
  ### include cells neighboring forested cells as having adequate forest coverage
  # if (for.neighbors==T){
  #   forested.cells <- unitData$cell_id[unitData$for_thresh==1]
  #   neighbors.to.forested.cells <- unique(adjacency[adjacency[,2] %in% forested.cells,][,1])
  #   forested.and.neighbors <- unique(c(forested.cells, neighbors.to.forested.cells))
  #   unitData$for_thresh[unitData$cell_id %in% forested.and.neighbors] = 1
  # }
  ## Calc observed sum stats ---------------------------------------------------------------
  ## melt SoI into df of unique outbreaks
  # setwd("/save_home/jslim/Mechanistic_modelling/Marryland/data/RDS")
  # obs.matsum <- readRDS("obs.soi-nobuff.rds") %>%
  #   reshape2::melt(varnames = c("unit", "week")) %>%
  #   dplyr::filter(value==6) %>%
  #   arrange(unit, week) %>%
  #   # mutate(week=week+(1/sigma.wb)) %>% # align observed data with model time frame
  #   group_by(unit) %>%
  #   complete(week = 0:end_week, fill = list(value=0)) %>%
  #   mutate(new.outbreak=case_when(value == 6 & value != lag(value) ~ 6)) %>%
  #   dplyr::select(-value) %>%
  #   filter(new.outbreak==6) %>% # df of new outbreaks
  #   dplyr::rename(value=new.outbreak) %>%
  #   filter(week!=0) # Due to the lag function in mutate (line 128), seed cases in week 1 was removed. Thus, the week starts at 0 and did muatate and removed week 0 data.
  # 
  # ## calculate sum stats function (if sourced from outside, objects in function not found, so must be duplicated here)
  # calc.sum.stats.ob <- function(matsum, n.breaks) {
  #   sum.stats <- matsum %>%
  #     filter(value == 6) %>%
  #     mutate(div = cut(
  #       week,
  #       breaks = (0:n.breaks) * ceiling(end_week / n.breaks),
  #       right = TRUE)
  #     ) %>%
  #     mutate(div.n = as.integer(div))
  #   
  #   # sum.stats <- sum.stats %>%
  #   #   left_join(unitData, by = c("unit" = "unit.id")) %>%
  #   #   # dplyr::select(unit, div, div.n, county, epi_unit) %>%
  #   #   # group_by(county, epi_unit, div.n) %>%
  #   #   dplyr::select(unit, div, div.n, county, host) %>%
  #   #   group_by(county, host, div.n) %>%
  #   #   count(name="inc") %>%
  #   #   ungroup() %>%
  #   #   complete(
  #   #     county = unique(unitData$county),
  #   #     host = unique(unitData$host),
  #   #     # epi_unit = unique(unitData$epi_unit),
  #   #     # div = factor(levels(sum.stats$div)),
  #   #     div.n = 1:n.breaks,
  #   #     fill = list(inc = 0)
  #   #   ) %>%
  #   #   mutate(div.n = factor(div.n)) %>%
  #   #   # arrange(epi_unit, county, div.n) %>%
  #   #   arrange(host, county, div.n) %>%
  #   #   group_by(host, county) %>%
  #   #   # group_by(epi_unit, county) %>%
  #   #   mutate(cum.inc = cumsum(inc)) %>%
  #   #   mutate(
  #   #     sum.stat = case_when(
  #   #       host == "domestic pig" ~ inc,  # to set if sum stat is inc or cum.inc for each epi_unit/host
  #   #       host == "wild boar" ~ inc
  #   #     )
  #   #   ) %>% 
  #   #   # mutate(
  #   #   #   sum.stat = case_when(
  #   #   #     epi_unit == "village" | epi_unit == "industrial" ~ inc,
  #   #   #     epi_unit == "cell" ~ inc
  #   #   #   )
  #   #   # ) %>%
  #   #   ungroup()
  #   
  #   sum.stats <- sum.stats %>%
  #     left_join(unitData, by = c("unit" = "unit.id")) %>%
  #     # dplyr::select(unit, div, div.n, county, epi_unit) %>%
  #     # group_by(county, epi_unit, div.n) %>%
  #     dplyr::select(unit, div.n, county) %>%
  #     group_by(county, div.n) %>%
  #     dplyr::count(name="inc") %>%
  #     ungroup() %>%
  #     complete(
  #       county = unique(unitData$county),
  #       # epi_unit = unique(unitData$epi_unit),
  #       # div = factor(levels(sum.stats$div)),
  #       # div.n = 0:n.breaks,
  #       div.n = 0:n.breaks,
  #       fill = list(inc = 0)
  #     ) %>%
  #     mutate(div.n = factor(div.n)) %>%
  #     # arrange(epi_unit, county, div.n) %>%
  #     arrange(county, div.n) %>%
  #     group_by(county) %>%
  #     # group_by(epi_unit, county) %>%
  #     mutate(cum.inc = cumsum(inc)) %>%
  #     mutate(
  #       sum.stat = inc
  #     ) %>% 
  #     # mutate(
  #     #   sum.stat = case_when(
  #     #     epi_unit == "village" | epi_unit == "industrial" ~ inc,
  #     #     epi_unit == "cell" ~ inc
  #     #   )
  #     # ) %>%
  #     ungroup()
  #   
  #   overall <- sum.stats %>%
  #     group_by(div.n) %>% 
  #     dplyr::summarise(dplyr::across(where(is.numeric), sum), .groups="drop") %>% 
  #     mutate(county="All") %>% 
  #     relocate(county, .before=div.n)
  #   
  #   sum.stats <- rbind(overall, sum.stats)
  #   return(sum.stats$sum.stat)
  #    
  # }
  #setwd("/save_home/jslim/Mechanistic_modelling/Marryland/code/parallel_computing/integrity_testing/")
  
  #sum.stats.obs <- readRDS("[20231012]artifical_data_model_2.RDS")
  # set I0 ------------------------------------------------------------------
  #I0 <- c(2989, 2990, 3068, 3028, 3186, 3068, 3107, 3030)[1:num_I0]
  I0 = I0_dtf
  # Fixed pars for sim ------------------------------------------------------
  # fixed.pars <- c(gamma.cell)
  
  ## List to send to global env ----------------------------------------------
  LIST = list(unitData=unitData, N=N, cell.ids=cell.ids, n.cells=n.cells,end_week = end_week,
              #forThresh=forThresh,  cellAdj=cellAdj, dt = dt,
              adjacency=adjacency, I0=I0, dt = dt)
  
  list2env(LIST, .GlobalEnv)
}